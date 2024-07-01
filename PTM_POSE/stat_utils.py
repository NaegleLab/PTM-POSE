import scipy.stats as stats
import numpy as np
import pandas as pd

def adjust_single_p(p, num_tests, rank, prev_p):
    """
    Adjusts a single p-value using the Benjamini-Hochberg method.

    Parameters
    ----------
    p : float
        p-value to be adjusted.
    num_tests : int
        Number of tests performed.
    rank : int
        Rank of the p-value.
    prev_p : float
        Previous p-value.
    
    Returns
    -------
    adj_p : float
        Adjusted p-value.
    """
    adj_p = p*(num_tests/rank)
    if adj_p < prev_p:
        adj_p = prev_p
    elif adj_p > 1:
        adj_p = 1
    else:
        prev_p = adj_p
    return adj_p, prev_p
    

def adjustP(sorted_p, method = 'BH'):
    """
    Adjusts a list of p-values using the Benjamini-Hochberg method.

    Parameters
    ----------
    sorted_p : list
        List of p-values to be adjusted.
    method : string, optional
        Method to use for adjustment. Currently only supports Benjamini-Hochberg (BH) and Bonferroni (Bonf).

    Returns
    -------
    adj_p_list : list
        List of adjusted p-values.
    
    """
    adj_p_list = []
    prev_p = 0
    for i in range(len(sorted_p)):
        if method == 'BH':
            adj_p, prev_p = adjust_single_p(sorted_p[i], len(sorted_p), i+1, prev_p)
            adj_p_list.append(adj_p)
        elif method == 'Bonf':
            adj_p = sorted_p[i]*len(sorted_p)
            if adj_p > 1:
                adj_p = 1
            adj_p_list.append(adj_p)
    return adj_p_list
    
def hypergeom(M, n, N, k):
    """
    Calculates the hypergeometric p-value.

    Parameters
    ----------
    M : int
        Total number of ptms.
    n : int
        Total number of instances in the population.
    N : int
        Total number of ptms in the subset.
    k : int
        Number of instances in the subset  

    Returns
    -------
    p : float
        Hypergeometric p-value.

    """
    p = stats.hypergeom(M=M, 
                n = n, 
                N=N).sf(k-1)
    return p

def convertToFishers(M, n, N, x):
    """
    Given hypergeometric parameters, convert to a fishers exact table

    Parameters
    ----------
    M : int
        Total number of ptms.
    n : int
        Total number of instances in the population.
    N : int
        Total number of ptms in the subset.
    x : int
        Number of instances in the subset

    Returns 
    -------
    table : list
        List of lists containing the fishers exact table.
    """
    table = [[x, n-x],
            [N-x, M- (n+N) + x]]
    return table

def constructPivotTable(annotated_ptms, reference_col, database = 'PhosphoSitePlus', collapse_on_similar = False, include_unknown = False):
    """
    Given a ptm dataframe and regulatory data from phosphositeplus, create a table with PTMs in the rows and annotations in the columns, with 1 indicating that the PTM has that annotation

    Parameters
    ----------
    ptms : pandas dataframe
        dataframe containing PTM data
    regulatory : pandas dataframe
        dataframe containing regulatory data from phosphositeplus
    reference_col : str, optional
        column in regulatory dataframe to use as annotations. The default is 'ON_FUNCTION'.
    collapse_on_similar : bool, optional
        whether to collapse similar annotations into one category. The default is False.

    Returns
    -------
    annotation : pandas dataframe
        dataframe with PTMs in the rows and annotations in the columns, with 1 indicating that the PTM has that annotation
    """
    #create matrix indicating function of each ptm: ptm in the rows, function in columns, and 1 indicating that the ptm has that function## create molecular function table, with
    annotation = annotated_ptms.copy()
    if include_unknown:
        annotation.loc[annotation[reference_col].isna(), reference_col] = 'unknown'
    annotation = annotation.dropna(subset = reference_col)
    annotation[reference_col] = annotation[reference_col].apply(lambda x: x.split(';') if x == x else x)
    annotation = annotation.explode(reference_col).reset_index()
    if collapse_on_similar:
        annotation[reference_col] = annotation[reference_col].apply(lambda x: x.split(',')[0].strip(' ') if x == x else x)
    else:
        annotation[reference_col] = annotation[reference_col].apply(lambda x: x.strip(' ') if x == x else x)
    annotation['value'] = 1
    annotation = annotation[['PTM',reference_col, 'value']].drop_duplicates()
    annotation = annotation.pivot(index = 'PTM', columns = reference_col, values = 'value')
    #remove any sites with no functions
    annotation = annotation.dropna(how = 'all')
    return annotation

def getEnrichment(M, n, N, k, fishers = True):
    """
    Given a list of PTMs and their annotations, calculate the enrichment of a given annotation for that subset of PTMs

    Parameters
    ----------
    function_class : str
        Annotation to be tested.
    all_data : pandas dataframe
        Dataframe containing all PTMs and their annotations.
    subset_list : list
        List of PTMs to be tested.
    fishers : bool, optional
        Whether to use fishers exact test or hypergeometric test. The default is True, which uses the fishers exact test.
    
    Returns
    -------
    n : int
        Total number of PTMs with the annotation.
    k : int
        Total number of PTMs with the annotation in the subset.
    p : float
        p-value of the enrichment.
    M : int
        Total number of PTMs.
    N : int
        Total number of PTMs in the subset.
    odds : float
        Odds ratio of the enrichment.


    """

    if fishers:
        table = convertToFishers(M, n, N, k)
        odds, p = stats.fisher_exact(table)
        return p, odds
    else:
        p = hypergeom(M, n, N, k)
        return p

def generate_site_enrichment(ptm_subset, reference_df, subset_name = 'Subset', type = 'Process', fishers = True): 
    """
    Given a dataframe of PTMs and a dataframe of annotations, calculate the enrichment of each annotation in the PTM dataframe using a fishers exact test

    Parameters
    ----------
    ptm_subset: list
        ptms belonging to the group of interest
    reference_df : pandas dataframe
        dataframe containing annotations from PhosphoSitePlus, with PTMs in the rows and annotations in the columns. Must have at least the PTMs found in df.
    type : str, optional
        Indicates the type of annotations found in reference df. The default is 'Process', which indicates that the annotations are biological processes.

    Returns
    -------
    results : pandas dataframe
        dataframe containing the results of the fishers exact test for each annotation. Data collected includes the annotation name, number across the entire reference, number across the subset, fraction of subset with the annotation (Fraction Conserved), p-value, odds ratio, and log2 odds ratio
    """
    #initialize dictionary
    if fishers:
        results = {type:[], 'Number Across All PTMs':[], 'Number Across PTM subset':[], 'Fraction in PTM subset': [], 'p':[], 'Odds Ratio':[]}
    else:
        results = {type:[], 'Number Across All PTMs':[], 'Number Across PTM subset':[], 'Fraction in PTM subset': [], 'p':[]}
    #iterate through each annotation, perform analysis and record data
    for ob in reference_df.columns:
        #perform enrichment analysis
        enrichment_result = getEnrichment(ob, reference_df, ptm_subset, fishers = fishers)
        #store in results dictionary
        results[type].append(ob)
        results['Number Across All PTMs'].append(enrichment_result[0])
        results['Number Across PTM subset'].append(enrichment_result[1])
        results['Fraction in PTM subset'].append(enrichment_result[1]/enrichment_result[0])
        results['p'].append(enrichment_result[2])
        if fishers:
            results['Odds Ratio'].append(enrichment_result[5])

    #convert dictionary to dataframe
    results = pd.DataFrame(results).sort_values(ascending = True, by = 'p')
    #log transform odds ratio
    if fishers:
        results['Log Odds'] = results['Odds Ratio'].apply(lambda x: np.log2(x) if x != 0 else 0)
    #adjust p-values using Benjamini-Hochberg method or Bonferonni method
    results['BH adjusted p'] = adjustP(results['p'].values, method = 'BH')
    results['Bonferonni adjusted p'] = adjustP(results['p'].values, method = 'Bonf')
    #add subset name to dataframe
    results.insert(0, 'Subset', subset_name)
    
    return results
