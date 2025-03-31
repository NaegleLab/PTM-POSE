import numpy as np
import pandas as pd

#plotting 
import matplotlib.pyplot as plt
import seaborn as sns



#custom stat functions
from ptm_pose import helpers

def combine_outputs(spliced_ptms, altered_flanks, report_removed_annotations = True,  include_stop_codon_introduction = False, remove_conflicting = True, **kwargs):
    """
    Given the spliced_ptms (differentially included) and altered_flanks (altered flanking sequences) dataframes obtained from project and flanking_sequences modules, combine the two into a single dataframe that categorizes each PTM by the impact on the PTM site

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    altered_flanks: pd.DataFrame
        Dataframe with PTMs associated with altered flanking sequences and with annotations appended from various databases
    include_stop_codon_introduction: bool
        Whether to include PTMs that introduce stop codons in the altered flanks. Default is False.
    remove_conflicting: bool
        Whether to remove PTMs that are both included and excluded across different splicing events. Default is True.
    kwargs: dict
        Additional keyword arguments to pass to the function, will be passed to `helpers.filter_ptms` if filtering is desired. Will automatically filter out insignificant events if not provided
    """
    #filter spliced_ptms and altered_flanks dataframes to remove insignificant events or PTMs with low evidence
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        helpers.check_filter_kwargs(filter_arguments)
        spliced_ptms = helpers.filter_ptms(spliced_ptms, **filter_arguments)
        filter_arguments['remove_novel'] = False  #keep novel PTMs in altered flanks, as these are not removed from isoform
        altered_flanks = helpers.filter_ptms(altered_flanks, **filter_arguments)

    #extract specific direction of splicing change and add to dataframe
    spliced_ptms['Impact'] = spliced_ptms['dPSI'].apply(lambda x: 'Included' if x > 0 else 'Excluded')

    #restrict altered flanks to those that are changed and are not disrupted by stop codons
    if altered_flanks['Stop Codon Introduced'].dtypes != bool:
        altered_flanks['Stop Codon Introduced'] = altered_flanks['Stop Codon Introduced'].astype(bool)
    if include_stop_codon_introduction:
        altered_flanks['Impact'] = altered_flanks['Stop Codon Introduced'].apply(lambda x: 'Stop Codon Introduced' if x else 'Altered Flank')
    else:
        altered_flanks = altered_flanks[~altered_flanks['Stop Codon Introduced']].copy()
        altered_flanks['Impact'] = 'Altered Flank'

    #identify annotations that are found in both datasets
    annotation_columns_in_spliced_ptms = [col for col in spliced_ptms.columns if ':' in col]
    annotation_columns_in_altered_flanks = [col for col in altered_flanks.columns if ':' in col]
    annotation_columns = list(set(annotation_columns_in_spliced_ptms).intersection(annotation_columns_in_altered_flanks))
    if len(annotation_columns) != annotation_columns_in_spliced_ptms and report_removed_annotations:
        annotation_columns_only_in_spliced = list(set(annotation_columns_in_spliced_ptms) - set(annotation_columns_in_altered_flanks))
        annotation_columns_only_in_altered = list(set(annotation_columns_in_altered_flanks) - set(annotation_columns_in_spliced_ptms))
        if len(annotation_columns_only_in_spliced) > 0:
            print(f'Warning: some annotations in spliced ptms dataframe not found in altered flanks dataframe: {", ".join(annotation_columns_only_in_spliced)}. These annotations will be ignored. To avoid this, make sure to add annotations to both dataframes, or annotate the combined dataframe.')
        if len(annotation_columns_only_in_altered) > 0:
            print(f'Warning: some annotations in altered flanks dataframe not found in spliced ptms dataframe: {", ".join(annotation_columns_only_in_altered)}. These annotations will be ignored. To avoid this, make sure to add annotations to both dataframes, or annotate the combined dataframe.')

    #check if dPSI or sig columns are in both dataframes
    sig_cols = []
    if 'dPSI' in spliced_ptms.columns and 'dPSI' in altered_flanks.columns:
        sig_cols.append('dPSI')
    if 'Significance' in spliced_ptms.columns and 'Significance' in altered_flanks.columns:
        sig_cols.append('Significance')

    shared_columns = ['Impact', 'Gene', 'UniProtKB Accession', 'Isoform ID', 'Isoform Type', 'Residue', 'PTM Position in Isoform', 'Modification Class'] + sig_cols + annotation_columns
    combined = pd.concat([spliced_ptms[shared_columns], altered_flanks[shared_columns]])
    combined = combined.groupby([col for col in combined.columns if col != 'Impact'], as_index = False, dropna = False)['Impact'].apply(lambda x: ';'.join(set(x)))

    #remove ptms that are both included and excluded across different events
    if remove_conflicting:
        combined = combined[~((combined['Impact'].str.contains('Included')) & (combined['Impact'].str.contains('Excluded')))]

    return combined

def get_modification_counts(ptms, **kwargs):
    """
    Given PTM data (either spliced ptms, altered flanks, or combined data), return the counts of each modification class

    Parameters
    ----------
    ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events or with altered flanking sequences

    Returns
    -------
    modification_counts: pd.Series
        Series with the counts of each modification class
    """
    #filter ptms based on kwargs if provided
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        helpers.check_filter_kwargs(**filter_arguments)
        ptms = helpers.filter_ptms(ptms, **filter_arguments)

    ptms['Modification Class'] = ptms['Modification Class'].apply(lambda x: x.split(';'))
    ptms = ptms.explode('Modification Class')
    modification_counts = ptms.groupby('Modification Class').size()
    modification_counts = modification_counts.sort_values(ascending = True)
    return modification_counts

def get_modification_class_data(ptms, mod_class):
    """
    Given ptm dataframe and a specific modification class, return a dataframe with only the PTMs of that class

    Parameters
    ----------
    ptms : pd.DataFrame
        Dataframe with ptm information, such as the spliced_ptms or altered_flanks dataframe obtained during projection
    mod_class : str
        
        The modification class to filter by, e.g. 'Phosphorylation', 'Acetylation', etc.
    """
    #check if specific modification class was provided and subset data by modification if so
    ptms_of_interest = ptms[ptms['Modification Class'] == mod_class].copy()
    if ptms_of_interest.empty:
        raise ValueError(f"No PTMs found for modification class '{mod_class}'. Please check the input data or choose a different modification class.")

    return ptms_of_interest



def plot_modification_breakdown(spliced_ptms = None, altered_flanks = None, colors = sns.color_palette('colorblind'), ax = None, **kwargs):
    """
    Plot the number of PTMs that are differentially included or have altered flanking sequences, separated by PTM type

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs that are differentially included
    altered_flanks: pd.DataFrame
        Dataframe with PTMs that have altered flanking sequences
    colors: list
        List of colors to use for the bar plot (first two will be used). Default is seaborn colorblind palette.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    kwargs: dict
        Additional keyword arguments to pass to the function, will be passed to `helpers.filter_ptms` if filtering is desired. Will automatically filter out insignificant events by min_dpsi and significance if the columns are present
    """
    if spliced_ptms is None and altered_flanks is None:
        raise ValueError('Either spliced_ptms or altered_flanks must be provided to plot modification breakdown. Both may be provided.')
    
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        helpers.check_filter_kwargs(**filter_arguments)
        if spliced_ptms is not None:
            spliced_ptms = helpers.filter_ptms(spliced_ptms, **filter_arguments)
        if altered_flanks is not None:
            altered_flanks = helpers.filter_ptms(altered_flanks, **filter_arguments)
    
    if ax is None:
        fig, ax = plt.subplots(figsize = (4,4))

    #separate rows into unique PTM types
    if spliced_ptms is not None and altered_flanks is not None:
        differentially_included_counts = get_modification_counts(spliced_ptms.copy())
        altered_flanks_counts = get_modification_counts(altered_flanks.copy())
        ax.barh(differentially_included_counts.index, differentially_included_counts.values, color = colors[0], label = 'Differentially Included PTMs')
        altered_flanks_counts = altered_flanks_counts.reindex(differentially_included_counts.index, fill_value = 0)
        ax.barh(altered_flanks_counts.index, altered_flanks_counts.values, left = differentially_included_counts.values, color = colors[1], label = 'PTMs with Altered Flank')
        ax.legend()

        #annotate with number of combined PTMs
        total_count = differentially_included_counts.add(altered_flanks_counts, fill_value = 0)
        for i, num_ptm in enumerate(total_count.values):
            ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')  

        ax.set_xlim([0, total_count.max()*1.1])

    elif spliced_ptms is not None:
        modification_counts = get_modification_counts(spliced_ptms)
        ax.barh(modification_counts.index, modification_counts.values, color = colors[0])

        #annotate with number of PTMs
        for i, num_ptm in enumerate(modification_counts.values):
            ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')
    elif altered_flanks is not None:
        modification_counts = get_modification_counts(altered_flanks)
        ax.barh(modification_counts.index, modification_counts.values, color = colors[0])

        #annotate with number of PTMs
        for i, num_ptm in enumerate(modification_counts.values):
            ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')

    ax.set_xlabel('Number of PTMs')

