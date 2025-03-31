import numpy as np
import pandas as pd

import os
import time

#plotting 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

#custom stat functions
from ptm_pose import stat_utils, pose_config, annotate, helpers
from ptm_pose.analyze import summarize

try:
    import gseapy as gp
except ImportError:
    gp = None


package_dir = os.path.dirname(os.path.abspath(__file__))





def get_available_annotations(ptms):
    """
    Given a PTM dataframe, indicate the annotations that are available for analysis and indicate whether they have already been appended to the PTM dataset

    Parameters
    ----------
    ptms : pd.DataFrame
        contains PTM information and may have annotations already appended, such as spliced_ptms and altered_flanks dataframes generated during projection
    
    Returns
    -------
    available_annots : pd.DataFrame  
        DataFrame indicating the available annotation types and their sources, as well as whether they have been appended to the PTM data.
    """
    available_gmt = annotate.get_available_gmt_annotations(format = 'dataframe')
    available_gmt['Appended to PTM data?'] = 'No'


    #check to see if any annotations have been added to the spliced_ptms dataframe
    annot_cols = [col for col in ptms.columns if ':' in col]
    database_list = []
    annot_type_list = []
    for a in annot_cols:
        database_list.append(a.split(':')[0])
        annot_type_list.append(a.split(':')[1])
    available_in_df = pd.DataFrame({'Database':database_list, 'Annotation Type':annot_type_list, 'Appended to PTM data?': 'Yes'})

    #combine the two dataframes and remove duplicates
    available_annots = pd.concat([available_in_df, available_gmt]).drop_duplicates(keep = 'first', subset = ['Database', 'Annotation Type']).reset_index(drop = True)
    available_annots.sort_values(by = 'Annotation Type')
    return available_annots



def get_ptm_annotations(ptms, annot_type = 'Function', database = 'PhosphoSitePlus', collapse_on_similar = False, min_dpsi = 0.1, alpha = 0.05,  **kwargs):
    """
    Given spliced ptm information obtained from project and annotate modules, grab PTMs in spliced ptms associated with specific PTM modules

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        PTMs projected onto splicing events and with annotations appended from various databases
    annot_type: str
        Type of annotation to pull from spliced_ptms dataframe. Available information depends on the selected database. Default is 'Function'.
    database: str
        database from which PTMs are pulled. Options include 'PhosphoSitePlus', 'ELM', or 'PTMInt'. ELM and PTMInt data will automatically be downloaded, but due to download restrictions, PhosphoSitePlus data must be manually downloaded and annotated in the spliced_ptms data using functions from the annotate module. Default is 'PhosphoSitePlus'.
    collapse_on_similar : bool
        Whether to collapse similar annotations (for example, "cell growth, increased" and "cell growth, decreased") into a single category. Default is False.
    min_dpsi: float
        Minimum change in PSI required to return a PTM as associated with the annotation. Default is 0.1. This can be used to filter out PTMs that are not significantly spliced.
    alpha : float
        Significance threshold to use to filter PTMs based on their significance. Default is 0.05. This can be used to filter out PTMs that are not significantly spliced.
    kwargs: additional keyword arguments
        Additional keyword arguments to pass to the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.

    Returns
    -------
    annotations : pd.DataFrame
        Individual PTM information of PTMs that have been associated with the requested annotation type. 
    annotation_counts : pd.Series or pd.DataFrame
        Number of PTMs associated with each annotation of the requested annotation type. If dPSI col is provided or impact col is present, will output annotation counts for each type of impact ('Included', 'Excluded', 'Altered Flank') separately. 
    """
    if database == 'Combined':
        if annot_type == 'Interactions':
            annot_col = 'Combined:Interactions'
            if annot_col not in ptms.columns:
                ptms = annotate.combine_interaction_data(ptms).copy()
            ptms = ptms.dropna(subset = annot_col)
        elif annot_type == 'Writer Enzyme' or annot_type == 'Writer_Enzyme':
            annot_col = 'Combined:Writer_Enzyme'
            if annot_col not in ptms.columns:
                ptms = annotate.combine_enzyme_data(ptms).copy()
            ptms = ptms.dropna(subset = annot_col)
        elif annot_type == 'Eraser Enzyme' or annot_type == 'Eraser_Enzyme':
            annot_col = 'Combined:Eraser_Enzyme'
            if annot_col not in ptms.columns:
                ptms = annotate.combine_enzyme_data(ptms).copy()
            ptms = ptms.dropna(subset = annot_col)
        elif annot_type == 'Enzyme':
            raise ValueError('Enzyme is not a valid annotation type. Please use "Writer Enzyme" or "Eraser Enzyme" instead.')
        else:
            raise ValueError(f'{annot_type} is not a valid annotation type for Combined database. Please use "Interactions", "Writer Enzyme", or "Eraser Enzyme" instead, or use a specific database.')
    else:
        if database == 'PTMsigDB':
            if annot_type == 'Perturbation':
                raise ValueError('PTMsigDB has multiple perturbation datasets. Please specify which dataset to use. Options are: "Perturbation-DIA","Perturbation-DIA2", "Perturbation-PRM".')
            elif annot_type == 'Pathway':
                raise ValueError('PTMsigDB has multiple pathway datasets. Please specify which dataset to use. Options are: "Pathway-WikiPathway", "Pathway-NetPath", "Pathway-BI"')
            else:
                annot_col = f'{database}:{annot_type}'
        elif database == 'OmniPath':
            annot_col = f'{database}:{annot_type.replace(' ', '_')}'
        else:
            annot_col = f'{database}:{annot_type}'
        
        #check to see if requested annotation information is annotated to the spliced_ptms dataframe, if not, append information
        if annot_col not in ptms.columns:
            ptms = annotate.append_from_gmt(ptms, database = database, annot_type = annot_type).copy()
        ptms = ptms.dropna(subset = annot_col)



    #check to make sure requested annotation is available
    filter_arguments = helpers.extract_filter_kwargs(min_dpsi = min_dpsi, alpha = alpha, **kwargs)
    helpers.check_filter_kwargs(filter_arguments)
    ptms_of_interest = helpers.filter_ptms(ptms, **filter_arguments)


    #extract relevant annotation and remove PTMs without an annotation
    if 'Impact' not in ptms_of_interest.columns and 'dPSI' in ptms_of_interest.columns:
        ptms_of_interest['Impact'] = ptms_of_interest['dPSI'].apply(lambda x: 'Included' if x > 0 else 'Excluded')

    optional_cols = [col for col in ptms_of_interest.columns if col in ['Impact', 'dPSI', 'Significance']]
    annotations = ptms_of_interest[['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class'] + [annot_col] + optional_cols].copy()
    annotations = annotations.dropna(subset = annot_col).drop_duplicates()

    if annotations.empty:
        print("No PTMs with associated annotation")
        return None, None
    
    #combine repeat entries for same PTM (with multiple impacts)
    annotations = annotations.groupby(['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'], as_index = False).agg(lambda x: ';'.join(set([str(i) for i in x if i == i])))

    #separate distinct modification annotations in unique rows
    annotations_exploded = annotations.copy()
    annotations_exploded[annot_col] = annotations_exploded[annot_col].apply(lambda x: x.split(';') if isinstance(x, str) else np.nan)
    annotations_exploded = annotations_exploded.explode(annot_col)
    annotations_exploded[annot_col] = annotations_exploded[annot_col].apply(lambda x: x.strip() if isinstance(x, str) else np.nan)
    
    #if desired collapse similar annotations (for example, same function but increasing or decreasing)
    if collapse_on_similar:
        annotations_exploded[annot_col] = annotate.collapse_annotations(annotations_exploded[annot_col].values, database = database, annot_type = annot_type)
        annotations_exploded.drop_duplicates(inplace = True)
        annotations = annotations_exploded.groupby([col for col in annotations_exploded.columns if col != annot_col], as_index = False, dropna = False)[annot_col].apply(lambda x: ';'.join(set(x)))
    
    #get the number of annotations found
    annotation_counts = annotations_exploded.drop_duplicates(subset = ['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'] + [annot_col])[annot_col].value_counts()


    #additional_counts
    sub_counts = []
    impact_type = []
    if 'Impact' in annotations_exploded.columns:
        for imp in ['Included', 'Excluded', 'Altered Flank']:
            tmp_annotations = annotations_exploded[annotations_exploded['Impact'] == imp].copy()
            if not tmp_annotations.empty:
                tmp_annotations = tmp_annotations.drop_duplicates(subset = ['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'] + [annot_col])
                sub_counts.append(tmp_annotations[annot_col].value_counts())
                impact_type.append(imp)
    
        annotation_counts = pd.concat([annotation_counts] + sub_counts, axis = 1)
        annotation_counts.columns = ['All Impacted'] + impact_type 
        annotation_counts = annotation_counts.replace(np.nan, 0)
    
    #combine repeat entries for same PTM (with multiple impacts)
    annotations = annotations.groupby(['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'], as_index = False).agg(lambda x: ';'.join(set([str(i) for i in x if i == i])))


    if database == 'PTMcode': #convert to readable gene name
        annotation_counts.index = [pose_config.uniprot_to_genename[i].split(' ')[0] if i in pose_config.uniprot_to_genename.keys() else i for i in annotation_counts.index]

    return annotations, annotation_counts


def get_background_annotation_counts(database = 'PhosphoSitePlus', annot_type = 'Function', **kwargs):
    """
    Given a database and annotation type, retrieve the counts of PTMs associated with the requested annotation type across all PTMs in the ptm_coordinates dataframe used for projection

    Parameters
    ----------
    database : str
        Source of annotation. Default is PhosphoSitePlus
    annot_type : str
        Type of annotation that can be found in indicated database. Default is 'Function'. Other options include 'Process', 'Disease', 'Enzyme', 'Interactions', etc.
    kwargs: additional keyword arguments
        Additional keyword arguments to pass to the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(report_removed = False, **kwargs)
        helpers.check_filter_kwargs(filter_arguments)
        background_ptms = helpers.filter_ptms(pose_config.ptm_coordinates, **filter_arguments)
        background_ptms = background_ptms.rename(columns = {'Gene name':'Gene'})
        _, background_annotation_count = get_ptm_annotations(background_ptms, annot_type = annot_type, database = database, report_removed = False, **kwargs)
        background_size = background_ptms.drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class']).shape[0]
    else:
        annotation_dict = annotate.process_database_annotations(database = database, annot_type=annot_type, key_type = 'annotation')
        background_annotation_count = pd.Series(annotation_dict).apply(len)
        background_size = pose_config.ptm_coordinates.drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class']).shape[0]
    return background_annotation_count, background_size


def get_enrichment_inputs(ptms,  annot_type = 'Function', database = 'PhosphoSitePlus', background_type = 'all', collapse_on_similar = False, mod_class = None, alpha = 0.05, min_dpsi = 0.1, **kwargs):
    """
    Given the spliced ptms, altered_flanks, or combined PTMs dataframe, identify the number of PTMs corresponding to specific annotations in the foreground (PTMs impacted by splicing) and the background (all PTMs in the proteome or all PTMs in dataset not impacted by splicing). This information can be used to calculate the enrichment of specific annotations among PTMs impacted by splicing. Several options are provided for constructing the background data: all (based on entire proteome in the ptm_coordinates dataframe) or significance (foreground PTMs are extracted from provided spliced PTMs based on significance and minimum delta PSI)

    Parameters
    ----------
    ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases. This can be either the spliced_ptms, altered_flanks, or combined dataframe.
    annot_type : str
        type of annotation to pull the annotations from. Default is 'Function'.
    database : str
        source of annotations. Default is 'PhosphoSitePlus'.
    background_type : str
        Type of background to construct. Options are either 'all' (all PTMs in proteome) or 'significance' (only PTMs in dataset). Note that significance option assumes that PTMs have not already been filtered for significance.
    collapse_on_similar : bool
        Whether to collapse similar annotations (for example, "cell growth, increased" and "cell growth, decreased") into a single category. Default is False.
    mod_class : str
        Type of modification to perform enrichment for
    min_dpsi: float
        Minimum change in PSI required to return a PTM as associated with the annotation. Default is 0.1. This can be used to filter out PTMs that are not significantly spliced.
    alpha : float
        Significance threshold to use to filter PTMs based on their significance. Default is 0.05. This can be used to filter out PTMs that are not significantly spliced.
    kwargs: additional keyword arguments
        Additional keyword arguments to pass to the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    if background_type == 'all':
        background_annotation_count, background_size = get_background_annotation_counts(database = database, annot_type = annot_type, **kwargs)

        background_size = pose_config.ptm_coordinates.drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class']).shape[0]

        #filter ptms
        filter_arguments = helpers.extract_filter_kwargs(min_dpsi = min_dpsi, alpha = alpha, modification_class = mod_class, report_removed = False, **kwargs)
        helpers.check_filter_kwargs(filter_arguments)
        ptms = helpers.filter_ptms(ptms, **filter_arguments)
    elif background_type == 'significance':
        if 'Significance' not in ptms.columns or 'dPSI' not in ptms.columns:
            raise ValueError('Significance and dPSI columns must be present in spliced_ptms dataframe to construct a background based on significance (these columns must be provided during projection).')
        #filter ptms if any kwargs are provided
        if kwargs:
            filter_arguments = helpers.extract_filter_kwargs(min_dpsi = 0, alpha = 1, modification_class = mod_class, report_removed = False, **kwargs)
            helpers.check_filter_kwargs(filter_arguments)
            ptms = helpers.filter_ptms(ptms, **filter_arguments)


        background = ptms.copy()
        #restrict sample to significantly spliced ptms
        ptms = ptms[(ptms['Significance'] <= alpha) & (ptms['dPSI'].abs() >= min_dpsi)].copy()


        #check to make sure there are significant PTMs in the data and that there is a difference in the number of significant and background PTMs
        if ptms.shape[0] == 0:
            raise ValueError('No significantly spliced PTMs found in the data')
        elif ptms.shape[0] == background.shape[0]:
            raise ValueError(f'The foreground and background PTM sets are the same size when considering significance. Please check data or use "all" for the background_type parameter.')
        else:
            #remove impact from background columns if present
            cols_to_remove = [col for col in ['dPSI', 'Significance', 'Impact'] if col in background.columns]
            if len(cols_to_remove) > 0:
                background = background.drop(columns = cols_to_remove)
        #get background counts
            background_size = background.drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class']).shape[0]
        #get background counts
        _, background_annotation_count = get_ptm_annotations(background, annot_type = annot_type, database = database, collapse_on_similar = collapse_on_similar)
    else:
        raise ValueError('Invalid background type. Must be all (default) or significance')

    #get foreground counts and size
    foreground_size = ptms.shape[0]
    annotation_details, foreground_annotation_count = get_ptm_annotations(ptms, annot_type = annot_type, database = database, collapse_on_similar = collapse_on_similar, min_dpsi = min_dpsi, alpha = alpha, report_removed = False, **kwargs)


    #process annotation details into usable format
    if annotation_details is None:
        print('No PTMs with requested annotation type, so could not perform enrichment analysis')
        return np.repeat(None, 5)
    else:
        annot_col = f'{database}:{annot_type}'
        annotation_details[annot_col] = annotation_details[annot_col].str.split(';')
        annotation_details = annotation_details.explode(annot_col)
        annotation_details[annot_col] = annotation_details[annot_col].str.strip()
        annotation_details['PTM'] = annotation_details['Gene'] + '_' + annotation_details['Residue'] + annotation_details['PTM Position in Isoform'].astype(int).astype(str)
        annotation_details = annotation_details.groupby(annot_col)['PTM'].agg(';'.join)
    
    return foreground_annotation_count, foreground_size, background_annotation_count, background_size, annotation_details


def annotation_enrichment(ptms, database = 'PhosphoSitePlus', annot_type = 'Function', background_type = 'all', collapse_on_similar = False, mod_class = None, alpha = 0.05, min_dpsi = 0.1, **kwargs):#
    """
    Given spliced ptm information (differential inclusion, altered flanking sequences, or both), calculate the enrichment of specific annotations in the dataset using a hypergeometric test. Background data can be provided/constructed in a few ways:

    1. Use annotations from the entire phosphoproteome (background_type = 'all')
    2. Use the alpha and min_dpsi parameter to construct a foreground that only includes significantly spliced PTMs, and use the entire provided spliced_ptms dataframe as the background. This will allow you to compare the enrichment of specific annotations in the significantly spliced PTMs compared to the entire dataset. Will do this automatically if alpha or min_dpsi is provided.

    Parameters
    ----------
    ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    database: str
        database from which PTMs are pulled. Options include 'PhosphoSitePlus', 'ELM', 'PTMInt', 'PTMcode', 'DEPOD', 'RegPhos', 'PTMsigDB'. Default is 'PhosphoSitePlus'.
    annot_type: str
        Type of annotation to pull from spliced_ptms dataframe. Available information depends on the selected database. Default is 'Function'.
    background_type: str
        how to construct the background data. Options include 'pregenerated' (default) and 'significance'. If 'significance' is selected, the alpha and min_dpsi parameters must be provided. Otherwise, will use whole proteome in the ptm_coordinates dataframe as the background.
    collapse_on_similar: bool
        Whether to collapse similar annotations (for example, increasing and decreasing functions) into a single category. Default is False.
    mod_class: str
        modification class to subset, if any
    alpha: float
        significance threshold to use to subset foreground PTMs. Default is None.
    min_dpsi: float
        minimum delta PSI value to use to subset foreground PTMs. Default is None.
    kwargs: additional keyword arguments
        Additional keyword arguments to pass to the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    foreground_annotation_count, foreground_size, background_annotations, background_size, annotation_details = get_enrichment_inputs(ptms, background_type = background_type, annot_type = annot_type, database = database, collapse_on_similar = collapse_on_similar, mod_class = mod_class, alpha = alpha, min_dpsi = min_dpsi)
    

    if foreground_annotation_count is not None:
        #iterate through all annotations and calculate enrichment with a hypergeometric test
        results = pd.DataFrame(columns = ['Fraction Impacted', 'p-value'], index = foreground_annotation_count.index)
        for i, n in background_annotations.items():
            #number of PTMs in the foreground with the annotation
            if i in foreground_annotation_count.index.values:
                if isinstance(foreground_annotation_count, pd.Series):
                    k = foreground_annotation_count[i]
                elif foreground_annotation_count.shape[1] == 1:
                    k = foreground_annotation_count.loc[i, 'count']
                elif foreground_annotation_count.shape[1] > 1:
                    k = foreground_annotation_count.loc[i, 'All Impacted']

                p = stat_utils.getEnrichment(background_size, n, foreground_size, k, fishers = False)
                results.loc[i, 'Fraction Impacted'] = f"{k}/{n}"
                results.loc[i, 'p-value'] = p

        results = results.sort_values('p-value')
        results['Adjusted p-value'] = stat_utils.adjustP(results['p-value'].values)
        results = pd.concat([results, annotation_details], axis = 1)
    else:
        results = None

    return results


def plot_available_annotations(ptms, only_annotations_in_data = False, show_all_ptm_count = False, ax = None, **kwargs):
    """
    Given a dataframe with ptm annotations added, show the number of PTMs associated with each annotation type

    Parameters
    ----------
    ptms: pd.DataFrame
        Dataframe with PTMs and annotations added
    only_annotations_in_data : bool
        Only plot annotations that are already appended to the dataset
    show_all_ptm_count: bool
        Whether to show the total number of PTMs in the dataset. Default is True.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    """
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(report_removed = False, **kwargs)
        helpers.check_filter_kwargs(filter_arguments)
        ptms = helpers.filter_ptms(ptms, **filter_arguments)
    else:
        ptms = ptms.copy()

    if ax is None:
        fig, ax = plt.subplots(figsize = (5,5))

    if show_all_ptm_count:
        num_ptms = [ptms.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class']).shape[0]]
        num_ptms_filters = ['All PTMs']
        filter_source = ['None']
    else:
        num_ptms = []
        num_ptms_filters = []
        filter_source = []

 

    available_annotations = get_available_annotations(ptms)
    available_annotations = available_annotations.sort_values(by = ['Database', 'Annotation Type'])
    if only_annotations_in_data:
        available_annotations = available_annotations[available_annotations['Appended to PTM data?'] == 'Yes']

    for _, row in available_annotations.iterrows():
        annot_col = f"{row['Database']}:{row['Annotation Type']}"
        if annot_col not in ptms.columns:
            ptms = annotate.append_from_gmt(ptms, database = row['Database'], annot_type = row['Annotation Type'])
        
        num_ptms.append(ptms.dropna(subset = annot_col).drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class']).shape[0])
        num_ptms_filters.append(f"{row['Annotation Type']}({row['Database']})")
        filter_source.append(row['Database'])

    plt_data = pd.DataFrame({'PTM Count': num_ptms, 'Annotation Type': num_ptms_filters, 'Source': filter_source})
    plt_data = plt_data.sort_values(by = 'PTM Count', ascending = True)
    #remove annotations with 0 PTMs
    plt_data = plt_data[plt_data['PTM Count'] > 0]

    
    #plot bar plot
    #color bars based on datasource
    databases = available_annotations['Database'].unique()
    #construct a color palette for the databases, with "None" as gray
    palette = {db: col for db, col in zip(databases, sns.color_palette("rainbow", len(databases)))}
    palette['None'] = 'gray' #set "None" to gray for all PTMs (whether or not associated with annotation)
    colors = []
    for source in plt_data['Source']:
        colors.append(palette[source])


    ax.barh(plt_data['Annotation Type'], plt_data['PTM Count'], color = colors)
    ax.set_xlabel('Number of PTMs with annotation')
    
    #annotate with number of PTMs
    for i, num_ptm in enumerate(plt_data['PTM Count']):
        ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')

    #create legend
    handles = [plt.Rectangle((0,0),1,1, color = color) for color in palette.values() if color != 'gray']
    labels = [source for source in palette.keys() if source != 'None']
    ax.legend(handles, labels, title = 'Annotation Source')
    plt.show()

def plot_annotation_counts(spliced_ptms = None, altered_flanks = None, database = 'PhosphoSitePlus', annot_type = 'Function', collapse_on_similar = True, colors = None, fontsize = 10, top_terms = 5, legend = False, legend_loc = (1.05,0.5), title = None, title_type = 'database', ax = None, **kwargs):
    """
    Given a dataframe with PTM annotations added, plot the top annotations associated with the PTMs

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with differentially included PTMs
    altered_flanks: pd.DataFrame
        Dataframe with PTMs associated with altered flanking sequences
    database: str
        Database to use for annotations. Default is 'PhosphoSitePlus'.
    annot_type: str
        Type of annotation to plot. Default is 'Function'.
    collapse_on_similar: bool
        Whether to collapse similar annotations into a single category. Default is True.
    colors: list
        List of colors to use for the bar plot. Default is None.
    top_terms: int
        Number of top terms to plot. Default is 5.
    legend: bool
        Whether to show the legend. Default is True.
    legend_loc: tuple
        Location of the legend. Default is None, which will place the legend in the upper right corner.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    title_type: str
        Type of title to use for the plot. Default is 'database'. Options include 'database' and 'detailed'.
    title: str
        Title to use for the plot. Default is None.
    fontsize: int
        Font size for the plot. Default is 10.
    legend_loc: tuple
        Location of the legend. Default is None, which will place the legend to the right of the plot.
    **kwargs: additional keyword arguments
        Additional keyword arguments, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    #check what was provided as input, combine if both are provided
    if spliced_ptms is None and altered_flanks is None:
        raise ValueError('Must provide either spliced_ptms or altered_flanks dataframe to plot annotation counts')
    elif spliced_ptms is not None and altered_flanks is not None:
        ptms = summarize.combine_outputs(spliced_ptms, altered_flanks, report_removed = False, report_removed_annotations = False,  **kwargs)
    elif spliced_ptms is not None:
        ptms = spliced_ptms.copy()
    else:
        ptms = altered_flanks.copy()
    
    #get number of PTMs associated with each annotation
    _, annotation_counts = get_ptm_annotations(ptms, annot_type = annot_type, database = database, collapse_on_similar = collapse_on_similar, report_removed = False, kwargs = kwargs)
    if annotation_counts is None:
        return None


    if ax is None:
        fig, ax = plt.subplots(figsize = (2,3))


    if colors is None:
        colors = ['lightgrey', 'gray', 'white']

    #grab counts for plot
    if isinstance(annotation_counts, pd.Series): #if only one impact
        annotation_counts = annotation_counts.head(top_terms).sort_values(ascending = True)
        if isinstance(colors, list) or isinstance(colors, np.ndarray):
            colors = colors[0]
        
        ax.barh(annotation_counts.index, annotation_counts.values, color = colors, edgecolor = 'black')
        legend = False
    else: # if multiple impact types (Included, Excluded, Altered Flank)
        legend_labels = ['All Impacted']
        annotation_counts = annotation_counts.head(top_terms).sort_values(by = 'All Impacted', ascending = True)
        left = 0
        if 'Excluded' in annotation_counts.columns:
            ax.barh(annotation_counts['Excluded'].index, annotation_counts['Excluded'].values, height = 1, edgecolor = 'black', color = colors[0])
            legend_labels.append('Excluded')
            left = annotation_counts['Excluded'].values
        if 'Included' in annotation_counts.columns:
            ax.barh(annotation_counts['Included'].index, annotation_counts['Included'].values, left = left, height = 1, color = colors[1], edgecolor = 'black')
            legend_labels.append('Included')
            left = left + annotation_counts['Included'].values
        if 'Altered Flank' in annotation_counts.columns:
            ax.barh(annotation_counts['Altered Flank'].index, annotation_counts['Altered Flank'].values, left = left, height = 1, color = colors[2], edgecolor = 'black')
            legend_labels.append('Altered Flank')
            left = left + annotation_counts['Altered Flank'].values
    #ax.set_xticks([0,50,100,150])
    ax.set_ylabel('', fontsize = fontsize)
    ax.set_xlabel('Number of PTMs', fontsize = fontsize)

    if title is not None:
        ax.set_title(title, fontsize = fontsize, weight = 'bold')
    elif title_type == 'detailed':
        ax.set_title(f'Top {top_terms} {database} {annot_type} Annotations', fontsize = fontsize, weight = 'bold')
    elif title_type == 'database':
        ax.set_title(f'{database}', fontsize = fontsize, weight = 'bold')


    x_label_dict = {'Function':'Number of PTMs\nassociated with Function', 'Process':'Number of PTMs\nassociated with Process', 'Disease':'Number of PTMs\nassociated with Disease', 'Kinase':'Number of Phosphosites\ntargeted by Kinase', 'Enzyme':'Number of PTMs\ntargeted by Enzyme', 'Interactions': 'Number of PTMs\nthat regulate interaction\n with protein','Motif Match':'Number of PTMs\nfound within a\nmotif instance', 'Intraprotein': 'Number of PTMs\nthat are important\for intraprotein\n interactions','Phosphatase':'Number of Phosphosites\ntargeted by Phosphatase', 'Perturbation (DIA2)': "Number of PTMs\nAffected by Perturbation\n(Measured by DIA)", 'Perturbation (PRM)': 'Number of PTMs\nAffected by Perturbation\n(Measured by PRM)', 'NetPath':'Number of PTMs/Genes\nassociated with NetPath', 'Perturbation':'Number of PTMs\nAffected by Perturbation'}
    if annot_type in x_label_dict.keys():
        ax.set_xlabel(x_label_dict[annot_type], fontsize = fontsize)
    else:
        ax.set_xlabel(f'Number of PTMs\nassociated with {annot_type}', fontsize = fontsize)

    
    #make a custom legend
    if legend:
        handles = []
        for i, impact in enumerate(['Included', 'Excluded', 'Altered Flank']):
            if impact in legend_labels:
                handles.append(mpatches.Patch(facecolor = colors[i], edgecolor = 'black', label = impact))


        ax.legend(handles = handles, ncol = 1, fontsize = fontsize - 2, title = 'Type of Impact', title_fontsize = fontsize - 1, loc = legend_loc)

    #return ax


def gene_set_enrichment(spliced_ptms = None, altered_flanks = None, sig_col = 'Significance', dpsi_col = 'dPSI', alpha = 0.05, min_dpsi = None, gene_sets = ['GO_Biological_Process_2023', 'Reactome_2022'], background = None, return_sig_only = True, max_retries = 5, delay = 10, **kwargs):
    """
    Given spliced_ptms and/or altered_flanks dataframes (or the dataframes combined from combine_outputs()), perform gene set enrichment analysis using the enrichr API

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with differentially included PTMs projected onto splicing events and with annotations appended from various databases. Default is None (will not be considered in analysis). If combined dataframe is provided, this dataframe will be ignored. 
    altered_flanks: pd.DataFrame
        Dataframe with PTMs associated with altered flanking sequences and with annotations appended from various databases. Default is None (will not be considered). If combined dataframe is provided, this dataframe will be ignored.
    combined: pd.DataFrame
        Combined dataframe with spliced_ptms and altered_flanks dataframes. Default is None. If provided, spliced_ptms and altered_flanks dataframes will be ignored.
    gene_sets: list
        List of gene sets to use in enrichment analysis. Default is ['KEGG_2021_Human', 'GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023','Reactome_2022']. Look at gseapy and enrichr documentation for other available gene sets
    background: list
        List of genes to use as background in enrichment analysis. Default is None (all genes in the gene set database will be used).
    return_sig_only: bool
        Whether to return only significantly enriched gene sets. Default is True.
    max_retries: int
        Number of times to retry downloading gene set enrichment data from enrichr API. Default is 5.
    delay: int
        Number of seconds to wait between retries. Default is 10.
    **kwargs: additional keyword arguments
        Additional keyword arguments to pass to the combine_outputs() function from the summarize module. These will be used to filter the spliced_ptms and altered_flanks dataframes before performing gene set enrichment analysis. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `combine_outputs()` function for more options.

    Returns
    -------
    results: pd.DataFrame
        Dataframe with gene set enrichment results from enrichr API

    """
    if gp is None:
        raise ImportError('gseapy package is required to perform gene set enrichment analysis. Please install it using "pip install gseapy"')
    
    if background is not None:
        raise ValueError('Background data not supported at this time, but will be added in future versions. Please set background = None, which will use all genes in the gene set database as the background.')
    
    #grab the genes associated with impacted PTMs
    if spliced_ptms is not None and altered_flanks is not None:
        #gene information (total and spliced genes)
        combined = summarize.combine_outputs(spliced_ptms, altered_flanks, report_removed_annotations=False, **kwargs)
        foreground = combined.copy()
        type = 'Differentially Included + Altered Flanking Sequences'

        #isolate the type of impact on the gene
        combined_on_gene = combined.groupby('Gene')['Impact'].apply(lambda x: ';'.join(set(x)))
        included = combined_on_gene.str.contains('Included')
        excluded = combined_on_gene.str.contains('Excluded')
        differential = included | excluded
        altered_flank = combined_on_gene.str.contains('Altered Flank')

        altered_flank_only = altered_flank & ~differential
        differential_only = differential & ~altered_flank
        both = differential & altered_flank

        altered_flank_only = combined_on_gene[altered_flank_only].index.tolist()
        differential_only = combined_on_gene[differential_only].index.tolist()
        both = combined_on_gene[both].index.tolist()
    elif spliced_ptms is not None:
        foreground = spliced_ptms.copy()
        combined = spliced_ptms.copy()
        type = 'Differentially Included'

        #isolate the type of impact on the gene
        altered_flank_only = []
        differential_only = spliced_ptms['Gene'].unique().tolist()
        both = []
    elif altered_flanks is not None:
        foreground = altered_flanks.copy()
        combined = altered_flanks.copy()
        type = 'Altered Flanking Sequences'

        #isolate the type of impact on the gene
        altered_flank_only = altered_flanks['Gene'].unique().tolist()
        differential_only = []
        both = []
    else:
        raise ValueError('No dataframes provided. Please provide spliced_ptms, altered_flanks, or the combined dataframe.')
    
    #restrict to significant ptms, if available
    if sig_col in foreground.columns and (min_dpsi is not None and dpsi_col in foreground.columns):
        foreground = foreground[foreground[sig_col] <= alpha].copy()
        foreground = foreground[foreground[dpsi_col].abs() >= min_dpsi]
    elif sig_col in foreground.columns:
        if min_dpsi is not None:
            print('`min_dpsi` provided but `dpsi_col` not found in dataframe. Ignoring `min_dpsi` parameter. Please indicate correct column name for dPSI values and rerun if desired.')
        foreground = foreground[foreground[sig_col] <= alpha].copy()
    elif min_dpsi is not None and 'dPSI' in combined.columns:
        print('`sig_col` not found in dataframe. Ignoring `alpha` parameter. Please indicate correct column name for Significance and rerun if desired.')
        foreground = combined[combined['dPSI'].abs() >= min_dpsi].copy()
    else:
        print('Significance column not found and min_dpsi not provided. All PTMs in dataframe will be considered as the foreground')

    foreground = foreground['Gene'].unique().tolist()   



    

    
    #perform gene set enrichment analysis and save data
    for i in range(max_retries):
        try:
            enr = gp.enrichr(foreground, gene_sets = gene_sets)
            break
        except Exception as e: 
            enrichr_error = e
            time.sleep(delay)
    else:
        raise Exception('Failed to run enrichr analysis after ' + str(max_retries) + ' attempts. Please try again later. Error given by EnrichR: ' + str(enrichr_error))
    
    results = enr.results.copy()
    results['Type'] = type

    #indicate the genes in each gene set associated with each type of impact
    results['Genes with Differentially Included PTMs only'] = results['Genes'].apply(lambda x: ';'.join(set(x.split(';')) & (set(differential_only))))
    results['Genes with PTM with Altered Flanking Sequence only'] = results['Genes'].apply(lambda x: ';'.join(set(x.split(';')) & (set(altered_flank_only))))
    results['Genes with Both'] = results['Genes'].apply(lambda x: ';'.join(set(x.split(';')) & (set(both))))


    if return_sig_only:
        return results[results['Adjusted P-value'] <= 0.05]
    else:
        return results
    

def draw_pie(dist, xpos, ypos, size,colors,edgecolor =None, type = 'donut', ax=None):
    """
    Draws pies individually, as if points on a scatter plot. This function was taken from this stack overflow post: https://stackoverflow.com/questions/56337732/how-to-plot-scatter-pie-chart-using-matplotlib
    
    Parameters
    ----------
    dist: list
        list of values to be represented as pie slices for a single point
    xpos: float
        x position of pie chart in the scatter plot
    ypos: float
        y position of pie chart in the scatter plot
    size: float
        size of pie chart
    colors: list
        list of colors to use for pie slices
    ax: matplotlib.Axes
        axis to plot on, if None, will create new figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))
    #remove slices with 0 size
    colors = [c for c, d in zip(colors, dist) if d != 0]
    dist = [d for d in dist if d != 0]
    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()

    num_colors = len(dist)
    for i, r1, r2 in zip(range(num_colors), pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])

        ax.scatter([xpos], [ypos], marker=xy, s=size, facecolor= colors[i], edgecolors=edgecolor, linewidth = 0.3)

        if type == 'donut': # add white circle in the middle
            donut_edgecolors = 'w' if edgecolor is None else edgecolor
            ax.scatter([xpos], [ypos], s=size/5, facecolor='w', edgecolors=donut_edgecolors, linewidth = 0.3)
    return ax



def plot_EnrichR_pies(enrichr_results, top_terms = None, terms_to_plot = None, colors = None, edgecolor = None, row_height = 0.3, type = 'circle', ax = None):
    """
    Given PTM-specific EnrichR results, plot EnrichR score for the provided terms, with each self point represented as a pie chart indicating the fraction of genes in the group with PTMs
    
    Parameters
    ----------
    ptm_results: pd.selfFrame
        selfFrame containing PTM-specific results from EnrichR analysis
    num_to_plot: int
        number of terms to plot, if None, will plot all terms. Ignored if specific terms are provided in terms to plot list
    terms_to_plot: list
        list of terms to plot
    colors: list
        list of colors to use for pie slices. Default is None, which will use seaborn colorblind palette
    edgecolor: str
        color to use for edge of pie slices. Default is None, which will use the same color as the pie slice
    row_height: float
        height of each row in the plot. Default is 0.3.
    type: str
        type of pie chart to plot. Default is 'circle'. Options include 'circle' and 'donut' (hole in center).
    ax: matplotlib.Axes
        axis to plot on, if None, will create new figure
    """
    if colors is None:
        colors = sns.color_palette('colorblind', n_colors = 3)


    plt_data = enrichr_results.copy()
    plt_data['Number with Differential Inclusion Only'] = plt_data['Genes with Differentially Included PTMs only'].apply(lambda x: len(x.split(';')))
    plt_data['Number with Altered Flank Only'] = plt_data['Genes with Differentially Included PTMs only'].apply(lambda x: len(x.split(';')))
    plt_data['Number with Both'] = plt_data['Genes with Both'].apply(lambda x: len(x.split(';')) if x != '' else 0)
    

    if terms_to_plot is None:
        plt_data = plt_data.sort_values(by = 'Combined Score')
        if top_terms is not None:
            plt_data = plt_data.iloc[-top_terms:] if top_terms < plt_data.shape[0] else plt_data
    else:
        plt_data = plt_data[plt_data['Term'].isin(terms_to_plot)].sort_values(by = 'Combined Score')
        if plt_data.shape[0] == 0:
            print('No significant terms found in EnrichR results. Please check the terms_to_plot list and try again.')
            return
        

    #remove gene ontology specific terms
    plt_data['Term'] = plt_data['Term'].apply(lambda x: x.split(' R-HSA')[0] +' (R)' if 'R-HSA' in x else x.split('(GO')[0]+' (GO)')
    #construct multiple piecharts for each term in 'Term' column, where location along x-axis is dictated by combined score and piechart is dictated by 'Fraction With PTMs'
    plt_data = plt_data.reset_index(drop = True)

    #set up figure
    if ax is None:
        figure_length = plt_data.shape[0]*row_height
        fig, ax = plt.subplots(figsize = (2, figure_length))
    
    #get non-inf max score and replace inf values with max score
    maxscore = np.nanmax(plt_data['Combined Score'][plt_data['Combined Score'] != np.inf])
    plt_data['Combined Score'] = plt_data['Combined Score'].replace([-np.inf, np.inf], maxscore)
    ax.set_xlim([maxscore*-0.05, maxscore*1.1])
    mult = 4
    ax.set_yticks(list(range(0,plt_data.shape[0]*mult,mult)))
    ax.set_yticklabels(plt_data['Term'].values)
    ax.set_ylim([-(mult/2), plt_data.shape[0]*mult-(mult/2)])
    type = 'circle'
    event_type = plt_data['Type'].values[0]
    for i, row in plt_data.iterrows():
        if event_type == 'Differentially Included + Altered Flanking Sequences':
            draw_pie([row['Number with Differential Inclusion Only'], row['Number with Altered Flank Only'], row['Number with Both']],xpos = row['Combined Score'], ypos = i*mult, colors = colors, edgecolor=edgecolor,ax = ax, type = type, size = 70)
        else:
            draw_pie([1],xpos = row['Combined Score'], ypos = i*mult, colors = colors, edgecolor=edgecolor,ax = ax, type = type, size = 70)
        
        ax.axhline(i*mult+(mult/2), c= 'k', lw = 0.5)
        ax.axhline(i*mult-(mult/2), c = 'k', lw = 0.5)
        #ax.tick_params(labelsize = )

    #make a custom legend
    if event_type == 'Differentially Included + Altered Flanking Sequences':
        handles = [mpatches.Patch(color = colors[2], label = 'Contains Both Events'), mpatches.Patch(color = colors[1], label = 'PTMs with Altered Flanking Sequence'), mpatches.Patch(color = colors[0], label = 'Differentially Included PTMs')]
        ax.legend(handles = handles, loc = 'upper center', borderaxespad = 0, bbox_to_anchor = (0.5, 1 + (1/figure_length)), ncol = 1, fontsize = 9)



    ax.set_xlabel('EnrichR Combined\nScore', fontsize = 11)
