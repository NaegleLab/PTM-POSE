import numpy as np
import pandas as pd
import networkx as nx

import os

#plotting 
import matplotlib.pyplot as plt

#analysis packages
import gseapy as gp

#custom stat functions
from ptm_pose import stat_utils, pose_config, annotate, helpers

def get_annotation_col(spliced_ptms, annotation_type = 'Function', database = 'PhosphoSitePlus'):
    """
    Given the database of interest and annotation type, return the annotation column that will be found in a annotated spliced_ptm dataframe

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTM annotations added from annotate module
    annotation_type: str
        Type of annotation to pull from spliced_ptms dataframe. Available information depends on the selected database. Default is 'Function'.
    database: str
        database from which PTMs are pulled. Options include 'PhosphoSitePlus', 'ELM', 'PTMInt', 'PTMcode', 'DEPOD', and 'RegPhos'. Default is 'PhosphoSitePlus'.

    Returns
    -------
    annotation_col: str
        Column name in spliced_ptms dataframe that contains the requested annotation
    """
    if database == 'Combined':
        if f'Combined:{annotation_type}' not in spliced_ptms.columns:
            raise ValueError(f'Requested annotation data has not yet been added to spliced_ptms dataframe. Please run the annotate.{pose_config.annotation_function_dict[database]} function to append this information.')
        return f'Combined:{annotation_type}'
    elif annotation_type in pose_config.annotation_col_dict[database].keys():
        annotation_col = pose_config.annotation_col_dict[database][annotation_type]
        if annotation_col not in spliced_ptms.columns:
            raise ValueError(f'Requested annotation data has not yet been added to spliced_ptms dataframe. Please run the annotate.{pose_config.annotation_function_dict[database]} function to append this information.')
        return annotation_col
    else:
        raise ValueError(f"Invalid annotation type for {database}. Available annotation data for {database} includes: {', '.join(pose_config.annotation_col_dict[database].keys())}")


def combine_outputs(spliced_ptms, altered_flanks, mod_class = None, include_stop_codon_introduction = False, remove_conflicting = True):
    """
    Given the spliced_ptms (differentially included) and altered_flanks (altered flanking sequences) dataframes obtained from project and flanking_sequences modules, combine the two into a single dataframe that categorizes each PTM by the impact on the PTM site

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    altered_flanks: pd.DataFrame
        Dataframe with PTMs associated with altered flanking sequences and with annotations appended from various databases
    mod_class: str
        modification class to subset, if any
    include_stop_codon_introduction: bool
        Whether to include PTMs that introduce stop codons in the altered flanks. Default is False.
    remove_conflicting: bool
        Whether to remove PTMs that are both included and excluded across different splicing events. Default is True.
    """
    #process differentially included PTMs and altered flanking sequences
    if mod_class is not None:
        spliced_ptms = get_modification_class_data(spliced_ptms, mod_class)
        altered_flanks = get_modification_class_data(altered_flanks, mod_class)

    #extract specific direction of splicing change and add to dataframe
    spliced_ptms['Impact'] = spliced_ptms['dPSI'].apply(lambda x: 'Included' if x > 0 else 'Excluded')

    #restrict altered flanks to those that are changed and are not disrupted by stop codons
    if include_stop_codon_introduction:
        altered_flanks = altered_flanks[~altered_flanks['Matched']].copy()
        altered_flanks['Impact'] = altered_flanks['Stop Codon Introduced'].apply(lambda x: 'Stop Codon Introduced' if x else 'Altered Flank')
    else:
        altered_flanks = altered_flanks[(~altered_flanks['Matched']) & (~altered_flanks['Stop Codon Introduced'])].copy()
        altered_flanks['Impact'] = 'Altered Flank'

    #identify annotations that are found in both datasets
    annotation_columns_in_spliced_ptms = [col for col in spliced_ptms.columns if any([annot in col for annot in ['PSP', 'ELM','PTMInt','PTMcode','RegPhos','Combined', 'DEPOD']])]
    annotation_columns_in_altered_flanks = [col for col in altered_flanks.columns if any([annot in col for annot in ['PSP', 'ELM','PTMInt','PTMcode','RegPhos','Combined', 'DEPOD']])]
    annotation_columns = list(set(annotation_columns_in_spliced_ptms).intersection(annotation_columns_in_altered_flanks))
    if len(annotation_columns) != annotation_columns_in_spliced_ptms:
        print('Some annotations in spliced ptms dataframe not found in altered flanks dataframe. These annotations will be ignored. To avoid this, make sure to add annotations to both dataframes, or annotate the combined dataframe.')

    #check if dPSI or sig columns are in both dataframes
    sig_cols = []
    if 'dPSI' in spliced_ptms.columns and 'dPSI' in altered_flanks.columns:
        sig_cols.append('dPSI')
    if 'Significance' in spliced_ptms.columns and 'Significance' in altered_flanks.columns:
        sig_cols.append('Significance')

    shared_columns = ['Impact', 'Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'] + sig_cols + annotation_columns
    combined = pd.concat([spliced_ptms[shared_columns], altered_flanks[shared_columns]])
    combined = combined.groupby([col for col in combined.columns if col != 'Impact'], as_index = False, dropna = False)['Impact'].apply(lambda x: ';'.join(set(x)))

    #remove ptms that are both included and excluded across different events
    if remove_conflicting:
        combined = combined[~((combined['Impact'].str.contains('Included')) & (combined['Impact'].str.contains('Excluded')))]

    return combined

def simplify_annotation(annotation, sep = ','):
    """
    Given an annotation, remove additional information such as whether or not a function is increasing or decreasing. For example, 'cell growth, induced' would be simplified to 'cell growth'

    Parameters
    ----------
    annotation: str
        Annotation to simplify
    sep: str
        Separator that splits the core annotation from additional detail. Default is ','. Assumes the first element is the core annotation.

    Returns
    -------
    annotation: str
        Simplified annotation
    """
    annotation = annotation.split(sep)[0].strip(' ') if annotation == annotation else annotation
    return annotation

def collapse_annotations(annotations, database = 'PhosphoSitePlus', annotation_type = 'Function'):
    sep_dict = {'PhosphoSitePlus':{'Function':',', 'Process':',','Interactions':'(', 'Disease':'->'}, 'ELM': {'Interactions': ' ', 'Motif Match': ' '}, 'PTMInt':{'Interactions':'->'}, 'PTMcode':{'Interactions':'_', 'Intraprotein':' '}, 'RegPhos':{'Kinase':' '}, 'DEPOD':{'Phosphatase':' '}, 'Combined':{'Kinase':' ', 'Interactions':'->'}}
                
    sep = sep_dict[database][annotation_type]
    collapsed = []
    for annot in annotations:
        if annot == annot:
            collapsed.append(simplify_annotation(annot, sep = sep))
        else:
            collapsed.append(annot)
    return collapsed


def get_modification_class_data(spliced_ptms, mod_class):
    #check if specific modification class was provided and subset data by modification if so
    if mod_class in spliced_ptms['Modification Class'].values:
        ptms_of_interest = spliced_ptms[spliced_ptms['Modification Class'].str.contains(mod_class)].copy()
    else:
        ptms_of_interest['Modification Class'] = ptms_of_interest['Modification Class'].apply(lambda x: x.split(';') if x == x else np.nan)
        ptms_of_interest = ptms_of_interest.explode('Modification Class').dropna(subset = 'Modification Class')
        available_ptms = ptms_of_interest['Modification Class'].unique()
        raise ValueError(f"Requested modification class not present in the data. The available modifications include {', '.join(available_ptms)}")

    return ptms_of_interest

def get_ptm_annotations(spliced_ptms, annotation_type = 'Function', database = 'PhosphoSitePlus', mod_class = None, collapse_on_similar = False):
    """
    Given spliced ptm information obtained from project and annotate modules, grab PTMs in spliced ptms associated with specific PTM modules

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        PTMs projected onto splicing events and with annotations appended from various databases
    annotation_type: str
        Type of annotation to pull from spliced_ptms dataframe. Available information depends on the selected database. Default is 'Function'.
    database: str
        database from which PTMs are pulled. Options include 'PhosphoSitePlus', 'ELM', or 'PTMInt'. ELM and PTMInt data will automatically be downloaded, but due to download restrictions, PhosphoSitePlus data must be manually downloaded and annotated in the spliced_ptms data using functions from the annotate module. Default is 'PhosphoSitePlus'.
    mod_class: str
        modification class to subset 
    """
    #check to make sure requested annotation is available
    if database != 'Combined':
        annotation_col = get_annotation_col(spliced_ptms, database = database, annotation_type = annotation_type)
    else:
        annotation_col = f'Combined:{annotation_type}'


    #check if specific modification class was provided and subset data by modification if so
    if mod_class is not None:
        ptms_of_interest = get_modification_class_data(spliced_ptms, mod_class)
    else:
        ptms_of_interest = spliced_ptms.copy()

    #extract relevant annotation and remove PTMs without an annotation
    optional_cols = [col for col in ptms_of_interest.columns if col in ['Impact', 'dPSI', 'Significance']]
    annotations = ptms_of_interest[['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'] + [annotation_col] + optional_cols].copy()
    annotations = annotations.dropna(subset = annotation_col).drop_duplicates()

    #separate distinct modification annotations in unique rows
    annotations_exploded = annotations.copy()
    annotations_exploded[annotation_col] = annotations_exploded[annotation_col].apply(lambda x: x.split(';') if isinstance(x, str) else np.nan)
    annotations_exploded = annotations_exploded.explode(annotation_col)
    annotations_exploded[annotation_col] = annotations_exploded[annotation_col].apply(lambda x: x.strip() if isinstance(x, str) else np.nan)
    
    #if desired collapse similar annotations (for example, same function but increasing or decreasing)
    if collapse_on_similar:
        annotations_exploded[annotation_col] = collapse_annotations(annotations_exploded[annotation_col].values, database = database, annotation_type = annotation_type)
        annotations_exploded.drop_duplicates(inplace = True)
        annotations = annotations_exploded.groupby([col for col in annotations_exploded.columns if col != annotation_col], as_index = False, dropna = False)[annotation_col].apply(lambda x: ';'.join(set(x)))
    
    #get the number of annotations found
    annotation_counts = annotations_exploded[annotation_col].value_counts()

    #additional_counts
    sub_counts = []
    if 'Impact' in annotations_exploded.columns:
        for imp in ['Included', 'Excluded', 'Altered Flank']:
            tmp_annotations = annotations_exploded[annotations_exploded['Impact'] == imp].copy()
            sub_counts.append(tmp_annotations[annotation_col].value_counts())
    
        annotation_counts = pd.concat([annotation_counts] + sub_counts, axis = 1)
        annotation_counts.columns = ['All Impacted', 'Included', 'Excluded', 'Altered Flank']
        annotation_counts = annotation_counts.replace(np.nan, 0)
    
    #annotation_counts = annotation_counts.set_index(annotation_col)

    return annotations, annotation_counts

def plot_annotations(spliced_ptms, database = 'PhosphoSitePlus', annot_type = 'Function', collapse_on_similar = True, colors = None, top_terms = 5, legend = True, leg_loc = (-2.5,-0.6), ax = None):
    _, annotation_counts = get_ptm_annotations(spliced_ptms, annotation_type = annot_type, database = database, collapse_on_similar = collapse_on_similar)
    if ax is None:
        fig, ax = plt.subplots(figsize = (2,3))

    if database == 'PTMcode': #convert to readable gene name
        annotation_counts.index = [pose_config.uniprot_to_genename[i].split(' ')[0] if i in pose_config.uniprot_to_genename.keys() else i for i in annotation_counts.index]

    if colors is None:
        colors = ['lightgrey', 'gray', 'white']
    if isinstance(annotation_counts, pd.Series):
        annotation_counts = annotation_counts.head(top_terms).sort_values(ascending = True)
        ax.barh(annotation_counts.index, annotation_counts.values, color = colors[0], edgecolor = 'black')
    else:
        annotation_counts = annotation_counts.head(top_terms).sort_values(by = 'All Impacted', ascending = True)
        ax.barh(annotation_counts['Excluded'].index, annotation_counts['Excluded'].values, height = 1, edgecolor = 'black', color = colors[0])
        ax.barh(annotation_counts['Included'].index, annotation_counts['Included'].values, left = annotation_counts['Excluded'].values, height = 1, color = colors[1], edgecolor = 'black')
        ax.barh(annotation_counts['Altered Flank'].index, annotation_counts['Altered Flank'].values, left = annotation_counts['Excluded'].values+annotation_counts['Included'].values, height = 1, color = colors[2], edgecolor = 'black')
    #ax.set_xticks([0,50,100,150])
    ax.set_ylabel('', fontsize = 10)
    ax.set_xlabel('Number of PTMs', fontsize = 10)

    ax.set_title(f'Top {top_terms} {database} {annot_type} Annotations', fontsize = 10, weight = 'bold')

    #label_dict = {'EXONT:Name':'Exon Ontology Term', 'PSP:ON_PROCESS':'Biological Process (PSP)', 'PSP:ON_FUNCTION':'Molecular Function (PSP)', 'Combined:Kinase':'Kinase'}
    #ax.text(-1*ax.get_xlim()[1]/10, top_terms-0.2, label_dict[term_to_plot], weight = 'bold', ha = 'right', fontsize = 8)
    x_label_dict = {'Function':'Number of PTMs\nassociated with Function', 'Process':'Number of PTMs\nassociated with Process', 'Disease':'Number of PTMs\nassociated with Disease', 'Kinase':'Number of Phosphosites\ntargeted by Kinase', 'Interactions': 'Number of PTMs\nthat regulate interaction\n with protein','Motif Match':'Number of PTMs\nfound within a\nmotif instance', 'Intraprotein': 'Number of PTMs\nthat are important\for intraprotein\n interactions','Phosphatase':'Number of Phosphosites\ntargeted by Phosphatase'}
    ax.set_xlabel(x_label_dict[annot_type], fontsize = 8)
    
    #make a custom legend
    if legend:
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(facecolor = colors[0], edgecolor = 'black', label = 'Excluded'), mpatches.Patch(facecolor = colors[1], edgecolor = 'black', label = 'Included'),mpatches.Patch(facecolor = colors[2], edgecolor = 'black', label = 'Altered Flank')]
        ax.legend(handles = handles, ncol = 1, fontsize = 7, title = 'Type of Impact', title_fontsize = 8)

def get_annotation_categories(spliced_ptms):
    """
    Given spliced ptm information, return the available annotation categories that have been appended to dataframe

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        PTMs projected onto splicing events and with annotations appended from various databases

    Returns
    -------
    annot_categories: pd.DataFrame
        Dataframe that indicates the available databases, annotations from each database, and column associated with that annotation
    """
    database_list = []
    type_list = []
    column_list = []
    #get available phosphositeplus annotations
    for col in spliced_ptms.columns:
        if ':' in col:
            database = col.split(':')[0] if 'PSP' not in col else 'PhosphoSitePlus'
            if database != 'Combined':
                col_dict = pose_config.annotation_col_dict[database]

                #flip through annotation types in col_dict and add the one that matches the column
                for key, value in col_dict.items():
                    if value == col:
                        type_list.append(key)
                        database_list.append(database)
                        column_list.append(col)
            else:
                type_list.append(col.split(':')[1])
                database_list.append('Combined')
                column_list.append(col)

    if len(type_list) > 0:
        annot_categories = pd.DataFrame({'database':database_list, 'annotation_type':type_list, 'column': column_list}).sort_values(by = 'database')
        return annot_categories
    else:
        print('No annotation information found. Please run functions from annotate module to append annotation information')
        return None
    

def construct_background(annotation_type = 'Function', database = 'PhosphoSitePlus', file = None):
    ptm_coordinates = pose_config.ptm_coordinates.copy()
    if database == 'PhosphoSitePlus':
        if file is None:
            raise ValueError('Please provide PhosphoSitePlus source file to construct the background dataframe')
        elif annotation_type in ['Function', 'Process', 'Interactions']:
            ptm_coordinates = annotate.add_PSP_regulatory_site_data(ptm_coordinates, file = file)
        elif annotation_type == 'Kinase':
            ptm_coordinates = annotate.add_PSP_kinase_substrate_data(ptm_coordinates, file = file)
        elif annotation_type == 'Disease':
            ptm_coordinates = annotate.add_PSP_disease_association(ptm_coordinates, file = file)
    if database == 'ELM':
        if annotation_type == 'Interactions':
            ptm_coordinates = annotate.add_ELM_interactions(ptm_coordinates, file = file)
        elif annotation_type == 'Motif Match':
            ptm_coordinates = annotate.add_ELM_motif_matches(ptm_coordinates, file = file)
    if database == 'PTMInt':
        ptm_coordinates = annotate.add_PTMInt_data(ptm_coordinates, file = file)
    if database == 'PTMcode':
        if annotation_type == 'Intraprotein':
            ptm_coordinates = annotate.add_PTMcode_intraprotein(ptm_coordinates, file = file)
        elif annotation_type == 'Interactions':
            ptm_coordinates = annotate.add_PTMcode_interprotein(ptm_coordinates, file = file)
    if database == 'RegPhos':
        ptm_coordinates = annotate.add_RegPhos_data(ptm_coordinates, file = file)
    if database == 'DEPOD':
        ptm_coordinates = annotate.add_DEPOD_data(ptm_coordinates)
    

    _, annotation_counts = get_ptm_annotations(ptm_coordinates, annotation_type = annotation_type, database = database, collapse_on_similar = False)
    if (database == 'PhosphoSitePlus' and annotation_type in ['Function', 'Process', 'Disease']) or annotation_type == 'Interactions':
        _, annotation_counts_collapsed = get_ptm_annotations(ptm_coordinates, annotation_type = annotation_type, database = database, collapse_on_similar = True)
        return annotation_counts, annotation_counts_collapsed
    else:
        return annotation_counts
    
def get_enrichment_inputs(spliced_ptms, background_type = 'pregenerated', background = None, annotation_type = 'Function', database = 'PhosphoSitePlus', collapse_on_similar = False, mod_class = None, alpha = 0.05, min_dPSI = None, annotation_file = None):
    if background is None:
        if background_type == 'pregenerated':
            if mod_class is None:
                print('Using pregenerated background information on all PTMs in the proteome.')
                background_annotation_count = pose_config.download_background(annotation_type = annotation_type, database = database, collapsed=collapse_on_similar)
                background_size = pose_config.ptm_coordinates.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).shape[0]
            else:
                raise ValueError('To do modification specific analysis, must use significance to construct the background.')
        elif background_type == 'significance':
            background = spliced_ptms.copy()
            #restrict sample to significantly spliced ptms
            spliced_ptms = spliced_ptms[spliced_ptms['Significance'] <= alpha].copy()

            #restrict by min_dPSI, if provided
            if min_dPSI is not None:
                spliced_ptms = spliced_ptms[spliced_ptms['dPSI'].abs() >= min_dPSI]

            #check to make sure there are significant PTMs in the data and that there is a difference in the number of significant and background PTMs
            if spliced_ptms.shape[0] == 0:
                raise ValueError('No significantly spliced PTMs found in the data')
            elif spliced_ptms.shape[0] == background.shape[0]:
                raise ValueError(f'The foreground and background PTM sets are the same size when considering significance. Please provide a different background set with the background_ptms parameter, or make sure spliced_ptms also includes non-significant PTMs. Instead using pregenerated background sets of the whole proteome.')
            else:
                if mod_class is not None:
                    background = get_modification_class_data(background, mod_class)

            #get background counts
                background_size = background.shape[0]
            _, background_annotation_count = get_ptm_annotations(background, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar)
    elif background is not None: #if custom background is provided
        print('Using the provided custom background')
        if isinstance(background, list) or isinstance(background, np.ndarray):
            #from list of PTM strings, separate into uniprot id, residue, and position
            uniprot_id = [ptm.split('_')[0] for ptm in background]
            residue = [ptm.split('_')[1][0] for ptm in background]
            position = [int(ptm.split('_')[1][1:]) for ptm in background]
            background = pd.DataFrame({'UniProtKB Accession':uniprot_id, 'Residue':residue, 'PTM Position in Canonical Isoform':position, 'Modification Class':mod_class})
        if isinstance(background, pd.DataFrame):
            #check to make sure ptm data has key columns to identify ptms
            if 'UniProtKB Accession' not in background.columns or 'Residue' not in background.columns or 'PTM Position in Canonical Isoform' not in background.columns or 'Modification Class' not in background.columns:
                raise ValueError('Background dataframe must have UniProtKB Accession, Residue, PTM Position in Canonical Isoform, and Modification Class columns to identify PTMs')
            
            #restrict to specific modification class
            if mod_class is not None and 'Modification Class' in background.columns:
                background = get_modification_class_data(background, mod_class)
            elif mod_class is not None:
                raise ValueError('Custom background dataframe must have a Modification Class column to subset by modification class.')
        else:
            raise ValueError('Custom backgrounds must be provided as a list/array of PTMs in the form of "P00533_Y1068" (Uniprot ID followed by site number) or as a custom background dataframe with UniProtKB Accession, Residue, PTM Position in Canonical Isoform, and Modification Class columns.')
        
        background = annotate.add_annotation(background, annotation_type = annotation_type, database = database, check_existing = True, file = annotation_file)    
        background_size = background.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).shape[0]

        #get background counts
        _, background_annotation_count = get_ptm_annotations(background, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar)
    elif background_type == 'custom':
        raise ValueError('Please provide a custom background dataframe or list of PTMs to use as the background if wanting to use custom background data.')    
    else:
        raise ValueError('Invalid background type. Must be pregenerated, significance, or custom')

    #get counts
    foreground_size = spliced_ptms.shape[0]
    annotation_details, foreground_annotation_count = get_ptm_annotations(spliced_ptms, annotation_type = annotation_type, database = database, collapse_on_similar=collapse_on_similar)

    #process annotation details into usable format
    if annotation_details.shape[0] == 0:
        print('No spliced PTMs with an annotation, so could not perform enrichment analysis')
    else:
        annotation_col = get_annotation_col(spliced_ptms, database = database, annotation_type = annotation_type)
        annotation_details[annotation_col] = annotation_details[annotation_col].str.split(';')
        annotation_details = annotation_details.explode(annotation_col)
        annotation_details[annotation_col] = annotation_details[annotation_col].str.strip()
        annotation_details['PTM'] = annotation_details['Gene'] + '_' + annotation_details['Residue'] + annotation_details['PTM Position in Canonical Isoform'].astype(int).astype(str)
        annotation_details = annotation_details.groupby(annotation_col)['PTM'].agg(';'.join)
    
    return foreground_annotation_count, foreground_size, background_annotation_count, background_size, annotation_details


def annotation_enrichment(spliced_ptms, background = None, background_type = 'pregenerated', database = 'PhosphoSitePlus', annotation_type = 'Function', collapse_on_similar = False, mod_class = None, alpha = None, min_dPSI = None, annotation_file = None):#
    """
    In progress, needs to be tested

    Given spliced ptm information (differential inclusion, altered flanking sequences, or both), calculate the enrichment of specific annotations in the dataset using a hypergeometric test. Background data can be provided/constructed in a few ways:

    1. Use preconstructed background data for the annotation of interest, which considers the entire proteome present in the ptm_coordinates dataframe. While this is the default, it may not be the most accurate representation of your data, so you may alternative decide to use the other options which will be more specific to your context.
    2. Use the alpha and min_dPSI parameter to construct a foreground that only includes significantly spliced PTMs, and use the entire provided spliced_ptms dataframe as the background. This will allow you to compare the enrichment of specific annotations in the significantly spliced PTMs compared to the entire dataset. Will do this automatically if alpha or min_dPSI is provided.
    3. Provide a custom background dataframe that contains at least 3 columns: UniProtKB Accession, Residue, and PTM Position in Canonical Isoform. This will be used to construct the background information for the annotation of interest. If the background dataframe does not contain the necessary annotation information, the function will attempt to add it using the annotation_file parameter. If the annotation_file is not provided, the function will attempt to download the necessary data from the internet, if possible. 
    4. Provide a custom list of PTMs in the form of "P00533_Y1068" (Uniprot ID followed by site number). This will be used to construct the background information for the annotation of interest. Annotation information will be found for the provided list using the annotation_file parameter, if provided. If the annotation_file is not provided, the function will attempt to download the necessary data from the internet, if possible. 

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    background: pd.DataFrame, list, or np.array
        Background data to use for the annotation enrichment analysis. If None, will use the entire proteome present in the ptm_coordinates dataframe or the. Default is None.
    """
    foreground_annotation_count, foreground_size, background_annotations, background_size, annotation_details = get_enrichment_inputs(spliced_ptms, background_type = background_type, background = background, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar, mod_class = mod_class, alpha = alpha, min_dPSI = min_dPSI, annotation_file = annotation_file)
    


    #iterate through all annotations and calculate enrichment with a hypergeometric test
    results = pd.DataFrame(columns = ['Fraction Impacted', 'p-value'], index = foreground_annotation_count.index)
    for i, row in background_annotations.iterrows():
        #number of PTMs in the background with the annotation
        n = row['count']
        #number of PTMs in the foreground with the annotation
        if i in foreground_annotation_count.index.values:
            if foreground_annotation_count.shape[1] == 1:
                k = foreground_annotation_count.loc[i, 'count']
            elif foreground_annotation_count.shape[1] > 1:
                k = foreground_annotation_count.loc[i, 'All Impacted']

            p = stat_utils.getEnrichment(background_size, n, foreground_size, k, fishers = False)
            results.loc[i, 'Fraction Impacted'] = f"{k}/{n}"
            results.loc[i, 'p-value'] = p

    results = results.sort_values('p-value')
    results['Adjusted p-value'] = stat_utils.adjustP(results['p-value'].values)
    results = pd.concat([results, annotation_details], axis = 1)

    return results


def gene_set_enrichment(spliced_ptms = None, altered_flanks = None, combined = None, alpha = 0.05, min_dPSI = None, gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023','Reactome_2022'], background = None, return_sig_only = True):
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

    Returns
    -------
    results: pd.DataFrame
        Dataframe with gene set enrichment results from enrichr API

    """
    if combined is not None:
        if spliced_ptms is not None or altered_flanks is not None:
            print('If combined dataframe is provided, you do not need to include spliced_ptms or altered_flanks dataframes. Ignoring these inputs.')

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
    elif spliced_ptms is not None and altered_flanks is not None:
        #gene information (total and spliced genes)
        combined = combine_outputs(spliced_ptms, altered_flanks)
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
        type = 'Differentially Included'

        #isolate the type of impact on the gene
        altered_flank_only = []
        differential_only = spliced_ptms['Gene'].unique().tolist()
        both = []
    elif altered_flanks is not None:
        foreground = altered_flanks.copy()
        type = 'Altered Flanking Sequences'

        #isolate the type of impact on the gene
        altered_flank_only = altered_flanks['Gene'].unique().tolist()
        differential_only = []
        both = []
    else:
        raise ValueError('No dataframes provided. Please provide spliced_ptms, altered_flanks, or the combined dataframe.')
    
    #restrict to significant ptms, if available
    if 'Significance' in combined.columns and (min_dPSI is not None and 'dPSI' in foreground.columns):
        foreground = combined[combined['Significance'] <= alpha].copy()
        foreground = foreground[foreground['dPSI'].abs() >= min_dPSI]
    elif 'Significance' in combined.columns:
        foreground = combined[combined['Significance'] <= alpha].copy()
    elif min_dPSI is not None and 'dPSI' in combined.columns:
        foreground = combined[combined['dPSI'].abs() >= min_dPSI].copy()
    else:
        print('Significance column not found and min_dPSI not provided. All PTMs in dataframe will be considered as the foreground')

    foreground = foreground['Gene'].unique().tolist()   

    #construct background
    if isinstance(background, list):
        pass
    elif isinstance(background, np.ndarray):
        background = list(background)
    elif background == 'Significance' and 'Significance' in foreground.columns:
        background = combined.copy()
        background = background['Gene'].unique().tolist()   
    

    
    #perform gene set enrichment analysis and save data
    enr = gp.enrichr(foreground, background = background, gene_sets = gene_sets, organism='human')
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
    


def kstar_enrichment(spliced_ptms, networks):
    pass

def get_interaction_network(spliced_ptms, node_type = 'Gene'):
    if node_type not in ['Gene', 'PTM']:
        raise ValueError("node_type parameter (which dictates whether to consider interactions at PTM or gene level) can be either Gene or PTM")
    
    #extract interaction information in provided data
    interactions = annotate.combine_interaction_data(spliced_ptms)
    interactions['Residue'] = interactions['Residue'] + interactions['PTM Position in Canonical Isoform'].astype(int).astype(str)
    interactions = interactions.drop(columns = ['PTM Position in Canonical Isoform'])

    #add regulation change information
    if 'dPSI' in spliced_ptms.columns:
        interactions['Regulation Change'] = interactions.apply(lambda x: '+' if x['Type'] != 'DISRUPTS' and x['dPSI'] > 0 else '+' if x['Type'] == 'DISRUPTS' and x['dPSI'] < 0 else '-', axis = 1)
        grouping_cols = ['Residue', 'Type', 'Source', 'dPSI', 'Regulation Change']
        interactions['dPSI'] = interactions['dPSI'].apply(str)
    else:
        grouping_cols = ['Residue', 'Type', 'Source']

    #extract gene_specific network information
    if node_type == 'Gene':
        network_data = interactions.groupby(['Modified Gene', 'Interacting Gene'], as_index = False)[grouping_cols].agg(helpers.join_unique_entries)
        #generate network with all possible PTM-associated interactions
        interaction_graph = nx.from_pandas_edgelist(network_data, source = 'Modified Gene', target = 'Interacting Gene')
    else:
        interactions['Spliced PTM'] = interactions['Modified Gene'] + '_' + interactions['Residue']
        network_data = interactions.groupby(['Spliced PTM', 'Interacting Gene'], as_index = False)[grouping_cols].agg(helpers.join_unique_entries)
        network_data = network_data.drop(columns = ['Residue'])
        
        #generate network with all possible PTM-associated interactions
        interaction_graph = nx.from_pandas_edgelist(network_data, source = 'Spliced PTM', target = 'Interacting Gene')


    return interaction_graph, network_data


def get_interaction_stats(interaction_graph):
    """
    Given the networkx interaction graph, calculate various network centrality measures to identify the most relevant PTMs or genes in the network
    """
    #calculate network centrality measures
    degree_centrality = nx.degree_centrality(interaction_graph)
    closeness_centrality = nx.closeness_centrality(interaction_graph)
    betweenness_centrality = nx.betweenness_centrality(interaction_graph)
    network_stats = pd.DataFrame({'Degree': dict(interaction_graph.degree()), 'Degree Centrality':degree_centrality, 'Closeness':closeness_centrality,'Betweenness':betweenness_centrality})
    return network_stats

def get_protein_interaction_network(protein, network_data):
    protein_network = network_data[network_data['Modified Gene'] == protein]
    protein_network = protein_network.drop(columns = ['Modified Gene'])
    protein_network = protein_network.rename(columns = {'Residue': 'Spliced PTMs facilitating Interacting'})
    return protein_network

def summarize_protein_network(protein, interaction_graph, network_data, network_stats = None):
    """
    Given a protein of interest, summarize the network data for that protein
    """
    protein_network = network_data[network_data['Modified Gene'] == protein]
    increased_interactions = protein_network.loc[protein_network['Regulation Change'] == '+', 'Interacting Gene'].values
    decreased_interactions = protein_network.loc[protein_network['Regulation Change'] == '-', 'Interacting Gene'].values
    ambiguous_interactions = protein_network.loc[protein_network['Regulation Change'].str.contains(';'), 'Interacting Gene'].values

    #print interactions
    if len(increased_interactions) > 0:
        print(f'Increased interaction likelihoods: {', '.join(increased_interactions)}')
    if len(decreased_interactions) > 0:
        print(f'Decreased interaction likelihoods: {', '.join(decreased_interactions)}')
    if len(ambiguous_interactions) > 0:
        print(f'Ambiguous interaction impact: {', '.join(ambiguous_interactions)}')

    if network_stats is None:
        network_stats = get_interaction_stats(interaction_graph)

    network_ranks = network_stats.rank(ascending = False).astype(int)
    print(f'Number of interactions: {network_stats.loc[protein, "Degree"]} (Rank: {network_ranks.loc[protein, "Degree"]})')
    print(f'Centrality measures - \t Degree = {network_stats.loc[protein, "Degree Centrality"]} (Rank: {network_ranks.loc[protein, "Degree Centrality"]})')
    print(f'                      \t Betweenness = {network_stats.loc[protein, "Betweenness"]} (Rank: {network_ranks.loc[protein, "Betweenness"]})')
    print(f'                      \t Eigenvector = {network_stats.loc[protein, "Eigenvector"]} (Rank: {network_ranks.loc[protein, "Eigenvector"]})')
    print(f'                      \t Closeness = {network_stats.loc[protein, "Closeness"]} (Rank: {network_ranks.loc[protein, "Closeness"]})')
    




def edit_sequence_for_kinase_library(seq):
    """
    Convert flanking sequence to version accepted by kinase library (modified residue denoted by asterick)
    """
    if seq == seq:
        seq = seq.replace('t','t*')
        seq = seq.replace('s','s*')
        seq = seq.replace('y','y*')
    else:
        return np.nan
    return seq

def process_data_for_kinase_library(altered_flanks, odir):
    """
    Extract flanking sequence information for 
    """
    #restrict to events with changed flanking sequences, no introduced stop codons, and phosphorylation modifications
    flank_data = altered_flanks[~altered_flanks['Matched']] 
    flank_data = flank_data[~flank_data['Stop Codon Introduced']]
    flank_data = flank_data[flank_data['Modification Class'].str.contains('Phosphorylation')]

    #generate files to input into Kinase Library (inclusion first then exclusion)
    inclusion_sequences = flank_data[['PTM', 'Inclusion Flanking Sequence']].drop_duplicates()
    inclusion_sequences['Inclusion Flanking Sequence'] = inclusion_sequences['Inclusion Flanking Sequence'].apply(edit_sequence_for_kinase_library)
    inclusion_sequences = inclusion_sequences.dropna(subset = 'Inclusion Flanking Sequence')
    #write sequences to text file
    with open(odir + 'inclusion_sequences_input.txt', 'w') as f:
        for _, row in inclusion_sequences.iterrows():
            f.write(row['Inclusion Flanking Sequence']+'\n')

    exclusion_sequences = flank_data[['PTM', 'Exclusion Flanking Sequence']].drop_duplicates()
    exclusion_sequences['Exclusion Flanking Sequence'] = exclusion_sequences['Exclusion Flanking Sequence'].apply(edit_sequence_for_kinase_library)
    exclusion_sequences = exclusion_sequences.dropna(subset = 'Exclusion Flanking Sequence')
    #write sequences to text file
    with open(odir + 'exclusion_sequences_input.txt', 'w') as f:
        for _, row in exclusion_sequences.iterrows():
            f.write(row['Exclusion Flanking Sequence']+'\n')


def process_kinase_library_output(altered_flanks, scores, flanking_sequence_col = 'Inclusion Flanking Sequence'):
    """
    Process output from Kinase Library to connect kinase library scores back to the PTMs in the altered flanks dataframe

    Parameters
    ----------
    altered_flanks: pd.DataFrame
        Dataframe with PTMs associated with altered flanking sequences
    scores: pd.DataFrame
        Dataframe with kinase library scores for flanking sequences (loaded from downloaded .tsv outputs from kinase library)
    flanking_sequence_col: str
        Column in altered_flanks dataframe that contains the flanking sequence to match with the kinase library scores. Default is 'Inclusion Flanking Sequence'. Can also be 'Exclusion Flanking Sequence'

    Returns
    -------
    percentiles_y: pd.DataFrame
        Dataframe with kinase library scores for tyrosine sites
    percentiles_st: pd.DataFrame
        Dataframe with kinase library scores for serine/threonine sites

    """
    #restrict to events with changed flanking sequences, no introduced stop codons, and phosphorylation modifications
    flank_data = altered_flanks[~altered_flanks['Matched']] 
    flank_data = flank_data[~flank_data['Stop Codon Introduced']]
    flank_data = flank_data[flank_data['Modification Class'].str.contains('Phosphorylation')]

    
    sequences = flank_data[['Region ID','PTM',  flanking_sequence_col]].drop_duplicates()
    sequences = sequences.dropna(subset =  flanking_sequence_col)
    sequences['Label'] = sequences['Region ID'] + ';' + sequences['PTM']
    sequences[flanking_sequence_col] = sequences[ flanking_sequence_col].apply(lambda x: x.upper().replace(' ', '_')+'_')

    sequences = sequences.merge(scores, left_on = 'Exclusion Flanking Sequence', right_on = 'sequence', how = 'left')
    #split info into tyrosine vs. serine/threonine
    sequences_y = sequences[sequences['Label'].str.contains('_Y')]
    sequences_st = sequences[(sequences['Label'].str.contains('_S')) | (sequences['Label'].str.contains('_T'))]

    #pivot table to get scores for each kinase
    percentiles_y = sequences_y.pivot_table(index = 'Label', columns = 'kinase', values = 'site_percentile')
    percentiles_st = sequences_st.pivot_table(index = 'Label', columns = 'kinase', values = 'site_percentile')

    return percentiles_y, percentiles_st

def get_kinase_library_differences(altered_flanks, inclusion_scores, exclusion_scores):
    """
    Given altered flanking sequences and kinase library scores for inclusion and exclusion flanking sequences, calculate the difference in kinase library site percentiles between the two

    Parameters
    ----------
    altered_flanks: pd.DataFrame
        Dataframe with PTMs associated with altered flanking sequences
    inclusion_scores: pd.DataFrame
        Dataframe with kinase library scores for inclusion flanking sequences (loaded from downloaded .tsv outputs from kinase library)
    exclusion_scores: pd.DataFrame
        Dataframe with kinase library scores for exclusion flanking sequences (loaded from downloaded .tsv outputs from kinase library)
    
    Returns
    -------
    percentiles_diff_y: pd.DataFrame
        Dataframe with the difference in kinase library scores for tyrosine sites
    percentiles_diff_st: pd.DataFrame
        Dataframe with the difference in kinase library scores for serine/threonine sites
    """
    inclusion_percentiles_y, inclusion_percentiles_st = process_kinase_library_output(altered_flanks, inclusion_scores, flanking_sequence_col = 'Inclusion Flanking Sequence')
    exclusion_percentiles_y, exclusion_percentiles_st = process_kinase_library_output(altered_flanks, exclusion_scores, flanking_sequence_col = 'Exclusion Flanking Sequence')

    #calculate the difference in percentiles
    labels= list(set(inclusion_percentiles_y.index).intersection(exclusion_percentiles_y.index))
    percentiles_diff_y = inclusion_percentiles_y.loc[labels].copy()
    percentiles_diff_y = percentiles_diff_y[exclusion_percentiles_y.columns]
    for i, row in percentiles_diff_y.iterrows():
        percentiles_diff_y.loc[i] = row - exclusion_percentiles_y.loc[i]

    labels= list(set(inclusion_percentiles_st.index).intersection(exclusion_percentiles_st.index))
    percentiles_diff_st = inclusion_percentiles_st.loc[labels].copy()
    percentiles_diff_st = percentiles_diff_st[exclusion_percentiles_st.columns]
    for i, row in percentiles_diff_st.iterrows():
        percentiles_diff_st.loc[i] = row - exclusion_percentiles_st.loc[i]

    return percentiles_diff_y, percentiles_diff_st
    
def process_data_for_exon_ontology(odir, spliced_ptms = None, altered_flanks = None):
    pass


class POSE_Analysis:
    def __init__(self, spliced_ptms = None, altered_flanks = None, alpha = 0.05, min_dPSI = 0):
        self.spliced_ptms = spliced_ptms
        self.altered_flanks = altered_flanks
        self.alpha = alpha
        self.min_dPSI = min_dPSI

    def load(self, idir):
        """
        Load data generated and saved by POSE_Project object which includes information about spliced PTMs and altered flanking sequences

        Parameters
        ----------
        idir: str
            Directory where data is saved
        """
        annotated_source_data = {}
        for _, _, file in os.walk(idir):
            if 'spliced_ptms.csv' in file:
                print('Spliced PTMs found in directory')
                self.spliced_ptms = pd.read_csv(file)
            elif 'altered_flanks.csv' in file:
                print('Altered flanking sequences found in directory')
                self.altered_flanks = pd.read_csv(file)
            elif 'source_data.csv' in file:
                file_label = file.split('/')[-1].split('_')[0]
                if file_label != 'annotated':
                    print(f'{file_label} source data found in directory')
                else:
                    print('Annotated source data found in directory')
                annotated_source_data[file_label] = pd.read_csv(file)
        if len(annotated_source_data) == 1:
            self.annotated_source_data = annotated_source_data[file_label]

        if self.spliced_ptms is not None and self.altered_flanks is not None:
            self.combined_ptms = combine_outputs(self.spliced_ptms, self.altered_flanks)

    def add_ptm_annotations(self, psp_regulatory_site_file = None, psp_ks_file = None, psp_disease_file = None, elm_interactions = False, elm_motifs = False, PTMint = False, PTMcode_intraprotein = False, PTMcode_interprotein = False, DEPOD = False, RegPhos = False, combine_similar = True):
        """
        Annotate spliced PTMs and altered flanking sequences with information from various databases

        Parameters
        ----------
        annotation_types: list
            List of databases to annotate data with. Default is ['PhosphoSitePlus', 'ELM', 'PTMInt', 'PTMcode', 'RegPhos', 'DEPOD']
        annotation_files: dict
            Dictionary with database names as keys and file paths as values to use for annotation. Default is None (will attempt to download data from the internet)
        """
        #add annotations to spliced_ptms
        if self.spliced_ptms is not None:
            self.spliced_ptms = annotate.annotate_ptms(self.spliced_ptms, psp_regulatory_site_file = psp_regulatory_site_file, psp_ks_file = psp_ks_file, psp_disease_file = psp_disease_file, elm_interactions = elm_interactions, elm_motifs = elm_motifs, PTMint = PTMint, PTMcode_intraprotein = PTMcode_intraprotein, PTMcode_interprotein = PTMcode_interprotein, DEPOD = DEPOD, RegPhos = RegPhos, combine_similar = combine_similar)
        #add annotations to altered_flanks
        if self.altered_flanks is not None:
            self.altered_flanks = annotate.annotate_ptms(self.altered_flanks, psp_regulatory_site_file = psp_regulatory_site_file, psp_ks_file = psp_ks_file, psp_disease_file = psp_disease_file, elm_interactions = elm_interactions, elm_motifs = elm_motifs, PTMint = PTMint, PTMcode_intraprotein = PTMcode_intraprotein, PTMcode_interprotein = PTMcode_interprotein, DEPOD = DEPOD, RegPhos = RegPhos, combine_similar = combine_similar)

        #recombine dataframes
        if self.spliced_ptms is not None and self.altered_flanks is not None:
            self.combined_ptms = combine_outputs(self.spliced_ptms, self.altered_flanks)

    def get_interaction_network(self, node_type = 'Gene'):
        if self.combined_ptms is not None:
            self.interaction_graph, self.network_data = get_interaction_network(self.combined_ptms, node_type = node_type)
        else:
            raise ValueError('Combined PTM data not found. Please combine spliced PTMs and altered flanking sequences before generating interaction network')


    


