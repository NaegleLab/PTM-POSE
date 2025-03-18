import numpy as np
import pandas as pd
import pickle

import os
import time

#plotting 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from ptm_pose import plots as pose_plots

#analysis packages
from Bio.Align import PairwiseAligner
import gseapy as gp
import networkx as nx
import re


#custom stat functions
from ptm_pose import stat_utils, pose_config, annotate, helpers
from ptm_pose.analyze import summarize


package_dir = os.path.dirname(os.path.abspath(__file__))



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
            raise ValueError(f'Requested annotation data has not yet been added to spliced_ptms dataframe. Please run the annotate.{pose_config.annotation_function_dict[database][annotation_type]} function to append this information.')
        return annotation_col
    else:
        raise ValueError(f"Invalid annotation type for {database}. Available annotation data for {database} includes: {', '.join(pose_config.annotation_col_dict[database].keys())}")


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
    sep_dict = {'PhosphoSitePlus':{'Function':',', 'Process':',','Interactions':'(', 'Disease':'->', 'Perturbation':'->'}, 'ELM': {'Interactions': ' ', 'Motif Match': ' '}, 'PTMInt':{'Interactions':'->'}, 'PTMcode':{'Interactions':'_', 'Intraprotein':' '}, 'RegPhos':{'Kinase':' '}, 'DEPOD':{'Phosphatase':' '}, 'Combined':{'Kinase':' ', 'Interactions':'->'}, 'PTMsigDB': {'WikiPathway':'->', 'NetPath':'->','mSigDB':'->', 'Perturbation (DIA2)':'->', 'Perturbation (DIA)': '->', 'Perturbation (PRM)':'->','Kinase':'->'}}
    
    if annotation_type == 'Kinase' and database != 'PTMsigDB':
        collapsed = annotations
    else:
        sep = sep_dict[database][annotation_type]
        collapsed = []
        for annot in annotations:
            if annot == annot:
                collapsed.append(simplify_annotation(annot, sep = sep))
            else:
                collapsed.append(annot)
    return collapsed


def get_ptm_annotations(ptms, annotation_type = 'Function', database = 'PhosphoSitePlus', collapse_on_similar = False, dPSI_col = None, sig_col = None, **kwargs):
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
        annotation_col = get_annotation_col(ptms, database = database, annotation_type = annotation_type)
    else:
        annotation_col = f'Combined:{annotation_type}'


    #check to make sure requested annotation is available
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        helpers.check_filter_kwargs(filter_arguments)
        ptms_of_interest = helpers.filter_ptms(ptms, **filter_arguments)
        

    #extract relevant annotation and remove PTMs without an annotation
    if 'Impact' not in ptms_of_interest.columns and 'dPSI' in ptms_of_interest.columns:
        ptms_of_interest['Impact'] = ptms_of_interest['dPSI'].apply(lambda x: 'Included' if x > 0 else 'Excluded')

    optional_cols = [col for col in ptms_of_interest.columns if col in ['Impact', 'dPSI', 'Significance'] or col == dPSI_col or col == sig_col ]
    annotations = ptms_of_interest[['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class'] + [annotation_col] + optional_cols].copy()
    annotations = annotations.dropna(subset = annotation_col).drop_duplicates()

    if annotations.empty:
        print("No PTMs with associated annotation")
        return None, None
    
    #combine repeat entries for same PTM (with multiple impacts)
    annotations = annotations.groupby(['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'], as_index = False).agg(lambda x: ';'.join(set([str(i) for i in x if i == i])))

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
    annotation_counts = annotations_exploded.drop_duplicates(subset = ['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'] + [annotation_col])[annotation_col].value_counts()

    #additional_counts
    sub_counts = []
    if 'Impact' in annotations_exploded.columns:
        for imp in ['Included', 'Excluded', 'Altered Flank']:
            tmp_annotations = annotations_exploded[annotations_exploded['Impact'] == imp].copy()
            tmp_annotations = tmp_annotations.drop_duplicates(subset = ['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'] + [annotation_col])
            sub_counts.append(tmp_annotations[annotation_col].value_counts())
    
        annotation_counts = pd.concat([annotation_counts] + sub_counts, axis = 1)
        annotation_counts.columns = ['All Impacted', 'Included', 'Excluded', 'Altered Flank']
        annotation_counts = annotation_counts.replace(np.nan, 0)
    
    #combine repeat entries for same PTM (with multiple impacts)
    annotations = annotations.groupby(['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform'], as_index = False).agg(lambda x: ';'.join(set([str(i) for i in x if i == i])))

    return annotations, annotation_counts

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
            if database != 'Combined' and database != 'Unnamed':
                col_dict = pose_config.annotation_col_dict[database]

                #flip through annotation types in col_dict and add the one that matches the column
                for key, value in col_dict.items():
                    if value == col:
                        type_list.append(key)
                        database_list.append(database)
                        column_list.append(col)
            elif database == 'Combined':
                type_list.append(col.split(':')[1])
                database_list.append('Combined')
                column_list.append(col)
            else:
                continue

    if len(type_list) > 0:
        annot_categories = pd.DataFrame({'database':database_list, 'annotation_type':type_list, 'column': column_list}).sort_values(by = 'database')
        return annot_categories
    else:
        print('No annotation information found. Please run functions from annotate module to append annotation information')
        return None


def construct_background(file = None, annotation_type = 'Function', database = 'PhosphoSitePlus', modification = None, collapse_on_similar = False, save = False):
    ptm_coordinates = pose_config.ptm_coordinates.copy()
    ptm_coordinates = ptm_coordinates.rename({'Gene name':'Gene'}, axis = 1)
    if modification is not None:
        ptm_coordinates = ptm_coordinates[ptm_coordinates['Modification Class'].str.contains(modification)].copy()
        if ptm_coordinates.empty:
            raise ValueError(f'No PTMs found with modification class {modification}. Please provide a valid modification class. Examples include Phosphorylation, Glycosylation, Ubiquitination, etc.')
    
        
    if database == 'PhosphoSitePlus':
        if file is None:
            raise ValueError('Please provide PhosphoSitePlus source file to construct the background dataframe')
        elif annotation_type in ['Function', 'Process', 'Interactions']:
            ptm_coordinates = annotate.add_PSP_regulatory_site_data(ptm_coordinates, file = file, report_success=False)
        elif annotation_type == 'Kinase':
            ptm_coordinates = annotate.add_PSP_kinase_substrate_data(ptm_coordinates, file = file, report_success=False)
        elif annotation_type == 'Disease':
            ptm_coordinates = annotate.add_PSP_disease_association(ptm_coordinates, file = file, report_success=False)
        elif annotation_type == 'Perturbation':
            ptm_coordinates = annotate.add_PTMsigDB_data(ptm_coordinates, file = file, report_success=False)
    if database == 'ELM':
        if annotation_type == 'Interactions':
            ptm_coordinates = annotate.add_ELM_interactions(ptm_coordinates, file = file, report_success = False)
        elif annotation_type == 'Motif Match':
            ptm_coordinates = annotate.add_ELM_matched_motifs(ptm_coordinates, file = file, report_success = False)
    if database == 'PTMInt':
        ptm_coordinates = annotate.add_PTMInt_data(ptm_coordinates, file = file, report_success=False)
    if database == 'PTMcode':
        if annotation_type == 'Intraprotein':
            ptm_coordinates = annotate.add_PTMcode_intraprotein(ptm_coordinates, file = file, report_success=False)
        elif annotation_type == 'Interactions':
            ptm_coordinates = annotate.add_PTMcode_interprotein(ptm_coordinates, file = file, report_success=False)
    if database == 'RegPhos':
        ptm_coordinates = annotate.add_RegPhos_data(ptm_coordinates, file = file, report_success=False)
    if database == 'DEPOD':
        ptm_coordinates = annotate.add_DEPOD_phosphatase_data(ptm_coordinates, report_success=False)
    if database == 'PTMsigDB':
        ptm_coordinates = annotate.add_PTMsigDB_data(ptm_coordinates, file = file, report_success=False)
    if database == 'Combined':
        raise ValueError('Combined information is not supported for constructing background data from entire proteome at this time. Please provide a specific database to construct background data.')
    

    _, annotation_counts = get_ptm_annotations(ptm_coordinates, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar)
    if save:
        package_dir = os.path.dirname(os.path.abspath(__file__))
        if collapse_on_similar and modification is not None:
            annotation_counts.to_csv(package_dir + f'/Resource_Files/background_annotations/{database}_{annotation_type}_{modification}_collapsed.csv')
        elif collapse_on_similar:
            annotation_counts.to_csv(package_dir + f'/Resource_Files/background_annotations/{database}_{annotation_type}_collapsed.csv')
        elif modification is not None:
            annotation_counts.to_csv(package_dir + f'/Resource_Files/background_annotations/{database}_{annotation_type}_{modification}.csv')
        else:
            annotation_counts.to_csv(package_dir + f'/Resource_Files/background_annotations/{database}_{annotation_type}.csv')

    return annotation_counts

def get_enrichment_inputs(ptms,  annotation_type = 'Function', database = 'PhosphoSitePlus', background_type = 'pregenerated', background = None, collapse_on_similar = False, mod_class = None, alpha = 0.05, min_dPSI = 0, annotation_file = None, save_background = False):
    """
    Given the spliced ptms, altered_flanks, or combined PTMs dataframe, identify the number of PTMs corresponding to specific annotations in the foreground (PTMs impacted by splicing) and the background (all PTMs in the proteome or all PTMs in dataset not impacted by splicing). This information can be used to calculate the enrichment of specific annotations among PTMs impacted by splicing. Several options are provided for constructing the background data: pregenerated (based on entire proteome in the ptm_coordinates dataframe) or significance (foreground PTMs are extracted from provided spliced PTMs based on significance and minimum delta PSI)

    Parameters
    ----------
    spliced_ptms: pd.DataFrame

    """
    if background_type == 'pregenerated':
        print('Using pregenerated background information on all PTMs in the proteome.')
        #first look for pregenerated background data
        try:
            background_annotation_count = pose_config.download_background(annotation_type = annotation_type, database = database, mod_class = mod_class, collapsed=collapse_on_similar)
        except:
            if annotation_file is None:
                print('Note: To avoid having to constructing background each time (which is slower), you can choose to set save_background = True to save the background data to Resource Files in package directory.')
            background_annotation_count = construct_background(file = annotation_file, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar, save = save_background)

        if mod_class is None:
            background_size = pose_config.ptm_coordinates.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0]
        else:
            background_size = pose_config.ptm_coordinates[pose_config.ptm_coordinates['Modification Class'].str.contains(mod_class)].drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0]

    elif background_type == 'significance':
        if 'Significance' not in spliced_ptms.columns or 'dPSI' not in spliced_ptms.columns:
            raise ValueError('Significance and dPSI columns must be present in spliced_ptms dataframe to construct a background based on significance (these columns must be provided during projection).')
        
        background = spliced_ptms.copy()
        #restrict sample to significantly spliced ptms
        spliced_ptms = spliced_ptms[(spliced_ptms['Significance'] <= alpha) & (spliced_ptms['dPSI'].abs() >= min_dPSI)].copy()


        #check to make sure there are significant PTMs in the data and that there is a difference in the number of significant and background PTMs
        if spliced_ptms.shape[0] == 0:
            raise ValueError('No significantly spliced PTMs found in the data')
        elif spliced_ptms.shape[0] == background.shape[0]:
            raise ValueError(f'The foreground and background PTM sets are the same size when considering significance. Please provide a different background set with the background_ptms parameter, or make sure spliced_ptms also includes non-significant PTMs. Instead using pregenerated background sets of the whole proteome.')
        else:
            if mod_class is not None:
                background = summarize.get_modification_class_data(background, mod_class)
            #remove impact from background columns if present
            cols_to_remove = [col for col in ['dPSI', 'Significance', 'Impact'] if col in background.columns]
            if len(cols_to_remove) > 0:
                background = background.drop(columns = cols_to_remove)
        #get background counts
            background_size = background.shape[0]
        _, background_annotation_count = get_ptm_annotations(background, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar)
    #elif background is not None: #if custom background is provided
    #    print('Using the provided custom background')
    #    if isinstance(background, list) or isinstance(background, np.ndarray):
    #        #from list of PTM strings, separate into uniprot id, residue, and position
    #        uniprot_id = [ptm.split('_')[0] for ptm in background]
    #        residue = [ptm.split('_')[1][0] for ptm in background]
    #        position = [int(ptm.split('_')[1][1:]) for ptm in background]
    #        background = pd.DataFrame({'UniProtKB Accession':uniprot_id, 'Residue':residue, 'PTM Position in Isoform':position, 'Modification Class':mod_class})
    #    if isinstance(background, pd.DataFrame):
    #        #check to make sure ptm data has key columns to identify ptms
    #        if 'UniProtKB Accession' not in background.columns or 'Residue' not in background.columns or 'PTM Position in Isoform' not in background.columns or #'Modification Class' not in background.columns:
    #            raise ValueError('Background dataframe must have UniProtKB Accession, Residue, PTM Position in Isoform, and Modification Class columns to identify PTMs')
            
    #        #restrict to specific modification class
    #        if mod_class is not None and 'Modification Class' in background.columns:
    #            background = get_modification_class_data(background, mod_class)
    #        elif mod_class is not None:
    #            raise ValueError('Custom background dataframe must have a Modification Class column to subset by modification class.')
    #    else:
    #        raise ValueError('Custom backgrounds must be provided as a list/array of PTMs in the form of "P00533_Y1068" (Uniprot ID followed by site number) or as a custom background dataframe with UniProtKB Accession, Residue, PTM Position in Isoform, and Modification Class columns.')
        
    #    background = annotate.add_annotation(background, annotation_type = annotation_type, database = database, check_existing = True, file = annotation_file)    
    #    background_size = background.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0]

        #get background counts
    #    _, background_annotation_count = get_ptm_annotations(background, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar)
    #elif background_type == 'custom':
    #    raise ValueError('Please provide a custom background dataframe or list of PTMs to use as the background if wanting to use custom background data.')    
    else:
        raise ValueError('Invalid background type. Must be pregenerated (default) or significance')

    #get counts
    foreground_size = spliced_ptms.shape[0]
    annotation_details, foreground_annotation_count = get_ptm_annotations(spliced_ptms, annotation_type = annotation_type, database = database, collapse_on_similar=collapse_on_similar)

    #process annotation details into usable format
    if annotation_details is None:
        print('No PTMs with requested annotation type, so could not perform enrichment analysis')
        return np.repeat(None, 5)
    else:
        annotation_col = get_annotation_col(spliced_ptms, database = database, annotation_type = annotation_type)
        annotation_details[annotation_col] = annotation_details[annotation_col].str.split(';')
        annotation_details = annotation_details.explode(annotation_col)
        annotation_details[annotation_col] = annotation_details[annotation_col].str.strip()
        annotation_details['PTM'] = annotation_details['Gene'] + '_' + annotation_details['Residue'] + annotation_details['PTM Position in Isoform'].astype(int).astype(str)
        annotation_details = annotation_details.groupby(annotation_col)['PTM'].agg(';'.join)
    
    return foreground_annotation_count, foreground_size, background_annotation_count, background_size, annotation_details


def annotation_enrichment(ptms, database = 'PhosphoSitePlus', annotation_type = 'Function', background_type = 'pregenerated', collapse_on_similar = False, mod_class = None, alpha = None, min_dPSI = None, annotation_file = None, save_background = False, **kwargs):#
    """
    In progress, needs to be tested

    Given spliced ptm information (differential inclusion, altered flanking sequences, or both), calculate the enrichment of specific annotations in the dataset using a hypergeometric test. Background data can be provided/constructed in a few ways:

    1. Use preconstructed background data for the annotation of interest, which considers the entire proteome present in the ptm_coordinates dataframe. While this is the default, it may not be the most accurate representation of your data, so you may alternative decide to use the other options which will be more specific to your context.
    2. Use the alpha and min_dPSI parameter to construct a foreground that only includes significantly spliced PTMs, and use the entire provided spliced_ptms dataframe as the background. This will allow you to compare the enrichment of specific annotations in the significantly spliced PTMs compared to the entire dataset. Will do this automatically if alpha or min_dPSI is provided.

    Parameters
    ----------
    ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    database: str
        database from which PTMs are pulled. Options include 'PhosphoSitePlus', 'ELM', 'PTMInt', 'PTMcode', 'DEPOD', 'RegPhos', 'PTMsigDB'. Default is 'PhosphoSitePlus'.
    annotation_type: str
        Type of annotation to pull from spliced_ptms dataframe. Available information depends on the selected database. Default is 'Function'.
    background_type: str
        how to construct the background data. Options include 'pregenerated' (default) and 'significance'. If 'significance' is selected, the alpha and min_dPSI parameters must be provided. Otherwise, will use whole proteome in the ptm_coordinates dataframe as the background.
    collapse_on_similar: bool
        Whether to collapse similar annotations (for example, increasing and decreasing functions) into a single category. Default is False.
    mod_class: str
        modification class to subset, if any
    alpha: float
        significance threshold to use to subset foreground PTMs. Default is None.
    min_dPSI: float
        minimum delta PSI value to use to subset foreground PTMs. Default is None.
    annotation_file: str
        file to use to annotate custom background data. Default is None.
    save_background: bool
        Whether to save the background data constructed from the ptm_coordinates dataframe into Resource_Files within package. Default is False.
    """
    
    foreground_annotation_count, foreground_size, background_annotations, background_size, annotation_details = get_enrichment_inputs(ptms, background_type = background_type, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar, mod_class = mod_class, alpha = alpha, min_dPSI = min_dPSI, annotation_file = annotation_file, save_background = save_background)
    

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


def plot_available_annotations(spliced_ptms, show_all_ptm_count = True, ax = None):
    """
    Given a dataframe with ptm annotations added, show the number of PTMs associated with each annotation type

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs and annotations added
    show_all_ptm_count: bool
        Whether to show the total number of PTMs in the dataset. Default is True.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize = (5,5))

    if show_all_ptm_count:
        num_ptms = [spliced_ptms.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0]]
        num_ptms_filters = ['All PTMs']
        filter_source = ['None']
    else:
        num_ptms = []
        num_ptms_filters = []
        filter_source = []

    #look for annotations and add counts to lists
    ylabel_dict = {'PSP:ON_PROCESS':'Biological Process (PSP)', 'PSP:ON_FUNCTION':'Molecular Function (PSP)', 'PSP:Kinase':'Kinase (PSP)', 'PSP:Disease_Association':'Disease Association (PSP)', 'PSP:ON_PROT_INTERACT':'Interactions (PSP)', 'PSP:ON_OTHER_INTERACT':'Nonprotein Interactions (PSP)', 'ELM:Interactions':'Interactions (ELM)', 'ELM:Motif Matches':'Motif Matches (ELM)', 'PTMInt:Interaction':'Interactions (PTMInt)', 'PTMcode:Intraprotein_Interactions':'Intraprotein (PTMcode)','PTMcode:Interprotein_Interactions':'Interactions (PTMcode)', 'DEPOD:Phosphatase':'Phosphatase (DEPOD)', 'RegPhos:Kinase':'Kinase (RegPhos)', 'Combined:Kinase':'Kinase (Combined)', 'Combined:Interactions':'Interactions (Combined)'}
    available_annotations = [col for col in spliced_ptms.columns if 'Combined' in col or 'PSP:' in col or 'ELM:Interactions' in col or 'PTMInt:' in col or 'PTMcode:' in col or 'DEPOD:' in col or 'RegPhos:' in col]
    for annotation in available_annotations:
        num_ptms.append(spliced_ptms.dropna(subset = annotation).drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0])
        num_ptms_filters.append(ylabel_dict[annotation])
        filter_source.append(annotation.split(':')[0])

    
    #plot bar plot
    #color bars based on datasource
    palette = {'None': 'gray', 'PSP': 'blue', 'ELM': 'green', 'PTMInt':'red', 'PTMcode':'purple', 'DEPOD':'orange', 'RegPhos':'gold', 'Combined':'black'}
    colors = []
    for source in filter_source:
        colors.append(palette[source])

    ax.barh(num_ptms_filters[::-1], num_ptms[::-1], color = colors[::-1])
    ax.set_xlabel('Number of PTMs with annotation')
    
    #annotate with number of PTMs
    for i, num_ptm in enumerate(num_ptms[::-1]):
        ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')

    #create legend
    handles = [plt.Rectangle((0,0),1,1, color = color) for color in palette.values() if color != 'gray']
    labels = [source for source in palette.keys() if source != 'None']
    ax.legend(handles, labels, title = 'Annotation Source')
    plt.show()

def plot_annotation_counts(spliced_ptms, database = 'PhosphoSitePlus', annot_type = 'Function', collapse_on_similar = True, colors = None, fontsize = 10, top_terms = 5, legend = True, legend_loc = None, title = None, title_type = 'database', ax = None, **kwargs):
    """
    Given a dataframe with PTM annotations added, plot the top annotations associated with the PTMs

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs and annotations added
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
        Location of the legend. Default is None, which will place the legend in the upper right corner.
    **kwargs: additional keyword arguments
        Additional keyword arguments, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    _, annotation_counts = get_ptm_annotations(spliced_ptms, annotation_type = annot_type, database = database, collapse_on_similar = collapse_on_similar, kwargs = kwargs)
    if annotation_counts is None:
        return None


    if ax is None:
        fig, ax = plt.subplots(figsize = (2,3))

    if database == 'PTMcode': #convert to readable gene name
        annotation_counts.index = [pose_config.uniprot_to_genename[i].split(' ')[0] if i in pose_config.uniprot_to_genename.keys() else i for i in annotation_counts.index]

    if colors is None:
        colors = ['lightgrey', 'gray', 'white']

    if isinstance(annotation_counts, pd.Series):
        annotation_counts = annotation_counts.head(top_terms).sort_values(ascending = True)
        if isinstance(colors, list) or isinstance(colors, np.ndarray):
            colors = colors[0]
        
        ax.barh(annotation_counts.index, annotation_counts.values, color = colors, edgecolor = 'black')
        legend = False
    else:
        annotation_counts = annotation_counts.head(top_terms).sort_values(by = 'All Impacted', ascending = True)
        ax.barh(annotation_counts['Excluded'].index, annotation_counts['Excluded'].values, height = 1, edgecolor = 'black', color = colors[0])
        ax.barh(annotation_counts['Included'].index, annotation_counts['Included'].values, left = annotation_counts['Excluded'].values, height = 1, color = colors[1], edgecolor = 'black')
        ax.barh(annotation_counts['Altered Flank'].index, annotation_counts['Altered Flank'].values, left = annotation_counts['Excluded'].values+annotation_counts['Included'].values, height = 1, color = colors[2], edgecolor = 'black')
    #ax.set_xticks([0,50,100,150])
    ax.set_ylabel('', fontsize = fontsize)
    ax.set_xlabel('Number of PTMs', fontsize = fontsize)

    if title is not None:
        ax.set_title(title, fontsize = fontsize, weight = 'bold')
    elif title_type == 'detailed':
        ax.set_title(f'Top {top_terms} {database} {annot_type} Annotations', fontsize = fontsize, weight = 'bold')
    elif title_type == 'database':
        ax.set_title(f'{database}', fontsize = fontsize, weight = 'bold')

    #label_dict = {'EXONT:Name':'Exon Ontology Term', 'PSP:ON_PROCESS':'Biological Process (PSP)', 'PSP:ON_FUNCTION':'Molecular Function (PSP)', 'Combined:Kinase':'Kinase'}
    #ax.text(-1*ax.get_xlim()[1]/10, top_terms-0.2, label_dict[term_to_plot], weight = 'bold', ha = 'right', fontsize = 8)
    x_label_dict = {'Function':'Number of PTMs\nassociated with Function', 'Process':'Number of PTMs\nassociated with Process', 'Disease':'Number of PTMs\nassociated with Disease', 'Kinase':'Number of Phosphosites\ntargeted by Kinase', 'Interactions': 'Number of PTMs\nthat regulate interaction\n with protein','Motif Match':'Number of PTMs\nfound within a\nmotif instance', 'Intraprotein': 'Number of PTMs\nthat are important\for intraprotein\n interactions','Phosphatase':'Number of Phosphosites\ntargeted by Phosphatase', 'Perturbation (DIA2)': "Number of PTMs\nAffected by Perturbation\n(Measured by DIA)", 'Perturbation (PRM)': 'Number of PTMs\nAffected by Perturbation\n(Measured by PRM)', 'NetPath':'Number of PTMs/Genes\nassociated with NetPath', 'Perturbation':'Number of PTMs\nAffected by Perturbation'}
    ax.set_xlabel(x_label_dict[annot_type], fontsize = fontsize)
    
    #make a custom legend
    if legend:
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(facecolor = colors[0], edgecolor = 'black', label = 'Excluded'), mpatches.Patch(facecolor = colors[1], edgecolor = 'black', label = 'Included'),mpatches.Patch(facecolor = colors[2], edgecolor = 'black', label = 'Altered Flank')]
        ax.legend(handles = handles, ncol = 1, fontsize = fontsize - 2, title = 'Type of Impact', title_fontsize = fontsize - 1, loc = 'upper right', bbox_to_anchor = legend_loc)

    return ax


def gene_set_enrichment(spliced_ptms = None, altered_flanks = None, sig_col = 'Significance', dpsi_col = 'dPSI', alpha = 0.05, min_dPSI = None, gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023','Reactome_2022'], background = None, return_sig_only = True, max_retries = 5, delay = 10, **kwargs):
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
    if background is not None:
        raise ValueError('Background data not supported at this time, but will be added in future versions. Please set background = None, which will use all genes in the gene set database as the background.')
    
    #if combined is not None:
    #    if spliced_ptms is not None or altered_flanks is not None:
    #        print('If combined dataframe is provided, you do not need to include spliced_ptms or altered_flanks dataframes. Ignoring these inputs.')

    #    foreground = combined.copy()
    #    type = 'Differentially Included + Altered Flanking Sequences'

        #isolate the type of impact on the gene
    #    combined_on_gene = combined.groupby('Gene')['Impact'].apply(lambda x: ';'.join(set(x)))
    #    included = combined_on_gene.str.contains('Included')
    #    excluded = combined_on_gene.str.contains('Excluded')
    #    differential = included | excluded
    #    altered_flank = combined_on_gene.str.contains('Altered Flank')

    #    altered_flank_only = altered_flank & ~differential
    #    differential_only = differential & ~altered_flank
    #    both = differential & altered_flank
        
    #    altered_flank_only = combined_on_gene[altered_flank_only].index.tolist()
    #    differential_only = combined_on_gene[differential_only].index.tolist()
    #    both = combined_on_gene[both].index.tolist()
    if spliced_ptms is not None and altered_flanks is not None:
        #gene information (total and spliced genes)
        combined = summarize.combine_outputs(spliced_ptms, altered_flanks, **kwargs)
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
    if sig_col in foreground.columns and (min_dPSI is not None and dpsi_col in foreground.columns):
        foreground = foreground[foreground[sig_col] <= alpha].copy()
        foreground = foreground[foreground[dpsi_col].abs() >= min_dPSI]
    elif sig_col in foreground.columns:
        if min_dPSI is not None:
            print('`min_dpsi` provided but `dpsi_col` not found in dataframe. Ignoring `min_dpsi` parameter. Please indicate correct column name for dPSI values and rerun if desired.')
        foreground = foreground[foreground[sig_col] <= alpha].copy()
    elif min_dPSI is not None and 'dPSI' in combined.columns:
        print('`sig_col` not found in dataframe. Ignoring `alpha` parameter. Please indicate correct column name for Significance and rerun if desired.')
        foreground = combined[combined['dPSI'].abs() >= min_dPSI].copy()
    else:
        print('Significance column not found and min_dPSI not provided. All PTMs in dataframe will be considered as the foreground')

    foreground = foreground['Gene'].unique().tolist()   

    #construct background
    #if isinstance(background, list):
    #    pass
    #elif isinstance(background, np.ndarray):
    #    background = list(background)
    #elif isinstance(background, pd.DataFrame):
    #    if 'Gene' not in background.columns:
    #        raise ValueError('Gene column not found in dataframe. Please provide a valid column name for genes.')
    #    background = background['Gene'].unique().tolist()
    #elif background == 'Significance':
    #    if sig_col not in foreground.columns:
    #        raise ValueError('Significance column not found in dataframe. Please provide a valid column name for significance.')
    #    background = combined.copy()
    #    background = background['Gene'].unique().tolist()  

    

    
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
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(color = colors[2], label = 'Contains Both Events'), mpatches.Patch(color = colors[1], label = 'PTMs with Altered Flanking Sequence'), mpatches.Patch(color = colors[0], label = 'Differentially Included PTMs')]
        ax.legend(handles = handles, loc = 'upper center', borderaxespad = 0, bbox_to_anchor = (0.5, 1 + (1/figure_length)), ncol = 1, fontsize = 9)



    ax.set_xlabel('EnrichR Combined\nScore', fontsize = 11)
