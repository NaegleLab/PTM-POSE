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

package_dir = os.path.dirname(os.path.abspath(__file__))

def get_modification_counts(ptms):
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
    ptms['Modification Class'] = ptms['Modification Class'].apply(lambda x: x.split(';'))
    ptms = ptms.explode('Modification Class')
    modification_counts = ptms.groupby('Modification Class').size()
    modification_counts = modification_counts.sort_values(ascending = True)
    return modification_counts

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
    if len(annotation_columns) != annotation_columns_in_spliced_ptms:
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

    shared_columns = ['Impact', 'Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class'] + sig_cols + annotation_columns
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

def get_ptm_annotations(spliced_ptms, annotation_type = 'Function', database = 'PhosphoSitePlus', mod_class = None, collapse_on_similar = False, dPSI_col = None, sig_col = None):
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

    

    
def get_enrichment_inputs(spliced_ptms,  annotation_type = 'Function', database = 'PhosphoSitePlus', background_type = 'pregenerated', background = None, collapse_on_similar = False, mod_class = None, alpha = 0.05, min_dPSI = 0, annotation_file = None, save_background = False):
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
                background = get_modification_class_data(background, mod_class)
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


def annotation_enrichment(spliced_ptms, database = 'PhosphoSitePlus', annotation_type = 'Function', background_type = 'pregenerated', collapse_on_similar = False, mod_class = None, alpha = None, min_dPSI = None, annotation_file = None, save_background = False):#
    """
    In progress, needs to be tested

    Given spliced ptm information (differential inclusion, altered flanking sequences, or both), calculate the enrichment of specific annotations in the dataset using a hypergeometric test. Background data can be provided/constructed in a few ways:

    1. Use preconstructed background data for the annotation of interest, which considers the entire proteome present in the ptm_coordinates dataframe. While this is the default, it may not be the most accurate representation of your data, so you may alternative decide to use the other options which will be more specific to your context.
    2. Use the alpha and min_dPSI parameter to construct a foreground that only includes significantly spliced PTMs, and use the entire provided spliced_ptms dataframe as the background. This will allow you to compare the enrichment of specific annotations in the significantly spliced PTMs compared to the entire dataset. Will do this automatically if alpha or min_dPSI is provided.

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
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
    foreground_annotation_count, foreground_size, background_annotations, background_size, annotation_details = get_enrichment_inputs(spliced_ptms, background_type = background_type, annotation_type = annotation_type, database = database, collapse_on_similar = collapse_on_similar, mod_class = mod_class, alpha = alpha, min_dPSI = min_dPSI, annotation_file = annotation_file, save_background = save_background)
    

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


def gene_set_enrichment(spliced_ptms = None, altered_flanks = None, combined = None, sig_col = 'Significance', dpsi_col = 'dPSI', alpha = 0.05, min_dPSI = None, gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023','Reactome_2022'], background = None, return_sig_only = True, max_retries = 5, delay = 10):
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

    Returns
    -------
    results: pd.DataFrame
        Dataframe with gene set enrichment results from enrichr API

    """
    if background is not None:
        raise ValueError('Background data not supported at this time, but will be added in future versions. Please set background = None, which will use all genes in the gene set database as the background.')
    
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
            enr = gp.enrichr(foreground, gene_sets = gene_sets, organism='human')
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
    
def compare_flanking_sequences(altered_flanks, flank_size = 5):
    sequence_identity_list = []
    altered_positions_list = []
    residue_change_list = []
    flank_side_list = []
    for i, row in altered_flanks.iterrows():
        #if there is sequence info for both and does not introduce stop codons, compare sequence identity
        if not row['Stop Codon Introduced'] and row['Inclusion Flanking Sequence'] == row['Inclusion Flanking Sequence'] and row['Exclusion Flanking Sequence'] == row['Exclusion Flanking Sequence']:
            #compare sequence identity
            sequence_identity = getSequenceIdentity(row['Inclusion Flanking Sequence'], row['Exclusion Flanking Sequence'])
            #identify where flanking sequence changes
            altered_positions, residue_change, flank_side = findAlteredPositions(row['Inclusion Flanking Sequence'], row['Exclusion Flanking Sequence'], flank_size = flank_size)
        else:
            sequence_identity = np.nan
            altered_positions = np.nan
            residue_change = np.nan
            flank_side = np.nan



        #add to lists
        sequence_identity_list.append(sequence_identity)
        altered_positions_list.append(altered_positions)
        residue_change_list.append(residue_change)
        flank_side_list.append(flank_side)

    altered_flanks['Sequence Identity'] = sequence_identity_list
    altered_flanks['Altered Positions'] = altered_positions_list
    altered_flanks['Residue Change'] = residue_change_list
    altered_flanks['Altered Flank Side'] = flank_side_list
    return altered_flanks



def compare_inclusion_motifs(flanking_sequences, elm_classes = None):
    """
    Given a DataFrame containing flanking sequences with changes and a DataFrame containing ELM class information, identify motifs that are found in the inclusion and exclusion events, identifying motifs unique to each case. This does not take into account the position of the motif in the sequence or additional information that might validate any potential interaction (i.e. structural information that would indicate whether the motif is accessible or not). ELM class information can be downloaded from the download page of elm (http://elm.eu.org/elms/elms_index.tsv).

    Parameters
    ----------
    flanking_sequences: pandas.DataFrame
        DataFrame containing flanking sequences with changes, obtained from get_flanking_changes_from_splice_data()
    elm_classes: pandas.DataFrame
        DataFrame containing ELM class information (ELMIdentifier, Regex, etc.), downloaded directly from ELM (http://elm.eu.org/elms/elms_index.tsv). Recommended to download this file and input it manually, but will download from ELM otherwise

    Returns
    -------
    flanking_sequences: pandas.DataFrame
        DataFrame containing flanking sequences with changes and motifs found in the inclusion and exclusion events

    """
    if elm_classes is None:
        elm_classes = pd.read_csv('http://elm.eu.org/elms/elms_index.tsv', sep = '\t', header = 5)

        

    only_in_inclusion = []
    only_in_exclusion = []

    for _, row in flanking_sequences.iterrows():
        #check if there is a stop codon introduced and both flanking sequences are present
        if not row['Stop Codon Introduced'] and row['Inclusion Flanking Sequence'] == row['Inclusion Flanking Sequence'] and row['Exclusion Flanking Sequence'] == row['Exclusion Flanking Sequence']:
            #get elm motifs that match inclusion or Exclusion Flanking Sequences
            inclusion_matches = find_motifs(row['Inclusion Flanking Sequence'], elm_classes)
            exclusion_matches = find_motifs(row['Exclusion Flanking Sequence'], elm_classes)

            #get motifs that are unique to each case
            only_in_inclusion.append(';'.join(set(inclusion_matches) - set(exclusion_matches)))
            only_in_exclusion.append(';'.join(set(exclusion_matches) - set(inclusion_matches)))
        else:
            only_in_inclusion.append(np.nan)
            only_in_exclusion.append(np.nan)

    #save data
    flanking_sequences["Motif only in Inclusion"] = only_in_inclusion
    flanking_sequences["Motif only in Exclusion"] = only_in_exclusion
    return flanking_sequences

def identify_change_to_specific_motif(altered_flanks, elm_motif_name, elm_classes = None, modification_class = None, residues = None, dPSI_col = None):
    if 'Altered Positions' not in altered_flanks.columns:
        altered_flanks = compare_flanking_sequences(altered_flanks)
    
    #grab elm motifs that match inclusion or Exclusion Flanking Sequences
    if 'Motif only in Inclusion' not in altered_flanks.columns:
        altered_flanks = compare_inclusion_motifs(altered_flanks, elm_classes = elm_classes)

    #grab only needed info
    motif_data = altered_flanks.dropna(subset = ['Inclusion Flanking Sequence', 'Exclusion Flanking Sequence'], how = 'all').copy()
    cols_to_keep = ['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class', 'Inclusion Flanking Sequence', 'Exclusion Flanking Sequence', 'Motif only in Inclusion', 'Motif only in Exclusion', 'Altered Positions', 'Residue Change']
    if dPSI_col is not None:
        cols_to_keep.append(dPSI_col)

    #go through motif data and identify motifs matching elm motif of interest
    motif_data = motif_data[cols_to_keep]
    for i, row in motif_data.iterrows():
        if row['Motif only in Inclusion'] == row['Motif only in Inclusion']:
            if elm_motif_name in row['Motif only in Inclusion']:
                motif_data.loc[i, 'Motif only in Inclusion'] = ';'.join([motif for motif in row['Motif only in Inclusion'].split(';') if elm_motif_name in motif])
            else:
                motif_data.loc[i, 'Motif only in Inclusion'] = np.nan

        if row['Motif only in Exclusion'] == row['Motif only in Exclusion']:
            if elm_motif_name in row['Motif only in Exclusion']:
                motif_data.loc[i, 'Motif only in Exclusion'] = ';'.join([motif for motif in row['Motif only in Exclusion'].split(';') if elm_motif_name in motif])
            else:
                motif_data.loc[i, 'Motif only in Exclusion'] = np.nan

    #restrict to events that are specific modification types or residues (for example, SH2 domain motifs should be phosphotyrosine)
    motif_data = motif_data.dropna(subset = ['Motif only in Inclusion', 'Motif only in Exclusion'], how = 'all')
    if modification_class is not None:
        motif_data = motif_data[motif_data['Modification Class'].str.contains(modification_class)]

    if residues is not None and isinstance(residues, str):
        motif_data = motif_data[motif_data['Residue'] == residues]
    elif residues is not None and isinstance(residues, list):
        motif_data = motif_data[motif_data['Residue'].isin(residues)]
    elif residues is not None:
        raise ValueError('residues parameter must be a string or list of strings')
    
    return motif_data

    



def findAlteredPositions(seq1, seq2, flank_size = 5):
    """
    Given two sequences, identify the location of positions that have changed

    Parameters
    ----------
    seq1, seq2: str
        sequences to compare (order does not matter)
    flank_size: int
        size of the flanking sequences (default is 5). This is used to make sure the provided sequences are the correct length
    
    Returns
    -------
    altered_positions: list
        list of positions that have changed
    residue_change: list
        list of residues that have changed associated with that position
    flank_side: str
        indicates which side of the flanking sequence the change has occurred (N-term, C-term, or Both)
    """
    desired_seq_size = flank_size*2+1
    altered_positions = []
    residue_change = []
    flank_side = []
    seq_size = len(seq1)
    flank_size = (seq_size -1)/2
    if seq_size == len(seq2) and seq_size == desired_seq_size:
        for i in range(seq_size):
            if seq1[i] != seq2[i]:
                altered_positions.append(i-(flank_size))
                residue_change.append(f'{seq1[i]}->{seq2[i]}')
        #check to see which side flanking sequence
        altered_positions = np.array(altered_positions)
        n_term = any(altered_positions < 0)
        c_term = any(altered_positions > 0)
        if n_term and c_term:
            flank_side = 'Both'
        elif n_term:
            flank_side = 'N-term only'
        elif c_term:
            flank_side = 'C-term only'
        else:
            flank_side = 'Unclear'
        return altered_positions, residue_change, flank_side
    else:
        return np.nan, np.nan, np.nan
        
def getSequenceIdentity(seq1, seq2):
    """
    Given two flanking sequences, calculate the sequence identity between them using Biopython and parameters definded by Pillman et al. BMC Bioinformatics 2011

    Parameters
    ----------
    seq1, seq2: str
        flanking sequence 

    Returns
    -------
    normalized_score: float
        normalized score of sequence similarity between flanking sequences (calculated similarity/max possible similarity)
    """
    #make pairwise aligner object
    aligner = PairwiseAligner()
    #set parameters, with match score of 10 and mismatch score of -2
    aligner.mode = 'global'
    aligner.match_score = 10
    aligner.mismatch_score = -2
    #calculate sequence alignment score between two sequences
    actual_similarity = aligner.align(seq1, seq2)[0].score
    #calculate sequence alignment score between the same sequence
    control_similarity = aligner.align(seq1, seq1)[0].score
    #normalize score
    normalized_score = actual_similarity/control_similarity
    return normalized_score

def find_motifs(seq, elm_classes):
    """
    Given a sequence and a dataframe containinn ELM class information, identify motifs that can be found in the provided sequence using the RegEx expression provided by ELM (PTMs not considered). This does not take into account the position of the motif in the sequence or additional information that might validate any potential interaction (i.e. structural information that would indicate whether the motif is accessible or not). ELM class information can be downloaded from the download page of elm (http://elm.eu.org/elms/elms_index.tsv).

    Parameters
    ----------
    seq: str
        sequence to search for motifs
    elm_classes: pandas.DataFrame
        DataFrame containing ELM class information (ELMIdentifier, Regex, etc.), downloaded directly from ELM (http://elm.eu.org/elms/elms_index.tsv)
    """
    matches = []
    for j, elm_row in elm_classes.iterrows():
        reg_ex = elm_row['Regex']
        if re.search(reg_ex, seq) is not None:
            matches.append(elm_row['ELMIdentifier'])

    return matches
    

class protein_interactions:
    """
    Class to assess interactions facilitated by PTMs in splicing network

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    include_enzyme_interactions: bool
        Whether to include interactions with enzymes in the network, such as kinase-substrate interactions. Default is True
    interaction_databases: list
        List of databases whose information to include in the network. Default is ['PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD', 'ELM']

    Attributes
    ----------
    interaction_graph: nx.Graph
        NetworkX graph object representing the interaction network, created from analyze.get_interaction_network
    network_data: pd.DataFrame
        Dataframe containing details about specifici protein interactions (including which protein contains the spliced PTMs)
    network_stats: pd.DataFrame
        Dataframe containing network statistics for each protein in the interaction network, obtained from analyze.get_interaction_stats(). Default is None, which will not label any proteins in the network.
    """
    def __init__(self, spliced_ptms, include_enzyme_interactions = True, interaction_databases = ['PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD', 'ELM']):
        self.spliced_ptms = spliced_ptms
        self.include_enzyme_interactions = include_enzyme_interactions
        self.interaction_databases = interaction_databases


    def get_interaction_network(self, node_type = 'Gene'):
        """
        Given the spliced PTM data, extract interaction information and construct a dataframe containing all possible interactions driven by PTMs, either centered around specific PTMs or specific genes.

        Parameters
        ----------
        node_type: str
            What to define interactions by. Can either be by 'PTM', which will consider each PTM as a separate node, or by 'Gene', which will aggregate information across all PTMs of a single gene into a single node. Default is 'Gene'
        """
        if node_type not in ['Gene', 'PTM']:
            raise ValueError("node_type parameter (which dictates whether to consider interactions at PTM or gene level) can be either Gene or PTM")
        
        #extract interaction information in provided data
        interactions = annotate.combine_interaction_data(self.spliced_ptms, include_enzyme_interactions=self.include_enzyme_interactions, interaction_databases=self.interaction_databases)
        if interactions.empty:
            raise ValueError('No interaction data found in spliced PTM data. Please provide spliced PTM data with interaction information to generate network. This can be done by using the annotate module')
        interactions['Residue'] = interactions['Residue'] + interactions['PTM Position in Isoform'].astype(int).astype(str)
        interactions = interactions.drop(columns = ['PTM Position in Isoform'])

        #add regulation change information
        if 'dPSI' in self.spliced_ptms.columns:
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

        self.network_data = network_data
        self.interaction_graph = interaction_graph


    def get_interaction_stats(self):
        """
        Given the networkx interaction graph, calculate various network centrality measures to identify the most relevant PTMs or genes in the network
        """
        #calculate network centrality measures
        degree_centrality = nx.degree_centrality(self.interaction_graph)
        closeness_centrality = nx.closeness_centrality(self.interaction_graph)
        betweenness_centrality = nx.betweenness_centrality(self.interaction_graph)
        network_stats = pd.DataFrame({'Degree': dict(self.interaction_graph.degree()), 'Degree Centrality':degree_centrality, 'Closeness':closeness_centrality,'Betweenness':betweenness_centrality})
        self.network_stats = network_stats

    def get_protein_interaction_network(self, protein):
        """
        Given a specific protein, return the network data for that protein

        Parameters
        ----------
        protein: str
            Gene name of the protein of interest

        Returns
        -------
        protein_network: pd.DataFrame
            Dataframe containing network data for the protein of interest
        """
        if not hasattr(self, 'network_data'):
            self.get_interaction_network()

        if protein not in self.network_data['Modified Gene'].unique():
            print(f'{protein} is not found in the network data. Please provide a valid gene name.')
            return None
        
        protein_network = self.network_data[self.network_data['Modified Gene'] == protein]
        protein_network = protein_network.drop(columns = ['Modified Gene'])
        protein_network = protein_network.rename(columns = {'Residue': 'Spliced PTMs facilitating Interacting'})
        return protein_network

    def summarize_protein_network(self, protein):
        """
        Given a protein of interest, summarize the network data for that protein

        Parameters
        ----------
        protein: str
            Gene name of the protein of interest
        """
        if not hasattr(self, 'network_data'):
            self.get_interaction_network()

        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        protein_network = self.network_data[self.network_data['Modified Gene'] == protein]
        increased_interactions = protein_network.loc[protein_network['Regulation Change'] == '+', 'Interacting Gene'].values
        decreased_interactions = protein_network.loc[protein_network['Regulation Change'] == '-', 'Interacting Gene'].values
        ambiguous_interactions = protein_network.loc[protein_network['Regulation Change'].str.contains(';'), 'Interacting Gene'].values

        #print interactions
        if len(increased_interactions) > 0:
            print(f"Increased interaction likelihoods: {', '.join(increased_interactions)}")
        if len(decreased_interactions) > 0:
            print(f"Decreased interaction likelihoods: {', '.join(decreased_interactions)}")
        if len(ambiguous_interactions) > 0:
            print(f"Ambiguous interaction impact: {', '.join(ambiguous_interactions)}")

        network_ranks = self.network_stats.rank(ascending = False).astype(int)
        print(f'Number of interactions: {self.network_stats.loc[protein, "Degree"]} (Rank: {network_ranks.loc[protein, "Degree"]})')
        print(f'Centrality measures - \t Degree = {self.network_stats.loc[protein, "Degree Centrality"]} (Rank: {network_ranks.loc[protein, "Degree Centrality"]})')
        print(f'                      \t Betweenness = {self.network_stats.loc[protein, "Betweenness"]} (Rank: {network_ranks.loc[protein, "Betweenness"]})')
        print(f'                      \t Closeness = {self.network_stats.loc[protein, "Closeness"]} (Rank: {network_ranks.loc[protein, "Closeness"]})')

    def plot_interaction_network(self, modified_color = 'red', modified_node_size = 10, interacting_color = 'lightblue', interacting_node_size = 5, defaultedgecolor = 'gray', color_edges_by = 'Same', seed = 200, ax = None, proteins_to_label = None, labelcolor = 'black'):
        """
        Given the interactiong graph and network data outputted from analyze.get_interaction_network, plot the interaction network, signifying which proteins or ptms are altered by splicing and the specific regulation change that occurs. by default, will only label proteins 

        Parameters
        ----------
        modified_color: str
            Color to use for nodes that are modified by splicing. Default is 'red'
        modified_node_size: int
            Size of nodes that are modified by splicing. Default is 10
        interacting_color: str
            Color to use for nodes that interact with modified nodes. Default is 'lightblue'
        interacting_node_size: int
            Size of nodes that interact with modified nodes. Default is 5
        edgecolor: str
            Color to use for edges in the network. Default is 'gray'. Can choose to color by database by providing a dictionary with database names as keys and colors as values or by specifying 'database' as color
        seed: int
            Seed to use for random number generator. Default is 200
        ax: matplotlib.pyplot.Axes
            Axes object to plot the network. Default is None
        """
        if not hasattr(self, 'interaction_graph'):
            self.get_interaction_network()

        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        pose_plots.plot_interaction_network(self.interaction_graph, self.network_data, self.network_stats, modified_color = modified_color, modified_node_size = modified_node_size, interacting_color = interacting_color, interacting_node_size = interacting_node_size, defaultedgecolor = defaultedgecolor, color_edges_by=color_edges_by, seed = seed, ax = ax, proteins_to_label = proteins_to_label, labelcolor = labelcolor)

    def plot_network_centrality(self,  centrality_measure = 'Degree', top_N = 10, modified_color = 'red', interacting_color = 'black', ax = None):
        if not hasattr(self, 'interaction_graph'):
            self.get_interaction_network()
        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        pose_plots.plot_network_centrality(self.network_stats, self.network_data, centrality_measure=centrality_measure,top_N = top_N, modified_color = modified_color, interacting_color = interacting_color, ax = ax)

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


class KL_flank_analysis:
    def __init__(self, altered_flanks, odir):
        self.altered_flanks = altered_flanks
        self.odir = odir

    def identify_sequences_of_interest(self):
        self.sequences_of_interest = self.altered_flanks[(~self.altered_flanks['Matched']) & (~self.altered_flanks['Stop Codon Introduced']) & (self.altered_flanks['Modification Class'].str.contains('Phosphorylation'))].copy() 

    
    def process_data_for_kinase_library(self):
        """
        Extract flanking sequence information for 
        """
        #restrict to events with changed flanking sequences, no introduced stop codons, and phosphorylation modifications
        if not hasattr(self, 'sequences_of_interest'):
            self.identify_sequences_of_interest()

        #generate files to input into Kinase Library (inclusion first then exclusion)
        inclusion_sequences = self.sequences_of_interest[['PTM', 'Inclusion Flanking Sequence']].drop_duplicates()
        inclusion_sequences['Inclusion Flanking Sequence'] = inclusion_sequences['Inclusion Flanking Sequence'].apply(edit_sequence_for_kinase_library)
        inclusion_sequences = inclusion_sequences.dropna(subset = 'Inclusion Flanking Sequence')
        #write sequences to text file
        with open(self.odir + 'inclusion_sequences_input.txt', 'w') as f:
            for _, row in inclusion_sequences.iterrows():
                f.write(row['Inclusion Flanking Sequence']+'\n')

        exclusion_sequences = self.sequences_of_interest[['PTM', 'Exclusion Flanking Sequence']].drop_duplicates()
        exclusion_sequences['Exclusion Flanking Sequence'] = exclusion_sequences['Exclusion Flanking Sequence'].apply(edit_sequence_for_kinase_library)
        exclusion_sequences = exclusion_sequences.dropna(subset = 'Exclusion Flanking Sequence')
        #write sequences to text file
        with open(self.odir + 'exclusion_sequences_input.txt', 'w') as f:
            for _, row in exclusion_sequences.iterrows():
                f.write(row['Exclusion Flanking Sequence']+'\n')

        print('Input files for Kinase Library generated. Please run upload the file to the "score sites" tab of Kinase Library (https://kinase-library.mit.edu/sites) and download the full results.')

    def format_sequences_to_match_output(self, sequence_type = 'Inclusion'):
        if not hasattr(self, 'sequences_of_interest'):
            self.identify_sequences_of_interest()

        sequences = self.sequences_of_interest[['Region ID','PTM', f'{sequence_type} Flanking Sequence']].drop_duplicates().copy()
        sequences = sequences.dropna(subset = 'Inclusion Flanking Sequence')
        sequences['Label'] = sequences['Region ID'] + ';' + sequences['PTM']
        sequences[f'{sequence_type} Flanking Sequence'] = sequences[f'{sequence_type} Flanking Sequence'].apply(lambda x: x.upper().replace(' ', '_')+'_')
        return sequences

    def process_kinase_library_output(self, scores, sequence_type = 'Inclusion'):
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
        if not hasattr(self, 'sequences_of_interest'):
            self.identify_sequences_of_interest()

        sequences = self.format_sequences_to_match_output(sequence_type = sequence_type)


        sequences = sequences.merge(scores, left_on = f'{sequence_type} Flanking Sequence', right_on = 'sequence', how = 'left')
        #split info into tyrosine vs. serine/threonine
        sequences_y = sequences[sequences['Label'].str.contains('_Y')]
        sequences_st = sequences[(sequences['Label'].str.contains('_S')) | (sequences['Label'].str.contains('_T'))]

        #pivot table to get scores for each kinase
        percentiles_y = sequences_y.pivot_table(index = 'Label', columns = 'kinase', values = 'site_percentile')
        percentiles_st = sequences_st.pivot_table(index = 'Label', columns = 'kinase', values = 'site_percentile')

        return percentiles_y, percentiles_st

    def get_kinase_library_differences(self, inclusion_scores_file, exclusion_scores_file):
        """
        Given altered flanking sequences and kinase library scores for inclusion and Exclusion Flanking Sequences, calculate the difference in kinase library site percentiles between the two

        Parameters
        ----------
        altered_flanks: pd.DataFrame
            Dataframe with PTMs associated with altered flanking sequences
        inclusion_scores: pd.DataFrame
            Dataframe with kinase library scores for Inclusion Flanking Sequences (loaded from downloaded .tsv outputs from kinase library)
        exclusion_scores: pd.DataFrame
            Dataframe with kinase library scores for Exclusion Flanking Sequences (loaded from downloaded .tsv outputs from kinase library)
        
        Returns
        -------
        percentiles_diff_y: pd.DataFrame
            Dataframe with the difference in kinase library scores for tyrosine sites
        percentiles_diff_st: pd.DataFrame
            Dataframe with the difference in kinase library scores for serine/threonine sites
        """
        inclusion_scores = pd.read_csv(inclusion_scores_file, sep = '\t')
        inclusion_percentiles_y, inclusion_percentiles_st = self.process_kinase_library_output(inclusion_scores, sequence_type = 'Inclusion')
        exclusion_scores = pd.read_csv(exclusion_scores_file, sep = '\t')
        exclusion_percentiles_y, exclusion_percentiles_st = self.process_kinase_library_output(exclusion_scores, sequence_type = 'Exclusion')

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

        #save all data
        self.inclusion_percentiles = {}
        self.inclusion_percentiles['Y'] = inclusion_percentiles_y
        self.inclusion_percentiles['ST'] = inclusion_percentiles_st

        self.exclusion_percentiles = {}
        self.exclusion_percentiles['Y'] = exclusion_percentiles_y
        self.exclusion_percentiles['ST'] = exclusion_percentiles_st

        self.percentile_difference = {}
        self.percentile_difference['Y'] = percentiles_diff_y
        self.percentile_difference['ST'] = percentiles_diff_st

    
#def process_data_for_exon_ontology(odir, spliced_ptms = None, altered_flanks = None):
#    pass



    

class kstar_enrichment:
    def __init__(self, significant_ptms, network_dir, background_ptms = None, phospho_type = 'Y'):
        """
        Given spliced ptm or PTMs with altered flanks and a single kstar network, get enrichment for each kinase in the network using a hypergeometric. Assumes the  data has already been reduced to the modification of interest (phosphotyrosine or phoshoserine/threonine)

        Parameters
        ----------
        network_dir : dict
            dictionary of networks with kinase-substrate information
        spliced_ptms : pandas dataframe
            all PTMs of interest
        background_ptms: pd.DataFrame
            PTMs to consider as the background for enrichment purposes, which should overlap with the spliced ptms information provided (an example might be all identified events, whether or not they are significant). If not provided, will use all ptms in the phosphoproteome.
        phospho_type : str 
            type of phosphorylation event to extract. Can either by phosphotyrosine ('Y') or phosphoserine/threonine ('ST'). Default is 'Y'.

        """
        #change phospho_type to list if not already
        if isinstance(phospho_type, list):
            pass
        elif phospho_type in ['Y', 'ST']:
            phospho_type = [phospho_type]
        elif phospho_type == 'All':
            phospho_type = ['Y', 'ST']
        else:
            raise ValueError('phospho_type must be either Y, ST, or All')
        
        self.significant_ptms = {}
        self.background_ptms = {}
        self.networks = {}
        for ptype in phospho_type:
            #process ptms to only include specific phosphorylation data needed
            self.significant_ptms[ptype] = self.process_ptms(significant_ptms, phospho_type = ptype)
            #load background, either from provided data or from all phosphoproteome data
            if background_ptms is not None:
                self.background_ptms[ptype] = self.process_ptms(background_ptms, phospho_type=ptype)
            else:
                self.background_ptms[ptype] = self.process_ptms(pose_config.ptm_coordinates.copy(), phospho_type = ptype)

            #load kstar networks
            self.networks[ptype] = self.load_networks(network_dir, phospho_type = ptype)


        #save info
        self.phospho_type = phospho_type
        self.median_enrichment = None

    def load_networks(self, network_dir, phospho_type = 'Y'):
        """
        Given network directory, load all kstar networks for specific phosphorylation type
        """
        #check if file exists and whether a pickle has been generated: if not, load each network file individually
        if not os.path.exists(network_dir):
            raise ValueError('Network directory not found')
        elif os.path.exists(f"{network_dir}/*.p"):
            networks = pickle.load(open(f"{network_dir}/network_{phospho_type}.p", "rb" ) )
        else:
            network_directory = network_dir + f'/{phospho_type}/INDIVIDUAL_NETWORKS/'
            networks = {}
            for file in os.listdir(network_directory):
                if file.endswith('.tsv'):
                    #get the value of the network number
                    file_noext = file.strip(".tsv").split('_')
                    key_name = 'nkin'+str(file_noext[1])
                    #print("Debug: key name is %s"%(key_name))
                    networks[key_name] = pd.read_csv(f"{network_directory}{file}", sep='\t')
        return networks


    def process_ptms(self, ptms, phospho_type = 'Y'):
        """
        Given ptm information, restrict data to include only the phosphorylation type of interest and add a PTM column for matching information from KSTAR

        Parameters
        ----------
        ptms: pd.DataFrame
            ptm information containing modification type and ptm locatin information, such as the output from projection or altered flanking sequence analysis
        phospho_type: str
            type of phosphorylation event to extract. Can either by phosphotyrosine ('Y') or phosphoserine/threonine ('ST')
        
        Returns
        ptms: pd.DataFrame
            trimmed dataframe containing only modifications of interest and new 'PTM' column
        """

        #restrict to ptms to phosphorylation type of interest
        if phospho_type == 'Y':
            ptms = ptms[ptms['Modification'].str.contains('Phosphotyrosine')].copy()
        elif phospho_type == 'ST':
            ptms = ptms[(ptms["Modification"].str.contains('Phosphoserine')) | (ptms['Modification'].str.contains('Phosphothreonine'))].copy()

        #construct PTM column that matches KSTAR information
        ptms['PTM'] = ptms['UniProtKB Accession'] + '_' + ptms['Residue'] + ptms['PTM Position in Isoform'].astype(int).astype(str)

        #filter out any PTMs that come from alternative isoforms
        ptms = ptms[~ptms['UniProtKB Accession'].str.contains('-')]
        return ptms

    
    def get_enrichment_single_network(self, network_key, phospho_type = 'Y'):
        """
        in progress
        """
        network = self.networks[phospho_type][network_key]
        network['PTM'] = network['KSTAR_ACCESSION'] + '_' + network['KSTAR_SITE']

        #add network information to all significant data
        sig_ptms = self.significant_ptms[phospho_type][['PTM']].drop_duplicates()
        sig_ptms_kstar = sig_ptms.merge(network[['KSTAR_KINASE','PTM']], on = 'PTM')

        #repeat for background data
        bg_ptms = self.background_ptms[phospho_type][['PTM']].drop_duplicates()
        bg_ptms_kstar = bg_ptms.merge(network[['KSTAR_KINASE','PTM']], on = 'PTM')

        results = pd.DataFrame(np.nan, index = sig_ptms_kstar['KSTAR_KINASE'].unique(), columns = ['k','n','M','N','p'])
        for kinase in sig_ptms_kstar['KSTAR_KINASE'].unique():
            #get numbers for a hypergeometric test to look for enrichment of kinase substrates
            k = sig_ptms_kstar.loc[sig_ptms_kstar['KSTAR_KINASE'] == kinase, 'PTM'].nunique()
            n = bg_ptms_kstar.loc[bg_ptms_kstar['KSTAR_KINASE'] == kinase, 'PTM'].nunique()
            M = bg_ptms['PTM'].nunique()
            N = sig_ptms_kstar['PTM'].nunique()

            #run hypergeometric test
            results.loc[kinase,'p'] = stat_utils.hypergeom(M,n,N,k)
            results.loc[kinase, 'M'] = M
            results.loc[kinase, 'N'] = N
            results.loc[kinase, 'k'] = k
            results.loc[kinase, 'n'] = n

        return results
    
    def get_enrichment_all_networks(self, phospho_type = 'Y'):
        """
        Given prostate data and a dictionary of kstar networks, get enrichment for each kinase in each network in the prostate data. Assumes the prostate data has already been reduced to the modification of interest (phosphotyrosine or phoshoserine/threonine)

        Parameters
        ----------
        networks : dict
            dictionary of kstar networks
        prostate : pandas dataframe
            all PTMs identified in tCGA prostate data, regardless of significance (reduced to only include mods of interest)
        sig_prostate : pandas dataframe
            significant PTMs identified in tCGA prostate data, p < 0.05 and effect size > 0.25 (reduced to only include mods of interest)
        """
        results = {}
        for network in self.networks[phospho_type]:
            results[network] = self.get_enrichment_single_network(network_key=network, phospho_type=phospho_type)
        return results

    def extract_enrichment(self, results):
        """
        Given a dictionary of results from get_enrichment_all_networks, extract the p-values for each network and kinase, and then calculate the median p-value across all networks for each kinase

        Parameters
        ----------
        results : dict
            dictionary of results from get_enrichment_all_networks
        """
        enrichment = pd.DataFrame(index = results['nkin0'].index, columns = results.keys())
        for network in results:
            enrichment[network] = results[network]['p']
        enrichment['median'] = enrichment.median(axis = 1)
        return enrichment
    
    def run_kstar_enrichment(self):
        """
        Run full kstar analysis to generate substrate enrichment across each of the 50 KSTAR networks and calculate the median p-value for each kinase across all networks
        """
        enrichment = {}
        median = {}
        for ptype in self.phospho_type:
            results = self.get_enrichment_all_networks(phospho_type=ptype)
            enrichment[ptype] = self.extract_enrichment(results)
            median[ptype] = enrichment[ptype]['median']
        self.enrichment_all = enrichment
        self.median_enrichment = median

    def return_enriched_kinases(self, alpha = 0.05):
        """
        Return kinases with a median p-value less than the provided alpha value (substrates are enriched among the significant PTMs)

        Parameters
        ----------
        alpha : float
            significance threshold to use to subset kinases. Default is 0.05.
        """
        if self.median_enrichment is None:
            self.run_kstar_enrichment()

        sig_enrichment = {}
        for ptype in self.phospho_type:
            sig_enrichment[ptype] = self.median_enrichment[ptype][self.median_enrichment[ptype] < alpha].index.values
        return sig_enrichment
    


