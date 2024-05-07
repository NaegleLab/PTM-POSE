import numpy as np
import pandas as pd

#stats packages
import scipy.stats as stats

import seaborn as sns

#plotting 
import matplotlib.pyplot as plt

#custom functions
import annotate, stat_utils



def show_available_annotations(spliced_ptms, show_all_ptm_count = True, figsize = (5, 5)):
    if show_all_ptm_count:
        num_ptms = [spliced_ptms.drop_duplicates(['UniProtKB Accession', 'Residue']).shape[0]]
        num_ptms_filters = ['All PTMs']
        filter_source = []
    else:
        num_ptms = []
        num_ptms_filters = []
        filter_source = []

    #look for annotations and add counts to lists
    if 'PSP:DOMAIN' in spliced_ptms.columns:
        #in domain or not
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:DOMAIN').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs In Domain')
        filter_source.append('PhosphoSitePlus')
       
        #with a biological process annotation
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:ON_PROCESS').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs associated with Biological Process')
        filter_source.append('PhosphoSitePlus')

        #with a molecular function annoation in PSP
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:ON_FUNCTION').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs associated with Molecular Function')
        filter_source.append('PhosphoSitePlus')
    if 'PSP:Kinase' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:Kinase').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs associated with Kinase')
        filter_source.append('PhosphoSitePlus')
    if 'PSP:Disease_Association' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:Disease_Association').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs with Disease Association')
        filter_source.append('PhosphoSitePlus')
    if 'ELM Interactions' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'ELM Interactions').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs with Known ELM Interactions')
        filter_source.append('ELM')
    if 'ELM Motif Matches' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms[spliced_ptms['ELM Motif Matches'] != ''].drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs with ELM Motif Matches')
        filter_source.append('ELM')

        #with ligand motifs
        num_ptms.append(spliced_ptms[spliced_ptms['ELM Motif Matches'].str.contains('LIG')].drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs with ELM Ligand Motif Matches')
        filter_source.append('ELM')

        #with phospho dependent motifs
        sh2_motifs = ['LIG_SH2_CRK', 'LIG_SH2_GRB2like', 'LIG_SH2_NCK_1', 'LIG_SH2_PTP2', 'LIG_SH2_SFK_2', 'LIG_SH2_SFK_CTail_3', 'LIG_SH2_STAP1', 'LIG_SH2_STAT3', 'LIG_SH2_STAT5', 'LIG_SH2_STAT6']
        ptb_motifs = ['LIG_PTB_Phospho_1']
        fha_motifs = ['LIG_FHA_1', 'LIG_FHA_2']
        other_motifs = ['LIG_TYR_ITAM', 'LIG_TYR_ITIM', 'LIG_TYR_ITSM', 'LIG_RBL1_LxSxE_2', 'LIG_BRCT_BRCA1_1', 'LIG_BRCT_BRCA1_2', 'LIG_BRCT_MDC1_1', 'LIG_DLG_GKlike_1']
        fourteen33_motifs = ['LIG_14-3-3_CanoR_1', 'LIG_14-3-3_CterR_2']
        phospho_dependent_motifs = sh2_motifs + ptb_motifs + fha_motifs + other_motifs + fourteen33_motifs

        num_ptms.append(spliced_ptms[(spliced_ptms['ELM Motif Matches'].str.contains('|'.join(phospho_dependent_motifs))) & (spliced_ptms['Modification Class'] == 'Phosphorylation')].shape[0])
        num_ptms_filters.append('PTMs with Phospho Dependent ELM Motif Matches')
        filter_source.append('ELM')
    if 'PTMInt:Interaction' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'PTMInt:Interaction').drop_duplicates(subset = ['UniProtKB Accession', 'Residue']).shape[0])
        num_ptms_filters.append('PTMs with Known PTMInt Interactions')
        filter_source.append('PTMInt')
    
    #plot bar plot
    #color bars based on datasource
    palette = {'None': 'gray', 'PhosphoSitePlus': 'blue', 'ELM': 'green', 'PTMInt':'red'}
    colors = []
    for source in filter_source:
        colors.append(palette[source])
    fig, ax = plt.subplots(figsize = figsize)
    ax.barh(num_ptms_filters[::-1], num_ptms[::-1], color = colors[::-1])
    ax.set_xlabel('Number of PTMs')
    ax.set_title('Number of PTMs in Dataset')

    #create legend
    handles = [plt.Rectangle((0,0),1,1, color = color) for color in palette.values() if color != 'gray']
    labels = [source for source in palette.keys() if source != 'None']
    ax.legend(handles, labels, title = 'Annotation Source')
    plt.show()

def density_plots(annotated_splice_data, mod_class = None, sig_col = None, alpha = 0.05, type = 'hist', perform_sig_tests = True, figsize = (4,3)):
    if sig_col is not None:
        significant_splice_data = annotated_splice_data[annotated_splice_data[sig_col] <= alpha].copy()

    #if requested perform significance tests
    if perform_sig_tests and sig_col is not None:
        #report direction of significance
        if np.mean(annotated_splice_data['PTM Density (PTMs/bp)']) < np.mean(significant_splice_data['PTM Density (PTMs/bp)']):
            print('Significant events have higher PTM density')
        else:
            print('Significant events have lower PTM density')

        f, p = stats.mannwhitneyu(annotated_splice_data['PTM Density (PTMs/bp)'].dropna().values, significant_splice_data['PTM Density (PTMs/bp)'].dropna().values)
        print('Mann Whitney U test: p =', p)
    elif sig_col is not None:
        print('Significance column not provided (sig_col = None), skipping significance testing')


    if type =='hist':
        fig, ax = plt.subplots(figsize=figsize)
        
        ax.hist(annotated_splice_data['PTM Density (PTMs/bp)'], density = True, bins = 20, alpha = 0.5, label = 'All Events')
        if sig_col is not None:
            ax.hist(significant_splice_data['PTM Density (PTMs/bp)'], density = True, bins = 20, alpha = 0.5, label = 'Significant Events')
            ax.legend()

        ax.set_xlabel('PTM Density (PTMs/bp)')
        ax.set_ylabel('Density')
    elif type == 'scatter':
        if sig_col is None:
            fig, ax = plt.subplots(1,1,figsize=figsize)
            ax.scatter(annotated_splice_data['Event Length'], annotated_splice_data['PTM Density (PTMs/bp)'], s = 1)
            ax.set_xlabel('PTM Density (PTM sites/bp)')
            ax.set_ylabel('Event Length (bp)')
        else:
            fig, ax = plt.subplots(1,2, figsize=figsize,sharey = True, sharex = True)
            #all ptms
            ax[0].scatter(annotated_splice_data['Event Length'], annotated_splice_data['PTM Density (PTMs/bp)'], s = 1)
            ax[0].set_ylabel('PTM Density (PTM sites/bp)')
            ax[0].set_xlabel('Exon Length (bp)')
            ax[0].set_title('All PTMs')

            #significant ptms
            ax[1].scatter(significant_splice_data['Event Length'], significant_splice_data['PTM Density (PTMs/bp)'], s = 1)
            ax[1].set_ylabel('PTM Density (PTMs/bp)')
            ax[1].set_xlabel('Exon Length (bp)')
            ax[1].set_title('Significant PTMs')
    elif type == 'boxplot':
        annotated_splice_data['Type of Event'] = 'All'
        if sig_col is None:
            plt_data = annotated_splice_data.copy()
        
        else:
            significant_splice_data['Type of Event'] = 'Significant'
            plt_data = pd.concat([annotated_splice_data, significant_splice_data]).copy()
        fig, ax = plt.subplots(figsize=figsize)
        sns.boxplot(data = plt_data, x = 'Type of Event', y = 'PTM Density (PTMs/bp)', ax = ax)
    else: 
        print("Plot type not recognized. Available plot types are 'hist', 'scatter', and 'boxplot'")
        ax = None


    return ax

def get_annotation_col(spliced_ptms, annotation_type = 'Function', database = 'PhosphoSitePlus'):
    """
    Given the database of interest and annotation type, return the annotation column that will be found in a annotated spliced_ptm dataframe
    """
    #check to make sure requested annotation is available
    if database == 'PhosphoSitePlus':
        col_dict = {'Function':'PSP:ON_FUNCTION', 'Process':'PSP:ON_PROCESS', 'Disease':"PSP:Disease_Associated", 'Kinase':'PSP:Kinase', 'Interactions': 'PSP:ON_INTERACT', 'Domain':'PSP:Domain'}
        if annotation_type in col_dict.keys():
            annotation_col = col_dict[annotation_type]
            if annotation_col not in spliced_ptms.columns:
                raise ValueError('Requested annotation data has not yet been added to spliced_ptms dataframe. Please run the appropriate annotate module function to append this information.')
        else:
            raise ValueError(f"Invalid annotation type for PhosphoSitePlus. Available annotation data for PhosphoSitePlus includes: {', '.join(col_dict.keys())}")
    elif database == 'ELM':
        col_dict = {'Interactions':'ELM:Interactions', 'Motif Match':'ELM Motif Matches'}
        if annotation_type in col_dict.keys():
            annotation_col = col_dict[annotation_type]
            if annotation_col not in spliced_ptms.columns:
                if annotation_type == 'Interactions':
                    raise ValueError('Requested annotation data has not yet been added to spliced_ptms dataframe. Please run the annotate.add_ELM_interactions() function to append this information.')
                elif annotation_type == 'Motif Match':
                    raise ValueError('Requested annotation data has not yet been added to spliced_ptms dataframe. Please run the annotate.add_ELM_matched_motifs() function to append this information.')
        else:
            raise ValueError(f"Invalid annotation type for ELM. Available annotation data for ELM includes: {', '.join(col_dict.keys())}")
    elif database == 'PTMInt':
        col_dict = {'Interactions':'PTMInt:Interactions'}
        if annotation_type in col_dict.keys():
            annotation_col = col_dict[annotation_type]
        else:
            raise ValueError(f"Invalid annotation type for PTMInt. Available annotation data for PhosphoSitePlus includes: {', '.join(col_dict.keys())}")
    else:
        raise ValueError("Invalid database. Available options include 'PhosphoSitePlus', 'ELM', and 'PTMInt'")

    return annotation_col

def get_modification_class_data(spliced_ptms, mod_class):
    #check if specific modification class was provided and subset data by modification if so
    if mod_class in spliced_ptms['Modification Class'].values:
        ptms_of_interest = spliced_ptms[spliced_ptms['Modification Class'] == mod_class].copy()
    else:
        raise ValueError(f"Requested modification class not present in the data. The available modifications include {', '.join(col_dict.keys())}")

    return ptms_of_interest

def get_ptm_annotations(spliced_ptms, annotation_type = 'Function', database = 'PhosphoSitePlus', mod_class = None):
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
    annotation_col = get_annotation_col(spliced_ptms, database = database, annotation_type = annotation_type)


    #check if specific modification class was provided and subset data by modification if so
    if mod_class is not None:
        ptms_of_interest = get_modification_class_data(spliced_ptms, mod_class)
    else:
        ptms_of_interest = spliced_ptms.copy()

    #extract relevant annotation and remove PTMs without an annotation
    annotations = ptms_of_interest[['UniProtKB Accession', 'Residue', 'Modification Class'] + [annotation_col]].copy()
    annotations = annotations.dropna(subset = annotation_col)

    #separate distinct modification annotations in unique rows
    annotations_exploded = annotations.copy()
    annotations_exploded[annotation_col] = annotations_exploded[annotation_col].apply(lambda x: x.split(';') if isinstance(x, str) else np.nan)
    annotations_exploded = annotations_exploded.explode(annotation_col)
    annotations_exploded[annotation_col] = annotations_exploded[annotation_col].apply(lambda x: x.strip() if isinstance(x, str) else np.nan)
    #get the number of annotations found
    annotation_counts = annotations_exploded[annotation_col].value_counts().reset_index()
    #annotation_counts[f'Number of Spliced PTMs Associated with {annotation_type}']
    return annotations, annotation_counts



def annotation_enrichment(spliced_ptms, sig_col = None, background_ptms = None, database = 'PhosphoSitePlus', annotation_type = 'Function', collapse_on_similar = False, mod_class = None, alpha = 0.05):
    """
    """
    #if background not provided, use background data
    if sig_col is None:
        background_ptms = spliced_ptms.copy()
        #restrict sample to significantly spliced ptms
        spliced_ptms = spliced_ptms[spliced_ptms[sig_col] <= alpha].copy()
    elif background_ptms is None:
        raise ValueError('General background not yet created')
    

    #check if specific modification class was provided and subset data by modification if so
    if mod_class is not None:
        spliced_ptms = get_modification_class_data(spliced_ptms, mod_class)
        background_ptms = get_modification_class_data(background_ptms, mod_class)

    #check to make sure requested annotation is available
    annotation_col = get_annotation_col(database = database, annotation_type = annotation_type)

    #get PTMs
    background_ptms['PTM'] = background_ptms['UniProtKB Accession'] + '_' + background_ptms['Residue']

    #grab all ptms spliced out for comparison to background
    spliced_ptm_list = np.unique(spliced_ptms['UniProtKB Accession'] + '_' + spliced_ptms['Residue'])

    #extract relevant annotation and remove PTMs without an annotation
    background_annotations = background_ptms[['PTM'] + [annotation_col]].copy()

    background_annotations[annotation_col] = background_annotations[annotation_col].apply(lambda x: x.split(';') if isinstance(x, str) else np.nan)
    background_annotations = background_annotations.explode(annotation_col)
    background_annotations[annotation_col] = background_annotations[annotation_col].apply(lambda x: x.strip() if isinstance(x, str) else np.nan)
    
    #if requested, group similar annotations (such as same function but increasing or decreasing)
    if collapse_on_similar:
        if database == 'PhosphoSitePlus' and annotation_type in ['Function', 'Process']:
            background_annotations[annotation_col] = background_annotations[annotation_col].apply(lambda x: x.split(',')[0].strip(' ') if x == x else x)
        else:
            background_annotations[annotation_col] = background_annotations[annotation_col].apply(lambda x: x.strip(' ') if x == x else x)
    else:
        background_annotations[annotation_col] = background_annotations[annotation_col].apply(lambda x: x.strip(' ') if x == x else x)
    
    #construct pivot table
    background_annotations['value'] = 1
    background_annotation = background_annotation[['PTM',annotation_col, 'value']].drop_duplicates()
    background_annotation = background_annotation.pivot(index = 'PTM', columns = annotation_col, values = 'value')

    #remove any sites with no annotations
    background_annotation = background_annotation.dropna(how = 'all')

    enrichment = stat_utils.get_site_enrichment(spliced_ptm_list, background_annotation, subset_name = 'Spliced', type = annotation_type, fishers = True)


    






