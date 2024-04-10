import numpy as np
import pandas as pd

#plotting 
import matplotlib.pyplot as plt

def enrichment_analysis(spliced_ptms, annotation):
    pass


def show_available_annotations(spliced_ptms, figsize = (5, 5)):
    num_ptms = [spliced_ptms.shape[0]]
    num_ptms_filters = ['All PTMs']
    filter_source = ['None']
    if 'PSP:DOMAIN' in spliced_ptms.columns:
        #in domain or not
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:DOMAIN').shape[0])
        num_ptms_filters.append('PTMs In Domain')
        filter_source.append('PhosphoSitePlus')
       
        #with a biological process annotation
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:ON_PROCESS').shape[0])
        num_ptms_filters.append('PTMs associated with Biological Process')
        filter_source.append('PhosphoSitePlus')

        #with a molecular function annoation in PSP
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:ON_FUNCTION').shape[0])
        num_ptms_filters.append('PTMs associated with Molecular Function')
        filter_source.append('PhosphoSitePlus')
    if 'PSP:Kinase' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:Kinase').shape[0])
        num_ptms_filters.append('PTMs associated with Kinase')
        filter_source.append('PhosphoSitePlus')
    if 'PSP:Disease_Association' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'PSP:Disease_Association').shape[0])
        num_ptms_filters.append('PTMs with Disease Association')
        filter_source.append('PhosphoSitePlus')
    if 'ELM Interactions' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms.dropna(subset = 'ELM Interactions').shape[0])
        num_ptms_filters.append('PTMs with Known ELM Interactions')
        filter_source.append('ELM')
    if 'ELM Motif Matches' in spliced_ptms.columns:
        num_ptms.append(spliced_ptms[spliced_ptms['ELM Motif Matches'] != ''].shape[0])
        num_ptms_filters.append('PTMs with ELM Motif Matches')
        filter_source.append('ELM')

        #with ligand motifs
        num_ptms.append(spliced_ptms[spliced_ptms['ELM Motif Matches'].str.contains('LIG')].shape[0])
        num_ptms_filters.append('PTMs with ELM Ligand Motif Matches')
        filter_source.append('ELM')

        #with phospho dependent motifs
        sh2_motifs = ['LIG_SH2_CRK', 'LIG_SH2_GRB2like', 'LIG_SH2_NCK_1', 'LIG_SH2_PTP2', 'LIG_SH2_SFK_2', 'LIG_SH2_SFK_CTail_3', 'LIG_SH2_STAP1', 'LIG_SH2_STAT3', 'LIG_SH2_STAT5', 'LIG_SH2_STAT6']
        ptb_motifs = ['LIG_PTB_Phospho_1']
        fha_motifs = ['LIG_FHA_1', 'LIG_FHA_2']
        other_motifs = ['LIG_TYR_ITAM', 'LIG_TYR_ITIM', 'LIG_TYR_ITSM', 'LIG_RBL1_LxSxE_2', 'LIG_BRCT_BRCA1_1', 'LIG_BRCT_BRCA1_2', 'LIG_BRCT_MDC1_1', 'LIG_DLG_GKlike_1']
        fourteen33_motifs = ['LIG_14-3-3_CanoR_1', 'LIG_14-3-3_CterR_2']
        phospho_dependent_motifs = sh2_motifs + ptb_motifs + fha_motifs + other_motifs + fourteen33_motifs

        num_ptms.append(spliced_ptms[(spliced_ptms['ELM Motif Matches'].str.contains('|'.join(phospho_dependent_motifs))) & (spliced_ptms['Modifications'].str.contains('Phospho'))].shape[0])
        num_ptms_filters.append('PTMs with Phospho Dependent ELM Motif Matches')
        filter_source.append('ELM')
    
    #plot bar plot
    #color bars based on datasource
    palette = {'None': 'gray', 'PhosphoSitePlus': 'blue', 'ELM': 'green'}
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



