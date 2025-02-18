from ptm_pose import pose_config, project, annotate, analyze
from ptm_pose import plots as pose_plots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


database_dir = '../../Database_Information/'
psp_regulatory_site_file = database_dir + 'PhosphoSitePlus/Regulatory_sites.gz'
psp_ks_file = database_dir + 'PhosphoSitePlus/Kinase_Substrate_Dataset.gz'
psp_disease_file = database_dir + 'PhosphoSitePlus/Disease-associated_sites.gz'
elm_interactions = database_dir + '/ELM/elm_interactions.tsv'
elm_motifs = database_dir + '/ELM/elm_classes.tsv'
PTMint = database_dir + '/PTMInt/PTM_experimental_evidence.csv'
PTMcode_interprotein = database_dir + '/PTMcode/PTMcode2_associations_between_proteins.txt.gz'
DEPOD = True
RegPhos = database_dir + '/RegPhos/RegPhos_Phos_human.txt'
ptmsigdb_file = '../../Database_Information/PTMsigDB/data_PTMsigDB_all_sites_v2.0.0.xlsx'

file_dict = {'PhosphoSitePlus': {'Function':psp_regulatory_site_file, 'Process':psp_regulatory_site_file, 'Disease':psp_disease_file, 'Kinase':psp_ks_file, 'Interactions': psp_regulatory_site_file, 'Perturbation': ptmsigdb_file},
                       'ELM': {'Interactions':elm_interactions, 'Motif Match':elm_motifs},
                       'PTMcode': {'Interactions':PTMcode_interprotein},
                       'PTMInt': {'Interactions':PTMint},
                       'RegPhos': {'Kinase':RegPhos},
                       'DEPOD': {'Phosphatase':DEPOD},
                       'PTMsigDB':{'WikiPathway':ptmsigdb_file, 'NetPath':ptmsigdb_file,'mSigDB':ptmsigdb_file, 'Perturbation (DIA2)':ptmsigdb_file, 'Perturbation (DIA)': ptmsigdb_file, 'Perturbation (PRM)':ptmsigdb_file,'Kinase':ptmsigdb_file}}

def configure_unit_test():
    pose_config.download_translator()
    #pose_config.download_ptm_coordinates()

def projection_test(dataset = 'Yang2016', check_type = 'no changes', splicegraph_file = None):
    """
    Test the projection of PTMs onto splicing events to identify differentially included PTMs and PTMs with altered flanking sequences using the project model. Compare to see if results changed from previous runs.

    Parameters
    ----------
    dataset : str
        Dataset to test. Options include 'Yang2016' and 'TCGA_PRAD' (default is 'Yang2016')
    check_type : str
        Type of check to perform. Options include 'no changes', which will raise an error if there are differences from previous results, and 'report changes', which will print PTMs that are different (default is 'no changes')
    splicegraph_file : str
        Path to splicegraph file for TCGA PRAD dataset, which is only needed for splicegraph (default is None)
    """
    #load dataset of interest and project PTMs
    if dataset == 'Yang2016':
        #load validation results for comparison
        spliced_ptms_true = pd.read_csv('./Test_Datasets/Yang2016/spliced_ptms.csv')
        altered_flanks_true = pd.read_csv('./Test_Datasets/Yang2016/altered_flanks.csv')
        #### Yang 2016 (hg19) ######
        siRNA_data = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 2)
        SE_events = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 2).iloc[0:179]
        MXE_events = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 183).iloc[0:37]
        fiveASS_events = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 222).iloc[0:6]
        threeASS_events = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 230).iloc[0:7]
        RI_events = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 239)
    
        #project ptms and identify altered flanking sequences using MATS function
        results = project.project_ptms_onto_MATS(SE_events = SE_events, MXE_events = MXE_events, RI_events = RI_events, threeASS_events = threeASS_events, fiveASS_events = fiveASS_events, coordinate_type = 'hg19', PROCESSES = 1, identify_flanking_sequences=True)
        splice_data, spliced_ptms, altered_flanks = results

        #check if TSC2 S981 is in the spliced PTMs
        spliced_ptms['Label'] = spliced_ptms['Gene'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(str)
        assert 'TSC2_S981' in spliced_ptms['Label'].values or 'TSC2_S981.0' in spliced_ptms['Label'].values

    elif dataset == 'TCGA_PRAD':
        if splicegraph_file is None:
            raise ValueError('splicegraph_dir must be specified for TCGA PRAD dataset')
        
        #load validation results for comparison
        spliced_ptms_true = pd.read_csv('./Test_Datasets/TCGA_PRAD/spliced_ptms.csv')
        altered_flanks_true = pd.read_csv('./Test_Datasets/TCGA_PRAD/altered_flanks.csv')

        #### TCGA PRAD (hg19) ######
        psi_data = pd.read_csv('./Test_datasets/TCGA_PRAD/sig_events_PRAD.csv')
        splicegraph = pd.read_csv(splicegraph_file, sep = '\t')

        results = project.project_ptms_onto_SpliceSeq(psi_data = psi_data, splicegraph = splicegraph, coordinate_type = 'hg19', PROCESSES = 1, identify_flanking_sequences=True, dPSI_col = 'deltaPSI_MW', sig_col = 'p-adj_MW', extra_cols = ['Effect Size_MW'])
        psi_data, spliced_ptms, altered_flanks = results

        
        #check if TSC2 S981 is in the spliced PTMs
        spliced_ptms['Label'] = spliced_ptms['Gene'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(str)
        assert 'TSC2_S981' in spliced_ptms['Label'].values or 'TSC2_S981.0' in spliced_ptms['Label'].values
    else:
        raise ValueError('Dataset not recognized. Options include Yang2016 and TCGA_PRAD')

    #compare results to expected values, either raising an error or printing differences
    if check_type == 'no changes':
        spliced_ptms_true = np.unique(spliced_ptms_true['UniProtKB Accession'] + '_' + spliced_ptms_true['Residue'] + spliced_ptms_true['PTM Position in Canonical Isoform'].astype(float).astype(int).astype(str))
        spliced_ptms_new = np.unique(spliced_ptms['UniProtKB Accession'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        intersection = set(spliced_ptms_true).intersection(set(spliced_ptms_new))
        assert len(intersection) == len(spliced_ptms_true)

        altered_flanks_true = np.unique(altered_flanks_true['UniProtKB Accession'] + '_' + altered_flanks_true['Residue'] + altered_flanks_true['PTM Position in Canonical Isoform'].astype(float).astype(int).astype(str))
        altered_flanks_new = np.unique(altered_flanks['UniProtKB Accession'] + '_' + altered_flanks['Residue'] + altered_flanks['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        intersection = set(altered_flanks_true).intersection(set(altered_flanks_new))
        assert len(intersection) == len(altered_flanks_true)
        print('No changes found! Projection test passed!')
    elif check_type == 'report changes':
        spliced_ptms_true = np.unique(spliced_ptms_true['UniProtKB Accession'] + '_' + spliced_ptms_true['Residue'] + spliced_ptms_true['PTM Position in Canonical Isoform'].astype(float).astype(int).astype(str))
        spliced_ptms_new = np.unique(spliced_ptms['UniProtKB Accession'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        only_in_reference = set(spliced_ptms_true).difference(set(spliced_ptms_new))
        only_in_new = set(spliced_ptms_new).difference(set(spliced_ptms_true))
        print('Differentially included PTMs:')
        if len(only_in_reference) > 0 or len(only_in_new) > 0:
            print('Only in reference: {}'.format(', '.join(only_in_reference)))
            print('Only in new: {}'.format(', '.join(only_in_new)))
            print('\n')
        else:
            print('No differences found')
            print('\n')

        altered_flanks_true = np.unique(altered_flanks_true['UniProtKB Accession'] + '_' + altered_flanks_true['Residue'] + altered_flanks_true['PTM Position in Canonical Isoform'].astype(float).astype(int).astype(str))
        altered_flanks_new = np.unique(altered_flanks['UniProtKB Accession'] + '_' + altered_flanks['Residue'] + altered_flanks['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        only_in_reference = set(altered_flanks_true).difference(set(altered_flanks_new))
        only_in_new = set(altered_flanks_new).difference(set(altered_flanks_true))
        print('Altered flanking sequences:')
        if len(only_in_reference) > 0 or len(only_in_new) > 0:
            print('Only in reference: {}'.format(', '.join(only_in_reference)))
            print('Only in new: {}'.format(', '.join(only_in_new)))
            print('\n')
        else:
            print('No differences found')
            print('\n')

    return spliced_ptms, altered_flanks

def add_annotations_test(spliced_ptms, altered_flanks):


    #save previous annotation data and then remove any previous annotations
    annotation_cols = [col for col in spliced_ptms.columns if ':' in col]
    spliced_annotation_counts = spliced_ptms[annotation_cols].count()
    spliced_ptms = spliced_ptms.drop(columns = annotation_cols)
    annotation_cols = [col for col in altered_flanks.columns if ':' in col]
    altered_annotation_counts = altered_flanks[annotation_cols].count()
    altered_flanks = altered_flanks.drop(columns = annotation_cols)

    print('Annotation differentially spliced PTMs')
    #differentially included ptms
    spliced_ptms = annotate.annotate_ptms(spliced_ptms, psp_regulatory_site_file=psp_regulatory_site_file, psp_ks_file=psp_ks_file, psp_disease_file=psp_disease_file, elm_interactions=elm_interactions, elm_motifs = elm_motifs, PTMint=PTMint,PTMcode_interprotein=PTMcode_interprotein, DEPOD=DEPOD, RegPhos=RegPhos, ptmsigdb_file=ptmsigdb_file)

    #check if annotations are the same, if not, print differences
    for annot in spliced_annotation_counts.index:
        if annot in spliced_ptms.columns:
            if spliced_annotation_counts[annot] != spliced_ptms[annot].count():
                print(f'Annotation {annot} has changed from {spliced_annotation_counts[annot]} to {spliced_ptms[annot].count()}')

    print('Annotation altered flanking sequences')
    #altered flanking sequences
    altered_flanks = annotate.annotate_ptms(altered_flanks, psp_regulatory_site_file=psp_regulatory_site_file, psp_ks_file=psp_ks_file, psp_disease_file=psp_disease_file, elm_interactions=elm_interactions, elm_motifs = elm_motifs, PTMint=PTMint, PTMcode_interprotein=PTMcode_interprotein, DEPOD=DEPOD, RegPhos=RegPhos, ptmsigdb_file=ptmsigdb_file)

    #check if annotations are the same, if not, print differences
    for annot in spliced_annotation_counts.index:
        if annot in spliced_ptms.columns:
            if spliced_annotation_counts[annot] != spliced_ptms[annot].count():
                print(f'Annotation {annot} has changed from {spliced_annotation_counts[annot]} to {spliced_ptms[annot].count()}')

    return spliced_ptms, altered_flanks

def get_file(database = 'PhosphoSitePlus', annotation_type = 'Function'):
    return file_dict[database][annotation_type]

def annotations_analysis_test(spliced_ptms, altered_flanks):
    print('Combining Interaction Information')
    combined = analyze.combine_outputs(spliced_ptms = spliced_ptms, altered_flanks = altered_flanks)
    print('Getting annotation summaries')
    print(analyze.get_annotation_categories(combined))

    pose_plots.show_available_annotations(spliced_ptms)
    plt.show()
    

    for database in pose_config.annotation_col_dict.keys():
        for annot_type in pose_config.annotation_col_dict[database]:
            if not (database == 'ELM' and annot_type == 'Motif Match') and not (database == 'PTMcode' and annot_type == 'Intraprotein'):
                print(f'Testing analysis on {database} {annot_type} annotation type')
                for collapse_annotations in [True, False]:
                    annotations, annotation_counts = analyze.get_ptm_annotations(combined, annotation_type = annot_type, database = database)
                    print('Calculating enrichment with precomputed background')
                    file = get_file(database = database, annotation_type = annot_type)
                    enrichment = analyze.annotation_enrichment(combined, annotation_type = annot_type, database = database, annotation_file = file, collapse_on_similar = collapse_annotations, save_background = True)

                    if annot_type == 'Function' and database == 'PhosphoSitePlus':
                        print('Calculating enrichment by significance')
                        enrichment = analyze.annotation_enrichment(combined, annotation_type = annot_type, database = database, collapse_on_similar = collapse_annotations, annotation_file=file, background_type = 'significance', alpha = 0.01, min_dPSI = 0.2)

                    if annot_type == 'Function' and database == 'PhosphoSitePlus':
                        print('Calculating enrichment for specific modifications')
                        enrichment = analyze.annotation_enrichment(combined, annotation_type = annot_type, database = database, collapse_on_similar = collapse_annotations, annotation_file = file, mod_class='Phosphorylation')
                    
                    pose_plots.plot_annotations(combined, database=database, annot_type=annot_type, collapse_on_similar = collapse_annotations)
                    plt.show()

    print('Annotation unit test passed!')

def analysis_unit_tests(spliced_ptms, altered_flanks, kstar_network_dir):
    print('Testing overview analysis')
    pose_plots.modification_breakdown(spliced_ptms = spliced_ptms, altered_flanks = altered_flanks)
    plt.show()

    print('Testing gene set enrichment')
    enr_results = analyze.gene_set_enrichment(spliced_ptms = spliced_ptms, altered_flanks=altered_flanks, min_dPSI = 0.1)
    pose_plots.plot_EnrichR_pies(enr_results)
    plt.show()

    print('Testing interaction network functions')
    protein_interactions = analyze.protein_interactions(spliced_ptms)
    protein_interactions.get_interaction_network()
    protein_interactions.get_interaction_stats()
    protein_interactions.plot_interaction_network()
    plt.show()
    protein_interactions.plot_network_centrality()
    plt.show()

    print('Testing altered flanking sequence analysis')
    altered_flanks = analyze.compare_flanking_sequences(altered_flanks)
    altered_flanks = analyze.compare_inclusion_motifs(altered_flanks)

    pose_plots.location_of_altered_flanking_residues(altered_flanks, modification_class = 'Phosphorylation')
    plt.show()

    pose_plots.alterations_matrix(altered_flanks.head(10))
    plt.show()

    print("Testing KSTAR analysis")
    kstar = analyze.kstar_enrichment(spliced_ptms, network_dir = kstar_network_dir, phospho_type = 'Y')
    print(kstar.return_enriched_kinases())
    kstar = analyze.kstar_enrichment(spliced_ptms, network_dir=kstar_network_dir, phospho_type = 'ST')
    print(kstar.return_enriched_kinases())

    print('Analysis unit test passed!')


def comprehensive_unit_test():
    pass