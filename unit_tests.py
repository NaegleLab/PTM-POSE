from ptm_pose import pose_config, project, annotate, helpers, nease_runner, flanking_sequences
from ptm_pose.analyze import interactions, enzyme, filter, flank_analysis, annotations, summarize
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


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

kstar_network_dir = '../../Database_Information/NETWORKS/NetworKIN/'



def configure_unit_test():
    pose_config.download_translator()

def projection_test(dataset = 'Yang2016_MATS', check_type = 'no changes', splicegraph_file = None):
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
    if dataset == 'Yang2016_MATS':
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

        #repeat to check multiprocessing for just SE and A3SS (also test filters)
        results2 = project.project_ptms_onto_MATS(SE_events = SE_events, fiveASS_events = fiveASS_events, coordinate_type = 'hg19', PROCESSES = 4, identify_flanking_sequences=True, min_studies = 2, remove_novel = True)

        #check if TSC2 S981 is in the spliced PTMs
        spliced_ptms['Label'] = spliced_ptms['Gene'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(str)
        assert 'TSC2_S981' in spliced_ptms['Label'].values or 'TSC2_S981.0' in spliced_ptms['Label'].values
    elif dataset == 'Yang_SE':
        #generic projection
        SE_events = pd.read_excel('./Test_datasets/Yang2016/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 2).iloc[0:179]

        #project ptms and identify altered flanking sequences using MATS function
        splice_data, spliced_ptms = project.project_ptms_onto_splice_events(splice_data = SE_events, gene_col = 'geneSymbol', strand_col = 'strand', chromosome_col='chr', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col = 'meanDeltaPSI', sig_col = 'FDR', coordinate_type = 'hg19', PROCESSES = 1)

        results2 = project.project_ptms_onto_splice_events(splice_data = SE_events, gene_col = 'geneSymbol', strand_col = 'strand', chromosome_col='chr', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col = 'meanDeltaPSI', sig_col = 'FDR', coordinate_type = 'hg19', PROCESSES = 4, min_MS_observations = 2, modification_class = 'Phosphorylation')
        
        
        altered_flanks = flanking_sequences.get_flanking_changes_from_splice_data(splice_data, chromosome_col='chr', strand_col = 'strand', spliced_region_end_col='exonEnd', spliced_region_start_col='exonStart_0base', first_flank_end_col='firstFlankingEE', first_flank_start_col='firstFlankingES', second_flank_start_col = 'secondFlankingES', second_flank_end_col='secondFlankingEE', gene_col='geneSymbol', dPSI_col = 'meanDeltaPSI', sig_col = 'FDR')
    elif dataset == 'TCGA_PRAD':
        if splicegraph_file is None:
            raise ValueError('splicegraph_dir must be specified for TCGA PRAD dataset')
        
        #load validation results for comparison
        spliced_ptms_true = pd.read_csv('./Test_Datasets/TCGA_PRAD/spliced_ptms.csv')
        altered_flanks_true = pd.read_csv('./Test_Datasets/TCGA_PRAD/altered_flanks.csv')

        #### TCGA PRAD (hg19) ######
        psi_data = pd.read_csv('./Test_datasets/TCGA_PRAD/sig_events_PRAD.csv')
        splicegraph = pd.read_csv(splicegraph_file, sep = '\t')

        results = project.project_ptms_onto_SpliceSeq(psi_data = psi_data, splicegraph = splicegraph, coordinate_type = 'hg19', PROCESSES = 1, identify_flanking_sequences=True, dPSI_col = 'deltaPSI_MW', sig_col = 'p-adj_MW', extra_cols = ['Effect Size_MW'], min_dpsi = 0.2)
        psi_data, spliced_ptms, altered_flanks = results

        #test multiprocessing and filtering
        #### TCGA PRAD (hg19) ######
        psi_data = pd.read_csv('./Test_datasets/TCGA_PRAD/sig_events_PRAD.csv')
        splicegraph = pd.read_csv(splicegraph_file, sep = '\t')
        results = project.project_ptms_onto_SpliceSeq(psi_data = psi_data, splicegraph = splicegraph, coordinate_type = 'hg19', PROCESSES = 4, identify_flanking_sequences=True, dPSI_col = 'deltaPSI_MW', sig_col = 'p-adj_MW', extra_cols = ['Effect Size_MW'], min_dpsi = 0.2, min_LTP_studies = 1)

        
        #check if TSC2 S981 is in the spliced PTMs
        spliced_ptms['Label'] = spliced_ptms['Gene'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(str)
        assert 'TSC2_S981' in spliced_ptms['Label'].values or 'TSC2_S981.0' in spliced_ptms['Label'].values
    else:
        raise ValueError('Dataset not recognized. Options include Yang2016 and TCGA_PRAD')

    #compare results to expected values, either raising an error or printing differences
    if check_type == 'no changes':
        spliced_ptms_true = np.unique(spliced_ptms_true['UniProtKB Accession'] + '_' + spliced_ptms_true['Residue'] + spliced_ptms_true['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        spliced_ptms_new = np.unique(spliced_ptms['UniProtKB Accession'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        intersection = set(spliced_ptms_true).intersection(set(spliced_ptms_new))
        assert len(intersection) == len(spliced_ptms_true)

        altered_flanks_true = np.unique(altered_flanks_true['UniProtKB Accession'] + '_' + altered_flanks_true['Residue'] + altered_flanks_true['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        altered_flanks_new = np.unique(altered_flanks['UniProtKB Accession'] + '_' + altered_flanks['Residue'] + altered_flanks['PTM Position in Isoform'].astype(float).astype(int).astype(str))
        intersection = set(altered_flanks_true).intersection(set(altered_flanks_new))
        assert len(intersection) == len(altered_flanks_true)
        print('No changes found! Projection test passed!')
    elif check_type == 'report changes':
        spliced_ptms_true = np.unique(spliced_ptms_true['UniProtKB Accession'] + '_' + spliced_ptms_true['Residue'] + spliced_ptms_true['PTM Position in Isoform'].astype(float).astype(int).astype(str))
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

        altered_flanks_true = np.unique(altered_flanks_true['UniProtKB Accession'] + '_' + altered_flanks_true['Residue'] + altered_flanks_true['PTM Position in Isoform'].astype(float).astype(int).astype(str))
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
    elif check_type == 'ignore':
        pass
    else:
        raise ValueError('check type must either be no changes, report changes, or ignore')

    return spliced_ptms, altered_flanks

def build_annotations_test(phosphositeplus = True, ptmsigdb = True, regphos = True, ptmcode = True,ptmint = True, depod = True, omnipath = True):
    if not os.path.exists('./Test/Test_Annotations/'):
        os.makedirs('./Test/Test_Annotations/')

    if phosphositeplus:
        try:
            annotate.construct_PhosphoSitePlus_gmt_files(regulatory_site_file=psp_regulatory_site_file, disease_association_file=psp_disease_file, kinase_substrate_file=psp_ks_file, odir = './Test/Test_Annotations/', overwrite = True)
        except Exception as e:
            
            print(f'Error building PhosphoSitePlus annotations: {e}')

    if ptmsigdb:
        try:
            annotate.construct_PTMsigDB_gmt_files(ptmsigdb_file, odir = './Test/Test_Annotations/', overwrite=True)
            print('PTMsigDB annotations built successfully!')
        except Exception as e:
            print(f'Error building PTMsigDB annotations: {e}')

    if ptmcode:
        try:
            annotate.construct_PTMcode_interprotein_gmt_file(odir = './Test/Test_Annotations/')
            print('PTMcode interprotein annotations built successfully!')
        except Exception as e:
            print(f'Error building PTMcode interprotein annotations: {e}')

    if depod:
        try:
            annotate.construct_DEPOD_gmt_file(odir = './Test/Test_Annotations/', overwrite=True)
            print('DEPOD annotations built successfully!')
        except Exception as e:
            print(f'Error building DEPOD annotations: {e}')
    
    if regphos:
        try:
            annotate.construct_RegPhos_gmt_file(odir = './Test/Test_Annotations/', file=RegPhos, overwrite=True)
            print('RegPhos annotations built successfully!')
        except Exception as e:
            print(f'Error building RegPhos annotations: {e}')

    if ptmint:
        try:
            annotate.construct_PTMInt_gmt_file(odir = './Test/Test_Annotations/')
            print('PTMint annotations built successfully!')
        except Exception as e:
            print(f'Error building PTMint annotations: {e}')

    if omnipath:
        try:
            
            annotate.construct_omnipath_gmt_file(odir = './Test/Test_Annotations/')
            print('Omnipath annotations built successfully!')
        except Exception as e:
            print(f'Error building Omnipath annotations: {e}')


def add_annotations_test(spliced_ptms, altered_flanks, check_type = 'report changes'):


    #save previous annotation data and then remove any previous annotations
    annotation_cols = [col for col in spliced_ptms.columns if ':' in col]
    spliced_annotation_counts = spliced_ptms[annotation_cols].count()
    spliced_ptms = spliced_ptms.drop(columns = annotation_cols)
    annotation_cols = [col for col in altered_flanks.columns if ':' in col]
    altered_annotation_counts = altered_flanks[annotation_cols].count()
    altered_flanks = altered_flanks.drop(columns = annotation_cols)

    print('Annotation differentially spliced PTMs')
    #differentially included ptms
    spliced_ptms = annotate.annotate_ptms(spliced_ptms, elm = True)

    #check if annotations are the same, if not, print differences
    if check_type == 'report changes':
        for annot in spliced_annotation_counts.index:
            if annot in spliced_ptms.columns:
                if spliced_annotation_counts[annot] != spliced_ptms[annot].count():
                    print(f'Annotation {annot} has changed from {spliced_annotation_counts[annot]} to {spliced_ptms[annot].count()}')
    elif check_type == 'ignore':
        pass
    else:
        raise ValueError('check_type must be either report changes or ignore')

    print('Annotation altered flanking sequences')
    #altered flanking sequences
    altered_flanks = annotate.annotate_ptms(altered_flanks)

    if check_type == 'report changes':
        #check if annotations are the same, if not, print differences
        for annot in spliced_annotation_counts.index:
            if annot in spliced_ptms.columns:
                if spliced_annotation_counts[annot] != spliced_ptms[annot].count():
                    print(f'Annotation {annot} has changed from {spliced_annotation_counts[annot]} to {spliced_ptms[annot].count()}')
    elif check_type == 'ignore':
        pass
    else:
        raise ValueError('check_type must be either report changes or ignore')

    return spliced_ptms, altered_flanks


def annotations_analysis_test(spliced_ptms, altered_flanks):
    print('Combining Interaction Information')
    combined = summarize.combine_outputs(spliced_ptms = spliced_ptms, altered_flanks = altered_flanks)
    print('Getting annotation summaries')
    available_annotations = annotations.get_available_annotations(combined)
    print(available_annotations)

    annotations.plot_available_annotations(spliced_ptms)
    plt.show()
    

    
    for i, row in available_annotations.iterrows():
        database = row['Database']
        annot_type = row['Annotation Type']
        print(f'Testing analysis on {database} {annot_type} annotation type')
        for collapse_annotations in [True, False]:
            print(annot_type, database)
            annot, annotation_counts = annotations.get_ptm_annotations(combined, annot_type = annot_type, database = database)
            print('Calculating enrichment with precomputed background')
            enrichment = annotations.annotation_enrichment(combined, annot_type = annot_type, database = database,  collapse_on_similar = collapse_annotations)

            if annot_type == 'Function' and database == 'PhosphoSitePlus':
                print('Calculating enrichment by significance')
                enrichment = annotations.annotation_enrichment(combined, annot_type = annot_type, database = database, collapse_on_similar = collapse_annotations, background_type = 'significance', alpha = 0.01, min_dPSI = 0.25)

            if annot_type == 'Function' and database == 'PhosphoSitePlus':
                print('Calculating enrichment for specific modifications')
                enrichment = annotations.annotation_enrichment(combined, annot_type = annot_type, database = database, collapse_on_similar = collapse_annotations, modification_class='Phosphorylation')
            
            annotations.plot_annotation_counts(combined, database=database, annot_type=annot_type, collapse_on_similar = collapse_annotations)
            plt.show()

    print('Annotation unit test passed!')

def misc_analysis_tests(spliced_ptms, altered_flanks):
    print('Testing overview analysis')
    summarize.plot_modification_breakdown(spliced_ptms = spliced_ptms, altered_flanks = altered_flanks)
    plt.show()

    print('Testing gene set enrichment')
    enr_results = annotations.gene_set_enrichment(spliced_ptms = spliced_ptms, altered_flanks=altered_flanks, min_dPSI = 0.1)
    annotations.plot_EnrichR_pies(enr_results)
    plt.show()

    print('Analysis unit test passed!')

def interactions_unit_test(spliced_ptms, splice_data):
    print('Testing interaction network functions')
    protein_interactions = interactions.protein_interactions(spliced_ptms, remove_novel = True)
    protein_interactions.get_interaction_stats()
    protein_interactions.summarize_protein_network('TSC2')
    protein_interactions.get_protein_interaction_network('TSC2')
    protein_interactions.plot_interaction_network()
    plt.show()
    protein_interactions.plot_interaction_network(color_edges_by='Database')
    protein_interactions.plot_network_centrality()
    plt.show()

    print('Testing NEASE runner')
    nease_output = nease_runner.run_nease(splice_data, region_start_col = 'exonStart_0base', gene_col = 'geneSymbol', chromosome_col = 'chr', strand_col = 'strand', gene_id_type='name', region_end_col='exonEnd', coordinate_type = 'hg19', dpsi_col = 'meanDeltaPSI', remove_non_in_frame = True, only_divisible_by_3 = True)

    nease_output = nease_runner.run_nease(splice_data, region_start_col = 'exonStart_0base', gene_col = 'GeneID', chromosome_col = 'chr', strand_col = 'strand', gene_id_type='Ensembl', region_end_col='exonEnd', coordinate_type = 'hg19', dpsi_col = 'meanDeltaPSI', remove_non_in_frame = False, only_divisible_by_3 = False)

    nease_runner.save_nease(nease_output, odir = './Test/NEASE/', file_type = 'excel')

    nease_interactions, nease_slims, nease_residues, nease_domains = nease_runner.load_nease(odir='./Test/NEASE/', file_type='excel')

    print('Testing NEASE compariosn')
    protein_interactions.compare_to_nease(nease_interactions)
    protein_interactions.plot_nease_comparison(nease_interactions)





def flank_analysis_test(altered_flanks):
    print('Testing altered flanking sequence analysis')
    altered_flanks = flank_analysis.compare_flanking_sequences(altered_flanks)
    altered_flanks = flank_analysis.compare_inclusion_motifs(altered_flanks)

    flank_analysis.plot_location_of_altered_flanking_residues(altered_flanks, modification_class = 'Phosphorylation', min_studies = 1)
    plt.show()

    flank_analysis.plot_alterations_matrix(altered_flanks.head(10), min_dpsi = 0)
    plt.show()

    test = flank_analysis.identify_change_to_specific_motif(altered_flanks, elm_motif_name = 'LIG_14-3-3_CterR_2', residues = 'S', modification_class='Phosphorylation', min_studies = 0)
    print(test)
    

    flank_ex = altered_flanks[(altered_flanks['Gene'] == 'TSC2') & (altered_flanks['PTM Position in Isoform'] == 946)].squeeze().dropna()
    flank_analysis.plot_sequence_differences(flank_ex['Inclusion Flanking Sequence'], flank_ex['Exclusion Flanking Sequence'], dpsi = flank_ex['dPSI'])
    plt.show()
    
def enzyme_analysis_test(spliced_ptms, altered_flanks, test_kstar = True, test_ksea = True, test_kinase_library = True):
    if test_kstar:
        print("Testing KSTAR analysis")
        kstar = enzyme.kstar_enrichment(spliced_ptms, network_dir = kstar_network_dir, phospho_type = ['Y'])
        kstar = enzyme.kstar_enrichment(spliced_ptms, network_dir = kstar_network_dir, phospho_type = ['Y', 'ST'], remove_novel = True, min_studies = 1)
        print(kstar.return_enriched_kinases())
        kstar.dotplot(ptype = 'Y', impact_types = ['Included', 'Excluded'], sig_kinases_only=True)

    if test_ksea:
        print('Testing KSEA analysis')
        ksea = enzyme.KSEA(spliced_ptms, database = 'PhosphoSitePlus')
        ksea.runKSEA()
        ksea = enzyme.KSEA(spliced_ptms, database = 'Combined Writer', remove_novel = True, min_studies = 1)
        ksea.runKSEA()
        ksea = enzyme.KSEA(spliced_ptms, database = 'OmniPath Writer', min_dpsi = 0.15)
        ksea.runKSEA()

        ksea.plot_results()
        plt.show()
    
    if test_kinase_library:
        print('Testing kinase library')
        KL = enzyme.KL_flank_analysis(altered_flanks, min_dpsi = 0.15, num_studies = 1)
        KL.analyze_single_ptm(gene = 'ARHGAP17', loc = 497)
        KL.analyze_all_ptms()

        KL.get_kinases_with_largest_changes(kinase_type = 'Y', difference_type = 'normal')
        KL.get_kinases_with_largest_changes(kinase_type = 'Y', difference_type = 'absolute')
        KL.get_kinases_with_largest_changes(kinase_type = 'ST', difference_type = 'relative')

        KL.plot_top_kinases(kinase_type = 'Y', difference_type = 'normal')
        KL.plot_top_kinases(kinase_type = 'ST', difference_type = 'relative')
        plt.show()

        KL.plot_top_changes()
        plt.show()
        KL.plot_top_changes(gene = 'TSC2')
        

def filter_analysis_test(spliced_ptms):
    min_dpsi = 0.1
    alpha = 0.05
    min_studies = 0
    min_MS_observations = 0
    min_LTP_studies = 0
    min_compendia = 2
    remove_novel = True

    print('testing plotting single filter impact')
    filter.plot_filter_impact(spliced_ptms, min_dpsi = min_dpsi, alpha = alpha, min_studies = min_studies, min_MS_observations = min_MS_observations, min_LTP_studies = min_LTP_studies, min_compendia = min_compendia, remove_novel = remove_novel)

    print('testing plotting filter range')
    filter.assess_filter_range(spliced_ptms, min_value = 0, max_value = 30, filter_type = 'min_studies')

    filter.assess_filter_range(spliced_ptms, min_value = 0, max_value = 30, filter_type = 'min_LTP', phospho_only_evidence_filter=True)

    filter.assess_filter_range(spliced_ptms, min_value = 0.05, max_value = 0.5, filter_type = 'min_MS')

    filter.assess_filter_range(spliced_ptms, min_value = 0.05, max_value = 0.5, filter_type = 'min_LTP')



def comprehensive_unit_test():
    pass