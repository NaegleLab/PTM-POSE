from ptm_pose import pose_config, project, annotate, analyze
from ptm_pose import plots as pose_plots
import matplotlib.pyplot as plt
import pandas as pd



psp_regulatory_site_file = '../../Database_Information/PhosphoSitePlus/Regulatory_sites.gz'
psp_ks_file = '../../Database_Information/PhosphoSitePlus/Kinase_Substrate_Dataset.gz'
psp_disease_file = '../../Database_Information/PhosphoSitePlus/Disease-associated_sites.gz'
elm_interactions = '../../Database_Information/ELM/elm_interactions.tsv'
elm_motifs = '../../Database_Information/ELM/elm_classes.tsv'
PTMint = True
PTMcode_interprotein = '../../Database_Information/PTMcode/PTMcode2_associations_between_proteins.txt.gz'
DEPOD = True
RegPhos = True

file_dict = {'PhosphoSitePlus': {'Function':psp_regulatory_site_file, 'Process':psp_regulatory_site_file, 'Disease':psp_disease_file, 'Kinase':psp_ks_file, 'Interactions': psp_regulatory_site_file},
                       'ELM': {'Interactions':elm_interactions, 'Motif Match':elm_motifs},
                       'PTMcode': {'Interactions':PTMcode_interprotein},
                       'PTMInt': {'Interactions':None},
                       'RegPhos': {'Kinase':None},
                       'DEPOD': {'Phosphatase':None}}


def short_unit_test(dataset = 'Yang2016', repeat_setup = False):
    if repeat_setup:
        pose_config.download_translator(save = True)
        pose_config.ptm_coordinates = pose_config.download_ptm_coordinates(save = True)

    if dataset == 'Yang2016':
        #### Yang 2016 (hg19) ######
        siRNA_data = pd.read_excel('./Test_datasets/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 2)
        SE_events = pd.read_excel('./Test_datasets/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 2).iloc[0:179]
        MXE_events = pd.read_excel('./Test_datasets/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 183).iloc[0:37]
        fiveASS_events = pd.read_excel('./Test_datasets/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 222).iloc[0:6]
        threeASS_events = pd.read_excel('./Test_datasets/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 230).iloc[0:7]
        RI_events = pd.read_excel('./Test_datasets/esrp1_knockdown_data_Yang2016.xlsx', sheet_name='rMATS ESRP KD', header = 239)
    

        results = project.project_PTMs_onto_MATS(SE_events = SE_events, MXE_events = MXE_events, RI_events = RI_events, threeASS_events = threeASS_events, fiveASS_events = fiveASS_events, coordinate_type = 'hg19', PROCESSES = 1, identify_flanking_sequences=True)
        splice_data, spliced_ptms, altered_flanks = results

        #check if TSC2 S981 is in the spliced PTMs
        spliced_ptms['Label'] = spliced_ptms['Gene'] + '_' + spliced_ptms['Residue'] + spliced_ptms['PTM Position in Canonical Isoform'].astype(str)
        assert 'TSC2_S981' in spliced_ptms['Label'].values or 'TSC2_S981.0' in spliced_ptms['Label'].values

        expected_spliced_ptm_size = 455
        expected_altered_flanks_size = 105

        #check dataframes are the right size
        assert len(splice_data) == expected_spliced_ptm_size
        assert len(spliced_ptms) == expected_altered_flanks_size


    spliced_ptms = annotate.annotate_ptms(spliced_ptms, psp_regulatory_site_file=psp_regulatory_site_file, psp_ks_file=psp_ks_file, psp_disease_file=psp_disease_file, elm_interactions=elm_interactions, elm_motifs = elm_motifs, PTMint=PTMint, PTMcode_intraprotein=False, PTMcode_interprotein=PTMcode_interprotein, DEPOD=DEPOD, RegPhos=RegPhos)
    altered_flanks = annotate.annotate_ptms(altered_flanks, psp_regulatory_site_file=psp_regulatory_site_file, psp_ks_file=psp_ks_file, psp_disease_file=psp_disease_file, elm_interactions=elm_interactions, elm_motifs = elm_motifs, PTMint=PTMint, PTMcode_intraprotein=False, PTMcode_interprotein=PTMcode_interprotein, DEPOD=DEPOD, RegPhos=RegPhos)

    #recheck dataframes are the right size
    assert len(splice_data) == expected_spliced_ptm_size
    assert len(spliced_ptms) == expected_altered_flanks_size

    combined = analyze.combine_outputs(spliced_ptms = spliced_ptms, altered_flanks = altered_flanks)
    print(analyze.get_annotation_categories(combined))

    pose_plots.show_available_annotations(spliced_ptms, figsize = (5,5))
    plt.show()
    

    for database in pose_config.annotation_col_dict.keys():
        print('Testing analysis on {} database'.format(database))
        for annot_type in pose_config.annotation_col_dict[database]:
            print('Testing analysis on {} annotation type'.format(annot_type))
            for collapse_annotations in [True, False]:
                annotations, annotation_counts = analyze.get_ptm_annotations(combined, annotation_type = annot_type, database = database)
                print('Calculating enrichment with precomputed background')
                enrichment = analyze.annotation_enrichment(combined, annotation_type = annot_type, database = database)
                
                #add spot to calculate enrichment from custom background and based on signficance
                ##### here #####

    #test plotting functions
    print('Testing annotation plotting functions')
    analyze.plot_annotations(spliced_ptms, annotation_type = 'Function', database = 'PhosphoSitePlus')
    plt.show()
    pose_plots.plot_annotations(spliced_ptms, annotation_type = 'Interactions', database = 'PTMcode')
    plt.show()

    print('Testing gene set enrichment')
    enr_results = analyze.gene_set_enrichment(combined = combined, min_dPSI = 0.1)
    pose_plots.plot_EnrichR_pies(enr_results)

    print('Testing interaction network functions')
    interaction_graph, network_data = analyze.create_interaction_network(combined, database = 'PTMcode')
    network_stats = analyze.get_network_stats(interaction_graph)
    pose_plots.plot_interaction_network(interaction_graph, network_data, network_stats)
    pose_plots.plot_network_centrality(network_stats, network_data = network_data)


def comprehensive_unit_test():
    pass