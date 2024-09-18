from ptm_pose import pose_config, project, annotate, analyze
from ptm_pose import plots as pose_plots
from ptm_pose import flanking_sequences as fs

import datetime
import os
import pandas as pd


print('Warning: this module is still in the development phase and is not fully functional. Please check documentation for recommended methods for projection.')

class POSE_Project:
    def __init__(self, splice_data, ptm_coordinates = None, splice_data_type = 'generic', chromosome_col = None, strand_col = None, region_start_col = None, region_end_col = None, first_flank_start_col = None, first_flank_end_col = None, second_flank_start_col = None, second_flank_end_col = None, gene_name_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None, identify_flanking_sequences = False, splicegraph = None):
        self.splice_data = splice_data
        self.splice_data_type = splice_data_type
        self.splicegraph = splicegraph

        #set or download ptm_coordinates
        if ptm_coordinates is None and pose_config.ptm_coordinates is None:
            self.ptm_coordinates = pose_config.download_ptm_coordinates()
        elif ptm_coordinates is None and pose_config.ptm_coordinates is not None:
            self.ptm_coordinates = pose_config.ptm_coordinates.copy()
        elif ptm_coordinates is not None:
            self.ptm_coordinates = ptm_coordinates

        if isinstance(splice_data, dict):
            for key in splice_data.keys():
                self.check_dataframe_size(splice_data = splice_data[key], splice_data_label = f'{key} event data')
        else:
            self.check_dataframe_size(splice_data = splice_data)


        self.get_column_names(splice_data_type = splice_data_type, chromosome_col = chromosome_col, strand_col = strand_col, region_start_col = region_start_col, region_end_col = region_end_col, first_flank_start_col = first_flank_start_col, first_flank_end_col = first_flank_end_col, second_flank_start_col = second_flank_start_col, second_flank_end_col = second_flank_end_col, gene_name_col = gene_name_col, dPSI_col = dPSI_col, sig_col = sig_col, event_id_col = event_id_col, extra_cols = extra_cols, identify_flanking_sequences = identify_flanking_sequences)

    def check_dataframe_size(self, splice_data, splice_data_label = 'splice data'):
        if len(splice_data.columns) < 4 and not self.identify_flanking_sequences:
            raise ValueError(f'The {splice_data_label} must have at least 4 columns, consisting of columns with information on the chromosome, DNA strand, region start, and region end')
        elif len(splice_data.columns) < 8 and self.identify_flanking_sequences:
            raise ValueError(f'The {splice_data_label} must have at least 8 columns, consisting of columns with information on the chromosome, DNA strand, region start, region end, first flank start, first flank end, second flank start, and second flank end. If flanking sequences are not desired and/or you do not have flanking region information, set identify_flanking_sequences = False.')

    def get_column_names(self, splice_data_type = None, chromosome_col = None, strand_col = None, region_start_col = None, region_end_col = None, first_flank_start_col = None, first_flank_end_col = None, second_flank_start_col = None, second_flank_end_col = None, gene_name_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None, identify_flanking_sequences = False):
        self.identify_flanking_sequences = identify_flanking_sequences
        if splice_data_type is None:
            splice_data_type = self.splice_data_type

        if splice_data_type == 'generic':
            self.chromosome_col = self.splice_data.columns[0] if chromosome_col is None else chromosome_col
            self.strand_col = self.splice_data.columns[1] if strand_col is None else strand_col
            self.region_start_col = self.splice_data.columns[2] if region_start_col is None else region_start_col
            self.region_end_col = self.splice_data.columns[3] if region_end_col is None else region_end_col
            if identify_flanking_sequences:
                self.first_flank_start_col = self.splice_data.columns[4] if first_flank_start_col is None else first_flank_start_col
                self.first_flank_end_col = self.splice_data.columns[5] if first_flank_end_col is None else first_flank_end_col
                self.second_flank_start_col = self.splice_data.columns[6] if second_flank_start_col is None else second_flank_start_col
                self.second_flank_end_col = self.splice_data.columns[7] if second_flank_end_col is None else second_flank_end_col
            else:
                self.first_flank_start_col = None
                self.first_flank_region_end_col = None
                self.second_flank_start_col = None
                self.second_flank_end_col = None

            self.gene_name_col = None if gene_name_col is None else gene_name_col
            self.dPSI_col = None if dPSI_col is None else dPSI_col
            self.sig_col = None if sig_col is None else sig_col
            self.event_id_col = None if event_id_col is None else event_id_col
            self.extra_cols = None if extra_cols is None else extra_cols

            #make sure all columns are the correct type
            #self.check_columns()
        elif splice_data_type == 'MATS':
            #check if data is provided as dictionary with each splice event file
            if not isinstance(self.splice_data, dict):
                raise TypeError('If providing data from MATS, please supply this as a dict object with each key being the name of the splicing event type (SE, MXE, A3SS, A5SS, RI) and the value being the corresponding file path')
            elif not all([x in self.splice_data.keys() for x in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']]):
                if not any([x in self.splice_data.keys() for x in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']]):
                    raise ValueError('No data found associated with typical MATS splice events (SE, MXE, A3SS, A5SS, or RI). Please check and fix keys to correspond to correct events')
                
                #identify unrecognized columns and raise warning
                unrecognized_cols = [x for x in self.splice_data.keys() if x not in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']]
                print('Warning: The following event keys were not recognized and will not be used: {}'.format(unrecognized_cols))


            self.chromosome_col = 'chr' if chromosome_col is None else chromosome_col
            self.strand_col = 'strand' if strand_col is None else strand_col

            #if region columns are provided, make sure they are provided in dict format
            self.region_start_col = {'SE': 'exonStart_0base', 'MXE':'ExonStart', 'A3SS':{'+':'longExonStart_0base', '-':'shortEE'},'A5SS':{'+':'shortEE','-':'longExonStart_0base'}, 'RI':'riExonStart_0base'}
            self.region_end_col = {'SE': 'exonEnd', 'MXE':'ExonEnd', 'A3SS':{'+':'shortES','-':'longExonEnd'},'A5SS':{'+':'longExonEnd','-':'shortES'},'RI':'riExonEnd'}
            self.gene_name_col = 'geneSymbol' if gene_name_col is None else gene_name_col
            self.dPSI_col = 'IncLevelDifference' if dPSI_col is None else dPSI_col
            self.sig_col = 'FDR' if sig_col is None else sig_col
            self.event_id_col = event_id_col
            self.extra_cols = extra_cols


        elif splice_data_type == 'SpliceSeq':
            #check if data is provided as dictionary with each splice event file
            if self.splicegraph is None:
                raise ValueError('SpliceSeq data requires splicegraph information to be provided as a dataframe. You can download the splicegraph from the TCGA SpliceSeq website (https://bioinformatics.mdanderson.org/TCGASpliceSeq/faq.jsp) or other sources')

            #make sure spliceseq data contains needed columns
            self.gene_name_col = 'symbol' if gene_name_col is None else gene_name_col
            if not all([x in self.splice_data.columns for x in ['exons', self.gene_name_col]]):
                raise ValueError('SpliceSeq data must contain columns "exons" and "{}" to add splicegraph information'.format(self.gene_name_col))

            self.chromosome_col = 'Chromosome' 
            self.strand_col = 'Strand'
            self.dPSI_col = dPSI_col
            self.sig_col = sig_col
            self.event_id_col =  event_id_col
            self.extra_cols = extra_cols
            


    def run(self, coordinate_type = 'hg38', separate_modification_types = False, PROCESSES = 1):
        self.coordinate_type = coordinate_type
        self.separate_modification_types = separate_modification_types

        #generic analysis for use with any data type
        if self.splice_data_type == 'generic':
            # get spliced ptms with differential inclusion levels
            self.splice_data, self.spliced_ptms = project.project_ptms_onto_splice_events(self.splice_data, coordinate_type = coordinate_type, chromosome_col = self.chromosome_col, strand_col = self.strand_col, region_start_col = self.region_start_col, region_end_col = self.region_end_col, gene_col = self.gene_name_col, dPSI_col = self.dPSI_col, sig_col = self.sig_col, event_id_col = self.event_id_col, extra_cols = self.extra_cols, annotate_original_df = True, separate_modification_types = separate_modification_types, PROCESSES = PROCESSES)

            # get altered flanking sequences
            if self.identify_flanking_sequences:
                self.altered_flanks = fs.get_flanking_changes_from_splice_data(self.splice_data, self.ptm_coordinates, chromosome_col = self.chromosome_col, strand_col = self.strand_col, spliced_region_start_col = self.region_start_col, spliced_region_end_col = self.region_end_col, first_flank_start_col = self.first_flank_start_col, first_flank_end_col = self.first_flank_end_col, second_flank_start_col = self.second_flank_start_col, second_flank_end_col = self.second_flank_end_col, gene_col = self.gene_name_col, dPSI_col = self.dPSI_col, sig_col = self.sig_col, event_id_col = self.event_id_col, coordinate_type = coordinate_type)
        #MATS specific analysis
        elif self.splice_data_type == 'MATS': #MATS specific analysis
            #extract splice data associated with each splice type provided
            SE_events = self.splice_data['SE'] if 'SE' in self.splice_data else None
            fiveASS_events = self.splice_data['A5SS'] if 'A5SS' in self.splice_data else None
            threeASS_events = self.splice_data['A3SS'] if 'A3SS' in self.splice_data else None
            MXE_events = self.splice_data['MXE'] if 'MXE' in self.splice_data else None
            RI_events = self.splice_data['RI'] if 'RI' in self.splice_data else None

            #get both spliced ptms and altered flanking sequences
            results = project.project_ptms_onto_MATS(ptm_coordinates = self.ptm_coordinates, SE_events = SE_events, fiveASS_events=fiveASS_events,threeASS_events=threeASS_events, MXE_events=MXE_events, RI_events=RI_events, coordinate_type = coordinate_type, separate_modification_types=separate_modification_types, identify_flanking_sequences=self.identify_flanking_sequences, PROCESSES = PROCESSES)

            #assign results to appropriate variables
            if self.identify_flanking_sequences:
                self.splice_data, self.spliced_ptms, self.altered_flanks = results
            else:
                self.splice_data, self.spliced_ptms = results
        #SpliceSeq specific analysis
        elif self.splice_data_type == "SpliceSeq":
            results = project.project_ptms_onto_SpliceSeq(self.splice_data, splicegraph = self.splicegraph, coordinate_type = coordinate_type, separate_modification_types = separate_modification_types, dPSI_col = self.dPSI_col, sig_col = self.sig_col, identify_flanking_sequences = self.identify_flanking_sequences, PROCESSES = PROCESSES)
            if self.identify_flanking_sequences:
                self.splice_data, self.spliced_ptms, self.altered_flanks = results
            else:
                self.splice_data, self.spliced_ptms = results


    #def summarize(self):
    #    num_spliced_ptms = len(self.spliced_ptms)
    #    fraction_events_spliced_ptms = len(self.splice_data)
    #    print(f"{self.spliced_ptms.drop_duplicates(subset = ['UniProtKB Accession','Residue', 'PTM Position in Canonical Isoform']}"))

    def save(self, odir):
        if isinstance(self.splice_data, dict):
            for key in self.splice_data.keys():
                self.splice_data[key].to_csv(odir + f'{key}_annotated_source_data.csv', index = False)
        else:
            self.splice_data.to_csv(odir + 'annotated_source_data.csv', index = False)

        self.spliced_ptms.to_csv(odir + 'spliced_ptms.csv', index = False)
        if self.identify_flanking_sequences:
            self.altered_flanks.to_csv(odir + 'altered_flanks.csv', index = False)

        #save additional parameters from the run
        with open(odir + 'run_parameters.txt', 'w') as f:
            f.write('Run date: {}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            f.write('Splice data type: {}\n'.format(self.splice_data_type))
            f.write('Coordinate type: {}\n'.format(self.coordinate_type))
            f.write('Separate modification types: {}\n'.format(self.separate_modification_types))

            #indicate which columns were used for identifying splice events
            f.write('Columns used for analysis:\n')
            f.write('Chromosome column: {}\n'.format(self.chromosome_col))
            f.write('Strand column: {}\n'.format(self.strand_col))
            f.write('Spliced region start column: {}\n'.format(self.region_start_col))
            f.write('Spliced end column: {}\n'.format(self.region_end_col))
            if self.identify_flanking_sequences:
                f.write('First flank start column: {}\n'.format(self.first_flank_start_col))
                f.write('First flank end column: {}\n'.format(self.first_flank_end_col))
                f.write('Second flank start column: {}\n'.format(self.second_flank_start_col))
                f.write('Second flank end column: {}\n'.format(self.second_flank_end_col))

            f.write('Gene name column: {}\n'.format(self.gene_name_col))
            f.write('dPSI column: {}\n'.format(self.dPSI_col))
            f.write('Significance column: {}\n'.format(self.sig_col))
            f.write('Event ID column: {}\n'.format(self.event_id_col))
            #write additional columns
            if self.extra_cols is not None:
                f.write('Extra columns: {}\n'.format(self.extra_cols))


class POSE_Analyze:
    def __init__(self, spliced_ptms = None, altered_flanks = None, alpha = 0.05, min_dPSI = 0, odir = None,  psp_regulatory_site_file = None, psp_ks_file = None, psp_disease_file = None, elm_interactions_file = None, elm_motifs_file = None, PTMint_file = None, PTMcode_intraprotein_file = None, PTMcode_interprotein_file = None, RegPhos_file = None):
        self.spliced_ptms = spliced_ptms
        self.altered_flanks = altered_flanks
        self.alpha = alpha
        self.min_dPSI = min_dPSI
        self.odir = odir

        #save annotation file information
        self.psp_regulatory_site_file = psp_regulatory_site_file
        self.psp_ks_file = psp_ks_file
        self.psp_disease_file = psp_disease_file
        self.elm_interactions_file = elm_interactions_file
        self.elm_motifs_file = elm_motifs_file
        self.PTMint_file = PTMint_file
        self.PTMcode_intraprotein_file = PTMcode_intraprotein_file
        self.PTMcode_interprotein_file = PTMcode_interprotein_file


    def load_from_folder(self, idir):
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
            self.combined_ptms = analyze.combine_outputs(self.spliced_ptms, self.altered_flanks)

    def update_annotation_file_details(self, psp_regulatory_site_file = None, psp_ks_file = None, psp_disease_file = None, elm_interactions = False, elm_motifs = False, PTMint = False, PTMcode_intraprotein = False, PTMcode_interprotein = False, DEPOD = False, RegPhos = False):
        """
        Given annotation file information
        """
        self.psp_regulatory_site_file = psp_regulatory_site_file if psp_regulatory_site_file is not None else self.psp_regulatory_site_file
        self.psp_ks_file = psp_ks_file if psp_ks_file is not None else self.psp_ks_file
        self.psp_disease_file = psp_disease_file if psp_disease_file is not None else self.psp_disease_file

        if elm_interactions:
            self.elm_interactions_file = elm_interactions
        if elm_motifs:
            self.elm_motifs_file = elm_motifs
        if PTMint:
            self.PTMint_file = PTMint
        if PTMcode_intraprotein:
            self.PTMcode_intraprotein_file = PTMcode_intraprotein
        if PTMcode_interprotein:
            self.PTMcode_interprotein_file = PTMcode_interprotein
        if RegPhos:
            self.RegPhos_file = RegPhos
        if DEPOD:
            self.DEPOD_file = DEPOD

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
        self.update_annotation_file_details(psp_regulatory_site_file=psp_regulatory_site_file, psp_ks_file=psp_ks_file, psp_disease_file=psp_disease_file, elm_interactions=elm_interactions, elm_motifs=elm_motifs, PTMint=PTMint, PTMcode_intraprotein=PTMcode_intraprotein, PTMcode_interprotein=PTMcode_interprotein, DEPOD=DEPOD, RegPhos=RegPhos)

        #add annotations to spliced_ptms
        if self.spliced_ptms is not None:
            print('Adding annotations to spliced PTMs results')
            self.spliced_ptms = analyze.annotate_ptms(self.spliced_ptms, psp_regulatory_site_file = self.psp_regulatory_site_file, psp_ks_file = self.psp_ks_file, psp_disease_file = self.psp_disease_file, elm_interactions = self.elm_interactions_file, elm_motifs = self.elm_motifs_file, PTMint = self.PTMint_file, PTMcode_intraprotein = self.PTMcode_intraprotein_file, PTMcode_interprotein = self.PTMcode_interprotein_file, DEPOD = self.DEPOD_file, RegPhos = self.RegPhos_file, combine_similar = combine_similar)
        #add annotations to altered_flanks
        if self.altered_flanks is not None:
            print('Adding annotations to altered flanking sequences results')
            self.altered_flanks = analyze.annotate_ptms(self.altered_flanks, psp_regulatory_site_file = self.psp_regulatory_site_file, psp_ks_file = self.psp_ks_file, psp_disease_file = self.psp_disease_file, elm_interactions = self.elm_interactions_file, elm_motifs = self.elm_motifs_file, PTMint = self.PTMint_file, PTMcode_intraprotein = self.PTMcode_intraprotein_file, PTMcode_interprotein = self.PTMcode_interprotein_file, DEPOD = self.DEPOD_file, RegPhos = self.RegPhos_file, combine_similar = combine_similar)

        #recombine dataframes
        if self.spliced_ptms is not None and self.altered_flanks is not None:
            self.combined_ptms = analyze.combine_outputs(self.spliced_ptms, self.altered_flanks)

        self.annotation_date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def show_available_annotations(self, type = 'all', figsize = (5,5)):
        if type == 'all':
            pose_plots.show_available_annotations(self.combined_ptms, figsize = figsize)
        elif type == 'inclusion':
            pose_plots.show_available_annotations(self.spliced_ptms, figsize = figsize)
        elif type == 'flanking':
            pose_plots.show_available_annotations(self.altered_flanks, figsize = figsize)
        


    def analyze_annotations(self, annotation_type = 'Function', database = 'PhosphoSitePlus'):
        pass


    
    def save(self):
        if self.spliced_ptms is not None:
            self.spliced_ptms.to_csv(self.odir + 'spliced_ptms.csv', index = False)
        if self.altered_flanks is not None:
            self.altered_flanks.to_csv(self.odir + 'altered_flanks.csv', index = False)
        if self.combined_ptms is not None:
            self.combined_ptms.to_csv(self.odir + 'combined_information.csv', index = False)




    
        
    
def check_columns(splice_data, expected_cols = [], expected_dtypes = []):
    #check to make sure columns exist in the dataframe
    if not all([x in splice_data.columns for x in expected_cols]):
        raise ValueError('Not all expected columns are present in the splice data. Please check the column names and provide the correct names for the following columns: {}'.format([x for x in expected_cols if x not in splice_data.columns]))
    
    #check to make sure columns are the correct data type
    for col, dtype in zip(expected_cols, expected_dtypes):
        if dtype is None:
            continue
        elif isinstance(dtype, list):
            if splice_data[col].dtype not in dtype:
                #try converting to the expected data type
                try:
                    splice_data[col] = splice_data[col].astype(dtype[0])
                except:
                    raise ValueError('Column {} is not the expected data type. Expected data type is one of {}, but found data type {}'.format(col, dtype, splice_data[col].dtype))
        else:
            if splice_data[col].dtype != dtype:
                #try converting to the expected data type
                splice_data[col] = splice_data[col].astype(dtype)

        

def load_POSE_run(spliced_ptms_fname, altered_flanks_fname = None):
    spliced_ptms = pd.read_csv(spliced_ptms_fname)
    if altered_flanks_fname is not None:
        altered_flanks = pd.read_csv(altered_flanks_fname, index_col = 0)
    else:
        altered_flanks = None

    return spliced_ptms, altered_flanks



                
