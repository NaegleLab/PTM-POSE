import project, annotate

import pandas as pd



class POSE:
    def __init__(self, splice_data, ptm_coordinates, splice_data_type = 'generic', chromosome_col = None, strand_col = None, region_start_col = None, region_end_col = None, gene_name_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None):
        self.splice_data = splice_data
        self.spliced_ptms = None
        self.splice_data_type = splice_data_type
        self.ptm_coordinates = ptm_coordinates

        if len(splice_data.columns) < 4:
            raise ValueError('The splice data must have at least 4 columns, consisting of columns with information on the chromosome, DNA strand, region start, and region end')

        if splice_data_type == 'generic':
            self.chromosome_col = splice_data.columns[0] if chromosome_col is None else chromosome_col
            self.strand_col = splice_data.columns[1] if strand_col is None else strand_col
            self.region_start_col = splice_data.columns[2] if region_start_col is None else region_start_col
            self.region_end_col = splice_data.columns[3] if region_end_col is None else region_end_col
            self.gene_name_col = None if gene_name_col is None else gene_name_col
            self.dPSI_col = None if dPSI_col is None else dPSI_col
            self.sig_col = None if sig_col is None else sig_col
            self.event_id_col = None if event_id_col is None else event_id_col
            self.extra_cols = None if extra_cols is None else extra_cols

            #to do: add checks for the different data types


        elif splice_data_type == 'MATS':
            #check if data is provided as dictionary with each splice event file
            if not isinstance(splice_data, dict):
                raise TypeError('If providing data from MATS, please supply this as a dict object with each key being the name of the splicing event type (SE, MXE, A3SS, A5SS, RI) and the value being the corresponding file path')
            elif not all([x in splice_data.keys() for x in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']]):
                if not any([x in splice_data.keys() for x in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']]):
                    raise ValueError('No data found associated with typical MATS splice events (SE, MXE, A3SS, A5SS, or RI). Please check and fix keys to correspond to correct events')
                
                #identify unrecognized columns and raise warning
                unrecognized_cols = [x for x in splice_data.keys() if x not in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']]
                print('Warning: The following event keys were not recognized and will not be used: {}'.format(unrecognized_cols))

            self.chromosome_col = 'chr' 
            self.strand_col = 'strand'
            self.region_start_col = {'SE': 'exonStart_0base', 'MXE':'ExonStart', 'A3SS':{'+':'longExonStart_0base', '-':'shortEE'},'A5SS':{'+':'shortEE','-':'longExonStart_0base'}, 'RI':'riExonStart_0base'}
            self.region_end_col = {'SE': 'exonEnd', 'MXE':'ExonEnd', 'A3SS':{'+':'shortES','-':'longExonEnd'},'A5SS':{'+':'longExonEnd','-':'shortES'},'RI':'riExonEnd'}
            self.gene_name_col = 'geneSymbol'
            self.dPSI_col = 'IncLevelDifference' if dPSI_col is None else dPSI_col
            self.sig_col = 'FDR' if sig_col is None else sig_col
            self.event_id_col = None if event_id_col is None else event_id_col
            self.extra_cols = None if extra_cols is None else extra_cols

    def run(self, coordinate_type = 'hg38', annotate_original_df = True, separate_modification_types = False, PROCESSES = 1):
        if self.splice_data_type == 'generic':
            self.splice_data, self.spliced_ptms = project.project_ptms_onto_splice_events(self.splice_data, coordinate_type = coordinate_type, chromosome_col = self.chromosome_col, strand_col = self.strand_col, region_start_col = self.region_start_col, region_end_col = self.region_end_col, gene_name_col = self.gene_name_col, dPSI_col = self.dPSI_col, sig_col = self.sig_col, event_id_col = self.event_id_col, extra_cols = self.extra_cols, annotate_original_df = annotate_original_df, separate_modification_types = separate_modification_types, PROCESSES = PROCESSES)

        elif self.splice_data_type == 'MATS':
            self.splice_data, self.spliced_ptms = project.project_ptms_onto_MATS(self.splice_data, coordinate_type = coordinate_type, annotate_original_df = annotate_original_df, separate_modification_types = separate_modification_types, PROCESSES = PROCESSES)
    
        
        


