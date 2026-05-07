from ptm_pose import helpers, pose_config, project, flanking_sequences, nease_runner
import pandas as pd
import json



class GenericDataset:
    def __init__(self, splice_data, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', first_flank_start_col = None, first_flank_end_col = None, second_flank_start_col = None, second_flank_end_col = None, min_dpsi = 0, alpha = 0.05, dpsi_col = 'IncLevelDifference', sig_col = 'FDR', coordinate_type = 'hg38', event_id_col = None, start_coordinate_system = '1-based', gene_col = None):

        #initialize variables
        if not isinstance(splice_data, pd.DataFrame) and not isinstance(splice_data, dict):
            raise ValueError("splice_data must be a pandas DataFrame or a dictionary of pandas DataFrames")
        
        self.splice_data = splice_data
        self.chromosome_col = chromosome_col
        self.strand_col = strand_col


        #save column information, first checking that region start and end columns are present in the data
        self.region_start_col = region_start_col
        self.region_end_col = region_end_col
        self.first_flank_start_col = first_flank_start_col
        self.first_flank_end_col = first_flank_end_col
        self.second_flank_start_col = second_flank_start_col
        self.second_flank_end_col = second_flank_end_col
        self.event_id_col = event_id_col
        self.gene_col = gene_col
        self.dpsi_col = dpsi_col
        self.sig_col = sig_col

        #save other attributes
        self.min_dpsi = min_dpsi
        self.alpha = alpha
        self.coordinate_type = coordinate_type
        self.start_coordinate_system = start_coordinate_system


    def project_ptms_generic(self, extra_cols = None, separate_modification_types = False, PROCESSES = 1, **kwargs):
        """
        Project PTMs onto the spliced regions of the MATS events for each event type, saving an annotated version of the MATS data and dataframe of PTMs impacted by the splice events

        Parameters
        ----------
        extra_cols : list, optional
            List of additional column names from the MATS data to include in the output ptms dataframe (default is None)
        separate_modification_types : bool, optional
            Whether to separate different types of PTMs into different rows in the output dataframe (default is False)
        PROCESSES : int, optional
            Number of processes to use for multiprocessing (default is 1). If the number of events is small, multiprocessing will be automatically disabled to avoid overhead.
        kwargs:
            Additional keyword arguments to pass to the project_ptms_onto_splice_events function, such as filtering parameters to filter PTMs with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the project_ptms_onto_splice_events function for more options.

        Postconditions
        --------------
        self.ptms : pandas.DataFrame
            DataFrame containing information about PTMs projected onto the spliced regions of the MATS events, including the type of event, the associated PTMs, and any additional columns specified in extra_cols
        self.annotated_events : dict or pandas.DataFrame
             Contains original event information annotated with PTM information.If the input splice_data is a dictionary of dataframes, this will be a dictionary with the same keys where the values are the annotated dataframes for each event type. If the input splice_data is a single dataframe, this will be a single dataframe containing the annotated events.
        """
        #check for any keyword arguments to use for filtering PTMs prior to projection
        if kwargs:
            filter_arguments = helpers.extract_filter_kwargs(**kwargs)
            #check any excess unused keyword arguments, report them
            helpers.check_filter_kwargs(filter_arguments)
            #filter ptm coordinates file to include only ptms with desired evidence
            ptm_coordinates = helpers.filter_ptms(pose_config.ptm_coordinates.copy(), **filter_arguments)

            #save filter arguments to class
            for farg in filter_arguments.keys():
                setattr(self, farg, filter_arguments[farg])
        else:
            ptm_coordinates = pose_config.ptm_coordinates.copy()

        self.ptms = []
        if isinstance(self.splice_data, dict):
            self.annotated_events = {}
            for event_type, data in self.splice_data.items():
                #check to make sure there is enough information to do multiprocessing if that is desired
                if PROCESSES*4 > data.shape[0]:
                    processes = 1
                else:
                    processes = PROCESSES

                # Project the spliced region onto the original MATS data
                annotated, ptms = project.project_ptms_onto_splice_events(data, ptm_coordinates=ptm_coordinates, chromosome_col = self.chromosome_col, strand_col = self.strand_col, region_start_col = self.region_start_col, region_end_col = self.region_end_col, event_id_col = self.event_id_col, dPSI_col=self.dpsi_col, sig_col = self.sig_col, gene_col = self.gene_col, coordinate_type=self.coordinate_type, start_coordinate_system=self.start_coordinate_system, extra_cols = extra_cols, taskbar_label = f'{event_type} Events', separate_modification_types=separate_modification_types, PROCESSES = processes, min_dpsi = self.min_dpsi, alpha = self.alpha)
                # Save PTMs for later use
                self.annotated_events[event_type] = annotated
                self.ptms.append(ptms)
            self.ptms = pd.concat(self.ptms, axis=0)
        else:
            #check to make sure there is enough information to do multiprocessing if that is desired
            if PROCESSES*4 > self.splice_data.shape[0]:
                processes = 1
            else:
                processes = PROCESSES

            # Project the spliced region onto the original MATS data
            self.annotated_events, self.ptms = project.project_ptms_onto_splice_events(self.splice_data, ptm_coordinates=ptm_coordinates, chromosome_col = self.chromosome_col, strand_col = self.strand_col, region_start_col = self.region_start_col, region_end_col = self.region_end_col, event_id_col = self.event_id_col, dPSI_col=self.dpsi_col, sig_col = self.sig_col,min_dpsi = self.min_dpsi, alpha = self.alpha, gene_col = self.gene_col, coordinate_type=self.coordinate_type, start_coordinate_system=self.start_coordinate_system, extra_cols = extra_cols, taskbar_label = f'Projecting onto events', separate_modification_types=separate_modification_types, PROCESSES = processes)

    def get_altered_flanks_generic(self, extra_cols = None, **kwargs):
        """
        Identify changes to flanking sequences around PTMs resulting from splicing events, saving a dataframe of the altered flanking sequences and their associated PTMs. This function requires that the flanking sequence information be included in the original MATS data (either as upstream/downstream or first/second flanking regions)

        Parameters
        ----------
        extra_cols : list, optional
            List of additional column names from the MATS data to include in the output dataframe (default is None)
        separate_modification_types : bool, optional
            Whether to separate different types of PTMs into different rows in the output dataframe (default is False)
        PROCESSES : int, optional
            Number of processes to use for multiprocessing (default is 1). If the number of events
            is small, multiprocessing will be automatically disabled to avoid overhead.
        kwargs:
            Additional keyword arguments to pass to the get_flanking_changes_from_splice_data function
        
        Postconditions
        --------------
        self.altered_flanks : pandas.DataFrame
            DataFrame containing information about changes to flanking sequences around PTMs resulting from splicing events, including the type of event, the associated PTMs, and any additional columns specified in extra_cols
        """
        #check for any keyword arguments to use for filtering PTMs prior to projection
        if kwargs:
            filter_arguments = helpers.extract_filter_kwargs(**kwargs)
            #check any excess unused keyword arguments, report them
            helpers.check_filter_kwargs(filter_arguments)
            #filter ptm coordinates file to include only ptms with desired evidence
            ptm_coordinates = helpers.filter_ptms(pose_config.ptm_coordinates.copy(), **filter_arguments)

            #save filter arguments to class
            for farg in filter_arguments.keys():
                setattr(self, farg, filter_arguments[farg])
        else:
            ptm_coordinates = pose_config.ptm_coordinates.copy()

        if isinstance(self.splice_data, dict):
            spliced_flanks = []
            for event_type, data in self.splice_data.items():
                event_flanks = flanking_sequences.get_flanking_changes_from_splice_data(data, ptm_coordinates = ptm_coordinates,chromosome_col = self.chromosome_col, strand_col = self.strand_col, spliced_region_start_col = self.region_start_col, spliced_region_end_col = self.region_end_col, first_flank_start_col = self.first_flank_start_col, first_flank_end_col = self.first_flank_end_col, second_flank_start_col = self.second_flank_start_col, second_flank_end_col = self.second_flank_end_col, dPSI_col=self.dpsi_col, sig_col = self.sig_col, min_dpsi = self.min_dpsi, alpha = self.alpha, gene_col = self.gene_col, event_id_col = self.event_id_col, extra_cols = extra_cols, coordinate_type=self.coordinate_type, start_coordinate_system=self.start_coordinate_system)
                event_flanks['Event Type'] = event_type
                spliced_flanks.append(event_flanks)
            self.altered_flanks = pd.concat(spliced_flanks, axis=0)
        else:
            event_flanks = flanking_sequences.get_flanking_changes_from_splice_data(data, ptm_coordinates = ptm_coordinates, chromosome_col = self.strand_col, strand_col = self.strand_col, spliced_region_start_col = self.region_start_col, spliced_region_end_col = self.region_end_col, first_flank_start_col = self.first_flank_start_col, first_flank_end_col = self.first_flank_end_col, second_flank_start_col = self.second_flank_start_col, second_flank_end_col = self.second_flank_end_col, dPSI_col=self.dpsi_col, sig_col = self.sig_col, min_dpsi = self.min_dpsi, alpha = self.alpha, gene_col = self.gene_col, event_id_col = self.event_id_col, extra_cols = extra_cols, coordinate_type=self.coordinate_type, start_coordinate_system=self.start_coordinate_system)
            self.altered_flanks = event_flanks

    def run_nease_generic(self, events_to_skip = []):
        """
        Run NEASE analysis on the spliced regions of the SpliceSeq exons, saving the results in the class for later retrieval

        Parameters
        ----------
        events_to_skip: list, optional
            List of event types to skip when running NEASE. This is only applicable if the input splice_data is a dictionary of dataframes, where each key is an event type. If the input splice_data is a single dataframe, this parameter will be ignored. Default is an empty list, meaning that NEASE will be run on all event types.
        """

        if isinstance(self.splice_data, dict):
            nease = {}
            for event_type, data in self.splice_data.items():
                if event_type not in events_to_skip:
                    # Run NEASE on the tmp DataFrame
                    nease_output = nease_runner.run_nease(data, region_start_col = 'event_start', region_end_col = 'event_end', gene_col = self.gene_col, dpsi_col = self.dpsi_col)
                    # Process the NEASE output as needed
                    nease[event_type] = nease_output

            #save in class
            self.nease = nease
        else:
            # Run NEASE on the tmp DataFrame
            nease_output = nease_runner.run_nease(self.splice_data, region_start_col = 'event_start', region_end_col = 'event_end', gene_col = self.gene_col, dpsi_col = self.dpsi_col)
            #save
            self.nease = nease_output

    def get_nease_edges(self):
        """
        Get the impacted protein-protein interactions from the NEASE analysis, combining results across event types into a single dataframe.
        """
        if not hasattr(self, 'nease'):
            self.run_nease()

        if isinstance(self.nease, dict):
            edge_list = []
            for event_type in self.nease.keys():
                # Process the NEASE output as needed
                if isinstance(self.nease[event_type].get_edges(), pd.DataFrame):

                    edges = self.nease[event_type].get_edges()
                    edges['Event Type'] = event_type
                    edge_list.append(edges)
            if len(edge_list) > 0:
                edge_df = pd.concat(edge_list, axis=0)
            else:
                edge_df = pd.DataFrame()
        else:
            edge_df = self.nease.get_edges()
        return edge_df
    
    def get_nease_domains(self):
        """
        Get the impacted protein domains from the NEASE analysis, combining results across event types into a single dataframe.
        """
        if not hasattr(self, 'nease'):
            self.run_nease()

        if isinstance(self.nease, dict):
            domain_list = []
            for event_type in self.nease.keys():
                # Process the NEASE output as needed
                if isinstance(self.nease[event_type].get_domains(), pd.DataFrame):

                    domains = self.nease[event_type].get_domains()
                    domains['Event Type'] = event_type
                    domain_list.append(domains)

            if len(domain_list) == 0:
                return pd.DataFrame()
            domain_df = pd.concat(domain_list, axis=0)
        else:
            domain_df = self.nease.get_domains()
        return domain_df
    
    def get_nease_motifs(self):
        """
        Get the impacted linear motifs from ELM from the NEASE analysis, combining results across event types into a single dataframe.
        """
        if not hasattr(self, 'nease'):
            self.run_nease()

        if isinstance(self.nease, dict):
            motif_list = []
            for event_type in self.nease.keys():
                # Process the NEASE output as needed
                if isinstance(self.nease[event_type].get_elm(), pd.DataFrame):

                    motifs = self.nease[event_type].get_elm()
                    motifs['Event Type'] = event_type
                    motif_list.append(motifs)
            if len(motif_list) == 0:
                return pd.DataFrame()
            motif_df = pd.concat(motif_list, axis=0)
        else:
            motif_df = self.nease.get_elm()
        return motif_df
    
    def save(self, odir):
        """
        Save the results of the MATS analysis to the specified output directory. This includes the projected PTMs, annotated MATS data for each event type, altered flanking sequences, and NEASE results (impacted domains, PPIs, and motifs).

        Parameters
        ----------
        odir : str
            Output directory where the results will be saved. The function will create CSV files for each type of result in this directory.

        """
        if hasattr(self, 'ptms'):
            self.ptms.to_csv(f'{odir}/spliced_ptms.csv', index=False)

        if hasattr(self, 'annotated_MATS'):
            for event_type, data in self.annotated_MATS.items():
                data.to_csv(f'{odir}/{event_type}_annotated.csv', index=False)

        if hasattr(self, 'altered_flanks'):
            self.altered_flanks.to_csv(f'{odir}/altered_flanking_sequences.csv', index=False)

        if hasattr(self, 'nease'):

            domains = self.get_nease_domains()
            domains.to_csv(f'{odir}/nease_impacted_domains.csv', index=False)

            ppi = self.get_nease_edges()
            ppi.to_csv(f'{odir}/nease_impacted_ppi.csv', index=False)

            motifs = self.get_nease_motifs()
            motifs.to_csv(f'{odir}/nease_impacted_motifs.csv', index=False)

        #save the attributes (non dataframe/list) of the class as a json file for later retrieval
        attributes = {attr: getattr(self, attr) for attr in self.__dict__.keys() if isinstance(getattr(self, attr), (str, int, float, bool))}
        with open(f'{odir}/run_attributes.json', 'w') as f:
            json.dump(attributes, f)



