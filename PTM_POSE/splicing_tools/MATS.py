from ptm_pose import helpers, project, nease_runner
import pandas as pd

def get_skipped_region(MATS_SE, include_flanking = False):
    """
    Given skipped exon data from MATS, identify the differentially spliced region and flanking sequences if desired

    Parameters
    ----------
    MATS_SE : pandas.DataFrame
        DataFrame containing skipped exon data from MATS
    include_flanking : bool, optional
        Whether to include flanking sequence columns (default is False)

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with event and flanking region information
    """
    MATS_SE['event_start'] = MATS_SE['exonStart_0base']
    MATS_SE['event_end'] = MATS_SE['exonEnd']

    if include_flanking:
        if 'upstreamES' in MATS_SE.columns:
            MATS_SE['first_flank_start'] = MATS_SE['upstreamES']
            MATS_SE['first_flank_end'] = MATS_SE['upstreamEE']
            MATS_SE['second_flank_start'] = MATS_SE['downstreamES']
            MATS_SE['second_flank_end'] = MATS_SE['downstreamEE']

        elif 'firstFlankingES' in MATS_SE.columns:
            MATS_SE['first_flank_start'] = MATS_SE['firstFlankingES']
            MATS_SE['first_flank_end'] = MATS_SE['firstFlankingEE']
            MATS_SE['second_flank_start'] = MATS_SE['secondFlankingES']
            MATS_SE['second_flank_end'] = MATS_SE['secondFlankingEE']
    return MATS_SE

def get_A3SS_region(MATS_A3SS, include_flanking = False):
    """
    Given alternative 3' splice site data from MATS, identify the differentially spliced region and flanking sequences if desired

    Parameters
    ----------
    MATS_A3SS : pandas.DataFrame
        DataFrame containing alternative 3' splice site data from MATS
    include_flanking : bool, optional
        Whether to include flanking sequence columns (default is False)

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with event and flanking region information
    """
    #set the relevent start and end regions of the spliced out region, which are different depending on the strand
    region_start = []
    region_end = []
    first_flank_start = []
    first_flank_end = []
    second_flank_start = []
    second_flank_end = []

    #iterate through events and identify spliced region, which depends on DNA strand
    for i, row in MATS_A3SS.iterrows():
        strand = row['strand'] #check strand
        if strand == '+':
            region_start.append(row['longExonStart_0base'])
            region_end.append(row['shortES'])
            if include_flanking:
                second_flank_start.append(row['flankingES'])
                second_flank_end.append(row['flankingEE'])
                first_flank_start.append(row['shortES'])
                first_flank_end.append(row['shortEE'])
        else:
            region_start.append(row['shortEE'])
            region_end.append(row['longExonEnd'])
            if include_flanking:
                second_flank_start.append(row['flankingES'])
                second_flank_end.append(row['flankingEE'])
                first_flank_start.append(row['shortES'])
                first_flank_end.append(row['shortEE'])

    #save region info
    MATS_A3SS['event_start'] = region_start
    MATS_A3SS['event_end'] = region_end
    if include_flanking:
        MATS_A3SS['first_flank_start'] = first_flank_start
        MATS_A3SS['first_flank_end'] = first_flank_end
        MATS_A3SS['second_flank_start'] = second_flank_start
        MATS_A3SS['second_flank_end'] = second_flank_end
    return MATS_A3SS

def get_A5SS_region(MATS_A5SS, include_flanking = False):
    """
    Given alternative 5' splice site data from MATS, identify the differentially spliced region and flanking sequences if desired

    Parameters
    ----------
    MATS_A5SS : pandas.DataFrame
        DataFrame containing alternative 5' splice site data from MATS
    include_flanking : bool, optional
        Whether to include flanking sequence columns (default is False)

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with event and flanking region information
    """
    #set the relevent start and end regions of the spliced out region, which are different depending on the strand
    region_start = []
    region_end = []
    first_flank_start = []
    first_flank_end = []
    second_flank_start = []
    second_flank_end = []

    for i, row in MATS_A5SS.iterrows():
        strand = row['strand']
        if strand == '+':
            region_start.append(row['shortEE'])
            region_end.append(row['longExonEnd'])
            if include_flanking:
                first_flank_start.append(row['shortES'])
                first_flank_end.append(row['shortEE'])
                second_flank_start.append(row['flankingES'])
                second_flank_end.append(row['flankingEE'])
        else:
            region_start.append(row['longExonStart_0base'])
            region_end.append(row['shortES'])
            if include_flanking:
                second_flank_start.append(row['shortES'])
                second_flank_end.append(row['shortEE'])
                first_flank_start.append(row['flankingES'])
                first_flank_end.append(row['flankingEE'])

    MATS_A5SS['event_start'] = region_start
    MATS_A5SS['event_end'] = region_end
    if include_flanking:
        MATS_A5SS['first_flank_start'] = first_flank_start
        MATS_A5SS['first_flank_end'] = first_flank_end
        MATS_A5SS['second_flank_start'] = second_flank_start
        MATS_A5SS['second_flank_end'] = second_flank_end

    return MATS_A5SS

def get_RI_region(MATS_RI, include_flanking = False):
    """
    Given retained intron data from MATS, identify the differentially spliced region and flanking sequences if desired

    Parameters
    ----------
    MATS_RI : pandas.DataFrame
        DataFrame containing retained intron splice site data from MATS
    include_flanking : bool, optional
        Whether to include flanking sequence columns (default is False)

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with event and flanking region information
    """

    MATS_RI['event_start'] = MATS_RI['upstreamES']
    MATS_RI['event_end'] = MATS_RI['downstreamEE']
    if include_flanking:
        MATS_RI['first_flank_start'] = MATS_RI['upstreamES']
        MATS_RI['first_flank_end'] = MATS_RI['upstreamEE']
        MATS_RI['second_flank_start'] = MATS_RI['downstreamES']
        MATS_RI['second_flank_end'] = MATS_RI['downstreamEE']
    return MATS_RI

def get_spliced_region(MATS, event_type, include_flanking = False):
    if event_type == 'SE':
        return get_skipped_region(MATS, include_flanking = include_flanking)
    elif event_type == 'A3SS':
        return get_A3SS_region(MATS, include_flanking = include_flanking)
    elif event_type == 'A5SS':
        return get_A5SS_region(MATS, include_flanking = include_flanking)
    elif event_type == 'RI':
        return get_RI_region(MATS, include_flanking = include_flanking)
    elif event_type == 'MXE':
        print('MXE does not currently work')
        return MATS
    else:
        raise ValueError("Invalid event type")

def process_MATS_data(MATS, event_type, min_junction_counts = None, sig_col = 'FDR', dPSI_col = 'IncLevelDifference', min_dpsi = 0, alpha = 0.05, include_flanking = True):

    #restrict to significant events if indicated
    if alpha:
        MATS = MATS[MATS[sig_col] <= alpha].copy()
    if min_dpsi:
        MATS = MATS[MATS[dPSI_col].abs() >= min_dpsi].copy()

    #filter by junction counts in experiment if provided
    if min_junction_counts is not None:
        print(f'Filtering {event_type} events based on minimum junction counts.')
        MATS = helpers.get_junction_counts(MATS, quant_type = 'MATS')
        MATS = MATS[(MATS['TJC_SAMPLE_1'] >= min_junction_counts) & (MATS['TJC_SAMPLE_2'] >= min_junction_counts)]


    if MATS['chr'].str.contains('chr').any():
        MATS['chr'] = MATS['chr'].apply(lambda x: x[3:])

    #add ID column
    MATS['AS ID'] = f"{event_type}_" + MATS.index.astype(str)

    #extract information about where splice region is
    MATS = get_spliced_region(MATS, event_type, include_flanking = include_flanking)
    return MATS

def project_on_MATS_data(processed_MATS, event_type, sig_col = 'FDR', dPSI_col = 'IncLevelDifference', coordinate_type = 'hg38', extra_cols = None, separate_modification_types = False, PROCESSES = 1, **kwargs):
    #check to make sure there is enough information to do multiprocessing if that is desired
    if PROCESSES*4 > processed_MATS.shape[0]:
        processes = 1
    else:
        processes = PROCESSES


    # Project the spliced region onto the original MATS data
    processed_MATS, ptms = project.project_ptms_onto_splice_events(processed_MATS,  chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', coordinate_type=coordinate_type, start_coordinate_system='0-based', extra_cols = extra_cols, taskbar_label = f'{event_type} Events', separate_modification_types=separate_modification_types, PROCESSES = processes, **kwargs)
    #record event type
    ptms['Event Type'] = event_type

    return processed_MATS, ptms

def run_nease_on_mats(processed_MATS, dpsi_col = 'IncLevelDifference', ):
    """
    Given MATS event data and type, run NEASE analysis, return NEASE object

    Parameters
    ----------
    MATS : pandas.DataFrame
        DataFrame containing MATS event data
    event_type : str
        Type of the event (e.g., 'SE', 'A3SS', 'A5SS', 'RI')
    dpsi_col : str, optional
        Column name for delta PSI values (default is 'IncLevelDifference')

    Returns
    -------
    nease_output : NEASE object
        The output of the NEASE analysis
    """
    # Run NEASE on the tmp DataFrame
    nease_output = nease_runner.run_nease(processed_MATS, region_start_col = 'event_start', region_end_col = 'event_end', gene_col = 'geneSymbol', dpsi_col = dpsi_col)
    return nease_output

class MATS:
    def __init__(self, SE = None, A3SS = None, A5SS = None, RI = None, MXE = None, min_dpsi = 0, alpha = 0.05, min_junction_counts = None,  dpsi_col = 'IncLevelDifference', sig_col = 'FDR', coordinate_type = 'hg38'):
        self.splice_data = {}
        for event_type, data in zip(['SE', 'A3SS', 'A5SS', 'RI', 'MXE'], [SE, A3SS, A5SS, RI, MXE]):
            if data is not None:
                #filter by significance
                data = process_MATS_data(data, event_type = event_type, min_junction_counts=min_junction_counts, min_dpsi = min_dpsi, alpha = alpha, sig_col = sig_col, dPSI_col = dpsi_col)

                self.splice_data[event_type] = data.copy()

        #save different parameters
        self.min_dpsi = min_dpsi
        self.alpha = alpha
        self.min_junction_counts = min_junction_counts
        self.dpsi_col = dpsi_col
        self.sig_col = sig_col
        self.coordinate_type = coordinate_type

    def project(self, extra_cols = None, separate_modification_types = False, PROCESSES = 1):
        self.ptms = []
        self.annotated_MATS = {}
        for event_type, data in self.splice_data.items():
            processed_MATS, ptms = project_on_MATS_data(data, event_type, sig_col=self.sig_col, dPSI_col=self.dpsi_col, coordinate_type=self.coordinate_type, extra_cols=extra_cols, separate_modification_types=separate_modification_types, PROCESSES=PROCESSES)
            self.splice_data[event_type] = processed_MATS
            # Save PTMs for later use
            self.annotated_MATS[event_type] = processed_MATS
            self.ptms.append(ptms)
        self.ptms = pd.concat(self.ptms, axis=0)
    
    def get_altered_flanks(self):
        pass

    def run_nease(self):

        nease = {}
        for event_type, data in self.splice_data.items():
            if event_type != 'MXE':

                # Run NEASE on the tmp DataFrame
                nease_output = nease_runner.run_nease(data, region_start_col = 'event_start', region_end_col = 'event_end', gene_col = 'geneSymbol', dpsi_col = self.dpsi_col)
                # Process the NEASE output as needed
                nease[event_type] = nease_output

        #save in class
        self.nease = nease

    def get_nease_edges(self):
        if not hasattr(self, 'nease'):
            self.run_nease()

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
        return edge_df
    
    def get_nease_domains(self):
        if not hasattr(self, 'nease'):
            self.run_nease()

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
        return domain_df
    
    def get_nease_motifs(self):
        if not hasattr(self, 'nease'):
            self.run_nease()

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
        return motif_df

