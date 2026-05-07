from ptm_pose import helpers, project, flanking_sequences, nease_runner, pose_config
from ptm_pose import database_interfacing as di
from ptm_pose.splicing_tools.base import GenericDataset
import pandas as pd
from tqdm import tqdm


def extract_region_from_splicegraph(splicegraph, region_id):
    """
    Given a region id and the splicegraph from SpliceSeq, extract the chromosome, strand, and start and stop locations of that exon. Start and stop are forced to be in ascending order, which is not necessarily true from the splice graph (i.e. start > stop for negative strand exons). This is done to make the region extraction consistent with the rest of the codebase.

    Parameters
    ----------
    spliceseq : pandas.DataFrame
        SpliceSeq splicegraph dataframe, with region_id as index
    region_id : str
        Region ID to extract information from, in the format of 'GeneName_ExonNumber'

    Returns
    -------
    list
        List containing the chromosome, strand (1 for forward, -1 for negative), start, and stop locations of the region
    """
    region_info = splicegraph.loc[region_id]

    #check to see how many regions correspond to id, if multiple, default to first entry
    if isinstance(region_info, pd.DataFrame):
        region_info = region_info.iloc[0]
        print(f'Warning: {region_id} has multiple entries in splicegraph. Defaulting to first entry.')
    
    strand = project.convert_strand_symbol(region_info['Strand'])
    if strand == 1:
        return [region_info['Chromosome'], strand,region_info['Chr_Start'], region_info['Chr_Stop']]
    else:
        return [region_info['Chromosome'], strand,region_info['Chr_Stop'], region_info['Chr_Start']]

    
def get_spliceseq_event_regions(gene_name, from_exon, spliced_exons, to_exon, splicegraph):
    """
    Given all exons associated with a splicegraph event, obtain the coordinates associated with the flanking exons and the spliced region. The spliced region is defined as the exons that are associated with psi values, while flanking regions include the "from" and "to" exons that indicate the adjacent, unspliced exons.

    Parameters
    ----------
    gene_name : str
        Gene name associated with the splice event
    from_exon : int
        Exon number associated with the first flanking exon
    spliced_exons : str
        Exon numbers associated with the spliced region, separated by colons for each unique exon
    to_exon : int
        Exon number associated with the second flanking exon
    splicegraph : pandas.DataFrame
        DataFrame containing information about individual exons and their coordinates

    Returns
    -------
    tuple
        Tuple containing the genomic coordinates of the first flanking region, spliced regions, and second flanking region
    """
    first_exon_region = extract_region_from_splicegraph(splicegraph, region_id = gene_name+'_'+str(from_exon))
    spliced_regions = [extract_region_from_splicegraph(splicegraph, gene_name+'_'+exon) if '.' in exon else extract_region_from_splicegraph(splicegraph, gene_name+'_'+exon+'.0') for exon in spliced_exons.split(':')]
    second_exon_region = extract_region_from_splicegraph(splicegraph, region_id = gene_name+'_'+str(to_exon))
    return first_exon_region, spliced_regions, second_exon_region


def get_spliceseq_flank_loc(ptm, strand, from_region_coords, to_region_coords, coordinate_type = 'hg19'):
    """
    Given ptm information for identifying flanking sequences from splicegraph information, extract the relative location of the ptm in the flanking region (where it is located in translation of the flanking region).

    Parameters
    ----------
    ptm : pandas.Series
        Series containing PTM information
    strand : int
        Strand associated with the splice event (1 for forward, -1 for negative)
    from_region_coords : list
        List containing the chromosome, strand, start, and stop locations of the first flanking region
    to_region_coords : list
        List containing the chromosome, strand, start, and stop locations of the second flanking region
    
    Returns
    -------
    int
        Relative location of the PTM in the flanking region
    """
    if strand == 1 and ptm['Which Flank'] == 'First':
        return ptm[f'Gene Location ({coordinate_type})']  - from_region_coords[-2]
    elif strand == 1 and ptm['Which Flank'] == 'Second':
        return ptm[f'Gene Location ({coordinate_type})']  - to_region_coords[-2]
    elif strand == -1 and ptm['Which Flank'] == 'First':
        return from_region_coords[-1] - ptm[f'Gene Location ({coordinate_type})']
    else:
        return to_region_coords[-1] - ptm[f'Gene Location ({coordinate_type})']

def get_ptms_in_splicegraph_flank(gene_name, chromosome, strand, flank_region_start, flank_region_end, coordinate_type = 'hg19', which_flank = 'First', flank_size = 5):
    """

    """
    #check for ptms in first flank region
    flank_ptms = project.find_ptms_in_region(ptm_coordinates = pose_config.ptm_coordinates, chromosome = chromosome, strand = strand, start = flank_region_start, end = flank_region_end, coordinate_type = coordinate_type, gene = gene_name)
    if not flank_ptms.empty and which_flank == 'First': #if ptms found region, grab those close enough to splice boundary to have impacted flanking sequence
        flank_ptms = flank_ptms[flank_ptms['Proximity to Region End (bp)'] < flank_size*3]
        flank_ptms['Which Flank'] = 'First'
    elif not flank_ptms.empty and which_flank == 'Second': #if ptms found region, grab those close enough to splice boundary to have impacted flanking sequence
        flank_ptms = flank_ptms[flank_ptms['Proximity to Region Start (bp)'] < flank_size*3]
        flank_ptms['Which Flank'] = 'Second'
    
    return flank_ptms





class SpliceSeq_Dataset(GenericDataset):
    def __init__(self, data, splicegraph, min_dpsi = 0, alpha = 0.05, dpsi_col = None, sig_col = None, coordinate_type = 'hg38'):
        print('Removing ME events from analysis')
        data = data.copy()
        data = data[data['splice_type'] != 'ME'].copy()

        #remove redundant columns, if any, that will be added from splicegraph
        overlapping_columns = set(data.columns).intersection({'Chromosome', 'Strand', 'Chr_Start', 'Chr_Stop'})
        if len(overlapping_columns) > 0:
            #drop columns that will be added from splicegraph
            data = data.drop(columns=overlapping_columns)

        #split exons into individual exons
        data['Individual exon'] = data['exons'].apply(lambda x: x.split(':'))
        data = data.explode('Individual exon').drop_duplicates()
        data['Individual exon'] = data['Individual exon'].astype(float)

        #add gene location information to psi data from spliceseq
        splicegraph['Region ID'] = splicegraph['Symbol'] + '_' + splicegraph['Exon'].astype(str)
        data = data.merge(splicegraph.copy(), left_on = ['symbol', 'Individual exon'], right_on = ['Symbol', 'Exon'], how = 'left')
        data = data.rename(columns = {'Chr_Start': 'event_start', 'Chr_Stop': 'event_end'})

        #initialize base class with spliceseq information
        super().__init__(splice_data=data, min_dpsi=min_dpsi, alpha=alpha, dpsi_col=dpsi_col, sig_col=sig_col, coordinate_type=coordinate_type, chromosome_col = 'Chromosome', strand_col = 'Strand', event_id_col = 'Region ID', gene_col = 'symbol', start_coordinate_system = '1-based')

        #save spliceseq specific data
        self.splicegraph = splicegraph
        self.original_data = data

    def get_flank_changes_from_splicegraph_single_event(self, event_row, extra_cols = None, flank_size = 5, coordinate_type = 'hg19'):
        region_id = event_row[self.event_id_col] if self.event_id_col is not None else None
        dPSI = event_row[self.dpsi_col] if self.dpsi_col is not None else None
        sig = event_row[self.sig_col] if self.sig_col is not None else None

        #get region info
        from_region_coords, spliced_region_coords, to_region_coords = get_spliceseq_event_regions(gene_name = event_row['symbol'], from_exon = event_row['from_exon'], spliced_exons = event_row['exons'], to_exon = event_row['to_exon'], splicegraph = self.splicegraph)
        chromosome = from_region_coords[0]
        strand = from_region_coords[1]

        from_flank_ptms = get_ptms_in_splicegraph_flank(event_row['symbol'], chromosome, strand, from_region_coords[-2], from_region_coords[-1], coordinate_type = coordinate_type, which_flank = 'First', flank_size = flank_size)
        to_flank_ptms = get_ptms_in_splicegraph_flank(event_row['symbol'], chromosome, strand, to_region_coords[-2], to_region_coords[-1], coordinate_type = coordinate_type, which_flank = 'Second', flank_size = flank_size)
        ptms_of_interest = pd.concat([from_flank_ptms, to_flank_ptms]).reset_index()


        #if any ptms found for event that could have altered flanking sequences extract sequence information
        if not ptms_of_interest.empty:
            #add additional context from splice data, if indicated
            if self.event_id_col is not None:
                ptms_of_interest['Region ID'] = region_id
                
            if self.dpsi_col is not None:
                ptms_of_interest['dPSI'] = dPSI
            
            if self.sig_col is not None:
                ptms_of_interest['Significance'] = sig
            
            if extra_cols is not None:
                for col in extra_cols:
                    ptms_of_interest[col] = event_row[col]


            region_list = [from_region_coords] + spliced_region_coords + [to_region_coords]
            seqs = di.get_region_sequences_from_list(region_list, coordinate_type = 'hg19')
            from_sequence = seqs[0]
            to_sequence = seqs[-1] 
            spliced_sequence = ''.join(seqs[1:-1]) #combine all sequences from spliced region (may be multiple exons)

            inclusion_sequence = seqs[0] + ''.join(seqs[1:-1]) + seqs[-1] #combine sequences if spliced region is included
            exclusion_sequence = seqs[0] + seqs[-1] #combine sequences if spliced region is excluded

            #initialize columns for flanking sequences
            ptms_of_interest['Inclusion Flanking Sequence'] = ''
            ptms_of_interest['Exclusion Flanking Sequence'] = ''
            for i, ptm in ptms_of_interest.iterrows():
                ptm_loc_in_flank = get_spliceseq_flank_loc(ptm, strand, from_region_coords, to_region_coords)
                #grab where ptm is located in both the inclusion and exclusion event
                inclusion_ptm_loc, exclusion_ptm_loc = flanking_sequences.get_ptm_locs_in_spliced_sequences(ptm_loc_in_flank, from_sequence, spliced_sequence, to_sequence, strand = strand, which_flank = ptm['Which Flank'], order_by = 'Translation')

                #extract expected flanking sequence based on location in sequence
                inclusion_flank = flanking_sequences.get_flanking_sequence(inclusion_ptm_loc, inclusion_sequence, ptm_residue = ptm['Residue'], flank_size = flank_size, full_flanking_seq = False)
                exclusion_flank = flanking_sequences.get_flanking_sequence(exclusion_ptm_loc, exclusion_sequence, ptm_residue = ptm['Residue'], flank_size = flank_size, full_flanking_seq = False)

                #add to dataframe
                ptms_of_interest.loc[i, 'Inclusion Flanking Sequence'] = inclusion_flank
                ptms_of_interest.loc[i, 'Exclusion Flanking Sequence'] = exclusion_flank

            #trim the expected flanking sequence
            #ptms_of_interest['Expected Flanking Sequence'] = ptms_of_interest['Expected Flanking Sequence'].apply(lambda x: x[int((len(x)-1)/2-flank_size):int((len(x)-1)/2+flank_size+1)] if x == x else np.nan)
            #find flanking sequences that have changed and only keep those
            ptms_of_interest['Matched'] = ptms_of_interest['Inclusion Flanking Sequence'] == ptms_of_interest['Exclusion Flanking Sequence']
            ptms_of_interest = ptms_of_interest[~ptms_of_interest['Matched']]
            ptms_of_interest = ptms_of_interest.drop(columns=['Matched'])
            ptms_of_interest['Stop Codon Introduced'] = (ptms_of_interest['Inclusion Flanking Sequence'].str.contains(r'\*')) | (ptms_of_interest['Exclusion Flanking Sequence'].str.contains(r'\*')) 

        
        return ptms_of_interest

    def get_flanking_changes_from_splicegraph(self, extra_cols = None, flank_size = 5, **kwargs):
        """
        Given a DataFrame containing information about splice events  obtained from SpliceSeq and the corresponding splicegraph, extract the flanking sequences of PTMs that are nearby the splice boundary (potential for flanking sequence to be altered). Coordinate information of individual exons should be found in splicegraph. You can also provide columns with specific psi or significance information. Extra cols not in these categories can be provided with extra_cols parameter.

        Parameters
        ----------
        psi_data : pandas.DataFrame
            DataFrame containing information about splice events obtained from SpliceSeq
        splicegraph : pandas.DataFrame
            DataFrame containing information about individual exons and their coordinates
        ptm_coordinates : pandas.DataFrame
            DataFrame containing PTM coordinate information for identify PTMs in the flanking regions
        dPSI_col : str, optional
            Column name indicating delta PSI value, by default None
        sig_col : str, optional
            Column name indicating significance of the event, by default None
        event_id_col : str, optional
            Column name indicating event ID, by default None
        extra_cols : list, optional
            List of column names for additional information to add to the results, by default None
        gene_col : str, optional
            Column name indicating gene symbol of spliced gene, by default 'symbol'
        flank_size : int, optional
            Number of amino acids to include flanking the PTM, by default 5
        coordinate_type : str, optional
            Coordinate system used for the regions, by default 'hg19'. Other options is hg38.
        **kwargs: additional keyword arguments
            Additional keyword arguments, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.

        Returns
        -------
        altered_flanks : pandas.DataFrame
            DataFrame containing the PTMs associated with the flanking regions that are altered, and the flanking sequences that arise depending on whether the flanking sequence is included or not
        """
        #check for any keyword arguments to use for filtering
        if kwargs:
            filter_arguments = helpers.extract_filter_arguments(**kwargs)
            #check any excess unused keyword arguments, report them
            helpers.check_filter_kwargs(**filter_arguments)
            #filter ptm coordinates file to include only ptms with desired evidence
            ptm_coordinates = helpers.filter_ptms_by_evidence(pose_config.ptm_coordinates, **filter_arguments)
        else:
            ptm_coordinates = pose_config.ptm_coordinates.copy()

        #load spliceseq
        self.splicegraph.index = self.splicegraph['Region ID'].values

        data_for_flanks = self.original_data.drop_duplicates().copy() 

        #extract relevant columns
        relevant_columns = ['as_id', 'splice_type', 'symbol', 'from_exon', 'exons', 'to_exon']
        if self.event_id_col is not None:
            relevant_columns.append(self.event_id_col)
        if self.dpsi_col is not None:
            relevant_columns.append(self.dpsi_col)
        if self.sig_col is not None:
            relevant_columns.append(self.sig_col)
        if extra_cols is not None:
            relevant_columns.extend(extra_cols)
        
        #get region ids for flanking regions to grab from splicegraph
        data_for_flanks = data_for_flanks[relevant_columns].drop_duplicates()
        data_for_flanks = data_for_flanks.dropna(subset = ['from_exon', 'to_exon'])
        data_for_flanks['from_region_id'] = data_for_flanks[self.gene_col]+'_'+data_for_flanks['from_exon'].astype(str)
        data_for_flanks['to_region_id'] = data_for_flanks['symbol']+'_'+data_for_flanks['to_exon'].astype(str)

        #get coordinates for the different regions
        altered_flanks = []
        for i, row in tqdm(data_for_flanks.iterrows(), total = data_for_flanks.shape[0], desc = 'Finding flanking changes for splicegraph events'):
            single_event_altered_flanks = self.get_flank_changes_from_splicegraph_single_event(row,  extra_cols = extra_cols, flank_size = flank_size)

            altered_flanks.append(single_event_altered_flanks)

        altered_flanks = pd.concat(altered_flanks)    
        self.altered_flanks = altered_flanks

    def run_pose(self, identify_altered_flanks = True, extra_cols = None, flank_size = 5, PROCESSES = 1, **kwargs):
        #check for any keyword arguments to use for filtering
        self.project_ptms_generic(extra_cols = extra_cols, PROCESSES = PROCESSES, **kwargs)
        if identify_altered_flanks:
            self.get_flanking_changes_from_splicegraph(extra_cols = extra_cols, flank_size = flank_size,**kwargs)

    def run_nease(self):
        self.run_nease_generic()


