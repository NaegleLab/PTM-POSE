import numpy as np
import pandas as pd

import multiprocessing
import datetime

from ptm_pose import pose_config, helpers
from ptm_pose import flanking_sequences as fs

from tqdm import tqdm

def find_ptms_in_region(ptm_coordinates, chromosome, strand, start, end, gene = None, coordinate_type = 'hg38'):
    """
    Given an genomic region in either hg38 or hg19 coordinates (such as the region encoding an exon of interest), identify PTMs that are mapped to that region. If so, return the exon number. If none are found, return np.nan.
    
    Parameters
    ----------
    chromosome: str
        chromosome where region is located
    strand: int
        DNA strand for region is found on (1 for forward, -1 for reverse)
    start: int
        start position of region on the chromosome/strand (should always be less than end)
    end: int
        end position of region on the chromosome/strand (should always be greater than start)
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.
    
    Returns
    -------
    ptms_in_region: pandas.DataFrame
        dataframe containing all PTMs found in the region. If no PTMs are found, returns np.nan.
        
    """
    #restrict to PTMs on the same chromosome and strand
    ptms_in_region = ptm_coordinates[(ptm_coordinates['Chromosome/scaffold name'] == chromosome) & (ptm_coordinates['Strand'] == strand)].copy()

    if coordinate_type in ['hg18', 'hg19','hg38']:
        loc_col = f'Gene Location ({coordinate_type})'
    else:
        raise ValueError('Coordinate type must be hg38 or hg19')

    #check to make sure the start value is less than the end coordinate. If it is not, treat the end coordinate as the start and the start coordinate as the end
    if start < end:
        ptms_in_region = ptms_in_region[(ptms_in_region[loc_col] >= start) & (ptms_in_region[loc_col] <= end)]
    else:
        ptms_in_region = ptms_in_region[(ptms_in_region[loc_col] <= start) & (ptms_in_region[loc_col] >= end)]
    

    #extract only PTM information from dataframe and return that and list (if not ptms, return empty dataframe)
    if not ptms_in_region.empty:
        #grab uniprot id and residue
        ptms_in_region = ptms_in_region[['Source of PTM', 'UniProtKB Accession','Isoform ID',
       'Isoform Type', 'Residue', 'PTM Position in Isoform', loc_col, 'Modification', 'Modification Class', 'Canonical Flanking Sequence', 'Constitutive', 'MS_LIT', 'MS_CST', 'LT_LIT', 'Compendia', 'Number of Compendia']]
        #check if ptm is associated with the same gene (if info is provided). if not, do not add
        if gene is not None:
            for i, row in ptms_in_region.iterrows():
                #if ';' in row['UniProtKB Accession']:
                #    uni_ids = row['UniProtKB Accession'].split(';')
                #    remove = True
                #    for uni in uni_ids:
                #        if row['UniProtKB Accession'] in pose_config.uniprot_to_genename:
                #            if gene in pose_config.uniprot_to_genename[uni.split('-')[0]].split(' '):
                #                remove = False
                #                break

                #    if remove:
                #        ptms_in_region.drop(i)
                if row['UniProtKB Accession'] in pose_config.uniprot_to_genename:
                    if gene not in pose_config.uniprot_to_genename[row['UniProtKB Accession']].split(' '):
                        ptms_in_region = ptms_in_region.drop(i)
                else:
                    ptms_in_region = ptms_in_region.drop(i)

            #make sure ptms still are present after filtering
            if ptms_in_region.empty:
                return pd.DataFrame()
            else:
                ptms_in_region.insert(0, 'Gene', gene)
        
        #calculate proximity to region start and end
        ptms_in_region['Proximity to Region Start (bp)'] = (ptms_in_region[loc_col] - start).abs()
        ptms_in_region['Proximity to Region End (bp)'] = (ptms_in_region[loc_col] - end).abs()
        ptms_in_region['Proximity to Splice Boundary (bp)'] = ptms_in_region.apply(lambda x: min(x['Proximity to Region Start (bp)'], x['Proximity to Region End (bp)']), axis = 1)


        return ptms_in_region
    else:
        return pd.DataFrame()
    
def convert_strand_symbol(strand):
    """
    Given DNA strand information, make sure the strand information is in integer format (1 for forward, -1 for reverse). This is intended to convert from string format ('+' or '-') to integer format (1 or -1), but will return the input if it is already in integer format.

    Parameters
    ----------
    strand: str or int
        DNA strand information, either as a string ('+' or '-') or an integer (1 or -1)

    Returns
    -------
    int
        DNA strand information as an integer (1 for forward, -1 for reverse)
    """
    if isinstance(strand, str):
        if strand == '+' or strand == '1':
            return 1
        elif strand == '-' or strand == '-1':
            return -1
    else:
        return strand

def find_ptms_in_many_regions(region_data, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', gene_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None, annotate_original_df = True, coordinate_type = 'hg38', start_coordinate_system = '1-based', end_coordinate_system = '1-based', separate_modification_types = False, taskbar_label = None):
    """
    Given a dataframe with a unique region in each row, project PTMs onto the regions. Assumes that the region data will have chromosome, strand, and genomic start/end positions, and each row corresponds to a unique region.

    Parameters
    ----------
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs
    region_data: pandas.DataFrame
        dataframe containing region information, including chromosome, strand, and genomic location of regions of interest
    chromosome_col: str
        column name in splice_data that contains chromosome information. Default is 'chr'. Expects it to be a str with only the chromosome number: 'Y', '1', '2', etc.
    strand_col: str
        column name in splice_data that contains strand information. Default is 'strand'. Expects it to be a str with '+' or '-', or integers as 1 or -1. Will convert to integers automatically if string format is provided.
    region_start_col: str
        column name in splice_data that contains the start position of the region of interest. Default is 'exonStart_0base'.
    region_end_col: str
        column name in splice_data that contains the end position of the region of interest. Default is 'exonEnd'.
    gene_col: str
        column name in splice_data that contains the gene name. If provided, will be used to make sure the projected PTMs stem from the same gene (some cases where genomic coordiantes overlap between distinct genes). Default is None.
    event_id_col: str
        column name in splice_data that contains the unique identifier for the splice event. If provided, will be used to annotate the ptm information with the specific splice event ID. Default is None.
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.
    separate_modification_types: bool
        Indicate whether to store PTM sites with  multiple modification types as multiple rows. For example, if a site at K100 was both an acetylation and methylation site, these will be separated into unique rows with the same site number but different modification types. Default is True.
    taskbar_label: str
        Label to display in the tqdm progress bar. Default is None, which will automatically state "Projecting PTMs onto regions using ----- coordinates".
    
    

    Returns
    -------
    spliced_ptm_info: pandas.DataFrame
        Contains the PTMs identified across the different splice events
    splice_data: pandas.DataFrame
        dataframe containing the original splice data with an additional column 'PTMs' that contains the PTMs found in the region of interest, in the format of 'SiteNumber(ModificationType)'. If no PTMs are found, the value will be np.nan.
    """
    if taskbar_label is None:
        taskbar_label = 'Projecting PTMs onto regions using ' + coordinate_type + ' coordinates.'

    if region_data[chromosome_col].str.contains('chr').any():
        region_data[chromosome_col] = region_data[chromosome_col].str.strip('chr')
    

    spliced_ptm_info = []
    spliced_ptms_list = []
    num_ptms_affected = []
    num_unique_ptm_sites = []

    #copy
    region_data = region_data.copy()

    #iterate through each row of the splice data and find PTMs in the region
    for index, row in tqdm(region_data.iterrows(), total = len(region_data), desc = taskbar_label):
        #grab region information from row
        chromosome = row[chromosome_col]
        strand = convert_strand_symbol(row[strand_col])
        start = row[region_start_col]
        end = row[region_end_col]
        #only provide these if column is given
        gene = row[gene_col] if gene_col is not None else None

        #adjust region coordinates if needed (make sure in 1-based coordinate system)
        if start_coordinate_system == '0-based':
            start += 1
        elif start_coordinate_system != '1-based':
            raise ValueError("Start coordinate system must be either '0-based' or '1-based'")
        
        if end_coordinate_system == '0-based':
            end += 1
        elif end_coordinate_system != '1-based':
            raise ValueError("End coordinate system must be either '0-based' or '1-based'")

        #project ptms onto region
        ptms_in_region = find_ptms_in_region(ptm_coordinates, chromosome, strand, start, end, gene = gene, coordinate_type = coordinate_type)
        
        extra_info = {}


        #add additional context from splice data, if indicated
        extra_info = {}
        if event_id_col is not None:
            extra_info['Region ID'] = row[event_id_col]
            
        if dPSI_col is not None:
            extra_info['dPSI'] = row[dPSI_col]
        
        if sig_col is not None:
            extra_info['Significance'] = row[sig_col]
        
        if extra_cols is not None:
            for col in extra_cols:
                extra_info[col] = row[col]

        #add extra info to ptms_in_region
        ptms_in_region = pd.concat([pd.DataFrame(extra_info, index = ptms_in_region.index), ptms_in_region], axis = 1)

        #if desired, add ptm information to the original splice event dataframe
        if annotate_original_df:
            if not ptms_in_region.empty:
            #split and separate unique modification types
                if separate_modification_types:
                    ptms_in_region['Modification Class'] = ptms_in_region['Modification Class'].str.split(';')
                    ptms_in_region = ptms_in_region.explode('Modification Class')

                ptms_info = ptms_in_region.apply(lambda x: x['UniProtKB Accession'] + '_' + x['Residue'] + str(x['PTM Position in Isoform']) + ' (' + x['Modification Class'] + ')', axis = 1)
                ptms_str = '/'.join(ptms_info.values)
                spliced_ptms_list.append(ptms_str)
                num_ptms_affected.append(ptms_in_region.shape[0])
                num_unique_ptm_sites.append(ptms_in_region.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).size().shape[0])
            else:
                spliced_ptms_list.append(np.nan)
                num_ptms_affected.append(0)
                num_unique_ptm_sites.append(0)

        spliced_ptm_info.append(ptms_in_region.copy())

    #combine all PTM information 
    spliced_ptm_info = pd.concat(spliced_ptm_info, ignore_index = True)

    #convert ptm position to float
    if spliced_ptm_info.shape[0] > 0:
        spliced_ptm_info['PTM Position in Isoform'] = spliced_ptm_info['PTM Position in Isoform'].astype(float)
            
    #add ptm info to original splice event dataframe
    if annotate_original_df:
        region_data['PTMs'] = spliced_ptms_list
        region_data['Number of PTMs Affected'] = num_ptms_affected
        region_data['Number of Unique PTM Sites by Position'] = num_unique_ptm_sites
        region_data['Event Length'] = (region_data[region_end_col] - region_data[region_start_col]).abs()
        region_data['PTM Density (PTMs/bp)'] = (region_data['Number of Unique PTM Sites by Position']*3)/region_data['Event Length'] #multiply by 3 to convert aa to bp (3 bp per codon)

    return region_data, spliced_ptm_info


    
def project_ptms_onto_splice_events(splice_data, ptm_coordinates = None, annotate_original_df = True, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col = None, sig_col = None, event_id_col = None, gene_col = None, extra_cols = None, separate_modification_types = False, coordinate_type = 'hg38', start_coordinate_system = '1-based', end_coordinate_system = '1-based', taskbar_label = None, PROCESSES = 1, **kwargs):
    """
    Given splice event quantification data, project PTMs onto the regions impacted by the splice events. Assumes that the splice event data will have chromosome, strand, and genomic start/end positions for the regions of interest, and each row of the splice_event_data corresponds to a unique region.

    
    Important note: PTM-POSE relies on Ensembl based coordinates (1-based), so if the coordinates are 0-based, make sure to indicate using the start_coordinate_system and end_coordinate_system parameters. For example, rMATS uses 0-based for the start coordinates, but 1-based for the end coordinates. In this case, set start_coordinate_system = '0-based' and end_coordinate_system = '1-based'. 

    Parameters
    ----------

    splice_data: pandas.DataFrame
        dataframe containing splice event information, including chromosome, strand, and genomic location of regions of interest
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs. If none, it will pull from the config file.
    chromosome_col: str
        column name in splice_data that contains chromosome information. Default is 'chr'. Expects it to be a str with only the chromosome number: 'Y', '1', '2', etc.
    strand_col: str
        column name in splice_data that contains strand information. Default is 'strand'. Expects it to be a str with '+' or '-', or integers as 1 or -1. Will convert to integers automatically if string format is provided.
    region_start_col: str
        column name in splice_data that contains the start position of the region of interest. Default is 'exonStart_0base'.
    region_end_col: str
        column name in splice_data that contains the end position of the region of interest. Default is 'exonEnd'.
    event_id_col: str
        column name in splice_data that contains the unique identifier for the splice event. If provided, will be used to annotate the ptm information with the specific splice event ID. Default is None.
    gene_col: str
        column name in splice_data that contains the gene name. If provided, will be used to make sure the projected PTMs stem from the same gene (some cases where genomic coordiantes overlap between distinct genes). Default is None.
    dPSI_col: str
        column name in splice_data that contains the delta PSI value for the splice event. Default is None, which will not include this information in the output
    sig_col: str
        column name in splice_data that contains the significance value for the splice event. Default is None, which will not include this information in the output.
    extra_cols: list
        list of additional columns to include in the output dataframe. Default is None, which will not include any additional columns.
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.
    start_coordinate_system: str
        indicates the coordinate system used for the start position. Either '0-based' or '1-based'. Default is '1-based'.
    end_coordinate_system: str
        indicates the coordinate system used for the end position. Either '0-based' or '1-based'. Default is '1-based'.
    separate_modification_types: bool
        Indicate whether to store PTM sites with  multiple modification types as multiple rows. For example, if a site at K100 was both an acetylation and methylation site, these will be separated into unique rows with the same site number but different modification types. Default is True.
    taskbar_label: str
        Label to display in the tqdm progress bar. Default is None, which will automatically state "Projecting PTMs onto regions using ----- coordinates".
    PROCESSES: int
        Number of processes to use for multiprocessing. Default is 1 (single processing)
    **kwargs: additional keyword arguments
        Additional keyword arguments to pass to the find_ptms_in_many_regions function, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.

    Returns
    -------
    spliced_ptm_info: pandas.DataFrame
        Contains the PTMs identified across the different splice events
    splice_data: pandas.DataFrame
        dataframe containing the original splice data with an additional column 'PTMs' that contains the PTMs found in the region of interest, in the format of 'SiteNumber(ModificationType)'. If no PTMs are found, the value will be np.nan.
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None:
        ptm_coordinates = pose_config.ptm_coordinates.copy()
    
    
    #check for any keyword arguments to use for filtering
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        #check any excess unused keyword arguments, report them
        ptm_coordinates = helpers.filter_ptms_by_evidence(ptm_coordinates, **filter_arguments)

    if taskbar_label is None:
        taskbar_label = 'Projecting PTMs onto splice events using ' + coordinate_type + ' coordinates.'

    #check for any keyword arguments to use for filtering
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        #check any excess unused keyword arguments, report them
        helpers.check_filter_kwargs(filter_arguments)
        #filter ptm coordinates file to include only ptms with desired evidence
        ptm_coordinates = helpers.filter_ptms(ptm_coordinates, **filter_arguments)

    #copy
    splice_data = splice_data.copy()

    #check columns to make sure they are present and correct data type
    check_columns(splice_data, chromosome_col=chromosome_col, strand_col=strand_col, region_start_col=region_start_col, region_end_col=region_end_col, dPSI_col=dPSI_col, sig_col=sig_col, event_id_col=event_id_col, gene_col=gene_col, extra_cols=extra_cols)

    if PROCESSES == 1:
        splice_data, spliced_ptm_info = find_ptms_in_many_regions(splice_data, ptm_coordinates, chromosome_col = chromosome_col, strand_col = strand_col, region_start_col = region_start_col, region_end_col = region_end_col, dPSI_col = dPSI_col, sig_col = sig_col, event_id_col = event_id_col, gene_col = gene_col, extra_cols = extra_cols, annotate_original_df = annotate_original_df, coordinate_type = coordinate_type,start_coordinate_system=start_coordinate_system, end_coordinate_system=end_coordinate_system, taskbar_label = taskbar_label, separate_modification_types=separate_modification_types)
    elif PROCESSES > 1:
        #check num_cpus available, if greater than number of cores - 1 (to avoid freezing machine), then set to PROCESSES to 1 less than total number of cores
        num_cores = multiprocessing.cpu_count()
        if PROCESSES > num_cores - 1:
            PROCESSES = num_cores - 1

        #split dataframe into chunks equal to PROCESSES
        splice_data_split = np.array_split(splice_data, PROCESSES)
        pool = multiprocessing.Pool(PROCESSES)
        #run with multiprocessing
        results = pool.starmap(find_ptms_in_many_regions, [(splice_data_split[i], ptm_coordinates, chromosome_col, strand_col, region_start_col, region_end_col, gene_col, dPSI_col, sig_col, event_id_col, extra_cols, annotate_original_df, coordinate_type, start_coordinate_system, end_coordinate_system, separate_modification_types, taskbar_label) for i in range(PROCESSES)])

        splice_data = pd.concat([res[0] for res in results])
        spliced_ptm_info = pd.concat([res[1] for res in results])

        #raise ValueError('Multiprocessing not yet functional. Please set PROCESSES = 1.')

    print(f'PTMs projection successful ({spliced_ptm_info.shape[0]} identified).\n')

    return splice_data, spliced_ptm_info



def project_ptms_onto_MATS(ptm_coordinates = None, SE_events = None, A5SS_events = None, A3SS_events = None, RI_events = None, MXE_events = None, coordinate_type = 'hg38', identify_flanking_sequences = False, dPSI_col = 'meanDeltaPSI', sig_col = 'FDR', extra_cols = None, separate_modification_types = False, PROCESSES = 1, **kwargs):
    """
    Given splice quantification from the MATS algorithm, annotate with PTMs that are found in the differentially included regions.

    Parameters
    ----------
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs
    SE_events: pandas.DataFrame
        dataframe containing skipped exon event information from MATS
    A5SS_events: pandas.DataFrame
        dataframe containing 5' alternative splice site event information from MATS
    A3SS_events: pandas.DataFrame
        dataframe containing 3' alternative splice site event information from MATS
    RI_events: pandas.DataFrame
        dataframe containing retained intron event information from MATS
    MXE_events: pandas.DataFrame
        dataframe containing mutually exclusive exon event information from MATS
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.
    dPSI_col: str
        Column name indicating delta PSI value. Default is 'meanDeltaPSI'.
    sig_col: str
        Column name indicating significance of the event. Default is 'FDR'.
    extra_cols: list
        List of column names for additional information to add to the results. Default is None.
    separate_modification_types: bool
        Indicate whether residues with multiple modifications (i.e. phosphorylation and acetylation) should be treated as separate PTMs and be placed in unique rows of the output dataframe. Default is False.
    PROCESSES: int
        Number of processes to use for multiprocessing. Default is 1.
    **kwargs: additional keyword arguments
        Additional keyword arguments to pass to the find_ptms_in_many_regions function, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None:
        ptm_coordinates = pose_config.ptm_coordinates.copy()
    

    #check for any keyword arguments to use for filtering
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        #check any excess unused keyword arguments, report them
        helpers.check_filter_kwargs(filter_arguments)
        #filter ptm coordinates file to include only ptms with desired evidence
        ptm_coordinates = helpers.filter_ptms(ptm_coordinates, **filter_arguments)


    print(f'Projecting PTMs onto MATS splice events using {coordinate_type} coordinates.')
    #reformat chromosome name format
    spliced_events = {}
    
    spliced_flanks = []
    spliced_ptms = []
    if SE_events is not None:   
        if SE_events['chr'].str.contains('chr').any():
            SE_events['chr'] = SE_events['chr'].apply(lambda x: x[3:]) 

        SE_events['AS ID'] = "SE_" + SE_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > SE_events.shape[0]:
            SE_processes = 1
        else:
            SE_processes = PROCESSES

        spliced_events['SE'], SE_ptms = project_ptms_onto_splice_events(SE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', event_id_col = 'AS ID', extra_cols = extra_cols, coordinate_type=coordinate_type, start_coordinate_system='0-based', taskbar_label = "Skipped Exon events", separate_modification_types=separate_modification_types, PROCESSES = SE_processes)
        SE_ptms['Event Type'] = 'SE'
        spliced_ptms.append(SE_ptms)

        if identify_flanking_sequences:
            print('Identifying flanking sequences for skipped exon events.')
            if 'upstreamES' in SE_events.columns:
                first_flank_start_col = 'upstreamES'
                first_flank_end_col = 'upstreamEE'
                second_flank_start_col = 'downstreamES'
                second_flank_end_col = 'downstreamEE'

            elif 'firstFlankingES' in SE_events.columns:
                first_flank_start_col = 'firstFlankingES'
                first_flank_end_col = 'firstFlankingEE'
                second_flank_start_col = 'secondFlankingES'
                second_flank_end_col = 'secondFlankingEE'
            else:
                raise ValueError('Could not find flanking sequence columns in skipped exon event data, based on what is typically outputted by MATS. Please check column names and provide the appropriate columns for the first and second flanking sequences')
            SE_flanks = fs.get_flanking_changes_from_splice_data(SE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'exonStart_0base', spliced_region_end_col = 'exonEnd', first_flank_start_col = first_flank_start_col, first_flank_end_col = first_flank_end_col, second_flank_start_col = second_flank_start_col, second_flank_end_col = second_flank_end_col, dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', event_id_col = 'AS ID', extra_cols = extra_cols, coordinate_type=coordinate_type, start_coordinate_system='0-based')
            SE_flanks['Event Type'] = 'SE'
            spliced_flanks.append(SE_flanks)

    else:
        print('Skipped exon event data (SE_events) not provided, skipping')
    
    if A5SS_events is not None:
        if A5SS_events['chr'].str.contains('chr').any():
            A5SS_events['chr'] = A5SS['chr'].apply(lambda x: x[3:])

        #set the relevent start and end regions of the spliced out region, which are different depending on the strand
        region_start = []
        region_end = []
        first_flank_start = []
        first_flank_end = []
        second_flank_end = []
        second_flank_start = []
        for i, row in A5SS_events.iterrows():
            strand = row['strand']
            if strand == '+':
                region_start.append(row['shortEE'])
                region_end.append(row['longExonEnd'])
                if identify_flanking_sequences:
                    first_flank_start.append(row['shortES'])
                    first_flank_end.append(row['shortEE'])
                    second_flank_start.append(row['flankingES'])
                    second_flank_end.append(row['flankingEE'])
            else:
                region_start.append(row['longExonStart_0base'])
                region_end.append(row['shortES'])
                if identify_flanking_sequences:
                    second_flank_start.append(row['shortES'])
                    second_flank_end.append(row['shortEE'])
                    first_flank_start.append(row['flankingES'])
                    first_flank_end.append(row['flankingEE'])

        A5SS_events['event_start'] = region_start
        A5SS_events['event_end'] = region_end
        if identify_flanking_sequences:
            A5SS_events['first_flank_start'] = first_flank_start
            A5SS_events['first_flank_end'] = first_flank_end
            A5SS_events['second_flank_start'] = second_flank_start
            A5SS_events['second_flank_end'] = second_flank_end
        

        #set specific as id

        A5SS_events['AS ID'] =  "5ASS_" + A5SS_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > A5SS_events.shape[0]:
            fiveASS_processes = 1
        else:
            fiveASS_processes = PROCESSES

        #identify PTMs found within spliced regions
        spliced_events['5ASS'], fiveASS_ptms = project_ptms_onto_splice_events(A5SS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', coordinate_type=coordinate_type, start_coordinate_system = '0-based', extra_cols = extra_cols, taskbar_label = "5' ASS events", separate_modification_types=separate_modification_types, PROCESSES = fiveASS_processes)
        fiveASS_ptms['Event Type'] = '5ASS'
        spliced_ptms.append(fiveASS_ptms)

        #identify ptms with altered flanking sequences
        if identify_flanking_sequences:
            print("Identifying flanking sequences for 5'ASS events.")
            fiveASS_flanks = fs.get_flanking_changes_from_splice_data(A5SS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'event_start', spliced_region_end_col = 'event_end', first_flank_start_col = 'first_flank_start', first_flank_end_col = 'first_flank_end', second_flank_start_col = 'second_flank_start', second_flank_end_col = 'second_flank_end',dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol',  event_id_col = 'AS ID', extra_cols = extra_cols, coordinate_type=coordinate_type, start_coordinate_system='0-based')
            fiveASS_flanks['Event Type'] = '5ASS'
            spliced_flanks.append(fiveASS_flanks)
    else:
        print("5' ASS event data (A5SS_events) not provided, skipping.")
    
    if A3SS_events is not None:

        if RI_events['chr'].str.contains('chr').any():
            RI_events['chr'] = RI_events['chr'].apply(lambda x: x[3:])

        if A3SS_events['chr'].str.contains('chr').any():
            A3SS_events['chr'] = A3SS_events['chr'].apply(lambda x: x[3:])

        #set the relevent start and end regions of the spliced out region, which are different depending on the strand
        region_start = []
        region_end = []
        first_flank_start = []
        first_flank_end = []
        second_flank_end = []
        second_flank_start = []
        for i, row in A3SS_events.iterrows():
            strand = row['strand']
            if strand == '+':
                region_start.append(row['longExonStart_0base'])
                region_end.append(row['shortES'])
                if identify_flanking_sequences:
                    second_flank_start.append(row['flankingES'])
                    second_flank_end.append(row['flankingEE'])
                    first_flank_start.append(row['shortES'])
                    first_flank_end.append(row['shortEE'])
            else:
                region_start.append(row['shortEE'])
                region_end.append(row['longExonEnd'])
                if identify_flanking_sequences:
                    second_flank_start.append(row['flankingES'])
                    second_flank_end.append(row['flankingEE'])
                    first_flank_start.append(row['shortES'])
                    first_flank_end.append(row['shortEE'])


        #save region info
        A3SS_events['event_start'] = region_start
        A3SS_events['event_end'] = region_end
        if identify_flanking_sequences:
            A3SS_events['first_flank_start'] = first_flank_start
            A3SS_events['first_flank_end'] = first_flank_end
            A3SS_events['second_flank_start'] = second_flank_start
            A3SS_events['second_flank_end'] = second_flank_end

        #add event ids
        A3SS_events['AS ID'] = "3ASS_" + A3SS_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > A3SS_events.shape[0]:
            threeASS_processes = 1
        else:
            threeASS_processes = PROCESSES

        spliced_events['3ASS'], threeASS_ptms = project_ptms_onto_splice_events(A3SS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', extra_cols = extra_cols, coordinate_type=coordinate_type, start_coordinate_system = '0-based', taskbar_label = "3' ASS events", separate_modification_types=separate_modification_types, PROCESSES = threeASS_processes)
        threeASS_ptms['Event Type'] = '3ASS'
        spliced_ptms.append(threeASS_ptms)

            #identify ptms with altered flanking sequences
        if identify_flanking_sequences:
            print("Identifying flanking sequences for 3' ASS events.")
            threeASS_flanks = fs.get_flanking_changes_from_splice_data(A3SS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'event_start', spliced_region_end_col = 'event_end', first_flank_start_col = 'first_flank_start', first_flank_end_col = 'first_flank_end', second_flank_start_col = 'second_flank_start', second_flank_end_col = 'second_flank_end', dPSI_col=dPSI_col, sig_col = dPSI_col, gene_col = 'geneSymbol',  event_id_col = 'AS ID', extra_cols = extra_cols, coordinate_type=coordinate_type, start_coordinate_system='0-based')
            threeASS_flanks['Event Type'] = '3ASS'
            spliced_flanks.append(threeASS_flanks)


    else:
        print("3' ASS event data (A3SS_events) not provided, skipping")

    if RI_events is not None:

        if RI_events['chr'].str.contains('chr').any():
            RI_events['chr'] = RI_events['chr'].apply(lambda x: x[3:])

        #add event id
        RI_events['AS ID'] = "RI_" + RI_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > RI_events.shape[0]:
            RI_processes = 1
        else:
            RI_processes = PROCESSES

        spliced_events['RI'], RI_ptms = project_ptms_onto_splice_events(RI_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'upstreamEE', region_end_col = 'downstreamES', event_id_col = 'AS ID', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', coordinate_type=coordinate_type, start_coordinate_system='0-based', extra_cols = extra_cols, taskbar_label = 'Retained Intron Events', separate_modification_types=separate_modification_types, PROCESSES = RI_processes)
        RI_ptms['Event Type'] = 'RI'
        spliced_ptms.append(RI_ptms)

         #identify ptms with altered flanking sequences
        if identify_flanking_sequences:
            print('Identifying flanking sequences for retained intron events.')
            RI_flanks = fs.get_flanking_changes_from_splice_data(RI_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'upstreamEE', spliced_region_end_col = 'downstreamES', first_flank_start_col = 'upstreamES', first_flank_end_col = 'upstreamEE', second_flank_start_col = 'downstreamES', second_flank_end_col = 'downstreamEE', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol',  event_id_col = 'AS ID', extra_cols = extra_cols, coordinate_type=coordinate_type, start_coordinate_system='0-based')
            RI_flanks['Event Type'] = 'RI'
            spliced_flanks.append(RI_flanks)

    if MXE_events is not None:
        if MXE_events['chr'].str.contains('chr').any():
            MXE_events['chr'] = MXE_events['chr'].apply(lambda x: x[3:])

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > MXE_events.shape[0]:
            MXE_processes = 1
        else:
            MXE_processes = PROCESSES

        #add AS ID
        MXE_events['AS ID'] = "MXE_" + MXE_events.index.astype(str)
        
        mxe_ptms = []
        #first mxe exon
        spliced_events['MXE_Exon1'], MXE_Exon1_ptms = project_ptms_onto_splice_events(MXE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = '1stExonStart_0base', region_end_col = '1stExonEnd', event_id_col = 'AS ID', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', coordinate_type=coordinate_type, start_coordinate_system = '0-based', taskbar_label = 'MXE, First Exon', extra_cols=extra_cols, separate_modification_types=separate_modification_types, PROCESSES = MXE_processes)
        MXE_Exon1_ptms['Event Type'] = 'MXE (First Exon)'
        mxe_ptms.append(MXE_Exon1_ptms)

        #second mxe exon
        spliced_events['MXE_Exon2'], MXE_Exon2_ptms = project_ptms_onto_splice_events(MXE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = '2ndExonStart_0base', region_end_col = '2ndExonEnd', event_id_col = 'AS ID', dPSI_col=dPSI_col, sig_col = sig_col, gene_col = 'geneSymbol', extra_cols=extra_cols, coordinate_type=coordinate_type, start_coordinate_system='0-based', taskbar_label = 'MXE, Second Exon', separate_modification_types=separate_modification_types, PROCESSES = MXE_processes)
        MXE_Exon2_ptms['Event Type'] = 'MXE (Second Exon)'
        mxe_ptms.append(MXE_Exon2_ptms)

        #combine mxe ptms, and then drop any PTMs that were found in both MXE's
        mxe_ptms = pd.concat([MXE_Exon1_ptms, MXE_Exon2_ptms])

        columns_to_check = ['UniProtKB Accession', 'Source of PTM', 'Residue', 'PTM Position in Isoform', 'Modification', 'Modification Class', 'Gene']
        if dPSI_col is not None:
            columns_to_check.append('dPSI')
        if sig_col is not None:
            columns_to_check.append('Significance')
        if extra_cols is not None:
            columns_to_check += extra_cols

        mxe_ptms = mxe_ptms.drop_duplicates(subset = columns_to_check, keep = False)
        #flip dPSI values for second exon
        if dPSI_col is not None:
            mxe_ptms['dPSI'] = mxe_ptms.apply(lambda x: x['dPSI']* -1 if x['Event Type'] == 'MXE (Second Exon)' else x['dPSI'], axis = 1)

        #add mxe ptms to spliced_ptms
        spliced_ptms.append(mxe_ptms)

    spliced_ptms = pd.concat(spliced_ptms)
    if identify_flanking_sequences:
        spliced_flanks = pd.concat(spliced_flanks)
        return spliced_events, spliced_ptms, spliced_flanks
    else:
        return spliced_events, spliced_ptms
    
#def project_ptms_onto_MAJIQ_dPSI(majiq_data, ptm_coordinates = None, coordinate_type = 'hg38', identify_flanking_sequences = False, dPSI_col = 'dPSI', sig_col = 'FDR', separate_modification_types = False, PROCESSES = 1):
#    print('in progress')
#    pass
    
def add_splicegraph_info(psi_data, splicegraph, purpose = 'inclusion'):
    psi_data = psi_data[psi_data['splice_type'] != 'ME'].copy()

    if purpose == 'inclusion':
        #split exons into individual exons
        psi_data['Individual exon'] = psi_data['exons'].apply(lambda x: x.split(':'))
        psi_data = psi_data.explode('Individual exon').drop_duplicates()
        psi_data['Individual exon'] = psi_data['Individual exon'].astype(float)

        #add gene location information to psi data from spliceseq
        psi_data = psi_data.merge(splicegraph, left_on = ['symbol', 'Individual exon'], right_on = ['Symbol', 'Exon'], how = 'left')
        psi_data = psi_data.rename(columns = {'Chr_Start': 'spliced_region_start', 'Chr_Stop': 'spliced_region_end'})
        return psi_data
    elif purpose == 'flanking':
        print('Not yet active. Please check back later.')
    else:
        raise ValueError('Purpose must be either inclusion or flanking. Please provide the correct purpose for the splicegraph information.')

def project_ptms_onto_SpliceSeq(psi_data, splicegraph, gene_col ='symbol', dPSI_col = None, sig_col = None, extra_cols = None, coordinate_type = 'hg19', separate_modification_types = False, identify_flanking_sequences = False, flank_size = 5, PROCESSES = 1, **kwargs):
    """
    Given splice event quantification from SpliceSeq (such as what can be downloaded from TCGASpliceSeq), annotate with PTMs that are found in the differentially included regions.

    Parameters
    ----------
    psi_data: pandas.DataFrame
        dataframe containing splice event quantification from SpliceSeq. Must contain the following columns: 'symbol', 'exons', 'splice_type'.
    splicegraph: pandas.DataFrame
        dataframe containing exon information from the splicegraph used during splice event quantification. Must contain the following columns: 'Symbol', 'Exon', 'Chr_Start', 'Chr_Stop'.
    gene_col: str
        column name in psi_data that contains the gene name. Default is 'symbol'.
    dPSI_col: str
        column name in psi_data that contains the delta PSI value for the splice event. Default is None, which will not include this information in the output.
    sig_col: str
        column name in psi_data that contains the significance value for the splice event. Default is None, which will not include this information in the output.
    extra_cols: list
        list of additional columns to include in the output dataframe. Default is None, which will not include any additional columns.
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg19'.
    separate_modification_types: bool
        Indicate whether to store PTM sites with multiple modification types as multiple rows. For example, if a site at K100 was both an acetylation and methylation site, these will be separated into unique rows with the same site number but different modification types. Default is True.
    identify_flanking_sequences: bool
        Indicate whether to identify and return the flanking sequences for the splice events. Default is False.
    flank_size: int
        Size of the flanking sequence to extract from the splice event. Default is 5, which will extract 5 bases upstream and downstream of the splice event. Only relevant if identify_flanking_sequences is True.
    PROCESSES: int
        Number of processes to use for multiprocessing. Default is 1 (single processing).
    **kwargs: additional keyword arguments
        Additional keyword arguments to pass to the find_ptms_in_many_regions function, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None:
        ptm_coordinates = pose_config.ptm_coordinates.copy()
    

    #check for any keyword arguments to use for filtering
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        #check any excess unused keyword arguments, report them
        helpers.check_filter_kwargs(filter_arguments)
        #filter ptm coordinates file to include only ptms with desired evidence
        ptm_coordinates = helpers.filter_ptms_by_evidence(ptm_coordinates, **filter_arguments)

    #remove ME events from this analysis
    print('Removing ME events from analysis')
    psi_data = psi_data[psi_data['splice_type'] != 'ME'].copy()

    #split exons into individual exons
    psi_data['Individual exon'] = psi_data['exons'].apply(lambda x: x.split(':'))
    psi_data = psi_data.explode('Individual exon').drop_duplicates()
    psi_data['Individual exon'] = psi_data['Individual exon'].astype(float)

    #add gene location information to psi data from spliceseq

    spliced_data = psi_data.merge(splicegraph, left_on = ['symbol', 'Individual exon'], right_on = ['Symbol', 'Exon'], how = 'left')
    spliced_data = spliced_data.rename(columns = {'Chr_Start': 'spliced_region_start', 'Chr_Stop': 'spliced_region_end'})

    print('Projecting PTMs onto SpliceSeq data')
    spliced_data, spliced_ptms = project_ptms_onto_splice_events(spliced_data, chromosome_col = 'Chromosome', strand_col = 'Strand', gene_col = 'symbol', region_start_col = 'spliced_region_start', region_end_col = 'spliced_region_end', event_id_col = 'as_id',dPSI_col = dPSI_col, sig_col = sig_col, extra_cols = extra_cols, separate_modification_types = separate_modification_types, coordinate_type = coordinate_type, PROCESSES = PROCESSES)

    ## add code for extracting flanking sequences (to do)
    if identify_flanking_sequences:
        altered_flanks = fs.get_flanking_changes_from_splicegraph(psi_data, splicegraph, dPSI_col = dPSI_col, sig_col = sig_col, extra_cols = extra_cols, gene_col = gene_col, coordinate_type=coordinate_type, flank_size = flank_size)

        return spliced_data, spliced_ptms, altered_flanks
    else:
        return spliced_data, spliced_ptms


#def project_ptms_onto_TCGA_SpliceSeq(tcga_cancer = 'PRAD'):
#    """
#    In progress. Will download and process TCGA SpliceSeq data for a specific cancer type, and project PTMs onto the spliced regions.
#    """
#    print('Not yet active. Please check back later.')
#    pass


def check_columns(splice_data, chromosome_col = None, strand_col = None, region_start_col = None, region_end_col = None, first_flank_start_col = None, first_flank_end_col = None, second_flank_start_col = None, second_flank_end_col = None, gene_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None):
    """
    Function to quickly check if the provided column names exist in the dataset and if they are the correct type of data
    """
    expected_cols = [chromosome_col, strand_col, region_start_col, region_end_col, first_flank_start_col, first_flank_end_col, second_flank_start_col, second_flank_end_col, gene_col, dPSI_col, sig_col, event_id_col]
    expected_dtypes = [[str, object], [str,int, object], [int,float], [int,float], [int,float], [int,float], [int,float], [int,float], [str, object], float, float, None]

    #remove cols with None and the corresponding dtype entry
    expected_dtypes = [dtype for col, dtype in zip(expected_cols, expected_dtypes) if col is not None]
    expected_cols = [col for col in expected_cols if col is not None]

    #add extra columns to the expected columns list
    if extra_cols is not None:
        expected_cols += extra_cols
        expected_dtypes += [None]*len(extra_cols)   #extra columns do not have dtype requirement


    #check to make sure columns exist in the dataframe
    if not all([x in splice_data.columns for x in expected_cols]):
        raise ValueError('Not all expected columns are present in the splice data. Please check the column names and provide the correct names for the following columns: {}'.format([x for x in expected_cols if x not in splice_data.columns]))
    
    #check to make sure columns are the correct data type
    for col, data_type in zip(expected_cols, expected_dtypes):
        if data_type is None:
            continue
        elif isinstance(data_type, list):
            if splice_data[col].dtype not in data_type:
                #try converting to the expected data type
                try:
                    splice_data[col] = splice_data[col].astype(data_type[0])
                except:
                    raise ValueError('Column {} is not the expected data type. Expected data type is one of {}, but found data type {}'.format(col, data_type, splice_data[col].dtype))
        else:
            if splice_data[col].dtype != data_type:
                #try converting to the expected data type
                try:
                    splice_data[col] = splice_data[col].astype(data_type)
                except:
                    raise ValueError('Column {} is not the expected data type. Expected data type is {}, but found data type {}'.format(col, data_type, splice_data[col].dtype))
                


            


