import numpy as np
import pandas as pd

import multiprocessing

from ptm_pose import pose_config as config
from ptm_pose import flanking_sequences as fs

from tqdm import tqdm

def find_PTMs_in_region(ptm_coordinates, chromosome, strand, start, end, gene = None, coordinate_type = 'hg38'):
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

    if coordinate_type in ['hg19','hg38']:
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
        ptms_in_region = ptms_in_region[['Source of PTM', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', loc_col, 'Modification', 'Modification Class']]
        #check if ptm is associated with the same gene (if info is provided). if not, do not add
        if gene is not None:
            for i, row in ptms_in_region.iterrows():
                if ';' in row['UniProtKB Accession']:
                    uni_ids = row['UniProtKB Accession'].split(';')
                    remove = True
                    for uni in uni_ids:
                        if gene in config.uniprot_to_genename[uni].split(' '):
                            remove = False
                            break

                    if remove:
                        ptms_in_region.drop(i)
                else:
                    if gene not in config.uniprot_to_genename[row['UniProtKB Accession']].split(' '):
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
        if strand == '+':
            return 1
        elif strand == '-':
            return -1
    else:
        return strand

def find_ptms_in_many_regions(region_data, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', gene_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None, annotate_original_df = True, coordinate_type = 'hg38', separate_modification_types = False, taskbar_label = None):
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
        region_data['chr'] = region_data['chr'].str.strip('chr')
    

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

        #project ptms onto region
        ptms_in_region = find_PTMs_in_region(ptm_coordinates, chromosome, strand, start, end, gene = gene, coordinate_type = coordinate_type)
        

        #add additional context from splice data, if indicated
        if event_id_col is not None:
            ptms_in_region['Region ID'] = row[event_id_col]
            
        if dPSI_col is not None:
            ptms_in_region['dPSI'] = row[dPSI_col]
        
        if sig_col is not None:
            ptms_in_region['Significance'] = row[sig_col]
        
        if extra_cols is not None:
            for col in extra_cols:
                ptms_in_region[col] = row[col]

        #if desired, add ptm information to the original splice event dataframe
        if annotate_original_df:
            if not ptms_in_region.empty:
            #split and separate unique modification types
                if separate_modification_types:
                    ptms_in_region['Modification Class'] = ptms_in_region['Modification Class'].str.split(';')
                    ptms_in_region = ptms_in_region.explode('Modification Class')

                ptms_info = ptms_in_region.apply(lambda x: x['UniProtKB Accession'] + '_' + x['Residue'] + str(x['PTM Position in Canonical Isoform']) + ' (' + x['Modification Class'] + ')', axis = 1)
                ptms_str = '/'.join(ptms_info.values)
                spliced_ptms_list.append(ptms_str)
                num_ptms_affected.append(ptms_in_region.shape[0])
                num_unique_ptm_sites.append(ptms_in_region.groupby(['UniProtKB Accession', 'Residue']).size().shape[0])
            else:
                spliced_ptms_list.append(np.nan)
                num_ptms_affected.append(0)
                num_unique_ptm_sites.append(0)

        spliced_ptm_info.append(ptms_in_region.copy())

    #combine all PTM information 
    spliced_ptm_info = pd.concat(spliced_ptm_info, ignore_index = True)

    #convert ptm position to float
    if spliced_ptm_info.shape[0] > 0:
        spliced_ptm_info['PTM Position in Canonical Isoform'] = spliced_ptm_info['PTM Position in Canonical Isoform'].astype(float)
            
    #add ptm info to original splice event dataframe
    if annotate_original_df:
        region_data['PTMs'] = spliced_ptms_list
        region_data['Number of PTMs Affected'] = num_ptms_affected
        region_data['Number of Unique PTM Sites by Position'] = num_unique_ptm_sites
        region_data['Event Length'] = (region_data[region_end_col] - region_data[region_start_col]).abs()
        region_data['PTM Density (PTMs/bp)'] = region_data['Number of Unique PTM Sites by Position']/(region_data[region_end_col] - region_data[region_start_col]).abs()

    return region_data, spliced_ptm_info
    
def project_ptms_onto_splice_events(splice_data, ptm_coordinates = None, annotate_original_df = True, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col = None, sig_col = None, event_id_col = None, gene_col = None, extra_cols = None, separate_modification_types = False, coordinate_type = 'hg38', taskbar_label = None, PROCESSES = 1):
    """
    Given splice event quantification data, project PTMs onto the regions impacted by the splice events. Assumes that the splice event data will have chromosome, strand, and genomic start/end positions for the regions of interest, and each row of the splice_event_data corresponds to a unique region.

    Parameters

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
    separate_modification_types: bool
        Indicate whether to store PTM sites with  multiple modification types as multiple rows. For example, if a site at K100 was both an acetylation and methylation site, these will be separated into unique rows with the same site number but different modification types. Default is True.
    taskbar_label: str
        Label to display in the tqdm progress bar. Default is None, which will automatically state "Projecting PTMs onto regions using ----- coordinates".
    PROCESSES: int
        Number of processes to use for multiprocessing. Default is 1 (single processing)

    Returns
    -------
    spliced_ptm_info: pandas.DataFrame
        Contains the PTMs identified across the different splice events
    splice_data: pandas.DataFrame
        dataframe containing the original splice data with an additional column 'PTMs' that contains the PTMs found in the region of interest, in the format of 'SiteNumber(ModificationType)'. If no PTMs are found, the value will be np.nan.
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None and config.ptm_coordinates is not None:
        ptm_coordinates = config.ptm_coordinates
    elif ptm_coordinates is None:
        raise ValueError('ptm_coordinates dataframe not provided and not found in the resource files. Please provide the ptm_coordinates dataframe with config.download_ptm_coordinates() or download the file manually. To avoid needing to download this file each time, run pose_config.download_ptm_coordinates(save = True) to save the file locally within the package directory (will take ~63MB of storage space)')

    if taskbar_label is None:
        taskbar_label = 'Projecting PTMs onto splice events using ' + coordinate_type + ' coordinates.'



    #copy
    splice_data = splice_data.copy()

    if PROCESSES == 1:
        splice_data, spliced_ptm_info = find_ptms_in_many_regions(splice_data, ptm_coordinates, chromosome_col = chromosome_col, strand_col = strand_col, region_start_col = region_start_col, region_end_col = region_end_col, dPSI_col = dPSI_col, sig_col = sig_col, event_id_col = event_id_col, gene_col = gene_col, extra_cols = extra_cols, annotate_original_df = annotate_original_df, coordinate_type = coordinate_type,taskbar_label = taskbar_label, separate_modification_types=separate_modification_types)
    elif PROCESSES > 1:
        #check num_cpus available, if greater than number of cores - 1 (to avoid freezing machine), then set to PROCESSES to 1 less than total number of cores
        num_cores = multiprocessing.cpu_count()
        if PROCESSES > num_cores - 1:
            PROCESSES = num_cores - 1

        #split dataframe into chunks equal to PROCESSES
        splice_data_split = np.array_split(splice_data, PROCESSES)
        pool = multiprocessing.Pool(PROCESSES)
        #run with multiprocessing
        results = pool.starmap(find_ptms_in_many_regions, [(splice_data_split[i], ptm_coordinates, chromosome_col, strand_col, region_start_col, region_end_col, gene_col, dPSI_col, sig_col, event_id_col, extra_cols, annotate_original_df, coordinate_type, separate_modification_types, taskbar_label) for i in range(PROCESSES)])

        splice_data = pd.concat([res[0] for res in results])
        spliced_ptm_info = pd.concat([res[1] for res in results])

        #raise ValueError('Multiprocessing not yet functional. Please set PROCESSES = 1.')

    print(f'PTMs projection successful ({spliced_ptm_info.shape[0]} identified).\n')

    return splice_data, spliced_ptm_info



def project_PTMs_onto_MATS(ptm_coordinates = None, SE_events = None, fiveASS_events = None, threeASS_events = None, RI_events = None, MXE_events = None, coordinate_type = 'hg38', identify_flanking_sequences = False, PROCESSES = 1):
    """
    Given splice quantification from the MATS algorithm, annotate with PTMs that are found in the differentially included regions.

    Parameters
    ----------
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs
    SE_events: pandas.DataFrame
        dataframe containing skipped exon event information from MATS
    fiveASS_events: pandas.DataFrame
        dataframe containing 5' alternative splice site event information from MATS
    threeASS_events: pandas.DataFrame
        dataframe containing 3' alternative splice site event information from MATS
    RI_events: pandas.DataFrame
        dataframe containing retained intron event information from MATS
    MXE_events: pandas.DataFrame
        dataframe containing mutually exclusive exon event information from MATS
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.
    identify_flanking_sequences: bool
        Indicate whether to look for altered flanking sequences from spliced events, in addition to those directly in the spliced region. Default is False. (not yet active)
    PROCESSES: int
        Number of processes to use for multiprocessing. Default is 1.
    """
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

        spliced_events['SE'], SE_ptms = project_ptms_onto_splice_events(SE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', event_id_col = 'AS ID', coordinate_type=coordinate_type, taskbar_label = "Skipped Exon events", PROCESSES = SE_processes)
        SE_ptms['Event Type'] = 'SE'
        spliced_ptms.append(SE_ptms)

        if identify_flanking_sequences:
            print('Identifying flanking sequences for skipped exon events.')
            SE_flanks = fs.get_flanking_changes_from_splice_data(SE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'exonStart_0base', spliced_region_end_col = 'exonEnd', first_flank_start_col = 'firstFlankingES', first_flank_end_col = 'firstFlankingEE', second_flank_start_col = 'secondFlankingES', second_flank_end_col = 'secondFlankingEE', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', event_id_col = 'AS ID', coordinate_type=coordinate_type)
            SE_flanks['Event Type'] = 'SE'
            spliced_flanks.append(SE_flanks)

    else:
        print('Skipped exon event data (SE_events) not provided, skipping')
    
    if fiveASS_events is not None:
        if fiveASS_events['chr'].str.contains('chr').any():
            fiveASS_events['chr'] = fiveASS_events['chr'].apply(lambda x: x[3:])

        #set the relevent start and end regions of the spliced out region, which are different depending on the strand
        region_start = []
        region_end = []
        first_flank_start = []
        first_flank_end = []
        second_flank_end = []
        second_flank_start = []
        for i, row in fiveASS_events.iterrows():
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

        fiveASS_events['event_start'] = region_start
        fiveASS_events['event_end'] = region_end
        if identify_flanking_sequences:
            fiveASS_events['first_flank_start'] = first_flank_start
            fiveASS_events['first_flank_end'] = first_flank_end
            fiveASS_events['second_flank_start'] = second_flank_start
            fiveASS_events['second_flank_end'] = second_flank_end
        
        #set specific as id
        if 'AS type' in fiveASS_events.columns:
            fiveASS_events['AS ID'] =  "5ASS_" + fiveASS_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > fiveASS_events.shape[0]:
            fiveASS_processes = 1
        else:
            fiveASS_processes = PROCESSES

        #identify PTMs found within spliced regions
        spliced_events['5ASS'], fiveASS_ptms = project_ptms_onto_splice_events(fiveASS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, taskbar_label = "5' ASS events", PROCESSES = fiveASS_processes)
        fiveASS_ptms['Event Type'] = '5ASS'
        spliced_ptms.append(fiveASS_ptms)

        #identify ptms with altered flanking sequences
        if identify_flanking_sequences:
            print("Identifying flanking sequences for 5'ASS events.")
            fiveASS_flanks = fs.get_flanking_changes_from_splice_data(fiveASS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'event_start', spliced_region_end_col = 'event_end', first_flank_start_col = 'first_flank_start', first_flank_end_col = 'first_flank_end', second_flank_start_col = 'second_flank_start', second_flank_end_col = 'second_flank_end',dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type)
            fiveASS_flanks['Event Type'] = '5ASS'
            spliced_flanks.append(fiveASS_flanks)
    else:
        print("5' ASS event data (fiveASS_events) not provided, skipping.")
    
    if threeASS_events is not None:
        if threeASS_events['chr'].str.contains('chr').any():
            threeASS_events['chr'] = threeASS_events['chr'].apply(lambda x: x[3:])

        #set the relevent start and end regions of the spliced out region, which are different depending on the strand
        region_start = []
        region_end = []
        first_flank_start = []
        first_flank_end = []
        second_flank_end = []
        second_flank_start = []
        for i, row in threeASS_events.iterrows():
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
        threeASS_events['event_start'] = region_start
        threeASS_events['event_end'] = region_end
        if identify_flanking_sequences:
            threeASS_events['first_flank_start'] = first_flank_start
            threeASS_events['first_flank_end'] = first_flank_end
            threeASS_events['second_flank_start'] = second_flank_start
            threeASS_events['second_flank_end'] = second_flank_end

        #add event ids
        threeASS_events['AS ID'] = "3ASS_" + threeASS_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > threeASS_events.shape[0]:
            threeASS_processes = 1
        else:
            threeASS_processes = PROCESSES

        spliced_events['3ASS'], threeASS_ptms = project_ptms_onto_splice_events(threeASS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, taskbar_label = "3' ASS events", PROCESSES = threeASS_processes)
        threeASS_ptms['Event Type'] = '3ASS'
        spliced_ptms.append(threeASS_ptms)

            #identify ptms with altered flanking sequences
        if identify_flanking_sequences:
            print("Identifying flanking sequences for 3' ASS events.")
            threeASS_flanks = fs.get_flanking_changes_from_splice_data(threeASS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'event_start', spliced_region_end_col = 'event_end', first_flank_start_col = 'first_flank_start', first_flank_end_col = 'first_flank_end', second_flank_start_col = 'second_flank_start', second_flank_end_col = 'second_flank_end', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type)
            threeASS_flanks['Event Type'] = '3ASS'
            spliced_flanks.append(threeASS_flanks)


    else:
        print("3' ASS event data (threeASS_events) not provided, skipping")

    if RI_events is not None:
        if RI_events['chr'].str.contains('chr').any():
            RI_events['chr'] = RI_events['chr'].apply(lambda x: x[3:])

        RI_events['AS ID'] = "RI_" + RI_events.index.astype(str)

        #check to make sure there is enough information to do multiprocessing if that is desired
        if PROCESSES*4 > RI_events.shape[0]:
            RI_processes = 1
        else:
            RI_processes = PROCESSES

        spliced_events['RI'], RI_ptms = project_ptms_onto_splice_events(RI_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'upstreamEE', region_end_col = 'downstreamES', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, taskbar_label = 'Retained Intron Events', PROCESSES = RI_processes)
        RI_ptms['Event Type'] = 'RI'
        spliced_ptms.append(RI_ptms)

         #identify ptms with altered flanking sequences
        if identify_flanking_sequences:
            print('Identifying flanking sequences for retained intron events.')
            RI_flanks = fs.get_flanking_changes_from_splice_data(RI_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', spliced_region_start_col = 'upstreamEE', spliced_region_end_col = 'downstreamES', first_flank_start_col = 'upstreamES', first_flank_end_col = 'upstreamEE', second_flank_start_col = 'downstreamES', second_flank_end_col = 'downstreamEE', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type)
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
        spliced_events['MXE_Exon1'], MXE_Exon1_ptms = project_ptms_onto_splice_events(MXE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = '1stExonStart_0base', region_end_col = '1stExonEnd', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, taskbar_label = 'MXE, First Exon', PROCESSES = MXE_processes)
        MXE_Exon1_ptms['Event Type'] = 'MXE (First Exon)'
        mxe_ptms.append(MXE_Exon1_ptms)

        #second mxe exon
        spliced_events['MXE_Exon2'], MXE_Exon2_ptms = project_ptms_onto_splice_events(MXE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = '2ndExonStart_0base', region_end_col = '2ndExonEnd', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, taskbar_label = 'MXE, Second Exon', PROCESSES = MXE_processes)
        MXE_Exon2_ptms['Event Type'] = 'MXE (Second Exon)'
        mxe_ptms.append(MXE_Exon2_ptms)

        #combine mxe ptms, and then drop any PTMs that were found in both MXE's
        mxe_ptms = pd.concat([MXE_Exon1_ptms, MXE_Exon2_ptms])
        mxe_ptms = mxe_ptms.drop_duplicates(subset = ['UniProtKB Accession', 'Source of PTM', 'Residue', 'PTM Position in Canonical Isoform', 'Modification', 'Modification Class', 'dPSI', 'Significance', 'Gene'], keep = False)
        mxe_ptms['dPSI'] = mxe_ptms.apply(lambda x: x['dPSI']* -1 if x['Event Type'] == 'MXE (Second Exon)' else x['dPSI'], axis = 1)

        #add mxe ptms to spliced_ptms
        spliced_ptms.append(mxe_ptms)

    spliced_ptms = pd.concat(spliced_ptms)
    if identify_flanking_sequences:
        spliced_flanks = pd.concat(spliced_flanks)
        return spliced_events, spliced_ptms, spliced_flanks
    else:
        spliced_ptms = pd.concat(spliced_ptms)
        return spliced_events, spliced_ptms