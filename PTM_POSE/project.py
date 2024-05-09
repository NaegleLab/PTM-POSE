import numpy as np
import pandas as pd

import multiprocessing

import POSE_config as config

from tqdm import tqdm

def find_PTMs_in_region(ptm_coordinates, chromosome, strand, start, end, gene = None, coordinate_type = 'hg38', consider_ragged = True):
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
    #if gene name is provided, make sure there is translator column to check gene name match
    #if gene is not None and config.translator is None:
    #    print('Gene name was provided, but a ID translator file was not, so cannot check to make sure projected PTMs belong to the correct gene. In general, this should not cause errors, but to be certain PTMs are properly projected it is recommended that ')

    #restrict to PTMs on the same chromosome and strand
    ptms_in_region = ptm_coordinates[(ptm_coordinates['Chromosome/scaffold name'] == chromosome) & (ptm_coordinates['Strand'] == strand)].copy()

    if coordinate_type in ['hg18','hg19','hg38']:
        loc_col = f'Gene Location ({coordinate_type})'
    else:
        raise ValueError('Coordinate type must be hg38, hg19, or hg18')

    #check to make sure the start value is less than the end coordinate. If it is not, treat the end coordinate as the start and the start coordinate as the end
    if start < end:
        ptms_in_region = ptms_in_region[(ptms_in_region[loc_col] >= start) & (ptms_in_region[loc_col] <= end)]
    else:
        ptms_in_region = ptms_in_region[(ptms_in_region[loc_col] <= start) & (ptms_in_region[loc_col] >= end)]
    

    #extract only PTM information from dataframe and return that and list (if not ptms, return empty dataframe)
    if not ptms_in_region.empty:
        #calculate proximity to region start and end
        ptms_in_region['Proximity to Region Start (bp)'] = (ptms_in_region[loc_col] - start).abs()
        ptms_in_region['Proximity to Region End (bp)'] = (ptms_in_region[loc_col] - end).abs()
        ptms_in_region['Proximity to Splice Boundary (bp)'] = ptms_in_region.apply(lambda x: min(x['Proximity to Region Start (bp)'], x['Proximity to Region End (bp)']), axis = 1)

        #grab uniprot id and residue
        ptms_in_region = ptms_in_region[['UniProtKB Accession', 'Source of PTM', 'Residue', 'PTM Position in Canonical Isoform', 'Modification', 'Modification Class', 'Proximity to Region Start (bp)', 'Proximity to Region End (bp)', 'Proximity to Splice Boundary (bp)']]

        #check if ptm is associated with the same gene (if info is provided). if not, do not add
        if gene is not None:
            for i, row in ptms_in_region.iterrows():
                if gene not in config.uniprot_to_gene[row['UniProtKB Accession']].split(' '):
                    ptms_in_region = ptms_in_region.drop(i)
            #make sure ptms still are present after filtering
            if ptms_in_region.empty:
                return pd.DataFrame()
            
        #add event id and dpsi column if provided
        #if event_id is not None:
        #    ptms_in_region['Region Name'] = event_id
        #if dPSI is not None: 
        #    ptms_in_region['dPSI'] = dPSI
        #if sig is not None:
        #    ptms_in_region['Significance'] = sig
        if gene is not None:
            ptms_in_region.insert(0, 'Gene', gene)

        return ptms_in_region
    else:
        return pd.DataFrame()
    
def convert_strand_symbol(strand):
    if isinstance(strand, str):
        if strand == '+':
            return 1
        elif strand == '-':
            return -1
    else:
        return strand

def find_ptms_in_many_regions(region_data, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', gene_col = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None, annotate_original_df = True, coordinate_type = 'hg38', separate_modification_types = False):
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

    Returns
    -------
    spliced_ptm_info: pandas.DataFrame
        Contains the PTMs identified across the different splice events
    splice_data: pandas.DataFrame
        dataframe containing the original splice data with an additional column 'PTMs' that contains the PTMs found in the region of interest, in the format of 'SiteNumber(ModificationType)'. If no PTMs are found, the value will be np.nan.
    """
    spliced_ptm_info = []
    spliced_ptms_list = []
    num_ptms_affected = []
    num_unique_ptm_sites = []

    #copy
    region_data = region_data.copy()

    #iterate through each row of the splice data and find PTMs in the region
    for index, row in tqdm(region_data.iterrows(), total = len(region_data)):
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
            
    #add ptm info to original splice event dataframe
    if annotate_original_df:
        region_data['PTMs'] = spliced_ptms_list
        region_data['Number of PTMs Affected'] = num_ptms_affected
        region_data['Number of Unique PTM Sites by Position'] = num_unique_ptm_sites
        region_data['Event Length'] = (region_data[region_end_col] - region_data[region_start_col]).abs()
        region_data['PTM Density (PTMs/bp)'] = region_data['Number of Unique PTM Sites by Position']/(region_data[region_end_col] - region_data[region_start_col]).abs()

    return region_data, spliced_ptm_info
    
def project_ptms_onto_splice_events(splice_data, ptm_coordinates, annotate_original_df = True, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col = None, sig_col = None, event_id_col = None, gene_col = None, extra_cols = None, separate_modification_types = False, coordinate_type = 'hg38', PROCESSES = 1):
    """
    Given splice event quantification data, project PTMs onto the regions impacted by the splice events. Assumes that the splice event data will have chromosome, strand, and genomic start/end positions for the regions of interest, and each row of the splice_event_data corresponds to a unique region.

    Parameters
    ----------
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs
    splice_data: pandas.DataFrame
        dataframe containing splice event information, including chromosome, strand, and genomic location of regions of interest
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
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.

    Returns
    -------
    spliced_ptm_info: pandas.DataFrame
        Contains the PTMs identified across the different splice events
    splice_data: pandas.DataFrame
        dataframe containing the original splice data with an additional column 'PTMs' that contains the PTMs found in the region of interest, in the format of 'SiteNumber(ModificationType)'. If no PTMs are found, the value will be np.nan.
    """
    #initialize lists to store spliced PTM information
    spliced_ptm_info = []
    spliced_ptms_list = []
    num_ptms_affected = []

    #copy
    splice_data = splice_data.copy()

    if PROCESSES == 1:
        splice_data, spliced_ptm_info = find_ptms_in_many_regions(splice_data, ptm_coordinates, chromosome_col = chromosome_col, strand_col = strand_col, region_start_col = region_start_col, region_end_col = region_end_col, dPSI_col = dPSI_col, sig_col = sig_col, event_id_col = event_id_col, gene_col = gene_col, extra_cols = extra_cols, annotate_original_df = annotate_original_df, coordinate_type = coordinate_type, separate_modification_types=separate_modification_types)
    elif PROCESSES > 1:
        #check num_cpus available, if greater than number of cores - 1 (to avoid freezing machine), then set to PROCESSES to 1 less than total number of cores
        num_cores = multiprocessing.cpu_count()
        if PROCESSES > num_cores - 1:
            PROCESSES = num_cores - 1

        #split dataframe into chunks equal to PROCESSES
        splice_data_split = np.array_split(splice_data, PROCESSES)
        pool = multiprocessing.Pool(PROCESSES)
        #run with multiprocessing
        results = pool.starmap(find_ptms_in_many_regions, [(splice_data_split[i], ptm_coordinates, chromosome_col, strand_col, region_start_col, region_end_col, gene_col, dPSI_col, sig_col, event_id_col, extra_cols, annotate_original_df, coordinate_type, separate_modification_types) for i in range(PROCESSES)])

        splice_data = pd.concat([res[0] for res in results])
        spliced_ptm_info = pd.concat([res[1] for res in results])

        #raise ValueError('Multiprocessing not yet functional. Please set PROCESSES = 1.')



    return splice_data, spliced_ptm_info



def project_PTMs_onto_MATS(ptm_coordinates, SE_events = None, fiveASS_events = None, threeASS_events = None, RI_events = None, MXE_events = None, coordinate_type = 'hg38', PROCESSES = 1):
    """
    Given splice quantification from the MATS algorithm, annotate with PTMs that are found in the differentially included regions.

    Parameters
    ----------
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs
    MATS_data: dict
        dictionary containing the differentially included region quantification data from the MATS algorithm. Keys should be the different types of splicing events (SE, A5SS, A3SS, RI, MXE), and the values should be the quantification data in pandas.DataFrame format.
    """
    print(f'Projecting PTMs onto MATS splice events using {coordinate_type} coordinates.')
    #reformat chromosome name format
    spliced_events = {}
    spliced_ptms = []
    if SE_events is not None:
        if SE_events['chr'].str.contains('chr').any():
            SE_events['chr'] = SE_events['chr'].apply(lambda x: x[3:]) 

        SE_events['AS ID'] = SE_events['AS type'] + "_" + SE_events.index.astype(str)
        spliced_events['SE'], SE_ptms = project_ptms_onto_splice_events(SE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, PROCESSES = PROCESSES)
        SE_ptms['Event Type'] = 'SE'
        spliced_ptms.append(SE_ptms)
    else:
        print('Skipped exon event data (SE_events) not provided, skipping')
    
    if fiveASS_events is not None:
        if fiveASS_events['chr'].str.contains('chr').any():
            fiveASS_events['chr'] = fiveASS_events['chr'].apply(lambda x: x[3:])

        #set the relevent start and end regions of the spliced out region, which are different depending on the strand
        region_start = []
        region_end = []
        for i, row in fiveASS_events.iterrows():
            strand = row['strand']
            if strand == '+':
                region_start.append(row['shortEE'])
                region_end.append(row['longExonEnd'])
            else:
                region_start.append(row['longExonStart_0base'])
                region_end.append(row['shortES'])
        fiveASS_events['event_start'] = region_start
        fiveASS_events['event_end'] = region_end

        fiveASS_events['AS ID'] = fiveASS_events['AS type'] + "_" + fiveASS_events.index.astype(str)
        spliced_events['5ASS'], fiveASS_ptms = project_ptms_onto_splice_events(fiveASS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, PROCESSES = PROCESSES)
        fiveASS_ptms['Event Type'] = '5ASS'
        spliced_ptms.append(fiveASS_ptms)
    else:
        print("5' ASS event data (fiveASS_events) not provided, skipping.")
    
    if threeASS_events is not None:
        if threeASS_events['chr'].str.contains('chr').any():
            threeASS_events['chr'] = threeASS_events['chr'].apply(lambda x: x[3:])

        #set the relevent start and end regions of the spliced out region, which are different depending on the strand
        region_start = []
        region_end = []
        for i, row in threeASS_events.iterrows():
            strand = row['strand']
            if strand == '+':
                region_start.append(row['longExonStart_0base'])
                region_end.append(row['shortES'])
            else:
                region_start.append(row['shortEE'])
                region_end.append(row['longExonEnd'])
        threeASS_events['event_start'] = region_start
        threeASS_events['event_end'] = region_end

        threeASS_events['AS ID'] = threeASS_events['AS type'] + "_" + threeASS_events.index.astype(str)
        spliced_events['3ASS'], threeASS_ptms = project_ptms_onto_splice_events(threeASS_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'event_start', region_end_col = 'event_end', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, PROCESSES = PROCESSES)
        threeASS_ptms['Event Type'] = '3ASS'
        spliced_ptms.append(threeASS_ptms)
    else:
        print("3' ASS event data (threeASS_events) not provided, skipping")

    if RI_events is not None:
        if RI_events['chr'].str.contains('chr').any():
            RI_events['chr'] = RI_events['chr'].apply(lambda x: x[3:])

        RI_events['AS ID'] = RI_events['AS type'] + "_" + RI_events.index.astype(str)
        spliced_events['RI'], RI_ptms = project_ptms_onto_splice_events(RI_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'riExonStart_0base', region_end_col = 'riExonEnd', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, PROCESSES = PROCESSES)
        RI_ptms['Event Type'] = 'RI'
        spliced_ptms.append(RI_ptms)

    if MXE_events is not None:
        if MXE_events['chr'].str.contains('chr').any():
            MXE_events['chr'] = MXE_events['chr'].apply(lambda x: x[3:])

        #add AS ID
        MXE_events['AS ID'] = MXE_events['AS type'] + "_" + MXE_events.index.astype(str)
        
        mxe_ptms = []
        #first mxe exon
        spliced_events['MXE_Exon1'], MXE_Exon1_ptms = project_ptms_onto_splice_events(MXE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = '1stExonStart_0base', region_end_col = '1stExonEnd', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, PROCESSES = PROCESSES)
        MXE_Exon1_ptms['Event Type'] = 'MXE (First Exon)'
        mxe_ptms.append(MXE_Exon1_ptms)

        #second mxe exon
        spliced_events['MXE_Exon2'], MXE_Exon2_ptms = project_ptms_onto_splice_events(MXE_events, ptm_coordinates, chromosome_col = 'chr', strand_col = 'strand', region_start_col = '2ndExonStart_0base', region_end_col = '2ndExonEnd', event_id_col = 'AS ID', dPSI_col='meanDeltaPSI', sig_col = 'FDR', gene_col = 'geneSymbol', coordinate_type=coordinate_type, PROCESSES = PROCESSES)
        MXE_Exon2_ptms['Event Type'] = 'MXE (Second Exon)'
        mxe_ptms.append(MXE_Exon2_ptms)

        #combine mxe ptms, and then drop any PTMs that were found in both MXE's
        mxe_ptms = pd.concat([MXE_Exon1_ptms, MXE_Exon2_ptms])
        mxe_ptms = mxe_ptms.drop_duplicates(subset = ['UniProtKB Accession', 'Source of PTM', 'Residue', 'PTM Position in Canonical Isoform', 'Modification', 'Modification Class', 'dPSI', 'Significance', 'Gene'], keep = False)
        mxe_ptms['dPSI'] = mxe_ptms.apply(lambda x: x['dPSI']* -1 if x['Event Type'] == 'MXE (Second Exon)' else x['dPSI'], axis = 1)

        #add mxe ptms to spliced_ptms
        spliced_ptms.append(mxe_ptms)

    spliced_ptms = pd.concat(spliced_ptms)

    return spliced_events, spliced_ptms