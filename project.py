import numpy as np
import pandas as pd

from tqdm import tqdm

def find_PTMs_in_region(ptm_coordinates, chromosome, strand, start, end, coordinate_type = 'hg38'):
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
    ptms_in_region = ptm_coordinates[(ptm_coordinates['Chromosome/scaffold name'] == chromosome) & (ptm_coordinates['Strand'] == strand)]

    #restrict to PTMs within the region of interest
    if coordinate_type == 'hg38':
        loc_col = 'Gene Location (NC)'
    elif coordinate_type == 'hg19':
        loc_col = 'HG19 Location'
    elif coordinate_type == 'hg18':
        loc_col = 'HG18 Location'
    else:
        raise ValueError('Coordinate type must be hg38, hg19, or hg18')
    ptms_in_region = ptms_in_region[(ptms_in_region[loc_col] >= start) & (ptms_in_region[loc_col] <= end)]
    #ptms_in_region['Proximity to Boundary (bp)'] = ptms_in_region.apply(lambda x: min(abs(x[loc_col] - start), abs(x[loc_col] - end)), axis = 1)



    #extract only PTM information from dataframe and return that and list (if not ptms, return empty dataframe)
    if not ptms_in_region.empty:
        ptms_in_region['UniProtKB Accession'] = ptms_in_region['Source of PTM'].apply(lambda x: x.split('_')[0])
        ptms_in_region['Residue'] = ptms_in_region['Source of PTM'].apply(lambda x: x.split('_')[1])
        ptms_in_region = ptms_in_region[['UniProtKB Accession', 'Residue', 'Modifications']]
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
    
def project_ptms_onto_splice_events(splice_data, ptm_coordinates, annotate_original_df = True, chromosome_col = 'chr', strand_col = 'strand', region_start_col = 'exonStart_0base', region_end_col = 'exonEnd', dPSI_col = None, event_id_col = None, coordinate_type = 'hg38'):
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

    #iterate through each row of the splice data and find PTMs in the region
    for index, row in tqdm(splice_data.iterrows(), total = len(splice_data)):
        #grab region information from row
        chromosome = row[chromosome_col]
        strand = convert_strand_symbol(row[strand_col])
        start = row[region_start_col]
        end = row[region_end_col]

        #project ptms onto region
        ptms_in_region = find_PTMs_in_region(ptm_coordinates, chromosome, strand, start, end, coordinate_type = coordinate_type)

        #if an event id and deltaPSI column are provided, add that to the PTM information

        if event_id_col is not None:
            ptms_in_region['Event ID'] = row[event_id_col]

        if dPSI_col is not None:
            ptms_in_region['dPSI'] = row[dPSI_col]

        #if desired, add ptm information to the original splice event dataframe
        if annotate_original_df:
            if not ptms_in_region.empty:
                ptms_in_region['PTM Info'] = ptms_in_region.apply(lambda x: x['UniProtKB Accession'] + '_' + x['Residue'] + ' (' + x['Modifications'] + ')', axis = 1)
                ptms_str = '/'.join(ptms_in_region['PTM Info'].values)
                spliced_ptms_list.append(ptms_str)
                num_ptms_affected.append(ptms_in_region.shape[0])
            else:
                spliced_ptms_list.append(np.nan)
                num_ptms_affected.append(0)

        spliced_ptm_info.append(ptms_in_region.copy())

    #combine all PTM information 
    spliced_ptm_info = pd.concat(spliced_ptm_info, ignore_index = True)
            
    #add ptm info to original splice event dataframe
    if annotate_original_df:
        splice_data['PTMs'] = spliced_ptms_list
        splice_data['Number of PTMs Affected'] = num_ptms_affected

    return splice_data, spliced_ptm_info


def annotate_spliced_ptms(spliced_ptm_info, ptm_annotations):
    pass


def project_PTMs_onto_MATS(ptm_coordinates, MATS_data = {'SE':None, 'A5SS':None, 'A3SS':None, 'RI':None, 'MXE':None}):
    """
    Given splice quantification from the MATS algorithm, annotate with PTMs that are found in the differentially included regions.

    Parameters
    ----------
    ptm_coordinates: pandas.DataFrame
        dataframe containing PTM information, including chromosome, strand, and genomic location of PTMs
    MATS_data: dict
        dictionary containing the differentially included region quantification data from the MATS algorithm. Keys should be the different types of splicing events (SE, A5SS, A3SS, RI, MXE), and the values should be the quantification data in pandas.DataFrame format.
    """
    pass
