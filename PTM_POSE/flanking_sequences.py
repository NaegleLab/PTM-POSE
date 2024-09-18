#biopython packages
from Bio.Data import CodonTable

#standard packages
import numpy as np
import pandas as pd
import re

import tqdm
import warnings

#PTM pose functions
from ptm_pose import database_interfacing as di
from ptm_pose import project
from ptm_pose import pose_config



# Get the standard codon table
codon_table = CodonTable.unambiguous_dna_by_name["Standard"]


def translate_flanking_sequence(seq, flank_size = 7, full_flanking_seq = True, lowercase_mod = True, first_flank_length = None, stop_codon_symbol = '*', unknown_codon_symbol = 'X'):
    """
    Given a DNA sequence, translate the sequence into an amino acid sequence. If the sequence is not the correct length, the function will attempt to extract the flanking sequence with spaces to account for missing parts if full_flanking_seq is not True. If the sequence is still not the correct length, the function will raise an error. Any unrecognized codons that are found in the sequence and are not in the standard codon table, including stop codons, will be translated as 'X' (unknown) or '*' (stop codon).

    Parameters
    ----------
    seq : str
        DNA sequence to translate
    flank_size : int, optional
        Number of amino acids to include flanking the PTM, by default 7
    full_flanking_seq : bool, optional
        Whether to require the flanking sequence to be the correct length, by default True
    lowercase_mod : bool, optional
        Whether to lowercase the amino acid associated with the PTM, by default True
    first_flank_length : int, optional
        Length of the flanking sequence in front of the PTM, by default None. If full_flanking_seq is False and sequence is not the correct length, this is required.
    stop_codon_symbol : str, optional
        Symbol to use for stop codons, by default '*'
    unknown_codon_symbol : str, optional
        Symbol to use for unknown codons, by default 'X'

    Returns
    -------
    str
        Amino acid sequence of the flanking sequence if translation was successful, otherwise np.nan
    """
    aa_seq = ''
    if len(seq) == flank_size*2*3+3:
        for i in range(0, len(seq), 3):
            if seq[i:i+3] in codon_table.forward_table.keys():
                aa = codon_table.forward_table[seq[i:i+3]]
            elif seq[i:i+3] in codon_table.stop_codons:
                aa = stop_codon_symbol
            else:
                aa = unknown_codon_symbol

            if i/3 == flank_size and lowercase_mod:
                aa = aa.lower()
            aa_seq += aa
    elif len(seq) % 3 == 0 and not full_flanking_seq:
        for i in range(0, len(seq), 3):
            if seq[i:i+3] in codon_table.forward_table.keys():
                aa = codon_table.forward_table[seq[i:i+3]]
            elif seq[i:i+3] in codon_table.stop_codons:
                aa = '*'
            else:
                aa = 'X'

            if lowercase_mod and i/3 == first_flank_length:
                aa = aa.lower()
            aa_seq += aa
    elif len(seq) % 3 == 0 and full_flanking_seq:
        raise ValueError('Provided sequence length does not match indicated flank size. Fix sequence or set full_flanking_seq = False, which requires indicating the length of the flanking sequence in front of the PTM.')
    elif len(seq) % 3 != 0:
        raise ValueError('Provided sequence is not a multiple of 3')
    else:
        raise ValueError('Unknown error with flanking sequence')
    return aa_seq

def get_ptm_locs_in_spliced_sequences(ptm_loc_in_flank, first_flank_seq, spliced_seq, second_flank_seq, strand, which_flank = 'First', order_by = 'Coordinates'):
    """
    Given the location of a PTM in a flanking sequence, extract the location of the PTM in the Inclusion Flanking Sequence and the Exclusion Flanking Sequence associated with a given splice event. Inclusion Flanking Sequence will include the skipped exon region, retained intron, or longer alternative splice site depending on event type. The PTM location should be associated with where the PTM is located relative to spliced region (before = 'First', after = 'Second').

    Parameters
    ----------
    ptm_loc_in_flank : int
        Location of the PTM in the flanking sequence it is found (either first or second)
    first_flank_seq : str
        Flanking exon sequence before the spliced region
    spliced_seq : str
        Spliced region sequence
    second_flank_seq : str
        Flanking exon sequence after the spliced region
    which_flank : str, optional
        Which flank the PTM is associated with, by default 'First'
    order_by : str, optional
        Whether the first, spliced and second regions are defined by their genomic coordinates (first has smallest coordinate, spliced next, then second), or if they are defined by their translation (first the first when translated, etc.)

    Returns
    -------
    tuple
        Tuple containing the PTM location in the Inclusion Flanking Sequence and the Exclusion Flanking Sequence
    """
    if order_by == 'Translation':
        if which_flank == 'First':
            inclusion_ptm_loc, exclusion_ptm_loc = ptm_loc_in_flank, ptm_loc_in_flank
        elif which_flank == 'Second':
            inclusion_ptm_loc = ptm_loc_in_flank+len(spliced_seq)+len(first_flank_seq)
            exclusion_ptm_loc = ptm_loc_in_flank+len(first_flank_seq)

    elif order_by == 'Coordinates':
        #grab codon associated with ptm in sequence
        if (which_flank == 'First' and strand == 1) or (which_flank == 'Second' and strand == -1):
            inclusion_ptm_loc, exclusion_ptm_loc = ptm_loc_in_flank, ptm_loc_in_flank
        elif (strand == -1 and which_flank == 'First'):
            inclusion_ptm_loc =  ptm_loc_in_flank+len(spliced_seq)+len(second_flank_seq)
            exclusion_ptm_loc =  ptm_loc_in_flank+len(second_flank_seq)
        elif (strand == 1 and which_flank == 'Second'):
            inclusion_ptm_loc =  ptm_loc_in_flank+len(spliced_seq)+len(first_flank_seq)
            exclusion_ptm_loc =  ptm_loc_in_flank+len(first_flank_seq)
    else:
        raise ValueError('Unknown order_by value, must be either Coordinates (first, spliced and second regions are determined by genomic coordinates) or Translation (first, spliced and second regions are determined by translation')

    return int(inclusion_ptm_loc), int(exclusion_ptm_loc)
    

def get_flanking_sequence(ptm_loc, seq, ptm_residue, flank_size = 5, lowercase_mod = True, full_flanking_seq = False):
    """
    Given a PTM location in a sequence of DNA, extract the flanking sequence around the PTM location and translate into the amino acid sequence. If the sequence is not the correct length, the function will attempt to extract the flanking sequence with spaces to account for missing parts if full_flanking_seq is not True. If the sequence is still not the correct length, the function will raise an error. Any unrecognized codons that are found in the sequence and are not in the standard codon table, including stop codons, will be translated as 'X' (unknown) or '*' (stop codon).

    Parameters
    ----------
    ptm_loc : int
        Location of the first base pair associated with PTM in the DNA sequence
    seq : str
        DNA sequence containing the PTM
    ptm_residue : str
        Amino acid residue associated with the PTM
    flank_size : int, optional
        Number of amino acids to include flanking the PTM, by default 5
    lowercase_mod : bool, optional
        Whether to lowercase the amino acid associated with the PTM, by default True
    full_flanking_seq : bool, optional
        Whether to require the flanking sequence to be the correct length, by default False

    Returns
    -------
    str
        Amino acid sequence of the flanking sequence around the PTM if translation was successful, otherwise np.nan
    """
    ptm_codon = seq[ptm_loc:ptm_loc+3]
    #check if ptm codon codes for amino acid and then extract flanking sequence
    if ptm_codon in codon_table.forward_table.keys():
        if codon_table.forward_table[ptm_codon] == ptm_residue:
            if len(seq) != 3*(flank_size*2+1):
                if full_flanking_seq:
                    raise ValueError('Flanking sequence is not the correct length, please fix or set full_flanking_seq to False')
                else:
                    #check where issue is, at start or end of sequence
                    enough_at_start = ptm_loc >= flank_size*3
                    enough_at_end = len(seq) - ptm_loc >= flank_size*3+3
                    #extract length with amino acids and add cushion for missing parts
                    front_length = flank_size*3 if enough_at_start else ptm_loc
                    start_cushion = (flank_size*3 - ptm_loc)*' ' if not enough_at_start else ''
                    end_length = flank_size*3 + 3 if enough_at_end else len(seq) - ptm_loc
                    end_cushion = (flank_size*3 - (len(seq) - ptm_loc))*' ' if not enough_at_end else ''
                    #reconstruct sequence with spaces to account for missing ends
                    flanking_seq_bp = start_cushion +  seq[ptm_loc-front_length:ptm_loc+end_length] + end_cushion
            else:
                flanking_seq_bp = seq[ptm_loc-(flank_size*3):ptm_loc+(flank_size*3)+3]
            flanking_seq_aa = translate_flanking_sequence(flanking_seq_bp, flank_size = flank_size, lowercase_mod=lowercase_mod, full_flanking_seq = full_flanking_seq)
        else:
            flanking_seq_aa = np.nan
    else:
        flanking_seq_aa = np.nan
    
    return flanking_seq_aa

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





def get_flanking_changes(ptm_coordinates, chromosome, strand, first_flank_region, spliced_region, second_flank_region, gene = None, dPSI = None, sig = None, event_id = None, flank_size = 5, coordinate_type = 'hg38', lowercase_mod = True, order_by = 'Coordinates'):
    """
    Currently has been tested with MATS splicing events.

    Given flanking and spliced regions associated with a splice event, identify PTMs that have potential to have an altered flanking sequence depending on whether spliced region is included or excluded (if PTM is close to splice boundary). For these PTMs, extract the flanking sequences associated with the inclusion and exclusion cases and translate into amino acid sequences. If the PTM is not associated with a codon that codes for the expected amino acid, the PTM will be excluded from the results. 

    Parameters
    ----------
    ptm_coordinates : pandas.DataFrame
        DataFrame containing PTM coordinate information for identify PTMs in the flanking regions
    chromosome : str
        Chromosome associated with the splice event
    strand : int
        Strand associated with the splice event (1 for forward, -1 for negative)
    first_flank_region : list
        List containing the start and stop locations of the first flanking region (first is currently defined based on location the genome not coding sequence)
    spliced_region : list
        List containing the start and stop locations of the spliced region
    second_flank_region : list
        List containing the start and stop locations of the second flanking region (second is currently defined based on location the genome not coding sequence)
    event_id : str, optional
        Event ID associated with the splice event, by default None
    flank_size : int, optional
        Number of amino acids to include flanking the PTM, by default 7
    coordinate_type : str, optional
        Coordinate system used for the regions, by default 'hg38'. Other options is hg19.
    lowercase_mod : bool, optional
        Whether to lowercase the amino acid associated with the PTM in returned flanking sequences, by default True
    order_by : str, optional
        Whether the first, spliced and second regions are defined by their genomic coordinates (first has smallest coordinate, spliced next, then second), or if they are defined by their translation (first the first when translated, etc.)
    

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the PTMs associated with the flanking regions and the amino acid sequences of the flanking regions in the inclusion and exclusion cases
    """
    strand = project.convert_strand_symbol(strand)
    #check first flank for ptms
    ptms_in_region_first_flank = project.find_ptms_in_region(ptm_coordinates, chromosome, strand, first_flank_region[0], first_flank_region[1], gene = gene, coordinate_type = coordinate_type)
    if not ptms_in_region_first_flank.empty:
        ptms_in_region_first_flank = ptms_in_region_first_flank[ptms_in_region_first_flank['Proximity to Region End (bp)'] < flank_size*3]
        ptms_in_region_first_flank['Region'] = 'First'
    #check second flank for ptms
    ptms_in_region_second_flank = project.find_ptms_in_region(ptm_coordinates, chromosome, strand, second_flank_region[0], second_flank_region[1], gene = gene, coordinate_type = coordinate_type)
    if not ptms_in_region_second_flank.empty:
        ptms_in_region_second_flank = ptms_in_region_second_flank[ptms_in_region_second_flank['Proximity to Region Start (bp)'] < flank_size*3]
        ptms_in_region_second_flank['Region'] = 'Second'

    #combine
    ptms_in_region = pd.concat([ptms_in_region_first_flank, ptms_in_region_second_flank])


    if ptms_in_region.empty:
        return pd.DataFrame()
    else:

        #add chromosome/strand info to region info for ensembl query
        first_flank_region_query = [chromosome, strand] + first_flank_region
        spliced_region_query = [chromosome, strand] + spliced_region
        second_flank_region_query = [chromosome, strand] + second_flank_region
        regions_list = [first_flank_region_query, spliced_region_query, second_flank_region_query]

        #get dna sequences associated with regions
        first_flank_seq, spliced_seq, second_flank_seq = di.get_region_sequences_from_list(regions_list, coordinate_type = coordinate_type)

        #combine sequences for inclusion and exclusion cases
        if strand == 1:
            inclusion_seq = first_flank_seq + spliced_seq + second_flank_seq
            exclusion_seq = first_flank_seq + second_flank_seq
        else:
            inclusion_seq = second_flank_seq + spliced_seq + first_flank_seq
            exclusion_seq = second_flank_seq + first_flank_seq

        #go through all ptms in region within range of splice boundary and grab flanking sequences
        translate_success_list = []
        inclusion_seq_list = []
        exclusion_seq_list = []
        flank_region_list = []
        for i, ptm in ptms_in_region.iterrows():
            ptm_loc = ptm[f'Gene Location ({coordinate_type})']
            flank_region = ptm['Region']
            flank_region_loc = ptm['Region']
            flank_region = first_flank_region if flank_region_loc == 'First' else second_flank_region
            #grab ptm loc based on which strand ptm is on
            if strand == 1:
                relative_ptm_loc = int(ptm_loc - flank_region[0])
            else:
                relative_ptm_loc = int(flank_region[1] - ptm_loc)


            #grab where ptm is located in both the inclusion and exclusion event
            inclusion_ptm_loc, exclusion_ptm_loc = get_ptm_locs_in_spliced_sequences(relative_ptm_loc, first_flank_seq, spliced_seq, second_flank_seq,strand = strand, which_flank = flank_region_loc, order_by = order_by)

            #grab codon associated with ptm in sequence 
            ptm_codon_inclusion = inclusion_seq[inclusion_ptm_loc:inclusion_ptm_loc+3]
            ptm_codon_exclusion = exclusion_seq[exclusion_ptm_loc:exclusion_ptm_loc+3]


            #check if ptm codon codes for amino acid and then extract flanking sequence
            correct_seq = False
            if ptm_codon_inclusion in codon_table.forward_table.keys() and ptm_codon_exclusion in codon_table.forward_table.keys():
                if codon_table.forward_table[ptm_codon_inclusion] == ptm['Residue'] and codon_table.forward_table[ptm_codon_exclusion] == ptm['Residue']  and exclusion_ptm_loc-(flank_size*3) >= 0 and len(exclusion_seq) >= exclusion_ptm_loc+(flank_size*3)+3:
                    inclusion_flanking_seq = inclusion_seq[inclusion_ptm_loc-(flank_size*3):inclusion_ptm_loc+(flank_size*3)+3]
                    exclusion_flanking_seq = exclusion_seq[exclusion_ptm_loc-(flank_size*3):exclusion_ptm_loc+(flank_size*3)+3]
                    correct_seq = True


            #check to make sure ptm matches expected residue
            if correct_seq:
                translate_success_list.append(True)

                #translate flanking sequences
                inclusion_aa = translate_flanking_sequence(inclusion_flanking_seq, flank_size = flank_size, lowercase_mod=lowercase_mod)
                exclusion_aa = translate_flanking_sequence(exclusion_flanking_seq, flank_size = flank_size, lowercase_mod=lowercase_mod)

                #append to lists
                inclusion_seq_list.append(inclusion_aa)
                exclusion_seq_list.append(exclusion_aa)
                flank_region_list.append(flank_region_loc)
            else:
                translate_success_list.append(False)
                inclusion_seq_list.append(np.nan)
                exclusion_seq_list.append(np.nan)
                flank_region_list.append(flank_region_loc)

        #grab useful info from ptm dataframe
        if gene is not None:
            ptms_in_region = ptms_in_region[['Source of PTM', 'Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']].reset_index(drop = True)
        else:
            ptms_in_region = ptms_in_region[['Source of PTM', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']].reset_index(drop = True)
        #add flanking sequence information to ptm dataframe
        ptms_in_region['Inclusion Flanking Sequence'] = inclusion_seq_list
        ptms_in_region['Exclusion Flanking Sequence'] = exclusion_seq_list
        ptms_in_region['Region'] = flank_region_list
        ptms_in_region['Translation Success'] = translate_success_list

        if event_id is not None:
            ptms_in_region.insert(0, 'Event ID', event_id)
        if dPSI is not None:
            ptms_in_region['dPSI'] = dPSI
        if sig is not None:
            ptms_in_region['Significance'] = sig

        return ptms_in_region


def get_flanking_changes_from_splice_data(splice_data, ptm_coordinates = None, chromosome_col = None, strand_col = None, first_flank_start_col = None, first_flank_end_col = None, spliced_region_start_col = None, spliced_region_end_col = None, second_flank_start_col = None, second_flank_end_col = None, dPSI_col = None,  sig_col = None, event_id_col = None, gene_col = None, flank_size = 5, coordinate_type = 'hg38', lowercase_mod = True):
    """
    Given a DataFrame containing information about splice events, extract the flanking sequences associated with the PTMs in the flanking regions if there is potential for this to be altered. The DataFrame should contain columns for the chromosome, strand, start and stop locations of the first flanking region, spliced region, and second flanking region. The DataFrame should also contain a column for the event ID associated with the splice event. If the DataFrame does not contain the necessary columns, the function will raise an error.

    Parameters
    ----------
    splice_data : pandas.DataFrame
        DataFrame containing information about splice events
    ptm_coordinates : pandas.DataFrame
        DataFrame containing PTM coordinate information for identify PTMs in the flanking regions
    chromosome_col : str, optional
        Column name indicating chromosome, by default None
    strand_col : str, optional
        Column name indicating strand, by default None
    first_flank_start_col : str, optional
        Column name indicating start location of the first flanking region, by default None
    first_flank_end_col : str, optional
        Column name indicating end location of the first flanking region, by default None
    spliced_region_start_col : str, optional
        Column name indicating start location of the spliced region, by default None
    spliced_region_end_col : str, optional
        Column name indicating end location of the spliced region, by default None
    second_flank_start_col : str, optional
        Column name indicating start location of the second flanking region, by default None
    second_flank_end_col : str, optional
        Column name indicating end location of the second flanking region, by default None
    event_id_col : str, optional
        Column name indicating event ID, by default None
    flank_size : int, optional
        Number of amino acids to include flanking the PTM, by default 7
    coordinate_type : str, optional
        Coordinate system used for the regions, by default 'hg38'. Other options is hg19.
    lowercase_mod : bool, optional
        Whether to lowercase the amino acid associated with the PTM in returned flanking sequences, by default True
    
    Returns
    -------
    list
        List containing DataFrames with the PTMs associated with the flanking regions and the amino acid sequences of the flanking regions in the inclusion and exclusion cases
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None and pose_config.ptm_coordinates is not None:
        ptm_coordinates = pose_config.ptm_coordinates
    elif ptm_coordinates is None:
        raise ValueError('ptm_coordinates dataframe not provided and not found in the resource files. Please provide the ptm_coordinates dataframe with config.download_ptm_coordinates() or download the file manually. To avoid needing to download this file each time, run pose_config.download_ptm_coordinates(save = True) to save the file locally within the package directory (will take ~63MB of storage space)')

    #check to make sure all required columns are provided
    if chromosome_col is None and strand_col is None and first_flank_start_col is None and first_flank_end_col is None and spliced_region_start_col is None and spliced_region_end_col is None and second_flank_start_col is None and second_flank_end_col is None:
        raise ValueError('Please provide column names for chromosome, strand, first flank start, first flank end, spliced region start, spliced region end, second flank start, and second flank end.')

    #if chromosome is labeled with 'chr', remove
    if splice_data[chromosome_col].str.contains('chr').any():
        splice_data['chr'] = splice_data['chr'].str.strip('chr')
    

    results = []
    for i, event in tqdm.tqdm(splice_data.iterrows(), total = splice_data.shape[0], desc = 'Finding flanking sequences for PTMs nearby splice junctions'):
        if event_id_col is None:
            event_id = i
        else:
            event_id = event[event_id_col]

        #get gene info
        chromosome = event[chromosome_col]
        strand = event[strand_col]
        gene = event[gene_col] if gene_col is not None else None
        dPSI = event[dPSI_col] if dPSI_col is not None else None
        sig = event[sig_col] if sig_col is not None else None

        #extract region inof
        first_flank_region = [event[first_flank_start_col], event[first_flank_end_col]]
        spliced_region = [event[spliced_region_start_col], event[spliced_region_end_col]]
        second_flank_region = [event[second_flank_start_col], event[second_flank_end_col]]

        #get flanking changes
        ptm_flanks = get_flanking_changes(ptm_coordinates, chromosome, strand, first_flank_region, spliced_region, second_flank_region, gene = gene, sig = sig, dPSI = dPSI, event_id = event_id, flank_size = flank_size, coordinate_type = coordinate_type, lowercase_mod=lowercase_mod)

        #append to results
        results.append(ptm_flanks)

    results = pd.concat(results)
    #combine and remove any failed translation attempts
    if not results.empty:
        results = results[results['Translation Success']]

        #do some quick comparison of flanking sequences
        if not results.empty:
            #find flanking sequences that have changed and only keep those
            results['Matched'] = results['Inclusion Flanking Sequence'] == results['Exclusion Flanking Sequence']
            results = results[~results['Matched']]
            results = results.drop(columns=['Matched'])
            results['Stop Codon Introduced'] = (results['Inclusion Flanking Sequence'].str.contains(r'\*')) | (results['Exclusion Flanking Sequence'].str.contains(r'\*')) 

        print(f'{results.shape[0]} PTMs found with potential for altered flanking sequences.')
    else:
        print('No PTMs found with potential for altered flanking sequences.')
    return results


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

def get_flank_changes_from_splicegraph_single_event(event_row, splicegraph, event_id_col = None, dPSI_col = None, sig_col = None, extra_cols = None, flank_size = 5, coordinate_type = 'hg19'):
    region_id = event_row[event_id_col] if event_id_col is not None else None
    dPSI = event_row[dPSI_col] if dPSI_col is not None else None
    sig = event_row[sig_col] if sig_col is not None else None

    #get region info
    from_region_coords, spliced_region_coords, to_region_coords = get_spliceseq_event_regions(gene_name = event_row['symbol'], from_exon = event_row['from_exon'], spliced_exons = event_row['exons'], to_exon = event_row['to_exon'], splicegraph = splicegraph)
    chromosome = from_region_coords[0]
    strand = from_region_coords[1]

    from_flank_ptms = get_ptms_in_splicegraph_flank(event_row['symbol'], chromosome, strand, from_region_coords[-2], from_region_coords[-1], coordinate_type = coordinate_type, which_flank = 'First', flank_size = flank_size)
    to_flank_ptms = get_ptms_in_splicegraph_flank(event_row['symbol'], chromosome, strand, to_region_coords[-2], to_region_coords[-1], coordinate_type = coordinate_type, which_flank = 'Second', flank_size = flank_size)
    ptms_of_interest = pd.concat([from_flank_ptms, to_flank_ptms]).reset_index()


    #if any ptms found for event that could have altered flanking sequences extract sequence information
    if not ptms_of_interest.empty:
        #add additional context from splice data, if indicated
        if event_id_col is not None:
            ptms_of_interest['Region ID'] = region_id
            
        if dPSI_col is not None:
            ptms_of_interest['dPSI'] = dPSI
        
        if sig_col is not None:
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
            inclusion_ptm_loc, exclusion_ptm_loc = get_ptm_locs_in_spliced_sequences(ptm_loc_in_flank, from_sequence, spliced_sequence, to_sequence, strand = strand, which_flank = ptm['Which Flank'], order_by = 'Translation')

            #extract expected flanking sequence based on location in sequence
            inclusion_flank = get_flanking_sequence(inclusion_ptm_loc, inclusion_sequence, ptm_residue = ptm['Residue'], flank_size = flank_size, full_flanking_seq = False)
            exclusion_flank = get_flanking_sequence(exclusion_ptm_loc, exclusion_sequence, ptm_residue = ptm['Residue'], flank_size = flank_size, full_flanking_seq = False)

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

def get_flanking_changes_from_splicegraph(psi_data, splicegraph, ptm_coordinates = None, dPSI_col = None, sig_col = None, event_id_col = None, extra_cols = None, gene_col = 'symbol', flank_size = 5, coordinate_type = 'hg19'):
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

    Returns
    -------
    altered_flanks : pandas.DataFrame
        DataFrame containing the PTMs associated with the flanking regions that are altered, and the flanking sequences that arise depending on whether the flanking sequence is included or not
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None and pose_config.ptm_coordinates is not None:
        ptm_coordinates = pose_config.ptm_coordinates.copy()
    elif ptm_coordinates is None:
        raise ValueError('ptm_coordinates dataframe not provided and not found in the resource files. Please provide the ptm_coordinates dataframe with config.download_ptm_coordinates() or download the file manually. To avoid needing to download this file each time, run pose_config.download_ptm_coordinates(save = True) to save the file locally within the package directory (will take ~63MB of storage space)')

    #load spliceseq
    splicegraph['Region ID'] = splicegraph['Symbol'] + '_' + splicegraph['Exon'].astype(str)
    splicegraph.index = splicegraph['Region ID'].values

    data_for_flanks = psi_data.drop_duplicates().copy() 

    #extract relevant columns
    relevant_columns = ['as_id', 'splice_type', 'symbol', 'from_exon', 'exons', 'to_exon']
    if event_id_col is not None:
        relevant_columns.append(event_id_col)
    if dPSI_col is not None:
        relevant_columns.append(dPSI_col)
    if sig_col is not None:
        relevant_columns.append(sig_col)
    if extra_cols is not None:
        relevant_columns.extend(extra_cols)

    data_for_flanks = data_for_flanks[relevant_columns].drop_duplicates()
    data_for_flanks = data_for_flanks.dropna(subset = ['from_exon', 'to_exon'])
    data_for_flanks['from_region_id'] = data_for_flanks[gene_col]+'_'+data_for_flanks['from_exon'].astype(str)
    data_for_flanks['to_region_id'] = data_for_flanks['symbol']+'_'+data_for_flanks['to_exon'].astype(str)

    #get coordinates for the different regions
    altered_flanks = []
    for i, row in tqdm.tqdm(data_for_flanks.iterrows(), total = data_for_flanks.shape[0], desc = 'Finding flanking changes for splicegraph events'):
        single_event_altered_flanks = get_flank_changes_from_splicegraph_single_event(row, splicegraph, event_id_col = event_id_col, dPSI_col = dPSI_col, sig_col = sig_col, extra_cols = extra_cols, flank_size = flank_size, coordinate_type = coordinate_type)

        altered_flanks.append(single_event_altered_flanks)

    altered_flanks = pd.concat(altered_flanks)    
    return altered_flanks






