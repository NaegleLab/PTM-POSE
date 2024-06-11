#biopython packages
from Bio.Data import CodonTable
from Bio import pairwise2

#standard packages
import numpy as np
import pandas as pd
import re

#PTM pose functions
from ptm_pose import database_interfacing as di
from ptm_pose import project
from ptm_pose import pose_config as config

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
    Given the location of a PTM in a flanking sequence, extract the location of the PTM in the inclusion sequence and the exclusion sequence associated with a given splice event. Inclusion sequence will include the skipped exon region, retained intron, or longer alternative splice site depending on event type. The PTM location should be associated with where the PTM is located relative to spliced region (before = 'First', after = 'Second').

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
        Tuple containing the PTM location in the inclusion sequence and the exclusion sequence
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

def extract_region_from_splicegraph(spliceseq, region_id):
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
    region_info = spliceseq.loc[region_id]
    strand = project.convert_strand_symbol(region_info['Strand'])
    if strand == 1:
        return [region_info['Chromosome'], strand,region_info['Chr_Start'], region_info['Chr_Stop']]
    else:
        return [region_info['Chromosome'], strand,region_info['Chr_Stop'], region_info['Chr_Start']]
    
def get_spliceseq_event_regions(spliceseq_event, splicegraph):
    first_exon_region = extract_region_from_splicegraph(splicegraph, region_id = spliceseq_event['symbol']+'_'+str(spliceseq_event['from_exon']))
    spliced_regions = [extract_region_from_splicegraph(splicegraph, spliceseq_event['symbol']+'_'+exon) if '.' in exon else extract_region_from_splicegraph(splicegraph, spliceseq_event['symbol']+'_'+exon+'.0') for exon in spliceseq_event['exons'].split(':')]
    second_exon_region = extract_region_from_splicegraph(splicegraph, region_id = spliceseq_event['symbol']+'_'+str(spliceseq_event['to_exon']))
    return first_exon_region, spliced_regions, second_exon_region





def get_flanking_changes(ptm_coordinates, chromosome, strand, first_flank_region, spliced_region, second_flank_region, gene = None, event_id = None, flank_size = 5, coordinate_type = 'hg38', lowercase_mod = True, order_by = 'Coordinates'):
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
    ptms_in_region_first_flank = project.find_PTMs_in_region(ptm_coordinates, chromosome, strand, first_flank_region[0], first_flank_region[1], gene = gene, coordinate_type = coordinate_type)
    if not ptms_in_region_first_flank.empty:
        ptms_in_region_first_flank = ptms_in_region_first_flank[ptms_in_region_first_flank['Proximity to Region End (bp)'] < flank_size*3]
        ptms_in_region_first_flank['Region'] = 'First'
    #check second flank for ptms
    ptms_in_region_second_flank = project.find_PTMs_in_region(ptm_coordinates, chromosome, strand, second_flank_region[0], second_flank_region[1], gene = gene, coordinate_type = coordinate_type)
    if not ptms_in_region_second_flank.empty:
        ptms_in_region_second_flank = ptms_in_region_second_flank[ptms_in_region_second_flank['Proximity to Region Start (bp)'] < flank_size*3]
        ptms_in_region_second_flank['Region'] = 'Second'

    #combine
    ptms_in_region = pd.concat([ptms_in_region_first_flank, ptms_in_region_second_flank])


    if ptms_in_region.empty:
        return pd.DataFrame()
    else:
        #restrict to ptms within boundary
        if ptms_in_region.empty:
            return pd.DataFrame()
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
                if codon_table.forward_table[ptm_codon_inclusion] == ptm['Residue'] and codon_table.forward_table[ptm_codon_exclusion] == ptm['Residue']:
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
            ptms_in_region = ptms_in_region[['Source of PTM', 'Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']].reset_index(drop = True)
        else:
            ptms_in_region = ptms_in_region[['Source of PTM', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']].reset_index(drop = True)
        #add flanking sequence information to ptm dataframe
        ptms_in_region['Inclusion Sequence'] = inclusion_seq_list
        ptms_in_region['Exclusion Sequence'] = exclusion_seq_list
        ptms_in_region['Region'] = flank_region_list
        ptms_in_region['Translation Success'] = translate_success_list

        if event_id is not None:
            ptms_in_region.insert(0, 'Event ID', event_id)

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
    if ptm_coordinates is None and config.ptm_coordinates is not None:
        ptm_coordinates = config.ptm_coordinates
    elif ptm_coordinates is None:
        raise ValueError('ptm_coordinates dataframe not provided and not found in the resource files. Please provide the ptm_coordinates dataframe with config.download_ptm_coordinates() or download the file manually. To avoid needing to download this file each time, run pose_config.download_ptm_coordinates(save = True) to save the file locally within the package directory (will take ~63MB of storage space)')

    #check to make sure all required columns are provided
    if chromosome_col is None and strand_col is None and first_flank_start_col is None and first_flank_end_col is None and spliced_region_start_col is None and spliced_region_end_col is None and second_flank_start_col is None and second_flank_end_col is None:
        raise ValueError('Please provide column names for chromosome, strand, first flank start, first flank end, spliced region start, spliced region end, second flank start, and second flank end.')

    #if chromosome is labeled with 'chr', remove
    if splice_data[chromosome_col].str.contains('chr').any():
        splice_data['chr'] = splice_data['chr'].str.strip('chr')
    

    results = []
    for i, event in splice_data.iterrows():
        if event_id_col is None:
            event_id = i
        else:
            event_id = event[event_id_col]

        #get gene info
        chromosome = event[chromosome_col]
        strand = event[strand_col]
        gene = event[gene_col] if gene_col is not None else None

        #extract region inof
        first_flank_region = [event[first_flank_start_col], event[first_flank_end_col]]
        spliced_region = [event[spliced_region_start_col], event[spliced_region_end_col]]
        second_flank_region = [event[second_flank_start_col], event[second_flank_end_col]]

        #get flanking changes
        ptm_flanks = get_flanking_changes(ptm_coordinates, chromosome, strand, first_flank_region, spliced_region, second_flank_region, event_id = event_id, flank_size = flank_size, coordinate_type = coordinate_type, lowercase_mod=lowercase_mod)

        #append to results
        results.append(ptm_flanks)

    #combine and remove any failed translation attempts
    results = pd.concat(results)
    results = results.dropna(subset = ['Inclusion Sequence', 'Exclusion Sequence'])

    #do some quick comparison of flanking sequences
    results['Matched'] = results['Inclusion Sequence'] == results['Exclusion Sequence']
    results['Stop Codon Introduced'] = (results['Inclusion Sequence'].str.contains(r'\*')) | (results['Exclusion Sequence'].str.contains(r'\*'))
    return results


def getSequenceIdentity(seq1, seq2):
    """
    Given two flanking sequences, calculate the sequence identity between them using Biopython and parameters definded by Pillman et al. BMC Bioinformatics 2011

    Parameters
    ----------
    seq1, seq2: str
        flanking sequence 

    Returns
    -------
    normalized_score: float
        normalized score of sequence similarity between flanking sequences (calculated similarity/max possible similarity)
    """
    #align canonical and alternative flanks, return only the score
    actual_similarity = pairwise2.align.globalxs(seq1, seq2, -10, -2, score_only = True)
    #aling the canonical flank to itself, return only the score
    control_similarity = pairwise2.align.globalxs(seq1, seq1, -10, -2, score_only = True)
    #normalize score
    normalized_score = actual_similarity/control_similarity
    return normalized_score

def findAlteredPositions(seq1, seq2, flank_size = 5):
    """
    Given two sequences, identify the location of positions that have changed

    Parameters
    ----------
    seq1, seq2: str
        sequences to compare (order does not matter)
    flank_size: int
        size of the flanking sequences (default is 5). This is used to make sure the provided sequences are the correct length
    
    Returns
    -------
    altered_positions: list
        list of positions that have changed
    residue_change: list
        list of residues that have changed associated with that position
    flank_side: str
        indicates which side of the flanking sequence the change has occurred (N-term, C-term, or Both)
    """
    desired_seq_size = flank_size*2+1
    altered_positions = []
    residue_change = []
    flank_side = []
    seq_size = len(seq1)
    flank_size = (seq_size -1)/2
    if seq_size == len(seq2) and seq_size == desired_seq_size:
        for i in range(seq_size):
            if seq1[i] != seq2[i]:
                altered_positions.append(i-(flank_size))
                residue_change.append(f'{seq1[i]}->{seq2[i]}')
        #check to see which side flanking sequence
        altered_positions = np.array(altered_positions)
        n_term = any(altered_positions < 0)
        c_term = any(altered_positions > 0)
        if n_term and c_term:
            flank_side = 'Both'
        elif n_term:
            flank_side = 'N-term only'
        elif c_term:
            flank_side = 'C-term only'
        else:
            flank_side = 'Unclear'
        return altered_positions, residue_change, flank_side
    else:
        return np.nan, np.nan, np.nan
    
def compare_flanking_sequences(flanking_sequences, flank_size = 5):
    sequence_identity_list = []
    altered_positions_list = []
    residue_change_list = []
    flank_side_list = []
    for i, row in flanking_sequences.iterrows():
        #if sequences are not the same and do not introduce stop codons, compare sequence identity
        if not row['Stop Codon Introduced'] and not row['Matched']:
            #compare sequence identity
            sequence_identity = getSequenceIdentity(row['Inclusion Sequence'], row['Exclusion Sequence'])
            #identify where flanking sequence changes
            altered_positions, residue_change, flank_side = findAlteredPositions(row['Inclusion Sequence'], row['Exclusion Sequence'], flank_size = flank_size)
        else:
            sequence_identity = np.nan
            altered_positions = np.nan
            residue_change = np.nan
            flank_side = np.nan

        #add to lists
        sequence_identity_list.append(sequence_identity)
        altered_positions_list.append(altered_positions)
        residue_change_list.append(residue_change)
        flank_side_list.append(flank_side)

    flanking_sequences['Sequence Identity'] = sequence_identity_list
    flanking_sequences['Altered Positions'] = altered_positions_list
    flanking_sequences['Residue Change'] = residue_change_list
    flanking_sequences['Altered Flank Side'] = flank_side_list
    return flanking_sequences

    
def find_motifs(seq, elm_classes):
    """
    Given a sequence and a dataframe containinn ELM class information, identify motifs that can be found in the provided sequence using the RegEx expression provided by ELM (PTMs not considered). This does not take into account the position of the motif in the sequence or additional information that might validate any potential interaction (i.e. structural information that would indicate whether the motif is accessible or not). ELM class information can be downloaded from the download page of elm (http://elm.eu.org/elms/elms_index.tsv).

    Parameters
    ----------
    seq: str
        sequence to search for motifs
    elm_classes: pandas.DataFrame
        DataFrame containing ELM class information (ELMIdentifier, Regex, etc.), downloaded directly from ELM (http://elm.eu.org/elms/elms_index.tsv)
    """
    matches = []
    for j, elm_row in elm_classes.iterrows():
        reg_ex = elm_row['Regex']
        if re.search(reg_ex, seq) is not None:
            matches.append(elm_row['ELMIdentifier'])

    return matches

def compare_inclusion_motifs(flanking_sequences, elm_classes = None):
    """
    Given a DataFrame containing flanking sequences with changes and a DataFrame containing ELM class information, identify motifs that are found in the inclusion and exclusion events, identifying motifs unique to each case. This does not take into account the position of the motif in the sequence or additional information that might validate any potential interaction (i.e. structural information that would indicate whether the motif is accessible or not). ELM class information can be downloaded from the download page of elm (http://elm.eu.org/elms/elms_index.tsv).

    Parameters
    ----------
    flanking_sequences: pandas.DataFrame
        DataFrame containing flanking sequences with changes, obtained from get_flanking_changes_from_splice_data()
    elm_classes: pandas.DataFrame
        DataFrame containing ELM class information (ELMIdentifier, Regex, etc.), downloaded directly from ELM (http://elm.eu.org/elms/elms_index.tsv). Recommended to download this file and input it manually, but will download from ELM otherwise

    Returns
    -------
    flanking_sequences: pandas.DataFrame
        DataFrame containing flanking sequences with changes and motifs found in the inclusion and exclusion events

    """
    if elm_classes is None:
        elm_classes = pd.read_csv('http://elm.eu.org/elms/elms_index.tsv', sep = '\t', header = 5)

    only_in_inclusion = []
    only_in_exclusion = []

    for i, row in flanking_sequences.iterrows():
        inclusion_matches = find_motifs(row['Inclusion Sequence'], elm_classes)
        exclusion_matches = find_motifs(row['Exclusion Sequence'], elm_classes)

        only_in_inclusion.append(';'.join(set(inclusion_matches) - set(exclusion_matches)))
        only_in_exclusion.append(';'.join(set(exclusion_matches) - set(inclusion_matches)))

    flanking_sequences["Motif only in Inclusion"] = only_in_inclusion
    flanking_sequences["Motif only in Exclusion"] = only_in_exclusion
    return flanking_sequences
