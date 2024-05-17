from Bio.Data import CodonTable

import numpy as np
import pandas as pd

import database_interfacing as di
import project


# Get the standard codon table
codon_table = CodonTable.unambiguous_dna_by_name["Standard"].forward_table.copy()


def translate_flanking_sequence(seq, flank_size = 7, full_flanking_seq = True, lowercase_mod = True, first_flank_length = None):
    aa_seq = ''
    if len(seq) == flank_size*2*3+3:
        for i in range(0, len(seq), 3):
            if seq[i:i+3] in codon_table.keys():
                aa = codon_table[seq[i:i+3]]
            else:
                aa = 'X'

            if i/3 == flank_size and lowercase_mod:
                aa = aa.lower()
            aa_seq += aa
    elif len(seq) % 3 == 0 and not full_flanking_seq:
        for i in range(0, len(seq), 3):
            if seq[i:i+3] in codon_table.keys():
                aa = codon_table[seq[i:i+3]]
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



def get_flanking_changes(ptm_coordinates, chromosome, strand, first_flank_region, spliced_region, second_flank_region, event_id = None, flank_size = 7, coordinate_type = 'hg38', lowercase_mod = True):
    strand = project.convert_strand_symbol(strand)
    #check first flank for ptms
    ptms_in_region_first_flank = project.find_PTMs_in_region(ptm_coordinates, chromosome, strand, first_flank_region[0], first_flank_region[1], coordinate_type = coordinate_type)
    if not ptms_in_region_first_flank.empty:
        ptms_in_region_first_flank = ptms_in_region_first_flank[ptms_in_region_first_flank['Proximity to Region End (bp)'] < flank_size*3]
        ptms_in_region_first_flank['Region'] = 'First Flank'
    #check second flank for ptms
    ptms_in_region_second_flank = project.find_PTMs_in_region(ptm_coordinates, chromosome, strand, second_flank_region[0], second_flank_region[1], coordinate_type = coordinate_type)
    if not ptms_in_region_second_flank.empty:
        ptms_in_region_second_flank = ptms_in_region_second_flank[ptms_in_region_second_flank['Proximity to Region Start (bp)'] < flank_size*3]
        ptms_in_region_second_flank['Region'] = 'Second Flank'

    #combine
    ptms_in_region = pd.concat([ptms_in_region_first_flank, ptms_in_region_second_flank])


    if ptms_in_region.empty:
        return pd.DataFrame()
    else:
        #restrict to ptms within boundary
        if ptms_in_region.empty:
            return pd.DataFrame()
        #get sequences associated with flanking regions and spliced regions
        first_flank_region_query = [chromosome, strand] + first_flank_region
        spliced_region_query = [chromosome, strand] + spliced_region
        second_flank_region_query = [chromosome, strand] + second_flank_region
        regions_list = [first_flank_region_query, spliced_region_query, second_flank_region_query]
        first_flank_seq, spliced_seq, second_flank_seq = di.get_region_sequences_from_list(regions_list, coordinate_type = coordinate_type)

        #combine sequences for inclusion and exclusion cases
        if strand == 1:
            inclusion_seq = first_flank_seq + spliced_seq + second_flank_seq
            exclusion_seq = first_flank_seq + second_flank_seq
        else:
            inclusion_seq = second_flank_seq + spliced_seq + first_flank_seq
            exclusion_seq = second_flank_seq + first_flank_seq

        translate_success_list = []
        inclusion_seq_list = []
        exclusion_seq_list = []
        flank_region_list = []
        for i, ptm in ptms_in_region.iterrows():
            ptm_loc = ptm['Gene Location (hg19)']
            flank_region = ptm['Region']
            flank_region_loc = ptm['Region']
            flank_region = first_flank_region if flank_region_loc == 'First Flank' else second_flank_region
            #grab ptm loc based on which strand ptm is on
            if strand == 1:
                relative_ptm_loc = int(ptm_loc - flank_region[0])
            else:
                relative_ptm_loc = int(flank_region[1] - ptm_loc)

            #grab codon associated with ptm in sequence
            if (flank_region_loc == 'First Flank' and strand == 1) or (flank_region_loc == 'Second Flank' and strand == -1):
                inclusion_ptm_loc, exclusion_ptm_loc = relative_ptm_loc, relative_ptm_loc
            elif (strand == -1 and flank_region_loc == 'First Flank'):
                inclusion_ptm_loc =  relative_ptm_loc+len(spliced_seq)+len(second_flank_seq)
                exclusion_ptm_loc =  relative_ptm_loc+len(second_flank_seq)
            elif (strand == 1 and flank_region_loc == 'Second Flank'):
                inclusion_ptm_loc =  relative_ptm_loc+len(spliced_seq)+len(first_flank_seq)
                exclusion_ptm_loc =  relative_ptm_loc+len(first_flank_seq)

            ptm_codon_inclusion = inclusion_seq[inclusion_ptm_loc:inclusion_ptm_loc+3]
            ptm_codon_exclusion = exclusion_seq[exclusion_ptm_loc:exclusion_ptm_loc+3]



            #check if ptm codon codes for amino acid and then extract flanking sequence
            correct_seq = False
            if ptm_codon_inclusion in codon_table.keys() and ptm_codon_exclusion in codon_table.keys():
                if codon_table[ptm_codon_inclusion] == ptm['Residue'] and codon_table[ptm_codon_exclusion] == ptm['Residue']:
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

        if event_id is not None:
            return pd.DataFrame({'Event ID':event_id, 'PTM':ptms_in_region['Source of PTM'],'Inclusion Sequence':inclusion_seq_list,'Exclusion Sequence':exclusion_seq_list,'Translation Success':translate_success_list, 'Region':flank_region_list})
        else:
            return pd.DataFrame({'PTM':ptms_in_region['Source of PTM'],'Inclusion Sequence':inclusion_seq_list,'Exclusion Sequence':exclusion_seq_list,'Translation Success':translate_success_list, 'Region':flank_region_list})


def get_flanking_changes_from_splice_data(splice_data, ptm_coordinates, chromosome_col = None, strand_col = None, first_flank_start_col = None, first_flank_end_col = None, spliced_region_start_col = None, spliced_region_end_col = None, second_flank_start_col = None, second_flank_end_col = None, event_id_col = None, flank_size = 7, coordinate_type = 'hg38', lowercase_mod = True):
    if chromosome_col is None and strand_col is None and first_flank_start_col is None and first_flank_end_col is None and spliced_region_start_col is None and spliced_region_end_col is None and second_flank_start_col is None and second_flank_end_col is None:
        raise ValueError('Please provide column names for chromosome, strand, first flank start, first flank end, spliced region start, spliced region end, second flank start, and second flank end.')


    results = []
    for i, event in splice_data.iterrows():
        if event_id_col is None:
            event_id = i
        else:
            event_id = event[event_id_col]
        event_id = event[event_id_col]
        chromosome = event[chromosome_col]
        strand = event[strand_col]
        first_flank_region = [event[first_flank_start_col], event[first_flank_end_col]]
        spliced_region = [event[spliced_region_start_col], event[spliced_region_end_col]]
        second_flank_region = [event[second_flank_start_col], event[second_flank_end_col]]
        ptm_flanks = get_flanking_changes(ptm_coordinates, chromosome, strand, first_flank_region, spliced_region, second_flank_region, event_id = event_id, flank_size = flank_size, coordinate_type = coordinate_type, lowercase_mod=lowercase_mod)
        results.append(ptm_flanks)

    return results