import pandas as pd
import numpy as np
from ptm_pose import pose_config

#try importing pyliftover (optional dependency) for genomic coordinate conversion
try:
    import pyliftover
except ImportError:
    pyliftover = None

def extract_filter_kwargs(**kwargs):
    """
    Given keyword arguments from another function, extract all arguments that correspond to filter argument
    """
    filter_kwargs = {'report_removed':kwargs.get('report_removed', True),
                     'alpha': kwargs.get('alpha', 0.05),
                     'min_dpsi':kwargs.get('min_dpsi', 0.2),
                     'modification_class': kwargs.get('modification_class', None),
                     'min_studies': kwargs.get('min_studies', 0),
                     'min_MS_observations': kwargs.get('min_MS_observations', 0), 
                     'min_LTP_studies': kwargs.get('min_LTP_studies', 0),
                     'min_compendia':kwargs.get('min_compendia', 0),
                     'phospho_only_evidence_filter': kwargs.get('phospho_only_evidence_filter', False),
                     'remove_novel': kwargs.get('remove_novel', False)}
    
    return filter_kwargs

def check_filter_kwargs(filter_kwargs):
    """
    Check that there are not any extra filter kwargs that are not used in the filter_ptms function and report
    """
    valid_keys = ['report_removed', 'alpha', 'min_dpsi', 'modification_class', 'min_studies', 'min_MS_observations', 'min_LTP_studies', 'min_compendia', 'phospho_only_evidence_filter', 'remove_novel']
    extra_keys = [key for key in filter_kwargs.keys() if key not in valid_keys]
    
    if len(extra_keys) > 0:
        print(f"Extra unrecognized kwargs found: {extra_keys}. Ignoring these, please fix to match arguments in `helper.filter_ptms()`.")
    
    return True

def quantile_filter(ptms, quantile = 0.9, filter_col = 'MS_LIT'):
    """
    Filter PTMs based on a quantile value for a given column. Will return PTMs that are in the top quantile for each modification class

    Parameters
    ----------
    ptms: pandas DataFrame
        DataFrame containing PTM information
    quantile: float
        Quantile value to filter PTMs on
    filter_col: str
        Column to filter PTMs on

    Returns
    -------
    filtered_ptms: pandas DataFrame
        DataFrame containing PTMs that are in the top quantile for each modification class
    """
    if isinstance(filter_col, str):
        quantiles = pose_config.ptm_coordinates.groupby('Modification Class')[filter_col].quantile(quantile)
        filtered_ptms = []
        for mod_class in ptms['Modification Class'].unique():
            tmp_ptms = ptms[(ptms['Modification Class'] == mod_class) & (ptms[filter_col] >= quantiles[mod_class])].copy()
            filtered_ptms.append(tmp_ptms)
    elif isinstance(filter_col, list):
        tmp = pose_config.ptm_coordinates.copy()
        tmp['filter data'] = tmp[filter_col].sum(axis = 1)
        quantiles = tmp.groupby('Modification Class')['filter data'].quantile(quantile)
        filtered_ptms = []
        for mod_class in ptms['Modification Class'].unique():
            tmp_ptms = ptms[ptms['Modification Class'] == mod_class].copy()
            tmp_ptms['filter data'] = tmp_ptms[filter_col].sum(axis = 1)
            tmp_ptms = tmp_ptms[tmp_ptms['filter data'] >= quantiles[mod_class]]
            filtered_ptms.append(tmp_ptms)     
    else:
        raise TypeError('filter_col must be a string or list of strings')
    return pd.concat(filtered_ptms)



def filter_ptms_by_evidence(ptms, report_removed = False, min_studies = 0, min_MS_observations = 0, min_LTP_studies = 0, min_compendia = 0, filter_phospho_only = True):
    """
    Filter PTMs on various criteria, including number of prior studies site has been observed in, and number of compendia site has been observed in. If indicated, only filter out phosphorylation sites based on limited evidence. If a fraction value, will filter PTMs not based on the strongest evidence sites for each modification class
    """
    if filter_phospho_only:
        filtered_ptms = ptms[ptms['Modification Class'] == 'Phosphorylation'].copy()
        other_ptms = ptms[ptms['Modification Class'] != 'Phosphorylation'].copy()
        mod_type = 'phosphorylation sites'
    else:
        filtered_ptms = ptms.copy()
        mod_type = 'PTMs'

    filtered_ptms = ptms.copy()
    original_shape = filtered_ptms.shape[0]
    current_shape = original_shape
    #restrict to PTMs with more evidence
    if 'LT_LIT' in ptms.columns and min_LTP_studies > 0:
        if min_LTP_studies < 1:
            filtered_ptms = quantile_filter(filtered_ptms, quantile = min_LTP_studies, filter_col = 'LT_LIT')
        else:
            filtered_ptms = filtered_ptms[filtered_ptms['LT_LIT'] >= min_LTP_studies].copy()

        if report_removed:
            num_removed = current_shape - filtered_ptms.shape[0]
            percent_removed = num_removed/original_shape*100
            if min_LTP_studies < 1:
                print(f'{num_removed:,} {mod_type} removed due to being in the bottom {min_LTP_studies*100}% of its PTM type with low-throughput observations in literature ({percent_removed:.2f}%)')
            else:
                print(f'{num_removed:,} {mod_type} removed due to be reported in fewer than {min_studies} low-throughput studies in literature ({percent_removed:.2f}%)')
            current_shape = filtered_ptms.shape[0]
    

    #restrict to PTMs with more evidence
    if 'MS_LIT' in ptms.columns and 'LT_LIT' in ptms.columns and min_studies > 0:
        if min_studies < 1:
            filtered_ptms = quantile_filter(filtered_ptms, quantile = min_studies, filter_col = ['MS_LIT', 'LT_LIT'])
        else:
            filtered_ptms = filtered_ptms[filtered_ptms[['MS_LIT', 'LT_LIT']].sum(axis = 1) >= min_studies]
        if report_removed:
            num_removed = current_shape - filtered_ptms.shape[0]
            percent_removed = num_removed/original_shape*100
            if min_studies < 1:
                print(f'{num_removed:,} {mod_type} removed due to being in the bottom {min_studies*100}% of its PTM type with MS and LT observations in literature ({percent_removed:.2f}%)')
            else:
                print(f'{num_removed:,} {mod_type} removed due to be reported in fewer than {min_studies} studies (MS or low throughput) in literature ({percent_removed:.2f}%)')
            current_shape = filtered_ptms.shape[0]

    if ("MS_CST" in ptms.columns or "MS_LIT" in ptms.columns) and min_MS_observations > 0:
        if min_MS_observations < 1:
            filtered_ptms = quantile_filter(filtered_ptms, quantile = min_MS_observations, filter_col = ['MS_CST', 'MS_LIT'])
        else:
            filtered_ptms = filtered_ptms[filtered_ptms[['MS_CST', 'MS_LIT']].sum(axis = 1) >= min_MS_observations]

        if report_removed:
            num_removed = current_shape - filtered_ptms.shape[0]
            percent_removed = num_removed/original_shape*100
            if min_MS_observations < 1:
                print(f'{num_removed:,} {mod_type} removed due to being in the bottom {min_MS_observations*100}% of its PTM type with MS observations (in literature and by Cell Signaling Technologies) ({percent_removed:.2f}%)')
            else:
                print(f'{num_removed:,} {mod_type} removed due to fewer than {min_MS_observations} MS observations (in literature and by Cell Signaling Technologies): ({percent_removed:.2f}%)')
            current_shape = filtered_ptms.shape[0]

    if 'Compendia' in ptms.columns and min_compendia > 0:
        if min_compendia < 1:
            filtered_ptms = quantile_filter(filtered_ptms, quantile = min_compendia, filter_col = 'Number of Compendia')
        else:
            filtered_ptms = filtered_ptms[filtered_ptms['Number of Compendia'] >= min_compendia].copy()

        if report_removed:
            num_removed = current_shape - filtered_ptms.shape[0]
            percent_removed = num_removed/original_shape*100
            if min_compendia < 1:
                print(f'{num_removed:,} {mod_type} removed due to being in the bottom {min_compendia*100}% of its PTM type with compendia observations ({percent_removed:.2f}%)')
            else:
                print(f'{num_removed:,} {mod_type} removed due to being recorded fewer than {min_compendia} compendia: ({percent_removed:.2f}%)')
            current_shape = filtered_ptms.shape[0]

    if filter_phospho_only:
        filtered_ptms = pd.concat([filtered_ptms, other_ptms])
    
    return filtered_ptms

def filter_ptms(ptms, report_removed = True, min_dpsi = 0.1, alpha = 0.05, modification_class = None, min_studies = 0, min_MS_observations = 0, min_LTP_studies = 0, min_compendia = 0, phospho_only_evidence_filter = False, remove_novel = False):
    """
    Filter PTMs on various criteria, including significance of splice event, number of prior studies site has been observed in, and number of compendia site has been observed in
    """
    #grab specific modification if desired
    if modification_class is not None:
        filtered_ptms = ptms[ptms['Modification Class'] == modification_class].copy()
    else:
        filtered_ptms = ptms.copy()

    original_shape = filtered_ptms.shape[0]
    current_shape = original_shape
    #restrict to PTMs with significant variation
    if 'Significance' in ptms.columns:
        filtered_ptms = filtered_ptms[filtered_ptms['Significance'] <= alpha]
    if 'dPSI' in ptms.columns:
        filtered_ptms = filtered_ptms[filtered_ptms['dPSI'].abs() >= min_dpsi]
    if report_removed and (current_shape != filtered_ptms.shape[0]):
        num_removed = current_shape - filtered_ptms.shape[0]
        percent_removed = num_removed/original_shape*100
        print(f'{num_removed:,} PTMs removed due to insignificant splice event (p < {alpha}, dpsi >= {min_dpsi}): ({percent_removed:.2f}%)')
        current_shape = filtered_ptms.shape[0]

    if remove_novel:
        filtered_ptms = filtered_ptms[~filtered_ptms['Constitutive']]
        if report_removed:
            num_removed = current_shape - filtered_ptms.shape[0]
            percent_removed = num_removed/original_shape*100
            print(f'{num_removed:,} PTMs removed due being novel splice events: ({percent_removed:.2f}%)')
            current_shape = filtered_ptms.shape[0]


    #restrict to PTMs with more evidence
    filtered_ptms = filter_ptms_by_evidence(filtered_ptms, report_removed = report_removed, min_studies = min_studies, min_MS_observations = min_MS_observations, min_LTP_studies = min_LTP_studies, min_compendia=min_compendia, filter_phospho_only = phospho_only_evidence_filter)

    #restrict to PTMs with more evidence (or just )

    if report_removed:
        print(f'Final number of PTMs to be assessed: {filtered_ptms.shape[0]:,}')

    return filtered_ptms

def get_ptm_label(ptm_row, id_type = 'uniprot', consider_isoforms = True):
    """
    Given a row of PTM information, return a string that uniquely identifies the PTM based on the UniProtKB accession, residue, and position in isoform
    """
    if id_type not in ['uniprot', 'gene']:
        raise ValueError("id_type must be either 'uniprot' or 'gene'")
    elif id_type == 'uniprot':
        if consider_isoforms:
            if ptm_row['Isoform Type'] == 'Canonical':
                return f"{ptm_row['UniProtKB Accession']}_{ptm_row['Residue']}{int(ptm_row['PTM Position in Isoform'])}"
            else: 
                return f"{ptm_row['Isoform ID']}_{ptm_row['Residue']}{int(ptm_row['PTM Position in Isoform'])}"
        else:
            return f"{ptm_row['UniProtKB Accession']}_{ptm_row['Residue']}{int(ptm_row['PTM Position in Isoform'])}"
    elif id_type == 'gene':
        return f"{ptm_row['Gene']}_{ptm_row['Residue']}{int(ptm_row['PTM Position in Isoform'])}"

def add_ptm_column(ptms, id_type = 'uniprot', consider_isoforms = True):
    """
    Given a dataframe of PTM information, add a column that uniquely identifies each PTM based on the UniProtKB accession, residue, and position in isoform
    """
    if id_type not in ['uniprot', 'gene']:
        raise ValueError("id_type must be either 'uniprot' or 'gene'")

    ptms['PTM'] = ptms.apply(lambda x: get_ptm_label(x, id_type = id_type, consider_isoforms = consider_isoforms), axis = 1)
    return ptms


def load_example_data(spliced_ptms = False, altered_flanks = False, annotated_data = False):
    """Download example data for PTM-POSE, which is generated from applying PTM-POSE pipeline to MATS data from Yang et al, 2016
    
    Parameters
    ----------
    spliced_ptms: bool
        If True, will return example data for differentially included PTMs
    altered_flanks: bool
        If True, will return example data for PTMs with altered flanking sequences
    
    Returns
    -------
    spliced_ptms: pandas DataFrame
        DataFrame containing example data for differentially included PTMs. Returns only if spliced_ptms is true
    altered_flanks: pandas DataFrame
        DataFrame containing example data for PTMs with altered flanking sequences. Returns only if altered_flanks is true
    
    """
    output_data = []
    if annotated_data:
        annotated_data = pd.read_csv(pose_config.resource_dir + '/Example/splice_data.csv')
        output_data.append(annotated_data)
    if spliced_ptms:
        spliced_ptms = pd.read_csv(pose_config.resource_dir + '/Example/spliced_ptms.csv')
        output_data.append(spliced_ptms)
    if altered_flanks:
        altered_flanks = pd.read_csv(pose_config.resource_dir + '/Example/altered_flanks.csv')
        output_data.append(altered_flanks)

    if len(output_data) == 1:
        return output_data[0]
    elif len(output_data) > 1:
        return output_data
    else:
        raise ValueError("No data was requested. Please set spliced_ptms, altered_flanks, or annotated_data to True to load example data.")
    
def convert_genomic_coordinates(loc, chromosome, strand, from_coord = 'hg19', to_coord = 'hg38', liftover_object = None):
    """
    Convert genomic coordinates from one genome assembly to another using pyliftover, if available. If not, will raise an error.

    Parameters
    ----------
    loc: int
        Genomic location to convert
    chromosome: str
        Chromosome of the genomic location to convert
    strand: str
        Strand of the genomic location to convert, either '+' or '-'
    """
    if liftover_object is None and pyliftover is not None:
        liftover_object = pyliftover.LiftOver(from_coord, to_coord)
    elif liftover_object is None:
        raise ValueError("Liftover object must be provided or pyliftover must be installed to convert genomic coordinates")
    
    #convert strand
    strand = convert_strand_symbol(strand, to = 'symbol')

    chromosome = f'chr{chromosome}'
    try:
        results = liftover_object.convert_coordinate(chromosome, loc - 1, strand)
    except:
        print('Error')
        print(chromosome, loc, strand)

    if len(results) > 0:
        new_chromosome = results[0][0]
        new_strand = results[0][2]
        if new_chromosome == chromosome and new_strand == strand:
            return int(results[0][1]) + 1
        else:
            return -1
    else:
        return np.nan
    
def convert_genomic_coordinates_df(df, from_coord = 'hg19', to_coord = 'hg38', loc_col = 'Genomic Location', chromosome_col = 'Chromosome', strand_col = 'Strand', output_col = None):
    """
    Convert genomic coordinates from one genome assembly to another using pyliftover, if available. If not, will raise an error.

    Parameters
    ----------
    df: pandas DataFrame
        DataFrame containing genomic coordinates to convert
    from_coord: str
        Genome assembly to convert from (default: 'hg19')
    to_coord: str
        Genome assembly to convert to (default: 'hg38')
    loc_col: str
        Column name for genomic location in the dataframe (default: 'Genomic Location')
    chromosome_col: str
        Column name for chromosome in the dataframe (default: 'Chromosome')
    strand_col: str
        Column name for strand in the dataframe (default: 'Strand')

    Returns
    -------
    pandas DataFrame
        DataFrame with converted genomic coordinates in a new column labeled 'Genomic Location (<to_coord>)'
    """
    if pyliftover is None:
        raise ImportError("pyliftover must be installed to convert genomic coordinates. Please install it using 'pip install pyliftover'.")
    
    liftover_object = pyliftover.LiftOver(from_coord, to_coord)
    
    if output_col is None:
        output_col = f'Genomic Location ({to_coord})'


    df = df.copy()
    df[output_col] = df.apply(lambda x: convert_genomic_coordinates(x[loc_col], x[chromosome_col], x[strand_col], from_coord=from_coord, to_coord=to_coord, liftover_object=liftover_object), axis=1)
    
    return df
    
def convert_strand_symbol(strand, to = 'int'):
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
    if to == 'int':
        if isinstance(strand, str):
            if strand == '+' or strand == '1':
                return 1
            elif strand == '-' or strand == '-1':
                return -1
        else:
            return strand
    elif to == 'symbol':
        if isinstance(strand, str):
            return strand
        elif strand == -1:
            return '-'
        elif strand == 1:
            return '+'
    else:
        raise ValueError("`to` must be either 'int' or 'symbol'")

def join_unique_entries(x, sep = ';'):
    """
    For use with groupby, combines all unique entries separated by ';', removing any NaN entries
    """
    #check if only nan entries
    if all(i != i for i in x):
        return np.nan
    if any(sep in i for i in x if i == i): #check if ';' already in entry, if so, split and remove any NaN entries
        split_list = [i.split(sep) for i in x if i == i]
        #split entries in list by ';' and flatten list
        flat_list = [item for sublist in split_list for item in sublist]
        return sep.join(set(flat_list))
    else:
        entry_list = [str(i) for i in x if i == i]
        return sep.join(set(entry_list))

def join_entries(x, sep = ';'):
    """
    For use with groupby, combines all entries separated by ';', removing any NaN entries
    """
    #check if only nan entries
    if all(i != i for i in x):
        return np.nan

    if any(sep in i for i in x if i == i): #check if ';' already in entry, if so, split and remove any NaN entries
        split_list = [i.split(sep) for i in x if i == i]
        #split entries in list by ';' and flatten list
        flat_list = [item for sublist in split_list for item in sublist]
        return sep.join(flat_list)

    else:
        entry_list = [str(i) for i in x if i == i]
        return sep.join(entry_list)

def join_except_self(df, group_col, value_col, new_col, sep = ';'):
    """
    For a given dataframe, combines all entries with the same information except for the current row, adds that to the new_col label, and returns the updated dataframe
    
    Parameters
    ----------
    df: pandas DataFrame
        The dataframe to be updated
    group_col: str
        The column to group the dataframe by
    value_col: str
        The column to be combined
    new_col: str
        The new column to be added to the dataframe with the grouped information (excluding the info from the current row)

    Returns
    -------
    df: pandas DataFrame
        updated dataframe with new col labeled with new_col value
    """
    df = df.copy()
    df[new_col] = df.groupby(group_col)[value_col].transform(join_unique_entries, sep)

    #go through each row and remove the value(s) in the new column that is in the value column
    new_values = []
    for i, row in df.iterrows():
        if row[new_col] == row[new_col] and row[value_col] == row[value_col]:
            new_values.append(';'.join([trans for trans in row[new_col].split(sep) if trans not in row[value_col].split(sep)]))
        elif row[value_col] != row[value_col]:
            new_values.append(row[new_col])
        else:
            new_values.append(np.nan)
    df[new_col] = new_values
    return df