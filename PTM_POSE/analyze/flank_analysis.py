
import numpy as np
import pandas as pd

from ptm_pose import helpers

#plotting 
import matplotlib.pyplot as plt
import seaborn as sns

#analysis packages
from Bio.Align import PairwiseAligner
import re



def compare_flanking_sequences(altered_flanks, flank_size = 5):
    """
    Given the altered_flanks dataframe, compare the flanking sequences in the inclusion and exclusion isoforms. This includes identifying sequence identity, altered positions, residue changes, and the side of the flanking sequence that is altered (N-term, C-term, or both).

    Parameters
    ----------
    altered_flanks: pandas.DataFrame
        DataFrame containing flanking sequences with changes, obtained from get_flanking_changes_from_splice_data(). It should contain the columns 'Inclusion Flanking Sequence' and 'Exclusion Flanking Sequence' for comparison.
    flank_size : int, optional
        Size of the flanking region to analyze. Default is 5. This should be less than half the length of the shortest sequence to ensure proper comparison.
    """
    sequence_identity_list = []
    altered_positions_list = []
    residue_change_list = []
    flank_side_list = []
    for i, row in altered_flanks.iterrows():
        #if there is sequence info for both and does not introduce stop codons, compare sequence identity
        if not row['Stop Codon Introduced'] and row['Inclusion Flanking Sequence'] == row['Inclusion Flanking Sequence'] and row['Exclusion Flanking Sequence'] == row['Exclusion Flanking Sequence']:
            #compare sequence identity
            sequence_identity = getSequenceIdentity(row['Inclusion Flanking Sequence'], row['Exclusion Flanking Sequence'])
            #identify where flanking sequence changes
            altered_positions, residue_change, flank_side = findAlteredPositions(row['Inclusion Flanking Sequence'], row['Exclusion Flanking Sequence'], flank_size = flank_size)
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

    altered_flanks['Sequence Identity'] = sequence_identity_list
    altered_flanks['Altered Positions'] = altered_positions_list
    altered_flanks['Residue Change'] = residue_change_list
    altered_flanks['Altered Flank Side'] = flank_side_list
    return altered_flanks



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

    for _, row in flanking_sequences.iterrows():
        #check if there is a stop codon introduced and both flanking sequences are present
        if not row['Stop Codon Introduced'] and row['Inclusion Flanking Sequence'] == row['Inclusion Flanking Sequence'] and row['Exclusion Flanking Sequence'] == row['Exclusion Flanking Sequence']:
            #get elm motifs that match inclusion or Exclusion Flanking Sequences
            inclusion_matches = find_motifs(row['Inclusion Flanking Sequence'], elm_classes)
            exclusion_matches = find_motifs(row['Exclusion Flanking Sequence'], elm_classes)

            #get motifs that are unique to each case
            only_in_inclusion.append(';'.join(set(inclusion_matches) - set(exclusion_matches)))
            only_in_exclusion.append(';'.join(set(exclusion_matches) - set(inclusion_matches)))
        else:
            only_in_inclusion.append(np.nan)
            only_in_exclusion.append(np.nan)

    #save data
    flanking_sequences["Motif only in Inclusion"] = only_in_inclusion
    flanking_sequences["Motif only in Exclusion"] = only_in_exclusion
    return flanking_sequences

def identify_change_to_specific_motif(altered_flanks, elm_motif_name, elm_classes = None, modification_class = None, residues = None, dPSI_col = None, **kwargs):
    """
    Given the altered_flanks dataframe, identify PTM flanking sequences that match different ELM motifs in the inclusion and exclusion isoform. This function allows for filtering based on specific ELM motifs, modification classes, and residues of interest. It also allows for additional filtering through kwargs.

    Parameters
    ----------
    altered_flanks: pandas.DataFrame
        DataFrame containing flanking sequences with changes, obtained from get_flanking_changes_from_splice_data() and analyzed with compare_flanking_sequences(). It should contain the columns 'Inclusion Flanking Sequence' and 'Exclusion Flanking Sequence' for comparison.
    elm_motif_name: str
        Name of ELM motif to look for in inclusion and exclusion isoforms
    elm_classes: pandas.DataFrame, optional
    
        DataFrame containing ELM class information (ELMIdentifier, Regex, etc.), downloaded directly from ELM (http://elm.eu.org/elms/elms_index.tsv). If not provided, the function will download it automatically.
    modification_class : str, optional
        
        Specific modification class to filter for. This is useful if you are looking for a specific type of PTM (e.g., phosphorylation, acetylation). Default is None, which will not filter based on modification class.
    residues : str
        Specific residue to focus on. If none will look at all residues
    dPSI_col : str, optional
        column containing deltaPSI information
    kwargs : additional keyword arguments
        additional keyword arguments to pass to the filter_ptms function
    """

    #filter ptms 
    filter_arguments = helpers.extract_filter_kwargs(**kwargs)
    helpers.check_filter_kwargs(filter_arguments)
    altered_flanks = helpers.filter_ptms(altered_flanks, **filter_arguments)

    #compare flanking sequences
    if 'Altered Positions' not in altered_flanks.columns:
        altered_flanks = compare_flanking_sequences(altered_flanks)
    
    #grab elm motifs that match inclusion or Exclusion Flanking Sequences
    if 'Motif only in Inclusion' not in altered_flanks.columns:
        altered_flanks = compare_inclusion_motifs(altered_flanks, elm_classes = elm_classes)

    #grab only needed info
    motif_data = altered_flanks.dropna(subset = ['Inclusion Flanking Sequence', 'Exclusion Flanking Sequence'], how = 'all').copy()
    cols_to_keep = ['Region ID', 'Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class', 'Inclusion Flanking Sequence', 'Exclusion Flanking Sequence', 'Motif only in Inclusion', 'Motif only in Exclusion', 'Altered Positions', 'Residue Change']
    if dPSI_col is not None:
        cols_to_keep.append(dPSI_col)

    #go through motif data and identify motifs matching elm motif of interest
    motif_data = motif_data[cols_to_keep]
    for i, row in motif_data.iterrows():
        if row['Motif only in Inclusion'] == row['Motif only in Inclusion']:
            if elm_motif_name in row['Motif only in Inclusion']:
                motif_data.loc[i, 'Motif only in Inclusion'] = ';'.join([motif for motif in row['Motif only in Inclusion'].split(';') if elm_motif_name in motif])
            else:
                motif_data.loc[i, 'Motif only in Inclusion'] = np.nan

        if row['Motif only in Exclusion'] == row['Motif only in Exclusion']:
            if elm_motif_name in row['Motif only in Exclusion']:
                motif_data.loc[i, 'Motif only in Exclusion'] = ';'.join([motif for motif in row['Motif only in Exclusion'].split(';') if elm_motif_name in motif])
            else:
                motif_data.loc[i, 'Motif only in Exclusion'] = np.nan

    #restrict to events that are specific modification types or residues (for example, SH2 domain motifs should be phosphotyrosine)
    motif_data = motif_data.dropna(subset = ['Motif only in Inclusion', 'Motif only in Exclusion'], how = 'all')
    if modification_class is not None:
        motif_data = motif_data[motif_data['Modification Class'].str.contains(modification_class)]

    if residues is not None and isinstance(residues, str):
        motif_data = motif_data[motif_data['Residue'] == residues]
    elif residues is not None and isinstance(residues, list):
        motif_data = motif_data[motif_data['Residue'].isin(residues)]
    elif residues is not None:
        raise ValueError('residues parameter must be a string or list of strings')
    
    return motif_data

    



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
    #make pairwise aligner object
    aligner = PairwiseAligner()
    #set parameters, with match score of 10 and mismatch score of -2
    aligner.mode = 'global'
    aligner.match_score = 10
    aligner.mismatch_score = -2
    #calculate sequence alignment score between two sequences
    actual_similarity = aligner.align(seq1, seq2)[0].score
    #calculate sequence alignment score between the same sequence
    control_similarity = aligner.align(seq1, seq1)[0].score
    #normalize score
    normalized_score = actual_similarity/control_similarity
    return normalized_score

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

def plot_sequence_differences(inclusion_seq, exclusion_seq, dpsi = None, flank_size = 5, figsize = (3,1)):
    """
    Given the flanking sequences for a PTM resulting from a specific splice event, plot the differences between the two sequences, coloring the changing residues in red. If dPSI is also provided, will add an arrow to the plot indicating the direction of change

    Parameters
    ----------
    inclusion_seq: str
        Sequence of the inclusion isoform (with spliced region included)
    exclusion_seq: str
        Sequence of the exclusion isoform (with spliced region excluded)
    dpsi: float
        Change in PSI for the specific splice event. Default is None, which will not add an arrow to the plot.
    flank_size: int
        Size of flanking region to plot. Default is 5. This must be less than half the length of the shortest sequence.
    figsize: tuple
        Size of the figure. Default is (3,1).
    """
    #convert sequence into list of residues
    inclusion_seq = list(inclusion_seq)
    exclusion_seq = list(exclusion_seq)

    if len(inclusion_seq) < flank_size*2 + 1 or len(exclusion_seq) < flank_size*2 + 1:
        min_seq_length = min(len(inclusion_seq), len(exclusion_seq))
        raise ValueError(f'Longest allowable flank size for provided sequences is {(min_seq_length -1)/2}, or {min_seq_length} residues. Please use smaller flank size or rerun flanking sequence analysis to get longer flank sequences')
    
    #if either provided seq is longer than flank size, trim to flank size
    if len(inclusion_seq) > flank_size*2 + 1:
        lowercase_pos = [i for i, x in enumerate(inclusion_seq) if x.islower()]
        inclusion_seq = inclusion_seq[int(lowercase_pos[0])-flank_size:int(lowercase_pos[0])+flank_size+1]
    if len(exclusion_seq) > flank_size*2 + 1:
        lowercase_pos = [i for i, x in enumerate(exclusion_seq) if x.islower()]
        exclusion_seq = exclusion_seq[int(lowercase_pos[0])-flank_size:int(lowercase_pos[0])+flank_size+1]

    #convert sequences to dataframe for plotting
    plt_data = pd.DataFrame({'Inclusion': inclusion_seq, 'Exclusion': exclusion_seq}, index = list(range(-flank_size, flank_size+1)))

    #note where inclusion and exclusion rows differ
    binary_map = plt_data['Inclusion'] != plt_data['Exclusion']
    binary_map = pd.concat([binary_map, binary_map], axis =1).astype(int).T
    binary_map.loc[0] = 0
    plt_data = plt_data.T

    fig, ax = plt.subplots(figsize = figsize)

    sns.heatmap(binary_map, cmap = 'coolwarm', cbar = False, ax = ax, annot = plt_data, fmt = '', annot_kws = {'size': 10}, vmax = 1.1, vmin =-0.1, linewidths = 0.6, ec = 'black', linecolor = 'black')
    ax.set_yticklabels(['Inclusion', 'Exclusion'], rotation = 0)

    #if dpsi is provided, add arrow to plot showing direction of change
    end_of_map = ax.get_xticks()[-1]+1
    if dpsi is not None:
        start_of_arrow = 1.8 if dpsi < 0 else 0.2
        end_of_arrow = 0.2 if dpsi < 0 else 1.8
        ax.annotate(f'', xy = (end_of_map, start_of_arrow), xytext = (end_of_map, end_of_arrow), arrowprops=dict(arrowstyle="->", lw = 1, color = 'black'), annotation_clip=False)

def plot_location_of_altered_flanking_residues(altered_flanks, figsize = (4,3), modification_class = None, residue = None, **kwargs):
    """
    Plot the number of PTMs with altered residues as specific positions relative to the PTM site. This includes the specific position of the residue (-5 to +5 from PTM site) and the specific side of the PTM site that is altered (N-term or C-term)

    Parameters
    ----------
    altered_flanks: pd.DataFrame
        Dataframe with altered flanking sequences, and annotated information added with analyze.compare_flanking_sequences
    figsize: tuple
        Size of the figure. Default is (4,3).
    modification_class: str
        Specific modification class to plot. Default is None, which will plot all modification classes.
    residue: str
        Specific residue to plot. Default is None, which will plot all residues.
    kwargs: additional keyword arguments
        Additional keyword arguments to pass to the filter_ptms function. This allows for filtering based on specific criteria (e.g., dPSI, evidence, etc.).

    """
    
    #filter ptms 
    filter_arguments = helpers.extract_filter_kwargs(**kwargs)
    helpers.check_filter_kwargs(filter_arguments)
    altered_flanks = helpers.filter_ptms(altered_flanks, **filter_arguments)

    
    fig, ax = plt.subplots(nrows = 2, figsize = figsize, height_ratios = [0.5,1])
    fig.subplots_adjust(hspace = 1)

    if modification_class is not None:
        altered_flanks = altered_flanks[altered_flanks['Modification Class'].str.contains(modification_class)].copy()
    
    if residue is not None:
        altered_flanks = altered_flanks[altered_flanks['Residue'] == residue].copy()

    #### plot of side of modification that flank is altered
    terminus = altered_flanks.groupby('Altered Flank Side').size()
    terminus = terminus[['N-term only', 'C-term only']] #focus on cases where only one side is altered for ease of plotting
    ax[0].bar(terminus.index, terminus.values, color = 'gray')
    ax[0].set_xlabel('Location of Altered Region', fontsize = 9)
    ax[0].set_xticklabels(['N-term\nonly', 'C-term\nonly'])
    ax[0].set_ylabel('# of PTMs', fontsize = 9)

    #### plot specific positions of altered residues relative to PTM
    position_breakdown = altered_flanks.explode(['Altered Positions', 'Residue Change']).copy()[['Gene', 'Residue', 'PTM Position in Isoform','Altered Positions', 'Residue Change']]
    position_breakdown = position_breakdown.groupby('Altered Positions').size()
    ax[1].bar(position_breakdown.index, position_breakdown.values, color = 'gray')
    ax[1].set_xlim([-5.5,5.5])
    ax[1].set_xlabel('Position Relative to PTM', fontsize = 9)
    ax[1].set_ylabel('# of Changed\nResidues', fontsize = 9)
    ax[1].set_xticks(np.arange(-5,6,1))


def plot_alterations_matrix(altered_flanks, modification_class = None, residue = None, title = '', ax = None, **kwargs):
    """
    Given the altered flanking sequences dataframe, plot a matrix showing the positions of altered residues for specific proteins, as well as the specific change

    Parameters
    ----------
    altered_flanks: pd.DataFrame
        Dataframe with altered flanking sequences, and annotated information added with analyze.compare_flanking_sequences

    modification_class: str
        Specific modification class to plot. Default is None, which will plot all modification classes.

    residue: str
        Specific residue to plot. Default is None, which will plot all residues.
    title: str
        Title of the plot. Default is '' (no title).
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    kwargs: additional keyword arguments
        Additional keyword arguments to pass to the filter_ptms function. This allows for filtering based on specific criteria (e.g., dPSI, evidence, etc.).
    """
    
    #filter ptms 
    filter_arguments = helpers.extract_filter_kwargs(**kwargs)
    helpers.check_filter_kwargs(filter_arguments)
    altered_flanks = helpers.filter_ptms(altered_flanks, **filter_arguments)

    #extract altered flanking sequences and make sure there is altered position data
    position_breakdown = altered_flanks.copy()
    if 'Altered Positions' not in position_breakdown.columns:
        position_breakdown = compare_flanking_sequences(position_breakdown)

    position_breakdown = position_breakdown.dropna(subset = ['Altered Positions', 'Residue Change'])

    #restrit to desired PTM types and residues
    if modification_class is not None:
        position_breakdown = position_breakdown[position_breakdown['Modification Class'].str.contains(modification_class)].copy()
    if residue is not None:
        position_breakdown = position_breakdown[position_breakdown['Residue'] == residue].copy()

    #add ptm column to position breakdown
    position_breakdown['PTM'] = position_breakdown['Gene'] + '_' + position_breakdown['Residue'] + position_breakdown['PTM Position in Isoform'].astype(str)

    #separate altered residue into individual rows
    position_breakdown = position_breakdown.explode(['Altered Positions', 'Residue Change']).copy()[['Gene', 'PTM','Altered Positions', 'Residue Change']]

    #convert altered positions to integers and remove duplicates
    position_breakdown['Altered Positions'] = position_breakdown['Altered Positions'].astype(int)
    position_breakdown = position_breakdown.drop_duplicates()

    #position_breakdown = position_breakdown.drop_duplicates(subset = ["PTM", "Altered_Positions"], keep = False)
    position_breakdown['PTM']
    position_matrix = position_breakdown.pivot(columns = 'Altered Positions', index = 'PTM', values= 'Residue Change')
    for i, pos in zip(range(11),range(-5, 6)):
        if pos not in position_matrix.columns:
            position_matrix.insert(i, pos, np.nan)


    #replace strings with 1 and nans with 0
    position_values = position_matrix.map(lambda x: 1 if x == x else np.nan).sort_values(by = 5)
    position_matrix = position_matrix.loc[position_values.index]
    #plot heatmap with black edges around cells and no colorbar, annotate with strings in position matrix

    if ax is None:
        fig, ax = plt.subplots(figsize = (4,3))
    sns.heatmap(position_values, cmap = 'Greens', vmin = 0, vmax = 2, ax = ax, cbar = False, linewidths = 0.5, linecolor = 'black', yticklabels=True )
    ax.set_facecolor("lightgrey")

    #annotate with strings in the position matrix
    for i in range(position_values.shape[0]):
        for j in range(position_values.shape[1]):
            if position_values.iloc[i,j] == 1:
                ax.text(j+0.5, i+0.5, position_matrix.iloc[i,j], ha = 'center', va = 'center', fontsize = 6)

    #adjust figure parameters
    ax.set_xticklabels(['-5','-4','-3','-2','-1','PTM','1','2','3','4','5'], fontsize = 8)
    ax.set_xlabel('')
    ax.tick_params(axis = 'y', labelsize = 8)
    ax.set_title(title, fontsize = 9)