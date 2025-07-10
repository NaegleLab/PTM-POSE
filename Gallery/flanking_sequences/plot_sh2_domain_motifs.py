r"""
Identify altered SH2 domain motifs
===========================================================

A potential consequence of altered flanking sequences is that the short linear motifs important for driving PTM-specific domain interactions are disrupted, such as for the phosphotyrosine binding domains SH2 or 14-3-3 proteins (which bind to phosphoserines and threonines). Using PTM-POSE and motif data from [ELM](http://elm.eu.org/searchdb.html), we can identify and visualize the altered 14-3-3 domain motifs due to splicing events.
"""

from ptm_pose import helpers
from ptm_pose.analyze import flank_analysis as afs

#load example altered flanking sequence data
altered_flanks = helpers.load_example_data(altered_flanks = True)

# %%
# First, we need to identify the linear motifs present for each altered flanking sequence event. We can do this with the `compare_inclusion_motifs` function, which will identify matching any matching motifs for both the inclusion and exclusion flanking sequences. 

altered_flanks = afs.compare_inclusion_motifs(altered_flanks)
altered_flanks[['Gene', 'Residue', 'PTM Position in Isoform', 'Motif only in Inclusion', 'Motif only in Exclusion']].head()

# %%
# We can then identify the instances in which 14-3-3 motifs are altered:

fourteen33_motifs = afs.identify_change_to_specific_motif(altered_flanks, elm_motif_name = '14-3-3', modification_class = 'Phosphorylation', residues = ['S','T'])
fourteen33_motifs[['Gene', 'Residue', 'PTM Position in Isoform', 'Motif only in Inclusion', 'Motif only in Exclusion']]
# %%
# And visualize the differences in sequence
afs.plot_alterations_matrix(fourteen33_motifs)
