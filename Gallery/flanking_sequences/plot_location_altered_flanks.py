r"""
Probing where and how PTM flanking sequences are altered
===============================================================

In order to understand how PTMs may be altered due to splicing events, it is useful to identify the flanking sequences of the PTMs and how they may be altered due to nearby splice events (as identified by flanking sequence module). Once we have, this information we can analyze and visualize where the alterations in the flanking sequences occur. First, we need to compare the flanking sequences of PTMs based on whether an exonic region is included or excluded using the `compare_flanking_sequences` function in PTM-POSE.
"""

from ptm_pose import helpers
from ptm_pose.analyze import flank_analysis as afs


#load example altered flanking sequence data
altered_flanks = helpers.load_example_data(altered_flanks = True)

altered_flanks = afs.compare_flanking_sequences(altered_flanks)
print('Comparison of flanking sequences:')
altered_flanks[['UniProtKB Accession', 'Residue', 'PTM Position in Isoform', 'Modification Class', 'Inclusion Flanking Sequence', 'Exclusion Flanking Sequence', 'Sequence Identity', 'Altered Positions', 'Residue Change', 'Altered Flank Side']].head()

# %%
# Note, we only calculate these metrics for cases where altered flanking sequences do not cause a stop codon to be introduced, as this is harder to interpret (such as for the first PTM in the list). The above table will indicate the positions in the flanking sequence that are altered, how similar the altered flanking sequence is to the original flanking sequence, and the specific residue change that takes place. We can also plot some of this information to get a better sense of the distribution of altered flanking sequences:

afs.plot_location_of_altered_flanking_residues(altered_flanks)

# %%
# We can even create the same plot for specific modification types or residues, as well as label the specific residue changes that occur:

afs.plot_location_of_altered_flanking_residues(altered_flanks, modification_class='Acetylation')