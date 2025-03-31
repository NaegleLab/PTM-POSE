r"""
Identify enriched gene sets associated with differentially spliced PTMs
=======================================================================

As is commonly done for exon-centric analyses, we have provided the ability to perform gene set enrichment analysis for gene associated with spliced PTMs, using the EnrichR API from the gseapy module. By default, we include gene ontology terms, KEGG pathways, and Reactome pathways, but you can also provide your own gene sets listed in EnrichR.
"""

from ptm_pose import helpers
from ptm_pose.analyze import annotations
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


# Load spliced ptm and altered flank data
spliced_ptms, altered_flanks = helpers.load_example_data(spliced_ptms = True, altered_flanks=True)


# %%
# Use the below function, we can identify enriched gene sets associated with spliced ptms, altered flanks, or both. We can also specify the gene sets to assess, alpha value for significance, the minimum change in PSI value to consider, and whether to return only significant gene sets.
genesets = annotations.gene_set_enrichment(spliced_ptms, altered_flanks, alpha = 0.05, min_dPSI = 0.1, gene_sets = ['GO_Biological_Process_2023','Reactome_2022'], return_sig_only = True)
genesets.head()

# %%
# You can then plot the enriched gene sets, including the proportion of genes associated with differentially included PTMs and those with altered flanking sequences. Here, let's restrict to looking at the top 5 enriched gene sets:

annotations.plot_EnrichR_pies(genesets, top_terms = 5)
plt.tight_layout()
plt.show()