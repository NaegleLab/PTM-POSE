r"""
Identify kinases with enriched substrates in differentially included exons, using an adapted version of KSTAR
=============================================================================================================

Given that phosphorlaiton are one of the most commonly impacted modifications, there is potential for kinases targeting these sites to be indirectly impacted by alternative splicing through changes in the availability of their substrates. While we provide functions for performing enrichment of known kinase substrates from databases like PhosphoSitePlus, RegPhos, and PTMsigDB, these resources are limited by the overall number of validated substrates (<5%). For this purpose, we have adapted a previously developed algorithm called KSTAR (Kinase Substrate to Activity Relationships) for use with spliced PTM data, which harnesses kinase-substrate predictions to expand the overall number of phosphorylation sites that can be used as evidence. This particularly important as you may find many of the spliced PTMs in your dataset are less well studied and may not have any annotated kinases.

In order to perform KSTAR analysis, you will first need to download KSTAR networks from the following [figshare](https://figshare.com/articles/dataset/NETWORKS/14944305?file=28768155).

Once you have downloaded the networks, all you need is your PTM data. You will need to run analysis for tyrosine kinases (Y) and serine/threonine kinases (ST)
"""

from ptm_pose import analyze
import pandas as pd

# Load spliced ptm and altered flank data
spliced_ptms = pd.read_csv('spliced_ptms.csv')

#perform kstar enrichment for tyrosine phosphorylation, denoted by "Y"
network_dir = './NetworKIN/'
kstar_enrichment = analyze.kstar_enrichment(spliced_ptms, network_dir = network_dir, phospho_type = 'Y')
kstar_enrichment.run_kstar_enrichment()
kstar_enrichment.return_enriched_kinases()

# %%
# You can also run the same analysis for serine/threonine kinases:
kstar_enrichment = analyze.kstar_enrichment(spliced_ptms, network_dir = network_dir, phospho_type = 'ST')
kstar_enrichment.run_kstar_enrichment()
kstar_enrichment.return_enriched_kinases()