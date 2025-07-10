r"""
Identify kinases with enriched substrates in differentially included exons, using an adapted version of KSTAR
=============================================================================================================

Given that phosphorlaiton are one of the most commonly impacted modifications, there is potential for kinases targeting these sites to be indirectly impacted by alternative splicing through changes in the availability of their substrates. While we provide functions for performing enrichment of known kinase substrates from databases like PhosphoSitePlus, RegPhos, and PTMsigDB, these resources are limited by the overall number of validated substrates (<5%). For this purpose, we have adapted a previously developed algorithm called KSTAR (Kinase Substrate to Activity Relationships) for use with spliced PTM data, which harnesses kinase-substrate predictions to expand the overall number of phosphorylation sites that can be used as evidence. This particularly important as you may find many of the spliced PTMs in your dataset are less well studied and may not have any annotated kinases.

.. note::
    In order to perform KSTAR analysis, you will first need to download KSTAR networks from the following [figshare](https://figshare.com/articles/dataset/NETWORKS/14944305?file=28768155). Once you have downloaded the networks, all you need is your PTM splicing data and the directory location of the networks to run the analysis.
"""

from ptm_pose import helpers
from ptm_pose.analyze import enzyme


#load example differential inclusion data
spliced_ptms = helpers.load_example_data(spliced_ptms = True)

#perform kstar enrichment for tyrosine phosphorylation, denoted by "Y"
network_dir = '../../Gallery/NetworKIN/'
kstar_enrichment = enzyme.kstar_enrichment(spliced_ptms, network_dir = network_dir, phospho_type = 'Y')
kstar_enrichment.run_kstar_enrichment()
kstar_enrichment.return_enriched_kinases()

# %%
# You can also run the same analysis for serine/threonine kinases:
kstar_enrichment = enzyme.kstar_enrichment(spliced_ptms, network_dir = network_dir, phospho_type = 'ST')
kstar_enrichment.run_kstar_enrichment()
kstar_enrichment.return_enriched_kinases()

# %% 
# Finally, you can visualize the results using a KSTAR dotplot, which will show the statistical strength of the relationship based on the size of the dot, with significant kinases colored in orange. 

kstar_enrichment.dotplot(ptype = 'ST', kinase_axis = 'x', impact_types = ['Excluded', 'Included'], size_legend = False, sig_kinases_only = True, figsize = (2,2))