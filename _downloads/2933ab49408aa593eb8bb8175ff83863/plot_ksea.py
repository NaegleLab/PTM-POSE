r"""
Kinase substrate enrichment analysis (KSEA)
============================================================================

In the realm of kinase activity analysis, Kinase Substrate Enrichment Analysis (KSEA) is a commonly used tool to identify kinases with substrates that are significantly changing between two conditions (based on phosphoproteomic data). See the publication `here <https://academic.oup.com/bioinformatics/article/33/21/3489/3892392>`_  for what the original algorithm looks like. We have adapted KSEA for use with splicing data, which uses the dPSI values to calculate a z-score to identify kinases with substrates that are either excluded or included from isoforms. This is a useful tool, but you should treat these results with caution as they will often be on only a couple kinase-substrates due to the limited kinase-substrate information available. First, let's load the example data and run KSEA on the differentially included PTMs dataframe

"""

from ptm_pose import helpers
from ptm_pose.analyze import enzyme

spliced_ptms = helpers.load_example_data(spliced_ptms=True)

# %%
# Now we can establish the KSEA class from the enzyme module. We will need to specify the database we want to use for KSEA analysis. This can be individual databases (PhosphoSitePlus, RegPhos, DEPOD, OmniPath, or iKiP) or the combined database (Combined). Note that some of these databases are specific to modification classes (e.g. PhosphoSitePlus is specific to phosphorylation). For this example, let's start with PhosphoSitePlus, a commonly used database for kinase-substrate relationships.

ksea = enzyme.KSEA(spliced_ptms, database = 'PhosphoSitePlus')
#run KSEA for all available kinases
ksea.runKSEA()

# %%
# Now we can visualize the results. The `plot_results()` function will plot the z-scores for each kinase, as well as whether that was found to be significant (based on FDR). You can optionally include notations indicating how many substrates were identified in the splicing data for each kinase.

ksea.plot_results(show_substrate_count = True)


