r"""
Kinase affinity due to PTM flanking sequence alterations
============================================================================

Kinases are known to recognize specific motifs surrounding Y, S, and T residues, and these motifs can be used to predict kinase interactions using a tool called Kinase Library. Here, we will use Kinase Library to identify kinases who may have a stronger preference for the flanking sequences of PTMs in one isoform vs. another. For this example, let's focus on phosphorylation sites that have larger changes >0.4 in the splicing data.

.. note::
    The optional dependency `kinase-library` is required to run this example. You can install it using `pip install kinase-library`.
"""

from ptm_pose import helpers
from ptm_pose.analyze import enzyme

altered_flanks = helpers.load_example_data(altered_flanks=True)

# initialize the kinase library class
klibrary = enzyme.KL_flank_analysis(altered_flanks, min_dpsi = 0.4)

# %%
# Now that we have the class initialized, we have two options. If there is a specific PTM of interest, we can run the analysis on that specific PTM. Otherwise, we can run the analysis on all PTMs in the dataframe. Let's analyzing across all PTMs first. This will return a dataframe containing scores of each kinase for every event in the altered_flanks dataframe for both the inclusion and exclusion isoforms. In addition, it will calculate the difference and the relative difference (normalized by dPSI) between the two isoforms:
# :


klibrary.analyze_all_ptms()


# %%
# With all the scores calculated, we can then identify the sites and events with the largest changes in affinity. To better represent how big of a change this can be, we can calculate a relative affinity change, which is the difference in percentile score change between inclusion and exclusion isoforms, normalized by the change in dPSI. This will give us a better idea of how much the affinity is changing relative to the change in splicing. We can then plot this data to visualize the results.

klibrary.plot_top_changes()

# %%
# Here, we can see that a site in CD44 has a shift towards interactions with CAMK2G upon the provided perturbation, and our US01 site has several different kinases which prefer the isoform expressed prior to perturbation. In addition to looking at individual sites, we can also see if there are any kinases with consistently large differences in scores/percentiles across the flanking sequence changes. This can be done using the `plot_top_kinases()` function, which will plot the top kinases with the largest differences in scores across all PTMs (based on median by default). This can be useful for identifying kinases that may be more broadly impacted by the splicing events.

klibrary.plot_top_kinases(top_n = 5)

# %%
# Shown are the top 5 kinases with the largest median difference across all assessed PTMs. You'll see that many of the PTMs don't result in a large change, but some kinases, like MEK5, have a large change in preference after perturbation for some sites.

# %%
# Now let's say we didn't want to focus on a specific PTM, such as the USO1 S486 site. We can do this by using the `analyze_single_ptm()` function, which will return a dataframe with the scores of each kinase for the inclusion and exclusion isoforms, as well as the difference and relative difference between the two isoforms for the specified PTM. First, let's take a look at what the flanking sequences look like for this PTM:


gene = 'USO1'
loc = 486
example = altered_flanks[(altered_flanks['Gene'] == gene) & (altered_flanks['PTM Position in Isoform'] == loc)].squeeze()

from ptm_pose.analyze import flank_analysis

flank_analysis.plot_sequence_differences(example['Inclusion Flanking Sequence'], example['Exclusion Flanking Sequence'])

# %%
# Now we can run the kinase library analysis on this specific PTM. This will return a dataframe with the kinases that have a higher affinity for the inclusion or exclusion isoform, as well as the p-value and fold change.

affinity_change = klibrary.analyze_single_ptm('ARHGAP17', 497)
print(affinity_change.head(10))


# %%
# From this analysis, we can start to identify events that might rewire kinase interactions. When paired with analysis of differentially included PTMs, we can start to identify kinases that are more or less likely to be influenced by changes to splicing patterns.



