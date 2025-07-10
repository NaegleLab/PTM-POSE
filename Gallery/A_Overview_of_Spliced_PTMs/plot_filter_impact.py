r"""
Impact of filtering PTMs
=============================================================================================================

In some cases, you may wish to filter out PTMs with weaker evidence (only identified in small number of experiments), splice events with small changes, or PTMs that have not previously been observed to be impacted by splicing. This can be done using the `filter_ptms` function in the `helpers` module. We can also visualize the impact of filtering using the `analyze.filter` module. 
"""
from ptm_pose import helpers
from ptm_pose.analyze import filter
import matplotlib.pyplot as plt


#load example differential inclusion and altered flank data
spliced_ptms = helpers.load_example_data(spliced_ptms = True)

# %%
# There are several ways to filter PTMs based on the parameter values below. Let's set these parameters so that we focus on PTMs that have been observed to impacted by splicing (remove_novel = True) and have been observed in at least 2 mass spectrometry experiments (min_experiments = 2). By default, the filter ptms will also remove any insignificant splice events with a dPSI value less than 0.1 (min_dpsi = 0.1, alpha = 0.05), but let's restrict that a little more to focus on events with changes of at least 20%.

remove_novel = True
min_MS_observations = 5
min_dpsi = 0.2

#other paramter values we could set to filter PTMs
alpha = 0.05 # p-value threshold for significance
min_studies = 0 # minimum number of literature publications that support the PTM (high and low throughput)


# %% 
# We can now assess the impact of these parameters on the number of PTMs, and type of PTMs, that are present in the dataset after filtering using the `filter.plot_filter_impact()` function

filter.plot_filter_impact(spliced_ptms, output_type = 'count', remove_novel = remove_novel, min_MS_observations = min_MS_observations, min_dpsi = min_dpsi, report_removed = False)

# %%
# Rather than the total number of PTMs, we can also assess how it impacts the proportion of PTMs that are present in the dataset. This is useful to see how the filtering parameters impact the type of PTMs that are present in the dataset, as more restrictive parameters may skew the dataset towards a particular type of PTM (e.g. phosphorylation).

filter.plot_filter_impact(spliced_ptms, output_type = 'fraction', remove_novel = remove_novel, min_MS_observations = min_MS_observations, min_dpsi = min_dpsi, report_removed = False)

# %%
# As you can see, this a pretty restrictive filter, and we are left with only a small number of PTMs (mostly phosphorylation). That may be useful for some analyses, but we may want to relax the parameters a little to include more PTMs. We could also choose to only filter out phosphoryation sites based on evidence, and keep all other PTMs regardless of the evidence (e.g. acetylation, methylation, etc.). This is not perfect, but as phosphorylation is the most commonly studied/measured PTM, it may be a good compromise to keep other PTMs that are not as well studied.

phospho_only_evidence_filter = True

filter.plot_filter_impact(spliced_ptms, output_type = 'count', remove_novel = remove_novel, min_MS_observations = min_MS_observations, min_dpsi = min_dpsi, report_removed = False, phospho_only_evidence_filter = phospho_only_evidence_filter)

# %% 
# Now you can see that there is the same number of phosphorylation sites, but now there are other modification types still present. 