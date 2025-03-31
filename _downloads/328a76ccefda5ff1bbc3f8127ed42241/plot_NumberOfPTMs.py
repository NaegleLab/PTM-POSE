r"""
Inspecting types of PTMs impacted by splicing
=============================================================================================================

One of the first things you will want to do when analyzing spliced PTMs is to check how many and what type of PTMs were identified as being differentially included or altered due to splicing events. This can be done easily using the plotting module in PTM-POSE.
"""
from ptm_pose import helpers
from ptm_pose.analyze import summarize
import matplotlib.pyplot as plt


#load example differential inclusion and altered flank data
spliced_ptms, altered_flanks = helpers.load_example_data(spliced_ptms = True, altered_flanks = True)

summarize.plot_modification_breakdown(spliced_ptms = spliced_ptms, altered_flanks = altered_flanks)
plt.tight_layout()
plt.show() 