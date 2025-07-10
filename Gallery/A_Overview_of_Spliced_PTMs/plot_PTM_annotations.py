r"""
Inspecting number of PTMs with annotation information available
=============================================================================================================

As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. However, PTM functional information is inherently sparse, and so most annotations will only provide information on a handful of PTMs. For this reason, it can be useful to probe how many PTMsTo better understand the types of annotations that are available, as well as the number of PTMs that have an annotation of that type. This can be done using the analyze function in PTM-POSE.

Note: This examples assumes that you have already run the PTM-POSE pipeline and have at annotated PTMs with at least one layer of information.
"""

from ptm_pose import helpers
from ptm_pose.analyze import annotations
import matplotlib.pyplot as plt

#load example differential inclusion data
spliced_ptms = helpers.load_example_data(spliced_ptms = True)

# %% You can first look at what annotations are available, first in a table format
available_annotations = annotations.get_available_annotations(spliced_ptms)
available_annotations

# %%
# You can also visualize the number of PTMs with annotation information available for each annotation type using a bar plot.

annotations.plot_available_annotations(spliced_ptms)
plt.tight_layout()
plt.show() 

# %%
# As you can, see there are only a few PTMs from each annotation that have available information, with the most being 9 PTMs out of the 184 differentially included sites having been associated with a biological process. While this this should be taken into consideration when analyzing these annotations, we can glean some useful information and identify potentially interesting proteins/sites to dig deeper into. Let’s look at the PTMs that have been associated with a biological process:

ptms_with_annotation, annotation_counts = annotations.get_ptm_annotations(spliced_ptms, database = "PhosphoSitePlus", annotation_type = 'Process')
print('Specific PTMs with annotation:')
ptms_with_annotation

# %%
# We can also look at the number of PTMs associated with each annotation:

print('Number of PTMs associated with each annotation:')
annotation_counts

# %%
# Note: you could also do this analysis for altered flanking sequences by replacing `spliced_ptms` with `altered_flanks` in the above code.
