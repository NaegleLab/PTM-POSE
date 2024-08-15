r"""
Inspecting number of PTMs with annotation information available
===============================================================

As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. However, PTM functional information is inherently sparse, and so most annotations will only provide information on a handful of PTMs. For this reason, it can be useful to probe how many PTMsTo better understand the types of annotations that are available, as well as the number of PTMs that have an annotation of that type. This can be done using the `analyze` function in PTM-POSE.

Note: This examples assumes that you have already run the PTM-POSE pipeline and have at annotated PTMs with at least one layer of information.
"""


from ptm_pose import analyze
import pandas as pd

# Load spliced ptm and altered flank data
spliced_ptms = pd.read_csv('spliced_ptms.csv')
altered_flanks = pd.read_csv('altered_flanks.csv')

analyze.show_available_annotations(spliced_ptms, figsize = (5, 3))




# %%
# As you can, see there are only a few PTMs from each annotation that have 
# available information, with the most being 9 PTMs out of the 184 differentially 
# included sites having been associated with a biological process. While this this 
# should be taken into consideration when analyzing these annotations, we can glean 
# some useful information and identify potentially interesting proteins/sites to dig 
# deeper into. Let's look at the PTMs that have been associated with a biological
# process: 
ptms_with_annotation, annotation_counts = analyze.get_ptm_annotations(spliced_ptms, database = "PhosphoSitePlus", annotation_type = 'Process')
print('Specific PTMs with annotation:')
ptms_with_annotation

# %%
# We can also look at the number of PTMs associated with each annotation:
print('Number of PTMs associated with each annotation:')
annotation_counts


