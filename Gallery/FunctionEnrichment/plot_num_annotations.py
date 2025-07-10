r"""
Identifying available PTM-specific annotation information
=====================================================

As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. However, PTM functional information is inherently sparse, and so most annotations will only provide information on a handful of PTMs. For this reason, it can be useful to probe how many PTMsTo better understand the types of annotations that are available, as well as the number of PTMs that have an annotation of that type. This can be done using the `analyze` function in PTM-POSE.

.. note::
    This example shows the annotations with pre-installed gmt files, but you can quickly create additional gmt files for any databases you might like. See the instructions for annotating PTMs for how to do so.
"""


from ptm_pose import helpers
from ptm_pose.analyze import annotations
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Load spliced ptm data
spliced_ptms = helpers.load_example_data(spliced_ptms = True)

annotations.plot_available_annotations(spliced_ptms)
plt.tight_layout()
plt.show()




# %%
# As you can, see there are only a few PTMs from each annotation that have 
# available information, with the most being 9 PTMs out of the 184 differentially 
# included sites having been associated with a biological process. While this this 
# should be taken into consideration when analyzing these annotations, we can glean 
# some useful information and identify potentially interesting proteins/sites to dig 
# deeper into. 



