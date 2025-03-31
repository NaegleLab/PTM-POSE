r"""
Protein interactions network
============================

Post translational modifications (PTMs) often facilitate protein interactions, either through direct binding of domains specific to that particular modification (e.g. SH2 domains binding to phosphorylated tyrosines) or through allosteric effects that change the conformation of the protein to either enhance or disrupt interactions. We provide functions to annotate spliced PTMs with relevant protein interactions and to identify key PTMs that may disrupt protein interaction networks.

Currently, we provide functions to process and analyze protein interaction data from PhosphoSitePlus, PTMInt, and PTMcode. We can also include enzyme-specific interactions (such as kinase substrate interactions through PhosphoSitePlus and RegPhos). First, we need to annotate the spliced PTMs with protein interactions (see rest of documentation for how to do this). Then, we can process the interactions across the different databases using the protein_interactions class to identify key PTMs that may disrupt protein interaction networks.
"""

from ptm_pose import helpers
from ptm_pose.analyze import interactions
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


#load example differential inclusion data
spliced_ptms = helpers.load_example_data(spliced_ptms = True)

interactions = interactions.protein_interactions(spliced_ptms)
interactions.get_interaction_network()

interactions.network_data.head()

# %%
# We can also calculate interaction stats to identify proteins that are most impacted or relevant to spliced PTMs and the protein interaction network
interactions.get_interaction_stats()

interactions.network_stats.head()

# %%
# If we want to focus on a specific protein, we can summarize information about a single protein in the network. In this case, let's look at TSC2, which loses pS981 upon ESRP1 knockdown

interactions.summarize_protein_network(protein = 'TSC2')

# %% 
# We can also visualize the network...

interactions.plot_interaction_network(interacting_node_size = 10)
plt.tight_layout()
plt.show() 

# %%
# ...and the centrality of proteins in the network

interactions.plot_network_centrality(centrality_measure='Degree')
plt.tight_layout()
plt.show() 