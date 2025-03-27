import numpy as np
import pandas as pd

#plotting 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

#analysis packages
import networkx as nx



#custom stat functions
from ptm_pose.analyze import annotations
from ptm_pose import annotate, helpers


def get_edge_colors(interaction_graph, network_data, defaultedgecolor = 'gray', color_edges_by = 'Database', database_color_dict = {'PSP/RegPhos':'red', 'PhosphoSitePlus':'green', 'PTMcode':'blue', 'PTMInt':'gold', 'Multiple':'purple'}):
    edge_colors = []
    legend_handles = {}
    for edge in interaction_graph.edges:
        trim = network_data[(network_data['Modified Gene'] == edge[0]) & (network_data['Interacting Gene'] == edge[1])]
        if trim.shape[0] == 0:
            trim = network_data[(network_data['Interacting Gene'] == edge[0]) & (network_data['Modified Gene'] == edge[1])]

        if color_edges_by == 'Database':
            if trim.shape[0] > 1: #specific color to indicate multiple databases
                edge_colors.append(database_color_dict['Multiple'])
                if 'Multiple Databases' not in legend_handles:
                    legend_handles['Multiple Databases'] = mlines.Line2D([0], [0], color = 'purple', label = 'Multiple Databases')
            elif ';' in trim['Source'].values[0]:
                edge_colors.append(database_color_dict['Multiple'])
                if 'Multiple Databases' not in legend_handles:
                    legend_handles['Multiple Databases'] = mlines.Line2D([0], [0], color = 'purple', label = 'Multiple Databases')
            elif trim["Source"].values[0] in database_color_dict:
                edge_colors.append(database_color_dict[trim["Source"].values[0]])
                if trim["Source"].values[0] not in legend_handles:
                    legend_handles[trim["Source"].values[0]] = mlines.Line2D([0], [0], color = database_color_dict[trim["Source"].values[0]], label = trim["Source"].values[0])
            else: #gray means not specific to database in list
                edge_colors.append(defaultedgecolor)
        else:
            edge_colors.append(defaultedgecolor)
            legend_handles = None

    if isinstance(legend_handles, dict):
        legend_handles = list(legend_handles.values())
        
    return edge_colors, legend_handles


def plot_interaction_network(interaction_graph, network_data, network_stats = None, modified_color = 'red', modified_node_size = 10, interacting_color = 'lightblue', interacting_node_size = 1, defaultedgecolor = 'gray', color_edges_by = 'Same', seed = 200, legend_fontsize = 8, ax = None, proteins_to_label = None, labelcolor = 'black', legend = True):
    """
    Given the interaction graph and network data outputted from analyze.protein_interactions, plot the interaction network, signifying which proteins or ptms are altered by splicing and the specific regulation change that occurs. by default, will only label proteins 

    Parameters
    ----------
    interaction_graph: nx.Graph
        NetworkX graph object representing the interaction network, created from analyze.get_interaction_network
    network_data: pd.DataFrame
        Dataframe containing details about specifici protein interactions (including which protein contains the spliced PTMs)
    network_stats: pd.DataFrame
        Dataframe containing network statistics for each protein in the interaction network, obtained from analyze.get_interaction_stats(). Default is None, which will not label any proteins in the network.
    modified_color: str
        Color to use for proteins that are spliced. Default is 'red'.
    modified_node_size: int
        Size of nodes that are spliced. Default is 10.
    interacting_color: str
        Color to use for proteins that are not spliced. Default is 'lightblue'.
    interacting_node_size: int
        Size of nodes that are not spliced. Default is 1.
    edgecolor: str
        Color to use for edges in the network. Default is 'gray'.
    seed: int
        Seed to use for spring layout of network. Default is 200.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    proteins_to_label: list, int, or str
        Specific proteins to label in the network. If list, will label all proteins in the list. If int, will label the top N proteins by degree centrality. If str, will label the specific protein. Default is None, which will not label any proteins in the network.
    labelcolor: str
        Color to use for labels. Default is 'black'.
    """
    node_colors = []
    node_sizes = []
    for node in interaction_graph.nodes:
        if node in network_data['Modified Gene'].unique():
            node_colors.append(modified_color)
            node_sizes.append(modified_node_size)
        else:
            node_colors.append(interacting_color)
            node_sizes.append(interacting_node_size)

    if 'Regulation Change' in network_data.columns:
        #adjust line style of edge depending on sign of deltaPSI_MW
        edge_style = []
        for edge in interaction_graph.edges:
            edge_data = network_data[((network_data['Modified Gene'] == edge[0]) & (network_data['Interacting Gene'] == edge[1])) | ((network_data['Modified Gene'] == edge[1]) & (network_data['Interacting Gene'] == edge[0]))]
            if '+' in edge_data['Regulation Change'].values[0] and '-' in edge_data['Regulation Change'].values[0]:
                edge_style.append('dashdot')
            elif '+' in edge_data['Regulation Change'].values[0]:
                edge_style.append('solid')
            else:
                edge_style.append('dotted')
    else:
        edge_style = 'solid'

    np.random.seed(seed)
    interaction_graph.pos = nx.spring_layout(interaction_graph, seed = seed)

    #set up subplot if not provied
    if ax is None:
        fig, ax = plt.subplots(figsize = (4,4), layout='constrained')

    edge_colors, edge_legend = get_edge_colors(interaction_graph, network_data, color_edges_by = color_edges_by, defaultedgecolor = defaultedgecolor)

    nx.draw(interaction_graph, node_size = node_sizes, node_color = node_colors, edge_color = edge_colors, style = edge_style, ax = ax)

    #add legend for colored nodes
    if legend:
        modified_node = mlines.Line2D([0], [0], color='w',marker = 'o', markersize=modified_node_size,linewidth = 0.2, markerfacecolor = modified_color, markeredgecolor=modified_color, label='Spliced Protein')
        interacting_node = mlines.Line2D([0], [0], color='w', markerfacecolor = interacting_color, markeredgecolor=interacting_color, marker = 'o', markersize=interacting_node_size, linewidth = 0.2, label='Interacting Protein')
        solid_line = mlines.Line2D([0], [0], color='gray', linestyle = 'solid', label = 'Interaction increases')
        dashdot_line = mlines.Line2D([0], [0], color='gray', linestyle = 'dashdot', label = 'Interaction impact unclear')
        dotted_line = mlines.Line2D([0], [0], color='gray', linestyle = 'dotted', label = 'Interaction decreases')
        handles = [solid_line,dashdot_line, dotted_line, modified_node, interacting_node]
        net_legend = ax.legend(handles = handles, loc = 'upper center', ncol = 2, fontsize = legend_fontsize, bbox_to_anchor = (0.5, 1.1))
        ax.add_artist(net_legend)

        if color_edges_by == 'Database':
            ax.legend(handles = edge_legend, loc = 'lower center', ncol = 1, fontsize = legend_fontsize, bbox_to_anchor = (0.5, -0.15), title = 'Interaction Source')


    #if requested, label specific proteins in the network
    if proteins_to_label is not None and isinstance(proteins_to_label, list):
        for protein in proteins_to_label:
            ax.text(interaction_graph.pos[protein][0], interaction_graph.pos[protein][1], protein, fontsize = 10, fontweight = 'bold', color = labelcolor)
    elif proteins_to_label is not None and isinstance(proteins_to_label, int):
        if network_stats is None:
            network_stats = get_interaction_stats(interaction_graph)
        
        network_stats = network_stats.sort_values(by = 'Degree', ascending = False).iloc[:proteins_to_label]
        for index, row in network_stats.iterrows():
            ax.text(interaction_graph.pos[index][0], interaction_graph.pos[index][1], index, fontsize = 10, fontweight = 'bold', color = labelcolor)
    elif proteins_to_label is not None and isinstance(proteins_to_label, str):
        ax.text(interaction_graph.pos[proteins_to_label][0], interaction_graph.pos[proteins_to_label][1], proteins_to_label, fontsize = 10, fontweight = 'bold', color = labelcolor)
    elif proteins_to_label is not None:
        print('Proteins to label must be a list of strings or a single string. Ignoring when plotting.')
    
def plot_network_centrality(network_stats, network_data = None, centrality_measure = 'Degree', top_N = 10, modified_color = 'red', interacting_color = 'black', ax = None):
    """
    Given the network statistics data obtained from analyze.get_interaction_stats(), plot the top N proteins in the protein interaction network based on centrality measure (Degree, Betweenness, or Closeness)

    Parameters
    ----------
    network_stats: pd.DataFrame
        Dataframe containing network statistics for each protein in the interaction network, obtained from analyze.get_interaction_stats()
    network_data: pd.DataFrame
        Dataframe containing information on which proteins are spliced and how they are altered. Default is None, which will plot all proteins the same color (interacting_color)
    centrality_measure: str
        Centrality measure to use for plotting. Default is 'Degree'. Options include 'Degree', 'Degree Centrality', 'Betweenness Centrality', 'Closeness Centrality'.
    top_N: int
        Number of top proteins to plot. Default is 10.
    modified_color: str
        Color to use for proteins that are spliced. Default is 'red'.
    interacting_color: str
        Color to use for proteins that are not spliced. Default is 'black'.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    
    Outputs
    -------
    bar plot showing the top N proteins in the interaction network based on centrality measure
    """
    if centrality_measure not in network_stats.columns:
        raise ValueError('Centrality measure not found in network_stats dataframe. Please check the inputted centrality_measure. Available measures include Degree, Degree Centrality, Betweenness Centrality, Closeness Centrality, and Eigenvector Centrality.')
    
    #get specific centrality measure and grab top N terms
    plt_data = network_stats.sort_values(by = centrality_measure, ascending = False).iloc[:top_N].sort_values(by = centrality_measure, ascending = True)
    
    #color bars based on whether protein is spliced or not
    if network_data is not None:
        colors = []
        for index, row in plt_data.iterrows():
            if index in network_data['Modified Gene'].unique():
                colors.append(modified_color)
            else:
                colors.append(interacting_color)
    else:
        colors = modified_color
    
    #establish figure
    if ax is None:
        fig, ax = plt.subplots(figsize = (3,3))

    #plot bar plot
    ax.barh(plt_data.index, plt_data[centrality_measure], color = colors)
    ax.set_xlabel(f'{centrality_measure}')


def get_interaction_stats(interaction_graph):
    """
    Given the networkx interaction graph, calculate various network centrality measures to identify the most relevant PTMs or genes in the network
    """
    #calculate network centrality measures
    degree_centrality = nx.degree_centrality(interaction_graph)
    closeness_centrality = nx.closeness_centrality(interaction_graph)
    betweenness_centrality = nx.betweenness_centrality(interaction_graph)
    network_stats = pd.DataFrame({'Degree': dict(interaction_graph.degree()), 'Degree Centrality':degree_centrality, 'Closeness':closeness_centrality,'Betweenness':betweenness_centrality})
    return network_stats

class protein_interactions:
    """
    Class to assess interactions facilitated by PTMs in splicing network

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs projected onto splicing events and with annotations appended from various databases
    include_enzyme_interactions: bool
        Whether to include interactions with enzymes in the network, such as kinase-substrate interactions. Default is True
    interaction_databases: list
        List of databases whose information to include in the network. Default is ['PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD', 'ELM']
    **kwargs: additional keyword arguments
        Additional keyword arguments, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.

    Attributes
    ----------
    interaction_graph: nx.Graph
        NetworkX graph object representing the interaction network, created from analyze.get_interaction_network
    network_data: pd.DataFrame
        Dataframe containing details about specifici protein interactions (including which protein contains the spliced PTMs)
    network_stats: pd.DataFrame
        Dataframe containing network statistics for each protein in the interaction network, obtained from analyze.get_interaction_stats(). Default is None, which will not label any proteins in the network.
    """
    def __init__(self, spliced_ptms, include_enzyme_interactions = True, interaction_databases = ['PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD', 'OmniPath'], **kwargs):
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        helpers.check_filter_kwargs(filter_arguments)
        spliced_ptms = helpers.filter_ptms(spliced_ptms, **filter_arguments)
        if spliced_ptms.empty:
            raise ValueError('No spliced PTMs found in provided data after filtering PTMs and events. Consider reducing filter stringency.')
        
        self.include_enzyme_interactions = include_enzyme_interactions
        self.interaction_databases = interaction_databases

        #go through and add interaction information for any database that are not appended to the spliced_ptms dataframe
        #available_annotations = annotations.get_available_annotations(spliced_ptms)
        #available_interactions = available_annotations[available_annotations['Annotation Type'] == 'Interactions']
        #if self.include_enzyme_interactions:
        #    enzyme_annotations = available_annotations[available_annotations['Annotation Type'] == 'Enzyme']


        #go through and add interaction information that is not in the spliced_ptms dataframe
        #for db in self.interaction_databases:
        #    if f'{db}:Interactions' not in spliced_ptms.columns and db in available_interactions['Database'].values:
        #        spliced_ptms = annotate.append_from_gmt(spliced_ptms, database = db, annot_type = 'Interactions')

        #    if f'{db}:Enzyme' not in spliced_ptms.columns and self.include_enzyme_interactions and db in enzyme_annotations['Database'].values:
        #        spliced_ptms = annotate.append_from_gmt(spliced_ptms, database = db, annot_type = 'Enzyme')

        self.spliced_ptms = spliced_ptms
        self.get_interaction_network()


    def get_interaction_network(self, node_type = 'Gene'):
        """
        Given the spliced PTM data, extract interaction information and construct a dataframe containing all possible interactions driven by PTMs, either centered around specific PTMs or specific genes.

        Parameters
        ----------
        node_type: str
            What to define interactions by. Can either be by 'PTM', which will consider each PTM as a separate node, or by 'Gene', which will aggregate information across all PTMs of a single gene into a single node. Default is 'Gene'
        """
        if node_type not in ['Gene', 'PTM']:
            raise ValueError("node_type parameter (which dictates whether to consider interactions at PTM or gene level) can be either Gene or PTM")
        
        #extract interaction information in provided data
        interactions = annotate.combine_interaction_data(self.spliced_ptms, include_enzyme_interactions=self.include_enzyme_interactions, interaction_databases=self.interaction_databases)
        if interactions.empty:
            raise ValueError('No interaction data found in spliced PTM data.')
        interactions['Residue'] = interactions['Residue'] + interactions['PTM Position in Isoform'].astype(int).astype(str)
        interactions = interactions.drop(columns = ['PTM Position in Isoform'])

        #add regulation change information
        if 'dPSI' in self.spliced_ptms.columns:
            interactions['Regulation Change'] = interactions.apply(lambda x: '+' if x['Type'] != 'DISRUPTS' and x['dPSI'] > 0 else '+' if x['Type'] == 'DISRUPTS' and x['dPSI'] < 0 else '-', axis = 1)
            grouping_cols = ['Residue', 'Type', 'Source', 'dPSI', 'Regulation Change']
            interactions['dPSI'] = interactions['dPSI'].apply(str)
        else:
            grouping_cols = ['Residue', 'Type', 'Source']

        #extract gene_specific network information
        if node_type == 'Gene':
            network_data = interactions.groupby(['Modified Gene', 'Interacting Gene'], as_index = False)[grouping_cols].agg(helpers.join_unique_entries)
            #generate network with all possible PTM-associated interactions
            interaction_graph = nx.from_pandas_edgelist(network_data, source = 'Modified Gene', target = 'Interacting Gene')
        else:
            interactions['Spliced PTM'] = interactions['Modified Gene'] + '_' + interactions['Residue']
            network_data = interactions.groupby(['Spliced PTM', 'Interacting Gene'], as_index = False)[grouping_cols].agg(helpers.join_unique_entries)
            network_data = network_data.drop(columns = ['Residue'])
            
            #generate network with all possible PTM-associated interactions
            interaction_graph = nx.from_pandas_edgelist(network_data, source = 'Spliced PTM', target = 'Interacting Gene')

        self.network_data = network_data
        self.interaction_graph = interaction_graph


    def get_interaction_stats(self):
        """
        Given the networkx interaction graph, calculate various network centrality measures to identify the most relevant PTMs or genes in the network
        """
        #calculate network centrality measures
        self.network_stats = get_interaction_stats(self.interaction_graph)

    def get_protein_interaction_network(self, protein):
        """
        Given a specific protein, return the network data for that protein

        Parameters
        ----------
        protein: str
            Gene name of the protein of interest

        Returns
        -------
        protein_network: pd.DataFrame
            Dataframe containing network data for the protein of interest
        """
        if not hasattr(self, 'network_data'):
            self.get_interaction_network()

        if protein not in self.network_data['Modified Gene'].unique():
            print(f'{protein} is not found in the network data. Please provide a valid gene name.')
            return None
        
        protein_network = self.network_data[self.network_data['Modified Gene'] == protein]
        protein_network = protein_network.drop(columns = ['Modified Gene'])
        protein_network = protein_network.rename(columns = {'Residue': 'Spliced PTMs facilitating Interacting'})
        return protein_network

    def summarize_protein_network(self, protein):
        """
        Given a protein of interest, summarize the network data for that protein

        Parameters
        ----------
        protein: str
            Gene name of the protein of interest
        """
        if not hasattr(self, 'network_data'):
            self.get_interaction_network()

        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        protein_network = self.network_data[self.network_data['Modified Gene'] == protein]
        increased_interactions = protein_network.loc[protein_network['Regulation Change'] == '+', 'Interacting Gene'].values
        decreased_interactions = protein_network.loc[protein_network['Regulation Change'] == '-', 'Interacting Gene'].values
        ambiguous_interactions = protein_network.loc[protein_network['Regulation Change'].str.contains(';'), 'Interacting Gene'].values

        #print interactions
        if len(increased_interactions) > 0:
            print(f"Increased interaction likelihoods: {', '.join(increased_interactions)}")
        if len(decreased_interactions) > 0:
            print(f"Decreased interaction likelihoods: {', '.join(decreased_interactions)}")
        if len(ambiguous_interactions) > 0:
            print(f"Ambiguous interaction impact: {', '.join(ambiguous_interactions)}")

        network_ranks = self.network_stats.rank(ascending = False).astype(int)
        print(f'Number of interactions: {self.network_stats.loc[protein, "Degree"]} (Rank: {network_ranks.loc[protein, "Degree"]})')
        print(f'Centrality measures - \t Degree = {self.network_stats.loc[protein, "Degree Centrality"]} (Rank: {network_ranks.loc[protein, "Degree Centrality"]})')
        print(f'                      \t Betweenness = {self.network_stats.loc[protein, "Betweenness"]} (Rank: {network_ranks.loc[protein, "Betweenness"]})')
        print(f'                      \t Closeness = {self.network_stats.loc[protein, "Closeness"]} (Rank: {network_ranks.loc[protein, "Closeness"]})')

    def compare_to_nease(self, nease_edges):
        """
        Given the network edges generated by NEASE, compare the network edges generated by NEASE to the network edges generated by the PTM-driven interactions
        """
        if not hasattr(self, 'interaction_graph'):
            self.get_interaction_network()

        nease_edges['Affected binding'] = nease_edges['Affected binding'].apply(lambda x: x.split(','))
        nease_edges = nease_edges.explode('Affected binding')

        nease_nw = nx.from_pandas_edgelist(nease_edges, source = 'Gene name', target = 'Affected binding')

        #construct graphs with only overlapping edges or unique edges
        common_graph = nx.intersection(nease_nw, self.interaction_graph)

        #nease only graph
        nease_only_graph = nease_nw.copy()
        nease_only_graph.remove_edges_from(e for e in nease_nw.edges() if e in self.interaction_graph.edges())

        #ptm only graph
        ptm_only_graph = self.interaction_graph.copy()
        ptm_only_graph.remove_edges_from(e for e in self.interaction_graph.edges() if e in self.nease_nw.edges())

        #convert ptm_only_graph edges to pandas dataframe
        ptm_edges = nx.to_pandas_edgelist(ptm_only_graph)
        ptm_edges['PTM-POSE'] = True
        ptm_edges['NEASE'] = False

        #convert nease_only_graph edges to pandas dataframe
        nease_edges = nx.to_pandas_edgelist(nease_only_graph)
        nease_edges['NEASE'] = True
        nease_edges['PTM-POSE'] = False

        #convert common_graph edges to pandas dataframe
        common_edges = nx.to_pandas_edgelist(common_graph)
        common_edges['NEASE'] = True
        common_edges['PTM-POSE'] = True

        #combine all edges into one dataframe
        combined_edges = pd.concat([ptm_edges, nease_edges, common_edges])
        self.nease_comparison = combined_edges

    def plot_nease_comparison(self, ax = None, **kwargs):
        """
        Given the comparison data generated by compare_to_nease, plot the number of edges identified by NEASE and PTM-POSE
        """
        #check if nease_comparison exists already, if not generate if data is provided
        if not hasattr(self, 'nease_comparison') and 'nease_edges' in kwargs:
            self.compare_to_nease(kwargs['nease_edges'])
        elif not hasattr(self, 'nease_comparison'):
            raise ValueError('No NEASE comparison data found. Please provide nease_edges to compare to PTM-POSE network or run `compare_to_nease()` first.')
        
        
        if ax is None:
            fig, ax = plt.subplots(figsize = (3,2))
        
        #extract counts
        common_edges = self.nease_comparison[self.nease_comparison['NEASE'] & self.nease_comparison['PTM-POSE']].shape[0]
        ptm_only_edges = self.nease_comparison[self.nease_comparison['PTM-POSE'] & ~self.nease_comparison['NEASE']].shape[0]
        nease_only_edges = self.nease_comparison[self.nease_comparison['NEASE'] & ~self.nease_comparison['PTM-POSE']].shape[0]

        #plot barplot
        ax.barh(['Identified with NEASE\nand PTM-POSE', 'PTM-POSE\n(PTM-Specific)','NEASE\n(Exon-Specific)'], [common_edges, ptm_only_edges, nease_only_edges], color = 'gray', edgecolor = 'black', height=0.9)
        ax.tick_params(labelsize = 9)
        ax.set_xlabel('Number of Interactions\nImpacted By Splicing', fontsize = 10)



    def plot_interaction_network(self, modified_color = 'red', modified_node_size = 10, interacting_color = 'lightblue', interacting_node_size = 5, defaultedgecolor = 'gray', color_edges_by = 'Same', seed = 200, ax = None, proteins_to_label = None, labelcolor = 'black', legend = True, include_nease_comparison = False):
        """
        Given the interactiong graph and network data outputted from analyze.get_interaction_network, plot the interaction network, signifying which proteins or ptms are altered by splicing and the specific regulation change that occurs. by default, will only label proteins 

        Parameters
        ----------
        modified_color: str
            Color to use for nodes that are modified by splicing. Default is 'red'
        modified_node_size: int
            Size of nodes that are modified by splicing. Default is 10
        interacting_color: str
            Color to use for nodes that interact with modified nodes. Default is 'lightblue'
        interacting_node_size: int
            Size of nodes that interact with modified nodes. Default is 5
        edgecolor: str
            Color to use for edges in the network. Default is 'gray'. Can choose to color by database by providing a dictionary with database names as keys and colors as values or by specifying 'database' as color
        seed: int
            Seed to use for random number generator. Default is 200
        ax: matplotlib.pyplot.Axes
            Axes object to plot the network. Default is None
        """
        if not hasattr(self, 'interaction_graph'):
            self.get_interaction_network()

        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        # if include nease comparison, combine graphs
        if nease_edges is not None:
            nease_edges['Affected binding'] = nease_edges['Affected binding'].apply(lambda x: x.split(','))
            nease_edges = nease_edges.explode('Affected binding')

            nease_nw = nx.from_pandas_edgelist(nease_edges, source = 'Gene name', target = 'Affected binding')

            #combine graphs
            plt_graph = nx.compose(self.interaction_graph, nease_nw)
        else:
            plt_graph = self.interaction_graph.copy()


        plot_interaction_network(plt_graph, self.network_data, self.network_stats, modified_color = modified_color, modified_node_size = modified_node_size, interacting_color = interacting_color, interacting_node_size = interacting_node_size, defaultedgecolor = defaultedgecolor, color_edges_by=color_edges_by, seed = seed, ax = ax, proteins_to_label = proteins_to_label, labelcolor = labelcolor, legend = legend)

    def plot_network_centrality(self,  centrality_measure = 'Degree', top_N = 10, modified_color = 'red', interacting_color = 'black', ax = None):
        if not hasattr(self, 'interaction_graph'):
            self.get_interaction_network()
        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        plot_network_centrality(self.network_stats, self.network_data, centrality_measure=centrality_measure,top_N = top_N, modified_color = modified_color, interacting_color = interacting_color, ax = ax)
