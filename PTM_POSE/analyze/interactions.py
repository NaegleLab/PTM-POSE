import numpy as np
import pandas as pd

#plotting 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

#analysis packages
import networkx as nx



#custom stat functions
from ptm_pose import annotate, helpers, pose_config




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

        Parameters
        ----------
        nease_edges : pd.DataFrame
            Interactions found by NEASE, which is the output of nease.get_edges() function after running NEASE (see nease_runner module for example)
        
        Returns
        -------
        nease_comparison : pd.DataFrame
            Dataframe containing the all edges found by NEASE and PTM-POSE. This will include edges that are unique to NEASE, unique to PTM-POSE, and common between the two.
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
        ptm_only_graph.remove_edges_from(e for e in self.interaction_graph.edges() if e in nease_nw.edges())

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

    def plot_nease_comparison(self, nease_edges = None, ax = None):
        """
        Given the comparison data generated by compare_to_nease, plot the number of edges identified by NEASE and PTM-POSE

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
            axes to plot on
        nease_edges : pd.DataFrame, optional
                Interactions found by NEASE, which is the output of nease.get_edges() function after running NEASE (see nease_runner module for example). Only needed if you have not run comparison_to_nease() previously
        """
        #check if nease_comparison exists already, if not generate if data is provided
        if not hasattr(self, 'nease_comparison') and nease_edges is not None:
            self.compare_to_nease(nease_edges)
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



    def plot_interaction_network(self, modified_color = 'red', modified_node_size = 10, interacting_color = 'lightblue', interacting_node_size = 5, defaultedgecolor = 'gray', color_edges_by = 'Same', seed = 200, ax = None, proteins_to_label = None, labelcolor = 'black', legend = True):
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
        defaultedgecolor: str
            Color to use for edges in the network. Default is 'gray'. Can choose to color by database by providing a dictionary with database names as keys and colors as values or by specifying 'database' as color
        color_edges_by: str
            
            How to color edges in the network. Default is 'Same', which will color all edges the same color. Can also specify 'Database' to color edges based on the database they are from. If using 'Database', please provide a dictionary with database names as keys and colors as values in defaultedgecolor parameter.
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
        #if nease_edges is not None:
        #    nease_edges['Affected binding'] = nease_edges['Affected binding'].apply(lambda x: x.split(','))
        #    nease_edges = nease_edges.explode('Affected binding')

        #    nease_nw = nx.from_pandas_edgelist(nease_edges, source = 'Gene name', target = 'Affected binding')

            #combine graphs
        #    plt_graph = nx.compose(self.interaction_graph, nease_nw)
        #else:
        #    plt_graph = self.interaction_graph.copy()
        plt_graph = self.interaction_graph.copy()


        plot_interaction_network(plt_graph, self.network_data, self.network_stats, modified_color = modified_color, modified_node_size = modified_node_size, interacting_color = interacting_color, interacting_node_size = interacting_node_size, defaultedgecolor = defaultedgecolor, color_edges_by=color_edges_by, seed = seed, ax = ax, proteins_to_label = proteins_to_label, labelcolor = labelcolor, legend = legend)

    def plot_network_centrality(self,  centrality_measure = 'Degree', top_N = 10, modified_color = 'red', interacting_color = 'black', ax = None):
        """
        Plot the centrality measure for the top N proteins in the interaction network based on a specified centrality measure. This will help identify the most relevant PTMs/genes in the network.

        Parameters
        ----------
        centrality_measure: str
            How to calculate centrality. Options include 'Degree', 'Degree Centrality', 'Betweenness Centrality', 'Closeness Centrality', and 'Eigenvector Centrality'. Default is 'Degree'.
        top_N: int
            Number of top proteins to plot based on the centrality measure. Default is 10.
        modified_color: str
            Color to use for proteins that are spliced. Default is 'red'.
        interacting_color: str
            Color to use for proteins that are not spliced. Default is 'black'.
        ax: matplotlib.pyplot.Axes
            Axes object to plot the centrality bar plot. If None, a new figure will be created. Default is None.
        """
        if not hasattr(self, 'interaction_graph'):
            self.get_interaction_network()
        if not hasattr(self, 'network_stats'):
            self.get_interaction_stats()

        plot_network_centrality(self.network_stats, self.network_data, centrality_measure=centrality_measure,top_N = top_N, modified_color = modified_color, interacting_color = interacting_color, ax = ax)


def get_edge_colors(interaction_graph, network_data, defaultedgecolor = 'gray', color_edges_by = 'Database', database_color_dict = {'PSP/RegPhos':'red', 'PhosphoSitePlus':'green', 'PTMcode':'blue', 'PTMInt':'gold', 'Multiple':'purple'}):
    """
    Get the edge colors to use for a provided networkx graph, either plotting all edges the same color or coloring them based on the database they are from.

    Parameters
    ----------
    interaction_graph: nx.Graph
        networkx graph containing the interaction network
    network_data: pd.DataFrame
        specific network edge data that contains information on which database the interaction is from and any other relevant information (such as regulation change)
    defaultedgecolor : 'str'
        Default color to use for edges if no specific database color is found. Default is 'gray'.
    color_edges_by : str
        How to color the edges. Options are 'Database' to color by the database they are from, or 'Same' to color all edges the same color. Default is 'Database'.
    database_color_dict : dict
        Colors to use for specific databases
    

    """
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

def reformat_pssm(pssm_longform, protein):
    """
    Reformat a long-form PSSM dataframe containing multiple PSSMs into a wide-form PSSM dataframe including only a single domain, with amino acids as rows and positions as columns. 
    
    Long-form PSSM: index is the protein name associated with PSSM data, columns indicate position and amino acid in the PSSM (e.g. '1A', '-1G', '2P', '-4K', etc.), and values are the scores for each amino acid at each position.
    
    Wide-form PSSM: index is the position, columns are the amino acid (for a single protein/domain)

    Parameters
    ----------
    pssm_longform: pd.DataFrame
        Dataframe containing the long-form PSSM data, with index as protein names and columns as positions and amino acids. This is the same format as published in KinaseLibrary publications for kinase and SH2 domain PSSMs
    protein: str
        Protein name to extract the PSSM for. This should match the index of the pssm_longform dataframe.

    Returns
    -------
    pssm_wideform: pd.DataFrame
        Dataframe containing the wide-form PSSM data, with index as positions and columns as amino acids. The values are the scores for each amino acid at each position.
    """
    #extract amino acids and positions from the columns of the PSSM
    amino_acid = []
    position = []
    for i in pssm_longform.columns:
        amino_acid.append(i[-1])
        position.append(int(i[:-1]))
    amino_acid = np.unique(amino_acid)
    position = np.unique(position)

    #go through and extract score for each position and amino acid, placing it in the correct place in the matrix
    pssm_matrix = pd.DataFrame(np.nan, index=position, columns=amino_acid)
    pssm_series = pssm_longform.loc[protein]
    for i in pssm_longform.columns:
        pos = int(i[:-1])
        aa = i[-1]
        pssm_matrix.loc[pos, aa] = pssm_series[i]
    return pssm_matrix

def get_percentile(score, background):
    """
    Given a score and a background distribution, calculate the percentile of the score in the background distribution.

    Parameters
    ----------
    score: float
        Score to calculate the percentile for.
    background: list or np.array
        Background distribution to calculate the percentile from.

    Returns
    -------
    percentile: float
        Percentile of the score in the background distribution.
    """
    if isinstance(background, list):
        background = np.array(background)
    
    if len(background) == 0:
        return np.nan

    return np.sum(background <= score) / len(background) * 100

### PSSMs
class PSSM:
    """
    Class to take a position-specific scoring matrix (PSSM) and score a peptide containing a PTM against the PSSM. 

    Code adapted from work by Gabriela Salazar Lopez
    
    Parameters
    ----------
    pssm: pd.DataFrame
        Dataframe containing the PSSM, with amino acids as rows and positions as columns. The values should be the scores for each amino acid at each position.
        
    """
    def __init__(self, pssm, modification = 'Phosphorylation', residue = 'Y', plus_residues = 4, minus_residues = 2):
        self.pssm = pssm
        self.residue = residue
        self.modification = modification
        self.plus_residues = plus_residues
        self.minus_residues = minus_residues


    def trim_pep (self, pep):
        #find where psite is in pep
        psite_loc = pep.find(self.residue.lower())
        #get start and end locations of peptide to trim to match pssm positions
        beg = psite_loc - self.minus_residues
        end = psite_loc + self.plus_residues + 1
        if beg < 0 or end > len(pep):
            return np.nan
        else:
            return pep[beg:end]
            
    def make_mask(self, pep):
        """
        Create a peptide mask matching PSSM positions
        """
        #strip 'y' from pep
        pep = pep.replace(self.residue.lower(), '')
        mask = np.zeros(self.pssm.shape)
        for i in range(self.pssm.shape[0]):
            for j in range(self.pssm.shape[1]):
                AA = self.pssm.columns[j]
                if pep[i] == AA:
                    mask[i, j] = 1
        return mask

    
    def score_SS(self, pep):
        """
        Use the scansite scoring method to score a peptide against the PSSM. This will first trim the peptide to match the PSSM positions (if not possible will return NaN), then calculate the raw score for the peptide, the optimal score for the PSSM given the peptide length, and finally calculate the final score as (optimal - raw)/optimal.
        """
        pep_trimmed = self.trim_pep(pep)
        if pep_trimmed == pep_trimmed:
            #get raw score for trimmed peptide
            mask = self.make_mask(pep_trimmed)
            masked_PSSM = self.pssm * mask
            bit_scores = np.log(masked_PSSM.sum(axis=1)+ 0.0000000001)/np.log(2) #score peptide and add small amount to all sums to avoid taking log of 0
            raw_score = bit_scores.sum()/(len(pep_trimmed)-1)

            ### Get optimal score for the PSSM given the peptide length
            opts = self.pssm.max(axis=1) + 0.0000000001 #take max of each row and add small amount to avoid taking log of 0
            opts = np.log(opts)/np.log(2) 
            opt_score = opts.sum()/(len(pep_trimmed)-1)

            ### Calculate final score
            score = (opt_score-raw_score)/opt_score
            return score
        else:
            return np.nan
    
    def score_peptides(self, peptides, method = 'scansite'):
        scores = []
        for pep in peptides:
            if method == 'scansite':
                scores.append(self.score_SS(pep))
            else:
                raise ValueError(f"Unknown scoring method: {method}. Available methods are: 'scansite'.")
        return scores
    
class ScorePSSMs:
    """
    Class to score a list of peptides against a PSSM. This will take a list of peptides and score each peptide against the PSSM, returning a dataframe with the scores.

    Parameters
    ----------
    pssm: dict of pd.DataFrame
        Dictionary with dataframes containing the PSSM, with amino acids as rows and positions as columns. The values should be the scores for each amino acid at each position.
    modification: str
        Type of modification to score against the PSSM. Default is 'Phosphorylation'.
    residue: str
        Residue to score against the PSSM. Default is 'Y'.
    plus_residues: int
        Number of residues to consider after the modified residue. Default is 4.
    minus_residues: int
        Number of residues to consider before the modified residue. Default is 2.
    """
    def __init__(self, pssms, spliced_ptms, background_scores = None, modification = 'Phosphorylation', residue = 'Y', plus_residues = 4, minus_residues = 2):
        #reduce spliced_ptms dataframe to only those with the specified modification and residue
        spliced_ptms = spliced_ptms[spliced_ptms['Modification Class'] == modification]
        spliced_ptms = spliced_ptms[spliced_ptms['Residue'] == residue]
        if spliced_ptms.empty:
            raise ValueError(f"No PTMs for associated modification and residue ({modification}, {residue}) found in spliced PTM data, so cannot perform requested analysis.")
        
        self.spliced_ptm = spliced_ptms
        self.modification = modification
        self.residue = residue
        self.plus_residues = plus_residues
        self.minus_residues = minus_residues
        self.background_scores = background_scores

        print(isinstance(pssms, pd.DataFrame))
        #instantiate PSSM objects for each PSSM in the dictionary
        if isinstance(pssms, pd.DataFrame):
            #check if in long form (e.g. KinaseLibrary PSSM format
            if '1A' in pssms.columns:
                self.pssms = {}
                for protein in pssms.index:
                    pssm = reformat_pssm(pssms, protein)
                    self.pssms[protein] = PSSM(pssm, modification = modification, residue = residue, plus_residues = plus_residues, minus_residues = minus_residues)
            elif 'A' in pssms.columns and 1 in pssms.index:
                #if in wide form, just create a single PSSM object
                self.pssms = {'PSSM': PSSM(pssms, modification = modification, residue = residue, plus_residues = plus_residues, minus_residues = minus_residues)}
            else:
                raise ValueError("Dataframe format not recognized. Please provide a dataframe with either long-form PSSM data for multiple proteins (protein in index and positions/amino acids as columns like '1A', '-1G', etc.) or wide-form PSSM data for a single domain/protein (positions as index and amino acids as columns like 'A', 'G', 'P', etc.).")
        elif isinstance(pssms, dict):
            self.pssms = {k: PSSM(v, modification = modification, residue = residue, plus_residues = plus_residues, minus_residues = minus_residues) for k, v in pssms.items()}
        else:
            raise ValueError('PSSMs must be a properly formatted dataframe or dictionary of dataframes')
        
    def print_available_pssms(self):
        """
        Print the available PSSMs in the ScorePSSMs object.
        """
        print("Available PSSMs:")
        for protein in self.pssms:
            print(f"- {protein}")
    

    def score_ptm(self, ptm):
        pep = ptm['Canonical Flanking Sequence']
        ptm_label = ptm['Isoform ID'] + '_' + ptm['Residue'] + str(int(ptm['PTM Position in Isoform']))
        scores = pd.Series(index = self.pssms.keys(), name = ptm_label, dtype = float)
        for protein in self.pssms:
            scores[protein] = self.pssms[protein].score_SS(pep)
        return scores
    
    def score_background(self, background_peptides = None):
        """
        Score a background set of peptides against the PSSMs and return a dataframe with the scores.
        
        Parameters
        ----------
        background: pd.DataFrame
            Dataframe containing the background peptides to score against the PSSMs. The index should be the peptide sequence and the columns should be the PSSM names.

        Returns
        -------
        pd.DataFrame
            Dataframe containing the scores for each peptide against each PSSM, with columns for the PSSM name and the score.
        """
        if background_peptides is None:
            background_peptides = pose_config.ptm_coordinates['Canonical Flanking Sequence'].unique()

        background_scores = []
        for pep in background_peptides:
            scores = self.score_ptm(pep)
            scores.name = pep
            background_scores.append(scores)

        self.background_scores = pd.concat(background_scores, axis=1)


    def score_all_ptms(self):
        """
        Score all PTMs in the spliced PTM dataframe against the PSSMs and return a dataframe with the scores.
        
        Returns
        -------
        pd.DataFrame
            Dataframe containing the scores for each PTM against each PSSM, with columns for the PSSM name and the score.
        """
        ptm_scores = []
        for index, ptm in self.spliced_ptm.iterrows():
            scores = self.score_ptm(ptm)
            ptm_scores.append(scores)

        ptm_scores = pd.concat(ptm_scores, axis = 1)
        self.ptm_scores = ptm_scores

    def get_percentiles(self, background_scores = None):
        """
        Given the background matrix
        """
        if not hasattr(self, 'background_scores') and background_scores is None:
            print('generating scores on entire background set. If you want a specific background set or have already generated this, please provide it in the background_scores parameter.')
            self.score_background()
            background_scores = self.background_scores
        elif background_scores is not None:
            self.background_scores = background_scores


        if not hasattr(self, 'ptm_scores'):
            self.score_all_ptms()

        percentiles = self.ptm_scores.copy()
        for domain in percentiles.columns:
            percentiles[domain] = percentiles[domain].apply(lambda x: get_percentile(x, background_scores[domain].dropna().values))
        self.ptm_percentiles = percentiles


    
    def get_top_domains(self, threshold = 2, threshold_type = 'scores', n = 8):
        """
        Get the top domains based on the specified threshold and number of top domains. Threshold can either based on scores or percentile
        """
        if not hasattr(self, 'ptm_scores'):
            self.score_all_ptms()
        
        if threshold_type == 'scores':
            top_domains = self.ptm_scores[self.ptm_scores > threshold].sort_values(ascending=False).head(n)
        top_domains = ';'.join(top_domains.index)
        return top_domains
    
    def append_top_domains(self, threshold = 2, threshold_type = 'scores', n = 8):
        """
        Append the top domains to the spliced PTM dataframe based on the specified threshold and number of top domains.
        
        Parameters
        ----------
        threshold: int
            Threshold for the scores to consider a domain as top. Default is 2.
        threshold_type: str
            Type of threshold to use. Default is 'scores', which will use the scores from the PSSMs.
        n: int
            Number of top domains to return. Default is 8.
        """
        if not hasattr(self, 'ptm_scores'):
            self.score_all_ptms()
        
        self.spliced_ptm['Top Domains'] = self.spliced_ptm.apply(lambda x: self.get_top_domains(threshold=threshold, threshold_type=threshold_type, n=n), axis=1)

    def construct_score_gmt_file(self):
        pass

    def construct_percentile_gmt_file(self):
        pass




