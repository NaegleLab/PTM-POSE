import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import networkx as nx

from ptm_pose import analyze


def show_available_annotations(spliced_ptms, show_all_ptm_count = True, figsize = (5, 5)):
    """
    Given a dataframe with ptm annotations added, show the number of PTMs associated with each annotation type

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs and annotations added
    show_all_ptm_count: bool
        Whether to show the total number of PTMs in the dataset. Default is True.
    figsize: tuple
        Size of the figure. Default is (5, 5).
    
    Outputs
    -------
    bar plot showing the number of PTMs associated with each annotation type
    """
    if show_all_ptm_count:
        num_ptms = [spliced_ptms.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).shape[0]]
        num_ptms_filters = ['All PTMs']
        filter_source = ['None']
    else:
        num_ptms = []
        num_ptms_filters = []
        filter_source = []

    #look for annotations and add counts to lists
    ylabel_dict = {'PSP:ON_PROCESS':'Biological Process (PSP)', 'PSP:ON_FUNCTION':'Molecular Function (PSP)', 'PSP:Kinase':'Kinase (PSP)', 'PSP:Disease_Association':'Disease Association (PSP)', 'PSP:ON_PROT_INTERACT':'Interactions (PSP)', 'PSP:ON_OTHER_INTERACT':'Nonprotein Interactions (PSP)', 'ELM:Interactions':'Interactions (ELM)', 'ELM:Motif Matches':'Motif Matches (ELM)', 'PTMInt:Interaction':'Interactions (PTMInt)', 'PTMcode:Intraprotein_Interactions':'Intraprotein (PTMcode)','PTMcode:Interprotein_Interactions':'Interactions (PTMcode)', 'DEPOD:Phosphatase':'Phosphatase (DEPOD)', 'RegPhos:Kinase':'Kinase (RegPhos)', 'Combined:Kinase':'Kinase (Combined)', 'Combined:Interactions':'Interactions (Combined)'}
    available_annotations = [col for col in spliced_ptms.columns if 'Combined' in col or 'PSP:' in col or 'ELM:Interactions' in col or 'PTMInt:' in col or 'PTMcode:' in col or 'DEPOD:' in col or 'RegPhos:' in col]
    for annotation in available_annotations:
        num_ptms.append(spliced_ptms.dropna(subset = annotation).drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).shape[0])
        num_ptms_filters.append(ylabel_dict[annotation])
        filter_source.append(annotation.split(':')[0])

    
    #plot bar plot
    #color bars based on datasource
    palette = {'None': 'gray', 'PSP': 'blue', 'ELM': 'green', 'PTMInt':'red', 'PTMcode':'purple', 'DEPOD':'orange', 'RegPhos':'gold', 'Combined':'black'}
    colors = []
    for source in filter_source:
        colors.append(palette[source])
    fig, ax = plt.subplots(figsize = figsize)
    ax.barh(num_ptms_filters[::-1], num_ptms[::-1], color = colors[::-1])
    ax.set_xlabel('Number of PTMs with annotation')
    
    #annotate with number of PTMs
    for i, num_ptm in enumerate(num_ptms[::-1]):
        ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')

    #create legend
    handles = [plt.Rectangle((0,0),1,1, color = color) for color in palette.values() if color != 'gray']
    labels = [source for source in palette.keys() if source != 'None']
    ax.legend(handles, labels, title = 'Annotation Source')
    plt.show()


def draw_pie(dist, xpos, ypos, size,colors,edgecolor =None, type = 'donut', ax=None):
    """
    Draws pies individually, as if points on a scatter plot. This function was taken from this stack overflow post: https://stackoverflow.com/questions/56337732/how-to-plot-scatter-pie-chart-using-matplotlib
    
    Parameters
    ----------
    dist: list
        list of values to be represented as pie slices for a single point
    xpos: float
        x position of pie chart in the scatter plot
    ypos: float
        y position of pie chart in the scatter plot
    size: float
        size of pie chart
    colors: list
        list of colors to use for pie slices
    ax: matplotlib.Axes
        axis to plot on, if None, will create new figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))
    #remove slices with 0 size
    colors = [c for c, d in zip(colors, dist) if d != 0]
    dist = [d for d in dist if d != 0]
    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()

    num_colors = len(dist)
    for i, r1, r2 in zip(range(num_colors), pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])

        ax.scatter([xpos], [ypos], marker=xy, s=size, facecolor= colors[i], edgecolors=edgecolor, linewidth = 0.3)

        if type == 'donut': # add white circle in the middle
            donut_edgecolors = 'w' if edgecolor is None else edgecolor
            ax.scatter([xpos], [ypos], s=size/5, facecolor='w', edgecolors=donut_edgecolors, linewidth = 0.3)
    return ax



def plot_EnrichR_pies(enrichr_results, terms_to_plot = None, colors = None, edgecolor = None, row_height = 0.3, leg_loc = (0.5,1.1), type = 'circle', ax = None):
    """
    Given PTM-specific EnrichR results, plot EnrichR score for the provided terms, with each self point represented as a pie chart indicating the fraction of genes in the group with PTMs
    
    Parameters
    ----------
    ptm_results: pd.selfFrame
        selfFrame containing PTM-specific results from EnrichR analysis
    terms_to_plot: list
        list of terms to plot
    ax: matplotlib.Axes
        axis to plot on, if None, will create new figure
    """
    if colors is None:
        colors = sns.color_palette('colorblind', n_colors = 3)


    plt_data = enrichr_results.copy()
    plt_data['Number with Differential Inclusion Only'] = plt_data['Genes with Differentially Included PTMs only'].apply(lambda x: len(x.split(';')))
    plt_data['Number with Altered Flank Only'] = plt_data['Genes with Differentially Included PTMs only'].apply(lambda x: len(x.split(';')))
    plt_data['Number with Both'] = plt_data['Genes with Both'].apply(lambda x: len(x.split(';')) if x != '' else 0)

    if terms_to_plot is None:
        plt_data = plt_data.sort_values(by = 'Combined Score')
    else:
        plt_data = plt_data[plt_data['Term'].isin(terms_to_plot)].sort_values(by = 'Combined Score')
        if plt_data.shape[0] == 0:
            print('No significant terms found in EnrichR results. Please check the terms_to_plot list and try again.')
            return

    #remove gene ontology specific terms
    plt_data['Term'] = plt_data['Term'].apply(lambda x: x.split(' R-HSA')[0] +' (R)' if 'R-HSA' in x else x.split('(GO')[0]+' (GO)')
    #construct multiple piecharts for each term in 'Term' column, where location along x-axis is dictated by combined score and piechart is dictated by 'Fraction With PTMs'
    plt_data = plt_data.reset_index(drop = True)

    #set up figure
    if ax is None:
        figure_length = plt_data.shape[0]*row_height
        fig, ax = plt.subplots(figsize = (2, figure_length))
    
    #get non-inf max score and replace inf values with max score
    maxscore = np.nanmax(plt_data['Combined Score'][plt_data['Combined Score'] != np.inf])
    plt_data['Combined Score'] = plt_data['Combined Score'].replace([-np.inf, np.inf], maxscore)
    ax.set_xlim([maxscore*-0.05, maxscore*1.1])
    mult = 4
    ax.set_yticks(list(range(0,plt_data.shape[0]*mult,mult)))
    ax.set_yticklabels(plt_data['Term'].values)
    ax.set_ylim([-(mult/2), plt_data.shape[0]*mult-(mult/2)])
    type = 'circle'
    event_type = plt_data['Type'].values[0]
    for i, row in plt_data.iterrows():
        if event_type == 'Differentially Included + Altered Flanking Sequences':
            draw_pie([row['Number with Differential Inclusion Only'], row['Number with Altered Flank Only'], row['Number with Both']],xpos = row['Combined Score'], ypos = i*mult, colors = colors, edgecolor=edgecolor,ax = ax, type = type, size = 70)
        else:
            draw_pie([1],xpos = row['Combined Score'], ypos = i*mult, colors = colors, edgecolor=edgecolor,ax = ax, type = type, size = 70)
        
        ax.axhline(i*mult+(mult/2), c= 'k', lw = 0.5)
        ax.axhline(i*mult-(mult/2), c = 'k', lw = 0.5)
        #ax.tick_params(labelsize = )

    #make a custom legend
    if event_type == 'Differentially Included + Altered Flanking Sequences':
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(color = colors[2], label = 'Contains Both Events'), mpatches.Patch(color = colors[1], label = 'PTMs with Altered Flanking Sequence'), mpatches.Patch(color = colors[0], label = 'Differentially Included PTMs')]
        ax.legend(handles = handles, loc = 'upper center', borderaxespad = 0, bbox_to_anchor = (0.5, 1 + (1/figure_length)), ncol = 1, fontsize = 9)



    ax.set_xlabel('EnrichR Combined Score', fontsize = 11)

def plot_interaction_network(interaction_graph, network_data, network_stats = None, modified_color = 'red', modified_node_size = 10, interacting_color = 'lightblue', interacting_node_size = 1, edgecolor = 'gray', seed = 200, ax = None, proteins_to_label = None, labelcolor = 'black'):
    """
    Given the interactiong graph and network data outputted from analyze.get_interaction_network, plot the interaction network, signifying which proteins or ptms are altered by splicing and the specific regulation change that occurs

    Parameters
    ----------

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
        fig, ax = plt.subplots(figsize = (6,6))

    nx.draw(interaction_graph, node_size = node_sizes, node_color = node_colors, edge_color = edgecolor, style = edge_style, ax = ax)

    #add legend for colored nodes
    modified_node = mlines.Line2D([0], [0], color='w',marker = 'o', markersize=modified_node_size,linewidth = 0.2, markerfacecolor = modified_color, markeredgecolor=modified_color, label='Spliced Protein')
    interacting_node = mlines.Line2D([0], [0], color='w', markerfacecolor = interacting_color, markeredgecolor=interacting_color, marker = 'o', markersize=interacting_node_size, linewidth = 0.2, label='Interacting Protein')
    solid_line = mlines.Line2D([0], [0], color='gray', linestyle = 'solid', label = 'Interaction increases')
    dashdot_line = mlines.Line2D([0], [0], color='gray', linestyle = 'dashdot', label = 'Interaction impact unclear')
    dotted_line = mlines.Line2D([0], [0], color='gray', linestyle = 'dotted', label = 'Interaction decreases')
    handles = [solid_line,dashdot_line, dotted_line, modified_node, interacting_node]
    ax.legend(handles = handles, loc = 'upper center', ncol = 2, fontsize = 6, bbox_to_anchor = (0.5, 1.1))

    #if requested, label specific proteins in the network
    if proteins_to_label is not None and isinstance(proteins_to_label, list):
        for protein in proteins_to_label:
            ax.text(interaction_graph.pos[protein][0], interaction_graph.pos[protein][1], protein, fontsize = 10, fontweight = 'bold', color = labelcolor)
    elif proteins_to_label is not None and isinstance(proteins_to_label, int):
        if network_stats is None:
            network_stats = analyze.get_interaction_stats(interaction_graph)
        
        network_stats = network_stats.sort_values(by = 'Degree', ascending = False).iloc[:proteins_to_label]
        for index, row in network_stats.iterrows():
            ax.text(interaction_graph.pos[index][0], interaction_graph.pos[index][1], index, fontsize = 10, fontweight = 'bold', color = labelcolor)
    elif proteins_to_label is not None and isinstance(proteins_to_label, str):
        ax.text(interaction_graph.pos[proteins_to_label][0], interaction_graph.pos[proteins_to_label][1], proteins_to_label, fontsize = 10, fontweight = 'bold', color = labelcolor)
    elif proteins_to_label is not None:
        print('Proteins to label must be a list of strings or a single string. Ignoring when plotting.')
    
def plot_network_centrality(network_stats, network_data = None, centrality_measure = 'Degree', top_N = 10, modified_color = 'red', interacting_color = 'black', ax = None):
    if centrality_measure not in network_stats.columns:
        raise ValueError('Centrality measure not found in network_stats dataframe. Please check the inputted centrality_measure. Available measures include Degree, Degree Centrality, Betweenness Centrality, Closeness Centrality, and Eigenvector Centrality.')
    
    #get specific centrality measure
    plt_data = network_stats.sort_values(by = centrality_measure, ascending = False).iloc[:top_N]
    
    if network_data is not None:
        colors = []
        for index, row in plt_data.iterrows():
            if index in network_data['Modified Gene'].unique():
                colors.append(modified_color)
            else:
                colors.append(interacting_color)
    else:
        colors = modified_color
    
    if ax is None:
        fig, ax = plt.subplots(figsize = (3,3))

    ax.barh(plt_data.index, plt_data[centrality_measure], color = colors)
    ax.set_xlabel(f'{centrality_measure}')
