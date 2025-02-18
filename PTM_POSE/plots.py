import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import networkx as nx

from ptm_pose import analyze, pose_config
from ptm_pose import flanking_sequences as fs


def modification_breakdown(spliced_ptms = None, altered_flanks = None, colors = sns.color_palette('colorblind'), ax = None):
    """
    Plot the number of PTMs that are differentially included or have altered flanking sequences, separated by PTM type

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs that are differentially included
    altered_flanks: pd.DataFrame
        Dataframe with PTMs that have altered flanking sequences
    colors: list
        List of colors to use for the bar plot (first two will be used). Default is seaborn colorblind palette.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    """
    if spliced_ptms is None and altered_flanks is None:
        raise ValueError('Either spliced_ptms or altered_flanks must be provided to plot modification breakdown. Both may be provided.')
    
    if ax is None:
        fig, ax = plt.subplots(figsize = (4,4))

    #separate rows into unique PTM types
    if spliced_ptms is not None and altered_flanks is not None:
        differentially_included_counts = analyze.get_modification_counts(spliced_ptms.copy())
        altered_flanks_counts = analyze.get_modification_counts(altered_flanks.copy())
        ax.barh(differentially_included_counts.index, differentially_included_counts.values, color = colors[0], label = 'Differentially Included PTMs')
        altered_flanks_counts = altered_flanks_counts.reindex(differentially_included_counts.index, fill_value = 0)
        ax.barh(altered_flanks_counts.index, altered_flanks_counts.values, left = differentially_included_counts.values, color = colors[1], label = 'PTMs with Altered Flank')
        ax.legend()

        #annotate with number of combined PTMs
        total_count = differentially_included_counts.add(altered_flanks_counts, fill_value = 0)
        for i, num_ptm in enumerate(total_count.values):
            ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')  

        ax.set_xlim([0, total_count.max()*1.1])

    elif spliced_ptms is not None:
        modification_counts = analyze.get_modification_counts(spliced_ptms)
        ax.barh(modification_counts.index, modification_counts.values, color = colors[0])

        #annotate with number of PTMs
        for i, num_ptm in enumerate(modification_counts.values):
            ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')
    elif altered_flanks is not None:
        modification_counts = analyze.get_modification_counts(altered_flanks)
        ax.barh(modification_counts.index, modification_counts.values, color = colors[0])

        #annotate with number of PTMs
        for i, num_ptm in enumerate(modification_counts.values):
            ax.text(num_ptm, i, str(num_ptm), ha = 'left', va = 'center')

    ax.set_xlabel('Number of PTMs')

def show_available_annotations(spliced_ptms, show_all_ptm_count = True, ax = None):
    """
    Given a dataframe with ptm annotations added, show the number of PTMs associated with each annotation type

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs and annotations added
    show_all_ptm_count: bool
        Whether to show the total number of PTMs in the dataset. Default is True.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize = (5,5))

    if show_all_ptm_count:
        num_ptms = [spliced_ptms.drop_duplicates(['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0]]
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
        num_ptms.append(spliced_ptms.dropna(subset = annotation).drop_duplicates(subset = ['UniProtKB Accession', 'Residue', 'PTM Position in Isoform']).shape[0])
        num_ptms_filters.append(ylabel_dict[annotation])
        filter_source.append(annotation.split(':')[0])

    
    #plot bar plot
    #color bars based on datasource
    palette = {'None': 'gray', 'PSP': 'blue', 'ELM': 'green', 'PTMInt':'red', 'PTMcode':'purple', 'DEPOD':'orange', 'RegPhos':'gold', 'Combined':'black'}
    colors = []
    for source in filter_source:
        colors.append(palette[source])

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

def plot_annotations(spliced_ptms, database = 'PhosphoSitePlus', annot_type = 'Function', collapse_on_similar = True, colors = None, fontsize = 10, top_terms = 5, legend = True, legend_loc = None, title = None, title_type = 'database', ax = None):
    """
    Given a dataframe with PTM annotations added, plot the top annotations associated with the PTMs

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe with PTMs and annotations added
    database: str
        Database to use for annotations. Default is 'PhosphoSitePlus'.
    annot_type: str
        Type of annotation to plot. Default is 'Function'.
    collapse_on_similar: bool
        Whether to collapse similar annotations into a single category. Default is True.
    colors: list
        List of colors to use for the bar plot. Default is None.
    top_terms: int
        Number of top terms to plot. Default is 5.
    legend: bool
        Whether to show the legend. Default is True.
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    title_type: str
        Type of title to use for the plot. Default is 'database'. Options include 'database' and 'detailed'.
    """
    _, annotation_counts = analyze.get_ptm_annotations(spliced_ptms, annotation_type = annot_type, database = database, collapse_on_similar = collapse_on_similar)
    if annotation_counts is None:
        return None


    if ax is None:
        fig, ax = plt.subplots(figsize = (2,3))

    if database == 'PTMcode': #convert to readable gene name
        annotation_counts.index = [pose_config.uniprot_to_genename[i].split(' ')[0] if i in pose_config.uniprot_to_genename.keys() else i for i in annotation_counts.index]

    if colors is None:
        colors = ['lightgrey', 'gray', 'white']

    if isinstance(annotation_counts, pd.Series):
        annotation_counts = annotation_counts.head(top_terms).sort_values(ascending = True)
        if isinstance(colors, list) or isinstance(colors, np.ndarray):
            colors = colors[0]
        
        ax.barh(annotation_counts.index, annotation_counts.values, color = colors, edgecolor = 'black')
        legend = False
    else:
        annotation_counts = annotation_counts.head(top_terms).sort_values(by = 'All Impacted', ascending = True)
        ax.barh(annotation_counts['Excluded'].index, annotation_counts['Excluded'].values, height = 1, edgecolor = 'black', color = colors[0])
        ax.barh(annotation_counts['Included'].index, annotation_counts['Included'].values, left = annotation_counts['Excluded'].values, height = 1, color = colors[1], edgecolor = 'black')
        ax.barh(annotation_counts['Altered Flank'].index, annotation_counts['Altered Flank'].values, left = annotation_counts['Excluded'].values+annotation_counts['Included'].values, height = 1, color = colors[2], edgecolor = 'black')
    #ax.set_xticks([0,50,100,150])
    ax.set_ylabel('', fontsize = fontsize)
    ax.set_xlabel('Number of PTMs', fontsize = fontsize)

    if title is not None:
        ax.set_title(title, fontsize = fontsize, weight = 'bold')
    elif title_type == 'detailed':
        ax.set_title(f'Top {top_terms} {database} {annot_type} Annotations', fontsize = fontsize, weight = 'bold')
    elif title_type == 'database':
        ax.set_title(f'{database}', fontsize = fontsize, weight = 'bold')

    #label_dict = {'EXONT:Name':'Exon Ontology Term', 'PSP:ON_PROCESS':'Biological Process (PSP)', 'PSP:ON_FUNCTION':'Molecular Function (PSP)', 'Combined:Kinase':'Kinase'}
    #ax.text(-1*ax.get_xlim()[1]/10, top_terms-0.2, label_dict[term_to_plot], weight = 'bold', ha = 'right', fontsize = 8)
    x_label_dict = {'Function':'Number of PTMs\nassociated with Function', 'Process':'Number of PTMs\nassociated with Process', 'Disease':'Number of PTMs\nassociated with Disease', 'Kinase':'Number of Phosphosites\ntargeted by Kinase', 'Interactions': 'Number of PTMs\nthat regulate interaction\n with protein','Motif Match':'Number of PTMs\nfound within a\nmotif instance', 'Intraprotein': 'Number of PTMs\nthat are important\for intraprotein\n interactions','Phosphatase':'Number of Phosphosites\ntargeted by Phosphatase', 'Perturbation (DIA2)': "Number of PTMs\nAffected by Perturbation\n(Measured by DIA)", 'Perturbation (PRM)': 'Number of PTMs\nAffected by Perturbation\n(Measured by PRM)', 'NetPath':'Number of PTMs/Genes\nassociated with NetPath', 'Perturbation':'Number of PTMs\nAffected by Perturbation'}
    ax.set_xlabel(x_label_dict[annot_type], fontsize = fontsize)
    
    #make a custom legend
    if legend:
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(facecolor = colors[0], edgecolor = 'black', label = 'Excluded'), mpatches.Patch(facecolor = colors[1], edgecolor = 'black', label = 'Included'),mpatches.Patch(facecolor = colors[2], edgecolor = 'black', label = 'Altered Flank')]
        ax.legend(handles = handles, ncol = 1, fontsize = fontsize - 2, title = 'Type of Impact', title_fontsize = fontsize - 1, loc = 'upper right', bbox_to_anchor = legend_loc)

    return ax



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



def plot_EnrichR_pies(enrichr_results, top_terms = None, terms_to_plot = None, colors = None, edgecolor = None, row_height = 0.3, type = 'circle', ax = None):
    """
    Given PTM-specific EnrichR results, plot EnrichR score for the provided terms, with each self point represented as a pie chart indicating the fraction of genes in the group with PTMs
    
    Parameters
    ----------
    ptm_results: pd.selfFrame
        selfFrame containing PTM-specific results from EnrichR analysis
    num_to_plot: int
        number of terms to plot, if None, will plot all terms. Ignored if specific terms are provided in terms to plot list
    terms_to_plot: list
        list of terms to plot
    colors: list
        list of colors to use for pie slices. Default is None, which will use seaborn colorblind palette
    edgecolor: str
        color to use for edge of pie slices. Default is None, which will use the same color as the pie slice
    row_height: float
        height of each row in the plot. Default is 0.3.
    type: str
        type of pie chart to plot. Default is 'circle'. Options include 'circle' and 'donut' (hole in center).
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
        if top_terms is not None:
            plt_data = plt_data.iloc[-top_terms:] if top_terms < plt_data.shape[0] else plt_data
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



    ax.set_xlabel('EnrichR Combined\nScore', fontsize = 11)

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


def plot_interaction_network(interaction_graph, network_data, network_stats = None, modified_color = 'red', modified_node_size = 10, interacting_color = 'lightblue', interacting_node_size = 1, defaultedgecolor = 'gray', color_edges_by = 'Same', seed = 200, legend_fontsize = 8, ax = None, proteins_to_label = None, labelcolor = 'black'):
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
            network_stats = analyze.get_interaction_stats(interaction_graph)
        
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

def location_of_altered_flanking_residues(altered_flanks, figsize = (4,3), modification_class = None, residue = None):
    """
    Plot the number of PTMs with altered residues as specific positions relative to the PTM site. This includes the specific position of the residue (-5 to +5 from PTM site) and the specific side of the PTM site that is altered (N-term or C-term)

    Parameters
    ----------
    altered_flanks: pd.DataFrame
        Dataframe with altered flanking sequences, and annotated information added with analyze.compare_flanking_sequences
    figsize: tuple
        Size of the figure. Default is (4,3).
    modification_class: str
        Specific modification class to plot. Default is None, which will plot all modification classes.
    residue: str
        Specific residue to plot. Default is None, which will plot all residues.

    """
    
    fig, ax = plt.subplots(nrows = 2, figsize = figsize, height_ratios = [0.5,1])
    fig.subplots_adjust(hspace = 1)

    if modification_class is not None:
        altered_flanks = altered_flanks[altered_flanks['Modification Class'].str.contains(modification_class)].copy()
    
    if residue is not None:
        altered_flanks = altered_flanks[altered_flanks['Residue'] == residue].copy()

    #### plot of side of modification that flank is altered
    terminus = altered_flanks.groupby('Altered Flank Side').size()
    terminus = terminus[['N-term only', 'C-term only']] #focus on cases where only one side is altered for ease of plotting
    ax[0].bar(terminus.index, terminus.values, color = 'gray')
    ax[0].set_xlabel('Location of Altered Region', fontsize = 9)
    ax[0].set_xticklabels(['N-term\nonly', 'C-term\nonly'])
    ax[0].set_ylabel('# of PTMs', fontsize = 9)

    #### plot specific positions of altered residues relative to PTM
    position_breakdown = altered_flanks.explode(['Altered Positions', 'Residue Change']).copy()[['Gene', 'Residue', 'PTM Position in Isoform','Altered Positions', 'Residue Change']]
    position_breakdown = position_breakdown.groupby('Altered Positions').size()
    ax[1].bar(position_breakdown.index, position_breakdown.values, color = 'gray')
    ax[1].set_xlim([-5.5,5.5])
    ax[1].set_xlabel('Position Relative to PTM', fontsize = 9)
    ax[1].set_ylabel('# of Changed\nResidues', fontsize = 9)
    ax[1].set_xticks(np.arange(-5,6,1))

def alterations_matrix(altered_flanks, modification_class = None, residue = None, title = '', ax = None):
    """
    Given the altered flanking sequences dataframe, plot a matrix showing the positions of altered residues for specific proteins, as well as the specific change

    Parameters
    ----------
    altered_flanks: pd.DataFrame
        Dataframe with altered flanking sequences, and annotated information added with analyze.compare_flanking_sequences

    modification_class: str
        Specific modification class to plot. Default is None, which will plot all modification classes.

    residue: str
        Specific residue to plot. Default is None, which will plot all residues.
    title: str
        Title of the plot. Default is '' (no title).
    ax: matplotlib.Axes
        Axis to plot on. If None, will create new figure. Default is None.
    """
    #extract altered flanking sequences and make sure there is altered position data
    position_breakdown = altered_flanks.copy()
    if 'Altered Positions' not in position_breakdown.columns:
        position_breakdown = fs.compare_flanking_sequences(position_breakdown)

    position_breakdown = position_breakdown.dropna(subset = ['Altered Positions', 'Residue Change'])

    #restrit to desired PTM types and residues
    if modification_class is not None:
        position_breakdown = position_breakdown[position_breakdown['Modification Class'].str.contains(modification_class)].copy()
    if residue is not None:
        position_breakdown = position_breakdown[position_breakdown['Residue'] == residue].copy()

    #add ptm column to position breakdown
    position_breakdown['PTM'] = position_breakdown['Gene'] + '_' + position_breakdown['Residue'] + position_breakdown['PTM Position in Isoform'].astype(str)

    #separate altered residue into individual rows
    position_breakdown = position_breakdown.explode(['Altered Positions', 'Residue Change']).copy()[['Gene', 'PTM','Altered Positions', 'Residue Change']]

    #convert altered positions to integers and remove duplicates
    position_breakdown['Altered Positions'] = position_breakdown['Altered Positions'].astype(int)
    position_breakdown = position_breakdown.drop_duplicates()

    #position_breakdown = position_breakdown.drop_duplicates(subset = ["PTM", "Altered_Positions"], keep = False)
    position_breakdown['PTM']
    position_matrix = position_breakdown.pivot(columns = 'Altered Positions', index = 'PTM', values= 'Residue Change')
    for i, pos in zip(range(11),range(-5, 6)):
        if pos not in position_matrix.columns:
            position_matrix.insert(i, pos, np.nan)


    #replace strings with 1 and nans with 0
    position_values = position_matrix.map(lambda x: 1 if x == x else np.nan).sort_values(by = 5)
    position_matrix = position_matrix.loc[position_values.index]
    #plot heatmap with black edges around cells and no colorbar, annotate with strings in position matrix

    if ax is None:
        fig, ax = plt.subplots(figsize = (4,3))
    sns.heatmap(position_values, cmap = 'Greens', vmin = 0, vmax = 2, ax = ax, cbar = False, linewidths = 0.5, linecolor = 'black', yticklabels=True )
    ax.set_facecolor("lightgrey")

    #annotate with strings in the position matrix
    for i in range(position_values.shape[0]):
        for j in range(position_values.shape[1]):
            if position_values.iloc[i,j] == 1:
                ax.text(j+0.5, i+0.5, position_matrix.iloc[i,j], ha = 'center', va = 'center', fontsize = 6)

    #adjust figure parameters
    ax.set_xticklabels(['-5','-4','-3','-2','-1','PTM','1','2','3','4','5'], fontsize = 8)
    ax.set_xlabel('')
    ax.tick_params(axis = 'y', labelsize = 8)
    ax.set_title(title, fontsize = 9)
