

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#ptm_pose imports
from ptm_pose import helpers


def plot_filter_impact(ptms, output_type = 'count', topn = 10, ax = None, **kwargs):
    """
    Given a dataframe of PTMs and a set of filter arguments to be passed to helpers.filter_ptms, this function will plot the number or fraction of PTMs that are retained after filtering for each modification type

    Parameters
    ----------
    ptms : pd.DataFrame
        Dataframe containing PTM data with a column 'Modification Class' that contains the type of modification (e.g. phosphorylation, acetylation, etc.)
    output_type : str, optional
        Type of output to plot, either 'count' or 'fraction'. The default is 'count'.
    topn : int, optional
        The number of top modification classes to plot. The default is 10.
    ax : matplotlib.axes.Axes, optional
        The axes to plot on. If None, a new figure and axes will be created. The default is None.
    **kwargs : keyword arguments
        Additional keyword arguments to be passed to the filter_ptms function (e.g. min_studies, min_compendia, etc.). These will be extracted and checked for validity.
    """
    filter_arguments = helpers.extract_filter_kwargs(**kwargs)
    helpers.check_filter_kwargs(filter_arguments)
    filtered_ptms=helpers.filter_ptms(ptms, **filter_arguments)

    original_mods = ptms['Modification Class'].value_counts()
    original_mods.name = 'Original'
    filtered_mods = filtered_ptms['Modification Class'].value_counts()
    filtered_mods.name = 'Filtered'
    #sort by top n values after filtering
    original_mods = original_mods.sort_values(ascending = False).head(topn)
    filtered_mods = filtered_mods[filtered_mods.index.isin(original_mods.index)]

    #grab y-axis label and labelpad based on output_type
    if output_type == 'fraction':
        #convert counts to fractions
        original_mods = original_mods / original_mods.sum()
        filtered_mods = filtered_mods / filtered_mods.sum()
        ylabel = 'Fraction of PTMs'
        labelpad = 0.02
    elif output_type == 'count':
        ylabel = 'Number of PTMs'
        labelpad = 10
    else:
        raise ValueError("output_type must be either 'count' or 'fraction'")
    
    plt_data = pd.concat([original_mods, filtered_mods], axis = 1)
    plt_data = plt_data.fillna(0)

    if ax is None:
        fig, ax = plt.subplots(figsize = (3,2))

    plt_data.plot(kind = 'bar', ax = ax)

    ax.set_ylabel(ylabel)

    #annotate tops of bars with numbers
    for p in ax.patches:
        if output_type == 'count':
            ax.annotate(str(int(p.get_height())), (p.get_x() * 1.005, p.get_height() + labelpad), fontsize = 8, rotation = 90, annotation_clip=False)
        elif output_type == 'fraction':
            #for fractions, convert to percentage and round to 2 decimal places
            ax.annotate(str(round(p.get_height() * 100, 2)) + '%', (p.get_x() * 1.005, p.get_height() + labelpad), fontsize = 8, rotation = 90, annotation_clip=False)

    #remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def assess_filter_range(ptms, min_value = 0, max_value = None, step = None, filter_type = 'min_studies', phospho_only_evidence_filter = True, ax = None, fontsize = 11):
    num_ptms = []
    frac_phospho = []

    #check filter_type is valid
    if filter_type not in ['min_studies', 'min_compendia', 'min_MS', 'min_LTP']:
        raise ValueError("filter_type must be one of ['min_studies', 'min_compendia', 'min_MS', 'min_LTP']")

    #grab max value if not provided
    if max_value is None:
        if filter_type == 'min_studies':
            max_value = int(ptms[['MS_LIT', 'LT_LIT']].sum(axis =1).max())
        elif filter_type == 'min_compendia':
            max_value = ptms['Number of Compendia'].max()
        elif filter_type == 'min_MS':
            max_value = ptms[['MS_LIT', 'MS_CST']].sum(axis = 1).max()
        elif filter_type == 'min_LTP':
            max_value = ptms['LT_LIT'].max()
    
    #if specific step value not provided, round to nearest 10% of max value
    if step is None:
        step = round(max_value/10)

    #filter PTMs using the indicated filter type method for value in range
    x = np.arange(min_value, int(max_value) + 1, step)
    for i in x:
        #filter PTMs
        if filter_type == 'min_studies':
            filtered_ptms = helpers.filter_ptms(ptms, report_removed = False, min_studies = i, phospho_only_evidence_filter = phospho_only_evidence_filter)
        elif filter_type == 'min_compendia':
            filtered_ptms = helpers.filter_ptms(ptms, report_removed = False, min_compendia = i, phospho_only_evidence_filter=phospho_only_evidence_filter)
        elif filter_type == 'min_MS':
            filtered_ptms = helpers.filter_ptms(ptms, report_removed = False, min_MS_observations = i, phospho_only_evidence_filter=phospho_only_evidence_filter)
        elif filter_type == 'min_LTP':
            filtered_ptms = helpers.filter_ptms(ptms, report_removed = False, min_LTP_studies = i, phospho_only_evidence_filter=phospho_only_evidence_filter)

        #save number of PTMs and the fraction that are phosphorylated
        num_ptms.append(filtered_ptms.shape[0])
        #fraction of PTMs that are phosphorylation sites
        if filtered_ptms.shape[0] > 0 and 'Phosphorylation' in filtered_ptms['Modification Class'].unique():  
            filtered_mods = filtered_ptms['Modification Class'].value_counts()
            phospho_fraction = filtered_mods['Phosphorylation']/filtered_mods.sum()
            frac_phospho.append(phospho_fraction)
        elif filtered_ptms.shape[0] == 0:
            frac_phospho.append(np.nan)
        else:
            frac_phospho.append(0)

    x_label_dict = {'min_studies': 'Minimum number of\nliterature reports', 'min_LTP': 'Minimum number of\nLow-throughput Studies', 'min_MS': 'Minimum number of\nMS Observations', 'min_compendia': 'Minimum number of\ncompendia'}

    if ax is None:
        fig, ax = plt.subplots(figsize = (3,3))

    ax.plot(x, num_ptms, color = 'blue')
    ax.set_ylabel('Number of PTMs', color = 'blue', fontsize = fontsize)
    #change color of tick labels
    ax.tick_params(axis='y', labelcolor='blue')
    ax.set_xlabel(x_label_dict[filter_type], fontsize = fontsize)
    ax2 = ax.twinx()
    ax2.plot(x, frac_phospho, color = 'red')
    ax2.set_ylabel('Phosphorylation\nFraction', color = 'red', fontsize = fontsize)
    ax2.tick_params(axis='y', labelcolor='red')