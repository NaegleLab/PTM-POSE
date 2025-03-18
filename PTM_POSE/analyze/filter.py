

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#ptm_pose imports
from ptm_pose import helpers


def plot_filter_impact(ptms, min_dpsi = 0.2, alpha = 0.05, min_studies = 0, min_MS_observations = 0, min_LTP_studies = 0, min_compendia = 0, remove_novel = False):
    filtered_ptms=helpers.filter_ptms(ptms, min_dpsi = min_dpsi, alpha = alpha, min_studies = min_studies, min_MS_observations = min_MS_observations, min_LTP_studies = min_LTP_studies, min_compendia = min_compendia, remove_novel = remove_novel)

def assess_filter_range(ptms, min_value = 0, max_value = None, step = None, filter_type = 'min_studies', phospho_only_evidence_filter = True, ax = None):
    num_ptms = []
    frac_phospho = []

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
    ax.set_ylabel('Number of PTMs', color = 'blue')
    #change color of tick labels
    ax.tick_params(axis='y', labelcolor='blue')
    ax.set_xlabel(x_label_dict[filter_type])
    ax2 = ax.twinx()
    ax2.plot(x, frac_phospho, color = 'red')
    ax2.set_ylabel('Fraction that are\nphosphorylation sites', color = 'red')
    ax2.tick_params(axis='y', labelcolor='red')