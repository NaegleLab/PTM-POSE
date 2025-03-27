import numpy as np
import pandas as pd
import pickle

import os
from tqdm import tqdm

#plotting functions
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns

#custom stat functions and other helper functions
from ptm_pose import stat_utils, pose_config, helpers, annotate
import scipy.stats as stats



#optional analysis packages
try:
    import kinase_library as kl
    kinase_library = True
except ImportError:
    kinase_library = False

package_dir = os.path.dirname(os.path.abspath(__file__))


def compare_KL_for_sequence(inclusion_seq, exclusion_seq, dpsi = None, comparison_type = 'percentile'):
    """
    Given two sequences, compare the kinase library scores, percentiles, or ranks for each sequence. Optionally, provide a dPSI value to calculate the relative change in preference for each kinase.

    Parameters
    ----------
    inclusion_seq : str
        sequence to score for inclusion preference, with modification lowercased
    exclusion_seq : str
        sequence to score for exclusion preference, with modification lowercased
    dpsi : float
        dPSI value for the PTM event, which will be used to calculate the relative change in preference for each kinase (score difference * dPSI). Default is None.
    comparison_type : str
        type of comparison to perform. Can be 'percentile', 'score', or 'rank'. Default is 'percentile'.
    """

    #get percentiles for inclusion sequence
    s_inc = kl.Substrate(inclusion_seq)
    inn_perc = s_inc.percentile(sort_by = 'name')

    #get percentiles for exclusion sequence
    s_exc = kl.Substrate(exclusion_seq)

    if comparison_type == 'percentile':
        inc_results = s_inc.percentile(sort_by = 'name')
        ex_results = s_exc.percentile(sort_by = 'name')
    elif comparison_type == 'score':
        inc_results = s_inc.score()
        ex_results = s_exc.score()
    elif comparison_type == 'rank':
        inc_results = s_inc.rank()
        ex_results = s_exc.rank()

    results = pd.concat([inc_results, ex_results], axis = 1)
    results.columns = [f'Inclusion {comparison_type}', f'Exclusion {comparison_type}']
    results['Difference'] = results[f'Inclusion {comparison_type}'] - results[f'Exclusion {comparison_type}']
    results['Absolute Difference'] = results['Difference'].abs()
    results = results.sort_values(by = 'Absolute Difference', ascending = False)
    
    if dpsi is not None:
        results['dPSI'] = dpsi
        results['Relative Change in Preference'] = results['Difference'] * dpsi

    return results

def get_all_KL_scores(seq_data, seq_col, kin_type = ['ser_thr', 'tyrosine'], score_type = 'percentiles'):
    """
    Given a dataset with flanking sequences, score each flanking sequence

    Parameters
    ----------
    seq_data : pandas dataframe
        processed dataframe containing flanking sequences to score
    seq_col : str
        column in seq_data containing the flanking sequences
    kin_type : list
        list of kinase types to score. Can be 'ser_thr' or 'ST' for serine/threonine kinases, or 'tyrosine' or 'Y' for tyrosine kinases. Default is ['ser_thr', 'tyrosine'].
    score_type : str
        type of score to calculate. Can be 'percentile', 'score', or 'rank'. Default is 'percentile'.

    Returns
    -------
    merged_data : dict
        dictionary containing the merged dataframes for each kinase type with KinaseLibrary scores
    
    
    """
    phospho_predict = kl.PhosphoProteomics(seq_data, seq_col = seq_col)
    if isinstance(kin_type, list):
        if len(kin_type) == 1:
            phospho_predict.predict(kin_type = kin_type)
            merged_data = phospho_predict.merge_data_scores(kin_type = kin_type, score_type = score_type)
        else:
            merged_data = {}
            for ktype in kin_type:
                if ktype in ['ser_thr', 'tyrosine', 'Y', 'ST']:
                    if ktype == 'Y':
                        ktype == 'tyrosine'
                    elif ktype == 'ST':
                        ktype = 'ser_thr'
   
                    phospho_predict.predict(kin_type = ktype)
                    merged_data[ktype] = phospho_predict.merge_data_scores(kin_type = ktype, score_type = score_type)
                else:
                    print(f"{ktype} not recognized. Must be 'ser_thr' or 'ST' for serine/threonine, or 'tyrosine' or 'Y' for tyrosine")
    elif isinstance(kin_type, str):   
        phospho_predict.predict(metric = score_type, kin_type = kin_type)
        merged_data = phospho_predict.merge_data_scores(kin_type = kin_type, score_type = score_type)
    return merged_data



class KL_flank_analysis:
    def __init__(self, altered_flanks, metric = 'percentile', alpha = 0.05, min_dpsi = 0.2, **kwargs):
        """
        Class for comparing kinase preferences for inclusion and exclusion sequences from altered flanking sequences, using motif scoring from Kinase Library.

        Parameters
        ----------
        altered_flanks : pandas dataframe
            dataframe containing altered flanking sequences (output during projection), including columns for 'Region ID', 'Gene', 'Residue', 'PTM Position in Isoform', 'Inclusion Flanking Sequence', 'Exclusion Flanking Sequence', 'Modification Class', and 'dPSI'.
        metric : str
            metric to use for scoring. Can be 'percentile', 'score', or 'rank'. Default is 'percentile'.
        alpha : float
            significance threshold for p-value. Default is 0.05.
        min_dpsi : float
            effect size threshold for dPSI. Default is 0.2.
        **kwargs : dict
            additional keyword arguments to pass to the filter_ptms function
        """

        #check to make sure kinase library is installed
        if not kinase_library:
            print('Kinase library package not installed. To use this functionality, please install the kinase library package by running `pip install kinase-library`')
            return None
        #reduce to phosphorylation sites
        altered_flanks = altered_flanks[altered_flanks['Modification Class'] == 'Phosphorylation'].copy()

        #restrict to only significant PTMs
        filter_arguments = helpers.extract_filter_kwargs(alpha = alpha, min_dpsi = min_dpsi, modification_class = 'Phosphorylation', **kwargs)
        altered_flanks = helpers.filter_ptms(altered_flanks, **filter_arguments)

        #restrict to only necessary columns and remove duplicate entries (from canonical and alternative isoforms)
        optional_cols = ['dPSI']
        cols = ['Region ID', 'Gene', 'Residue', 'PTM Position in Isoform', 'Inclusion Flanking Sequence', 'Exclusion Flanking Sequence'] + [c for c in optional_cols if c in altered_flanks.columns]
        altered_flanks = altered_flanks.sort_values(by = 'Isoform Type', ascending = False)
        altered_flanks = altered_flanks[cols].drop_duplicates(keep = 'first')
        altered_flanks.reset_index(drop = True, inplace = True)

        #store result
        self.altered_flanks = altered_flanks.copy()
        self.metric = metric
        if 'dPSI' in altered_flanks.columns:
            self.optional_cols = ['dPSI']


        

    def analyze_single_ptm(self, gene, loc):
        """
        Score a single PTM for its flanking sequence in the inclusion and exclusion isoforms using Kinase-Library
        """
        ptm_data = self.altered_flanks[(self.altered_flanks['Gene'] == gene) & (self.altered_flanks['PTM Position in Isoform'] == loc)].copy()
        #remove duplicate entries with same splice event and .drop_duplicates(subset = ['Inclusion Flanking Sequence', 'Exclusion Flanking Sequence', 'Region ID'])
        if ptm_data.shape[0] == 0:
            raise ValueError(f"No flanking sequence data found for gene {gene} at position {loc}")
        elif ptm_data.shape[0] == 1:
            # If there's only one row, squeeze it to remove the extra dimension
            ptm_data = ptm_data.squeeze()
            dpsi = ptm_data['dPSI'] if 'dPSI' in ptm_data else None
            results = compare_KL_for_sequence(ptm_data['Inclusion Flanking Sequence'], ptm_data['Exclusion Flanking Sequence'], dpsi = dpsi)
            return results
        elif ptm_data.shape[0] > 1:
            # if multiple entries, perform analysis on all of them
            print(f"Multiple splice events impacting the flanking sequence of {ptm_data['Residue'].values[0]}{loc} in {gene} found. Performing analysis on all of them.")
            results = {}
            for i, row in ptm_data.iterrows():
                dpsi = row['dPSI'] if 'dPSI' in row else None
                results[row['Region ID']] = compare_KL_for_sequence(row['Inclusion Flanking Sequence'], row['Exclusion Flanking Sequence'], dpsi = dpsi)
        return results
    
    def score_all_ptms(self):
        """
        Score the flanking sequences of PTMs for both the inclusion and exclusion isoforms using Kinase-Library
        """
        inclusion_sequences = self.altered_flanks[[col for col in self.altered_flanks.columns if 'Exclusion' not in col]].copy()
        exclusion_sequences = self.altered_flanks[[col for col in self.altered_flanks.columns if 'Inclusion' not in col]].copy()

        print('Score sequences for inclusion isoforms')
        inclusion_percentile = get_all_KL_scores(inclusion_sequences, 'Inclusion Flanking Sequence')

        print('\nScoring sequences for exclusion isoforms')
        exclusion_percentile = get_all_KL_scores(exclusion_sequences, 'Exclusion Flanking Sequence')

        self.inclusion_percentile = inclusion_percentile
        self.exclusion_percentile = exclusion_percentile
            
    def melt_score_df(self, kin_type = 'ser_thr', flank_type = 'Inclusion'):
        """
        Melt the percentile dataframe to long format for easier plotting
        """
        if not hasattr(self, 'inclusion_percentile') or not hasattr(self, 'exclusion_percentile'):
            self.score_all_ptms()

        #check to make sure kin_type is valid
        if kin_type not in ['ser_thr', 'tyrosine']:
            raise ValueError(f"kin_type must be either 'ser_thr' or 'ST' for serine/threonine, or 'tyrosine'/'Y' for tyrosine, but got {kin_type}")
        elif kin_type == 'Y':
            kin_type = 'tyrosine'
        elif kin_type == 'ST':
            kin_type = 'ser_thr'

        #grab percentile data
        seq_col = f'{flank_type} Flanking Sequence'
        if flank_type == 'Inclusion':
            percentile_df = self.inclusion_percentile[kin_type].copy()
        else:
            percentile_df = self.exclusion_percentile[kin_type].copy()
        #remove unneeded columns
        percentile_df = percentile_df.drop(columns = ['phos_res', 'Sequence'])

        #melt data
        percentile_df = percentile_df.melt(id_vars = ['Region ID', 'Gene', 'Residue', 'PTM Position in Isoform', seq_col] + self.optional_cols, var_name = 'Kinase', value_name = self.metric.capitalize())
        return percentile_df

    def combine_score_dfs(self):
        """
        Combine the inclusion and exclusion score dataframes into a single dataframe that calculates differences in predicted preference between inclusion and exclusion isoforms
        """
        if not hasattr(self, 'inclusion_percentile') or not hasattr(self, 'exclusion_percentile'):
            self.score_all_ptms()

        melted_inclusion = pd.concat(
            [self.melt_percentile_df(flank_type = 'Inclusion', kin_type = 'ser_thr'), 
            self.melt_percentile_df(flank_type = 'Inclusion', kin_type = 'tyrosine')]
        )
        melted_inclusion = melted_inclusion.rename(columns = {'Percentile': f'Inclusion {self.score_type}'})

        melted_exclusion = pd.concat(
            [self.melt_percentile_df(flank_type = 'Exclusion', kin_type = 'ser_thr'), 
            self.melt_percentile_df(flank_type = 'Inclusion', kin_type = 'tyrosine')]
        )
        melted_exclusion = melted_exclusion.rename(columns = {'Percentile': f'Exclusion {self.score_type}'})
        combined = melted_inclusion.merge(melted_exclusion, on = ['Region ID', 'Gene', 'Residue', 'PTM Position in Isoform', 'Kinase'] + self.optional_cols, how = 'inner')
        
        #calculate the difference
        combined['Difference'] = combined[f'Inclusion {self.score_type}'] - combined[f'Exclusion {self.score_type}']
        combined['Absolute Difference'] = combined['Difference'].abs()
        combined = combined.sort_values(by = 'Absolute Difference', ascending = False)

        if 'dPSI' in combined.columns:
            combined['Relative Change in Preference'] = combined['Difference'] * combined['dPSI']
        
        self.summary_df = combined
        
    def analyze_all_ptms(self):
        """
        Perform the full analysis on all PTMs, including scoring and melting the dataframes
        """
        self.score_all_ptms()
        self.combine_score_dfs()

    def get_kinases_with_largest_changes(self, kinase_type = 'ST', top_n = 5, difference_type = 'relative', agg = 'median'):
        """
        Get the top n kinases with the largest changes in preference for inclusion or exclusion isoforms. Difference can be calculated as the normal difference in preference, the absolute difference in preference, or the relative change in preference (normalized by magnitude of splicing change, or dPSI).

        Parameters
        ----------
        kinase_type : str
            type of kinase to restrict analysis to. Can be 'ST' or 'ser_thr' for serine/threonine kinases, or 'Y' or 'tyrosine' for tyrosine kinases. Default is 'ST'.
        top_n : int
            number of top kinases to return. Default is 5.
        difference_type : str
            type of difference to calculate. Can be 'normal' (difference in preference), 'absolute' (absolute difference in preference), or 'relative' (relative change in preference normalized by magnitude of splicing change, or dPSI). Default is 'relative'.
        agg : str
            aggregation function to use across all interactions for calculating top kinases. Default is 'median'.
        """
        #grab kinase of interest
        if kinase_type == 'ST' or kinase_type == 'ser_thr':
            kinase_df = self.summary_df[self.summary_df['Residue'].str.contains('S|T')].copy()
        elif kinase_type == 'Y' or kinase_type == 'tyrosine':
            kinase_df = self.summary_df[self.summary_df['Residue'].str.contains('Y')].copy()

        if difference_type == 'normal':
            comp_col = 'Absolute Difference'
        elif difference_type == 'absolute':
            comp_col == 'Absolute Difference'
        elif difference_type == 'relative':
            comp_col = 'Relative Change in Preference'
            kinase_df['Absolute Relative'] = kinase_df['Relative Change in Preference'].abs()
            comp_col = 'Absolute Relative'

        else:
            raise ValueError('difference_type must be either "normal", "absolute", or "relative"')

        #check to make sure provided kinase type is valid
        if kinase_type not in ['ST','Y','ser_thr','tyrosine']:
            raise ValueError('kinase_type must be either "ST" or "ser_thr" for serine/threonine kinases or "Y" or "tyrosine" for tyrosine kinases')



        top_kinases = list(kinase_df.groupby('Kinase')[comp_col].agg(agg).sort_values(ascending = False).head(top_n).index)
        kinase_df = kinase_df[kinase_df['Kinase'].isin(top_kinases)]
        return kinase_df
    
    def plot_top_kinases(self, top_n = 5, kinase_type = 'ST', difference_type = 'relative', agg = 'median', metric_type = 'percentile', ax = None):
        """
        Plot a swarm plot showing the top n kinases with the largest changes in preference for inclusion or exclusion isoforms"
        """
        top_kinase_df = self.get_kinases_with_largest_changes(kinase_type = kinase_type, top_n = top_n, difference_type = difference_type, agg = agg)

        if ax is None:
            fig, ax = plt.subplots(figsize = (5, 4))

        #plot swarm plot
        y_data = 'Relative Change in Preference' if difference_type == 'relative' else 'Absolute Difference' if difference_type == 'absolute' else 'Difference'
        sns.swarmplot(data = top_kinase_df, x = 'Kinase', y = y_data, dodge = True, ax = ax)
        ax.axhline(0, color = 'black', linestyle = '--', linewidth = 0.5)
        if difference_type == 'relative':
            ax.set_ylabel('Relative Change\nin Preference\n(Change*dPSI)')
        elif difference_type == 'absolute':
            ax.set_ylabel(f'Absolute Difference\nin {metric_type.capitalize()}')
        elif difference_type == 'normal':
            ax.set_ylabel(f'Difference\nin {metric_type.capitalize()}')

    def plot_top_changes(self, gene = None, loc = None, top_n = 10, difference_type = 'relative', metric_type = 'percentile', ax = None):
        """
        Plot the top n changes in preference after an altered flanking sequence event, either for a all events, a specific gene, or a specific PTM. Top changes can be calculated by percentile change, absolute percentile change, or relative change in preference (percentile change in affinity*dPSI)

        Parameters
        ----------
        gene : str
            gene to restrict analysis to. Default is None, which will analyze all PTMs.
        loc : int
            location of PTM to restrict analysis to. Default is None, which will analyze all PTMs for the given gene or for all genes (if gene = None).
        top_n : int
            number of top changes to plot. Default is 10.
        difference_type : str
            type of difference to plot. Can be 'normal' (difference in preference), 'absolute' (absolute difference in preference), or 'relative' (relative change in preference normalized by magnitude of splicing change, or dPSI). Default is 'relative'.
        metric_type : str
            type of metric to plot. Default is 'percentile'.
        ax : matplotlib axis
            axis to plot on. Default is None, which will create a new figure and axis.
        """
        if not hasattr(self, 'summary_df'):
            self.analyze_all_ptms()


        if gene is None and loc is None:
            ptm = self.summary_df.copy()
            ptm['Kinase'] = ptm['Kinase'] + '->' + ptm['Gene'] + ' ' + ptm['Residue'] + ptm['PTM Position in Isoform'].astype(str)
        elif gene is not None:
            if loc is None:
                ptm = self.summary_df[self.summary_df['Gene'] == gene].copy()
                ptm['Kinase'] = ptm['Kinase'] + '->' + ptm['Gene'] + ' ' + ptm['Residue'] + ptm['PTM Position in Isoform'].astype(str)
            else:
                ptm = self.summary_df[(self.summary_df['Gene'] == gene) & (self.summary_df['PTM Position in Isoform'] == loc)].copy()

        if difference_type == 'normal':
            ptm = ptm.sort_values(by = 'Absolute Difference', ascending = False)
            diff_col = 'Difference'
            y_label = f'Difference\nin {metric_type.capitalize()}'
        elif difference_type == 'absolute':
            ptm = ptm.sort_values(by = 'Absolute Difference', ascending = False)
            diff_col = 'Absolute Difference'
            y_label = f'Absolute Difference\nin {metric_type.capitalize()}'
        elif difference_type == 'relative':
            ptm['Relative Absolute'] = ptm['Relative Change in Preference'].abs()
            ptm = ptm.sort_values(by = 'Relative Absolute', ascending = False)
            diff_col = 'Relative Change in Preference'
            y_label = 'Relative Change\nin Preference\n(Change*dPSI)'
        else:
            raise ValueError('difference_type must be either "normal", "absolute", or "relative"')

        plt_data = ptm.head(top_n).sort_values(by = diff_col, ascending = True)

        #color bars by sign
        colors = ['coral' if x > 0 else 'lightblue' for x in plt_data[diff_col]]
        if ax is None:
            fig, ax = plt.subplots(figsize = (3, 4))

        ax.barh(plt_data['Kinase'], plt_data[diff_col], color = colors, edgecolor = 'black', height=1)
        ax.set_xlabel(y_label)

        if gene is not None and loc is not None:
            ax.set_title(f'{gene} {loc} - Top {top_n} Changes')
        elif gene is not None:
            ax.set_title(f'{gene} - Top {top_n} Changes')
        else:
            ax.set_title(f'Top Changes')

        #construct legend
        if difference_type != 'absolute':
            handles = [plt.Rectangle((0,0),1,1, facecolor='coral', edgecolor='black'), plt.Rectangle((0,0),1,1, facecolor='lightblue', edgecolor='black')]
            labels = ['Increase in Preference\nupon Perturbation', 'Decrease in Preference \nupon Perturbation'] if difference_type == 'relative' else ['Prefers Inclusion Isoform', 'Prefers Exclusion Isoform']
            ax.legend(handles, labels, loc=(1.05,0.5), fontsize=8)


class KSEA:
    def __init__(self, ptms, database = 'OmniPath Writer', modification = None, min_dpsi = None, alpha = None, **kwargs):
        """
        Adapted version of the Kinase Substrate Enrichment Algorithm (KSEA) from Casado et al. 2013. This version is designed to work with the PTM-POSE pipeline and uses the OmniPath database as default for enzyme-substrate information.

        Parameters
        ----------
        ptms : pandas dataframe
            dataframe containing PTM information, including columns for 'PTM', 'Modification Class', 'p-value', and 'dPSI'. This dataframe should be the output from the PTM-POSE pipeline.
        alpha : float
            significance threshold for p-value. Default is 0.05.
        dpsi : float
            effect size threshold for dPSI. Default is 0.2.
        database : str
            database to use for enzyme-substrate information. Default is 'OmniPath Writer'. Other options include 'OmniPath Eraser', 'PhosphoSitePlus', or 'DEPOD'.
        modification : str
            modification type to restrict analysis to. Default is None, which includes all modifications. If provided, should be one of the modification classes in the 'Modification Class' column of the ptms dataframe.
        
        """

        #filter ptms if desired
        if min_dpsi is None:
            min_dpsi = 0
        if alpha is None:
            alpha = 1.1
        #filter ptms
        ptms = helpers.filter_ptms(ptms, modification_class = modification, min_dpsi = min_dpsi, alpha = alpha, **kwargs)

        #if omnipath database is used, make sure it was specified to either look at reader or writer enzymes
        if 'OmniPath' in database and database not in ['OmniPath Writer', 'OmniPath Eraser']:
            raise ValueError('OmniPath database must be either "OmniPath Writer" or "OmniPath Eraser"')
        elif database == 'OmniPath Writer':
            self.annot_col = 'OmniPath:Writer Enzyme'  #use writer enzyme annotations
        elif database == 'OmniPath Eraser':
            self.annot_col = 'OmniPath:Eraser Enzyme'
        elif database not in ['PhosphoSitePlus', 'DEPOD', 'RegPhos', 'OmniPath Writer', 'OmniPath Eraser']:
            raise ValueError('database must be either "OmniPath Writer", "OmniPath Eraser", "PhosphoSitePlus", "DEPOD", or "RegPhos"')
        else:
            self.annot_col = f'{database}:Enzyme'  #use enzyme annotations for other databases
        
        #load annotations
        if database == 'OmniPath Writer':
            self.annotations = annotate.process_database_annotations(database = 'OmniPath', annot_type = 'Writer Enzyme')
        elif database == 'OmniPath Eraser':
            self.annotations = annotate.process_database_annotations(database = 'OmniPath', annot_type = 'Eraser Enzyme')
        else:
            self.annotations = annotate.process_database_annotations(database = database, annot_type = 'Enzyme')

        self.enzymes = list(self.annotations.keys())

        #add labels to the ptms dataframe
        ptms = helpers.add_ptm_column(ptms)
        ptms['Label'] = ptms['PTM'] + '-' + ptms['Modification Class']
        self.ptms = ptms

    def runKSEA_singleenzyme(self, enzyme):
        """
        Given a single enzyme, calculate the zscore for the enzyme based on the dPSI values of the PTMs in the dataframe.

        Parameters
        ----------
        enzyme : str
            enzyme to calculate zscore for

        Returns
        -------
        z : float
            zscore for the enzyme
        m : int
            number of identified substrates for the enzyme in the dataset
        """
        annot = self.annotations[enzyme]
        #check the possible residues and modifications for the current enzyme (for establishing background)
        residues = set([x.split('_')[1][0] for x in annot])
        mod = set([x.split('-')[1] for x in annot])

        #filter the ptms dataframe to only include those that match the current enzyme
        tmp_ptms = self.ptms[(self.ptms['Residue'].isin(residues)) & (self.ptms['Modification Class'].isin(mod))]
        if tmp_ptms.empty:
            return np.nan, np.nan
        
        #get number of sites associated with enzyme
        m = tmp_ptms[tmp_ptms['Label'].isin(annot)].shape[0] 
        if m == 0:
            return np.nan, np.nan

        
        #calculate the mean dPSI for all residues
        p_mean = tmp_ptms['dPSI'].mean()
        #calculate the mean dPSI for all residues modified by the current enzyme
        s_mean = tmp_ptms[tmp_ptms['Label'].isin(annot)]['dPSI'].mean()
        #calculate the stdev of dPSI for all residues
        p_std = tmp_ptms['dPSI'].std()
        #calculate the zscore for the current enzyme
        z = ((s_mean - p_mean) * (m**0.5)) / (p_std )
        return z, m



    def runKSEA(self):
        """
        Perform KSEA-style analysis for each kinase with substrates identified in experiment using dPSI values as quantitative measure.

            
        Returns
        -------
        all_results: pandas dataframe
            melted dataframe including all KSEA results, including zscore, number of identified substrates (m), and false discovery rate
        zscore: pandas dataframe
            dataframe which contains zscores from all data columns, updated with new zscores calculated for data_col
        pvals: pandas dataframe
            dataframe which contains p-values from all data columns, updated with new p-values calculated for data_col
        FDR: pandas dataframe
            dataframe which contains false discover rate from all data columns, updated with new false discovery rate calculated for data_col
        m: pandas dataframe
            dataframe which contains number of identified substrates for each kianse from all data columns, updated with new m calculated for data_col
        """
        results = pd.DataFrame(np.nan, columns = ['z', 'p', 'FDR', 'm'], index = self.enzymes)
        #iterate through each possible enzyme
        for enzyme in self.annotations:
            #get the zscore and number of substrates for the current enzyme
            z, m = self.runKSEA_singleenzyme(enzyme)
            results.loc[enzyme, 'm'] = m
            results.loc[enzyme, 'z'] = z

        #calculate p-values and FDR for the results
        results = results.dropna(how = 'all')
        results['p'] = stats.norm.cdf(-abs(results['z'].values))
        results = results.sort_values(by = 'p', ascending = True)
        results['FDR'] = stat_utils.adjustP(results['p'].values, method = 'BH')
        self.results = results

    def plot_results(self, show_substrate_count = True, ax = None):
        """
        Create a horizontal barplot of KSEA results, coloring based on significance and annotating with the number of substrates (if desired)

        Parameters
        ----------
        show_substrate_count : bool
            whether to show the number of substrates on the right side of the plot. Default is True.
        ax : matplotlib axis
            axis to plot on. Default is None, which will create a new figure and axis.
        """
        #create horizontal barplot of KSEA results, coloring based on significance


        plt_data = self.results.copy()
        #sort results by zscore
        plt_data =  plt_data.sort_values('z')

        #plot results
        if ax is None:
            fig_height = len(plt_data)*0.2
            fig, ax = plt.subplots(figsize = (3.5, fig_height))

        colors = ['grey' if x > 0.05 else 'red' for x in plt_data['FDR']]
        ax.barh(plt_data.index, plt_data['z'], color = colors, edgecolor = 'black')
        ax.axvline(0, color = 'black')

        #set xlim to be equal on either side of 0
        xlim = round(plt_data['z'].abs().max())
        ax.set_xlim(-xlim, xlim)

        ax.set_xlabel('z-score')

        #add legend for significance
        legend_elements = [Patch(facecolor='grey', edgecolor='black', label='FDR > 0.05'), Patch(facecolor='red', edgecolor='black', label=r'$FDR \leq 0.05$')]
        ax.legend(handles = legend_elements)

        #annotate with number of substrates on the right side of plot
        if show_substrate_count:
            ax.text(xlim+0.1, len(plt_data)+1, 'm', va = 'center', ha = 'center', fontweight = 'bold')
            for i, v in enumerate(plt_data['m']):
                ax.text(xlim+0.1, i, str(int(v)), va = 'center', ha = 'center')
        
    def save_results(self, fname = 'ksea_results', odir = ''):
        """
        Saves zscores, FDR, and m for all data columns in seperate .tsv files
        """
        self.results.to_csv(f'{odir}{fname}.tsv', sep = '\t')


class kstar_enrichment:
    def __init__(self, spliced_ptms, network_dir, background_ptms = None, impact_type = ['All', 'Included', 'Excluded'], phospho_type = 'Y', **kwargs):
        """
        Given spliced ptm or PTMs with altered flanks and a single kstar network, get enrichment for each kinase in the network using a hypergeometric. Assumes the  data has already been reduced to the modification of interest (phosphotyrosine or phoshoserine/threonine)

        Parameters
        ----------
        network_dir : dict
            dictionary of networks with kinase-substrate information
        ptms : pandas dataframe
            differentially included ptms
        background_ptms: pd.DataFrame
            PTMs to consider as the background for enrichment purposes, which should overlap with the spliced ptms information provided (an example might be all identified events, whether or not they are significant). If not provided, will use all ptms in the phosphoproteome.
        impact_type: list of str or str
            type of impacts to perform enrichment analysis for. Can be 'All' (all significantly differentially included sites), 'Included' (sites with dPSI > 0), and/or 'Excluded' (sites with dPSI < 0). Default is to perform enrichment on all three groups.
        phospho_type : str 
            type of phosphorylation event to extract. Can either by phosphotyrosine ('Y') or phosphoserine/threonine ('ST'). Default is 'Y'.

        """
        #change phospho_type to list if not already
        if isinstance(phospho_type, list):
            pass
        elif phospho_type in ['Y', 'ST']:
            phospho_type = [phospho_type]
        elif phospho_type == 'All':
            phospho_type = ['Y', 'ST']
        else:
            raise ValueError('phospho_type must be either Y, ST, or All')

        if isinstance(impact_type, list):
            #if impact type includes inclusion and/or exclusion, check if dPSI is present and assign impact type accordingly
            if 'Included' in impact_type or 'Excluded' in impact_type:
                if 'dPSI' in spliced_ptms.columns:
                    spliced_ptms['Impact'] = spliced_ptms['dPSI'].apply(lambda x: 'Included' if x > 0 else 'Excluded')
                else:
                    raise ValueError('spliced_ptms must contain dPSI column to determine impact type, please either provide dPSI or restrict impact_type to "All"')
            self.impact_type = impact_type
        elif isinstance(impact_type, str):
            if impact_type == 'Included' or impact_type == 'Excluded':
                if 'dPSI' not in spliced_ptms.columns:
                    raise ValueError('spliced_ptms must contain dPSI column to determine impact type, please either provide dPSI or restrict impact_type to "All"')
                else:
                    spliced_ptms['Impact'] = spliced_ptms['dPSI'].apply(lambda x: 'Included' if x > 0 else 'Excluded')
                    self.impact_type = [impact_type]
            elif impact_type == 'All':
                self.impact_type = ['All']
            else:
                raise ValueError('impact_type must be at least one of either "All", "Included", or "Excluded"')
        
        self.significant_ptms = {}
        self.background_ptms = {}
        self.networks = {}
        for ptype in phospho_type:
            if ptype == 'Y':
                print('\nProcessing differentially included phosphotyrosine data')
            elif ptype == 'ST':
                print('\nProcessing differentially included phosphoserine/threonine data')
            #process ptms to only include specific phosphorylation data needed
            self.significant_ptms[ptype] = self.process_ptms(spliced_ptms, phospho_type = ptype, **kwargs)

            #load background, either from provided data or from all phosphoproteome data
            if ptype == 'Y':
                print('\nProcessing background phosphotyrosine data')
            elif ptype == 'ST':
                print('\nProcessing background phosphoserine/threonine data')
            if background_ptms is not None:
                self.background_ptms[ptype] = self.process_ptms(background_ptms, phospho_type=ptype, min_dpsi = 0, alpha = 1.1, **kwargs)
            else:
                self.background_ptms[ptype] = self.process_ptms(pose_config.ptm_coordinates.copy(), phospho_type = ptype, min_dpsi = 0, alpha = 1.1, **kwargs)

            #load kstar networks
            self.networks[ptype] = self.load_networks(network_dir, phospho_type = ptype)


        #save info
        self.phospho_type = phospho_type
        self.median_enrichment = None

    def load_networks(self, network_dir, phospho_type = 'Y'):
        """
        Given network directory, load all kstar networks for specific phosphorylation type
        """
        #check if file exists and whether a pickle has been generated: if not, load each network file individually
        if not os.path.exists(network_dir):
            raise ValueError('Network directory not found')
        elif os.path.exists(f"{network_dir}/*.p"):
            networks = pickle.load(open(f"{network_dir}/network_{phospho_type}.p", "rb" ) )
        else:
            network_directory = network_dir + f'/{phospho_type}/INDIVIDUAL_NETWORKS/'
            networks = {}
            for file in os.listdir(network_directory):
                if file.endswith('.tsv'):
                    #get the value of the network number
                    file_noext = file.strip(".tsv").split('_')
                    key_name = 'nkin'+str(file_noext[1])
                    #print("Debug: key name is %s"%(key_name))
                    networks[key_name] = pd.read_csv(f"{network_directory}{file}", sep='\t')
        return networks


    def process_ptms(self, ptms, phospho_type = 'Y', **kwargs):
        """
        Given ptm information, restrict data to include only the phosphorylation type of interest and add a PTM column for matching information from KSTAR

        Parameters
        ----------
        ptms: pd.DataFrame
            ptm information containing modification type and ptm locatin information, such as the output from projection or altered flanking sequence analysis
        phospho_type: str
            type of phosphorylation event to extract. Can either by phosphotyrosine ('Y') or phosphoserine/threonine ('ST')
        
        Returns
        ptms: pd.DataFrame
            trimmed dataframe containing only modifications of interest and new 'PTM' column
        """

        #restrict to ptms to phosphorylation type of interest
        ptms = ptms[(ptms['Modification Class'] == 'Phosphorylation')].copy()
        if phospho_type == 'Y':
            ptms = ptms[ptms['Residue'] == 'Y'].copy()
        elif phospho_type == 'ST':
            ptms = ptms[ptms['Residue'].isin(['S', 'T'])].copy()
        else:
            raise ValueError('phospho_type must be either Y or ST')
        
        ptms = helpers.filter_ptms(ptms, **kwargs)

        #construct PTM column that matches KSTAR information
        ptms['PTM'] = ptms['UniProtKB Accession'] + '_' + ptms['Residue'] + ptms['PTM Position in Isoform'].astype(int).astype(str)
        return ptms

    
    def get_enrichment_single_network(self, network_key, impact_type = 'All', phospho_type = 'Y'):
        """
        in progress
        """
        network = self.networks[phospho_type][network_key]
        network['PTM'] = network['KSTAR_ACCESSION'] + '_' + network['KSTAR_SITE']

        #add network information to all significant data
        if impact_type == 'All':
            sig_ptms = self.significant_ptms[phospho_type][['PTM']].drop_duplicates()
        elif impact_type == 'Included':
            sig_ptms = self.significant_ptms[phospho_type][self.significant_ptms[phospho_type]['Impact'] == 'Included'][['PTM']].drop_duplicates()
        elif impact_type == 'Excluded':
            sig_ptms = self.significant_ptms[phospho_type][self.significant_ptms[phospho_type]['Impact'] == 'Excluded'][['PTM']].drop_duplicates()
        else:
            raise ValueError('impact_type must be either "All", "Included", or "Excluded"')
        
        sig_ptms_kstar = sig_ptms.merge(network[['KSTAR_KINASE','PTM']], on = 'PTM')

        #repeat for background data
        bg_ptms = self.background_ptms[phospho_type][['PTM']].drop_duplicates()
        bg_ptms_kstar = bg_ptms.merge(network[['KSTAR_KINASE','PTM']], on = 'PTM')

        results = pd.DataFrame(np.nan, index = sig_ptms_kstar['KSTAR_KINASE'].unique(), columns = ['k','n','M','N','p'])
        for kinase in sig_ptms_kstar['KSTAR_KINASE'].unique():
            #get numbers for a hypergeometric test to look for enrichment of kinase substrates
            k = sig_ptms_kstar.loc[sig_ptms_kstar['KSTAR_KINASE'] == kinase, 'PTM'].nunique()
            n = bg_ptms_kstar.loc[bg_ptms_kstar['KSTAR_KINASE'] == kinase, 'PTM'].nunique()
            M = bg_ptms['PTM'].nunique()
            N = sig_ptms_kstar['PTM'].nunique()

            #run hypergeometric test
            results.loc[kinase,'p'] = stat_utils.hypergeom(M,n,N,k)
            results.loc[kinase, 'M'] = M
            results.loc[kinase, 'N'] = N
            results.loc[kinase, 'k'] = k
            results.loc[kinase, 'n'] = n

        return results
    
    def get_enrichment_all_networks(self, impact_type = 'All', phospho_type = 'Y'):
        """
        Given prostate data and a dictionary of kstar networks, get enrichment for each kinase in each network in the prostate data. Assumes the prostate data has already been reduced to the modification of interest (phosphotyrosine or phoshoserine/threonine)

        Parameters
        ----------
        networks : dict
            dictionary of kstar networks
        prostate : pandas dataframe
            all PTMs identified in tCGA prostate data, regardless of significance (reduced to only include mods of interest)
        sig_prostate : pandas dataframe
            significant PTMs identified in tCGA prostate data, p < 0.05 and effect size > 0.25 (reduced to only include mods of interest)
        """
        results = {}
        for network in self.networks[phospho_type]:
            results[network] = self.get_enrichment_single_network(network_key=network, impact_type=impact_type, phospho_type=phospho_type)
        return results

    def extract_enrichment(self, results):
        """
        Given a dictionary of results from get_enrichment_all_networks, extract the p-values for each network and kinase, and then calculate the median p-value across all networks for each kinase

        Parameters
        ----------
        results : dict
            dictionary of results from get_enrichment_all_networks
        """
        enrichment = pd.DataFrame(index = results['nkin0'].index, columns = results.keys())
        for network in results:
            enrichment[network] = results[network]['p']
        enrichment['median'] = enrichment.median(axis = 1)
        return enrichment
    
    def run_kstar_enrichment(self):
        """
        Run full kstar analysis to generate substrate enrichment across each of the 50 KSTAR networks and calculate the median p-value for each kinase across all networks
        """
        enrichment = {}
        median = {}
        for ptype in self.phospho_type:
            enrichment[ptype] = {}
            median_impacts = {}
            for impact in tqdm(self.impact_type, desc = f'Running enrichment for {ptype} data'):
                results = self.get_enrichment_all_networks(phospho_type=ptype, impact_type = impact)
                enrichment[ptype][impact] = self.extract_enrichment(results)
                median_impacts[impact] = enrichment[ptype][impact]['median']
            
            median[ptype] = pd.DataFrame(median_impacts)

        self.enrichment_all = enrichment
        self.median_enrichment = median

    def return_enriched_kinases(self, impact_type = 'All', alpha = 0.05):
        """
        Return kinases with a median p-value less than the provided alpha value (substrates are enriched among the significant PTMs)

        Parameters
        ----------
        alpha : float
            significance threshold to use to subset kinases. Default is 0.05.
        """
        if self.median_enrichment is None:
            self.run_kstar_enrichment()

        sig_enrichment = {}
        for ptype in self.phospho_type:
            tmp_data = self.median_enrichment[ptype][impact_type].copy()
            sig_enrichment[ptype] = tmp_data[tmp_data < alpha].index.values
        return sig_enrichment
    
    def dotplot(self, ptype = 'Y', impact_types = ['All', 'Included', 'Excluded'], kinase_axis = 'x', ax = None, facecolor = 'white', title = '', size_legend = False, color_legend = True, max_size = None, sig_kinases_only = True, alpha = 0.05, dotsize = 5, 
                 colormap={0: '#6b838f', 1: sns.color_palette('colorblind')[1]},
                 labelmap = {0: 'FPR > %0.2f'%(0.05), 1:'FPR <= %0.2f'%(0.05)},
                 legend_title = 'p-value', size_number = 5, size_color = 'gray', 
                 color_title = 'Significant', markersize = 10, 
                 legend_distance = 1.0, figsize = (5,5)):
        """
        Generates the dotplot plot, where size is determined by values dataframe and color is determined by significant dataframe. This is a stripped down version of the code used in KSTAR to generate the dotplot for the kinase activities
        
        Parameters
        -----------
        ax : matplotlib Axes instance, optional
            axes dotplot will be plotted on. If None then new plot generated
        """ 
        multiplier = 10
        offset = 5
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        ax.set_facecolor(facecolor)
        ax.set_title(title)

        plt_data = self.median_enrichment[ptype][impact_types]
        if isinstance(plt_data, pd.Series):
            plt_data = pd.DataFrame(plt_data)
        
        if sig_kinases_only:
            plt_data = plt_data[(plt_data < alpha).any(axis = 1)]
            if plt_data.empty:
                raise ValueError('No significant kinases found with p-value < %0.2f'%(alpha))
        
        if kinase_axis == 'x':
            plt_data = plt_data.T
        elif kinase_axis == 'y':
            plt_data = plt_data.copy()
        else:
            raise ValueError('kinase_axis must be either "x" or "y"')
        
        # Transform Data
        columns = list(plt_data.columns)
        values = -np.log10(plt_data).copy()
        colors = ((plt_data <= alpha) *1).copy()
        values['row_index'] = np.arange(len(values)) * multiplier + offset
        colors['row_index'] = np.arange(len(colors)) * multiplier + offset



        melt = values.melt(id_vars = 'row_index')
        values.drop(columns = ['row_index'], inplace = True)
        melt['var'] = melt.apply(lambda row : columns.index(row.iloc[1]) * multiplier + offset, axis = 1)
        
        melt_color = colors.melt(id_vars = 'row_index')
        melt_color['var'] = melt_color.apply(lambda row : columns.index(row.iloc[1]) * multiplier + offset, axis = 1)
        colors.drop(columns = ['row_index'], inplace = True)

        # Plot Data
        x = melt['var']
        y = melt['row_index'][::-1]    #needs to be done in reverse order to maintain order in the dataframe
        
        
        s = melt.value * dotsize
        
        #check to see if more than 2 values are given (fprs). Otherwise get color based on binary significance

        #get color for each datapoint based on significance
        melt_color['color'] = [colormap.get(l,'black') for l in melt_color.value]

            
        c = melt_color['color']
        scatter = ax.scatter(x, y, c=c, s=s)
        
        # Add Color Legend
        if color_legend:
            #create the legend
            color_legend = []
            for color_key in colormap.keys():
                color_legend.append(
                    Line2D([0], [0], marker='o', color='w', label=labelmap[color_key],
                            markerfacecolor= colormap[color_key], markersize=markersize),
                )     
            legend1 = ax.legend(handles=color_legend, loc=f'upper right', bbox_to_anchor=(legend_distance,1), title = color_title)  

            legend1.set_clip_on(False)
            ax.add_artist(legend1)
            


        # Add Size Legend
        if size_legend:
            #check to see if max pval parameter was given: if so, use to create custom legend
            if max_size is not None:
                s_label = np.arange(max_size/size_number,max_size+1,max_size/size_number).astype(int)
                dsize = [s*dotsize for s in s_label]
                legend_elements = []
                for element, s in zip(s_label, dsize):
                    legend_elements.append(Line2D([0],[0], marker='o', color = 'w', markersize = s**0.5, markerfacecolor = size_color, label = element))
                legend2 = ax.legend(handles = legend_elements, loc = f'lower right', title = legend_title, bbox_to_anchor=(legend_distance,0))        
            else:
                kw = dict(prop="sizes", num=size_number, color=size_color, func=lambda s: s/dotsize) 
                legend2 = ax.legend(*scatter.legend_elements(**kw),
                        loc=f'lower right', title=legend_title, bbox_to_anchor=(legend_distance,0)) 
            ax.add_artist(legend2)

        
        # Add Additional Plotting Information
        ax.tick_params(axis = 'x', rotation = 90)
        ax.yaxis.set_ticks(np.arange(len(values)) * multiplier + offset)
        ax.xaxis.set_ticks(np.arange(len(columns)) * multiplier + offset)


        ax.set_xticklabels(plt_data.columns)
        ax.set_yticklabels(plt_data.index[::-1])  #reverse order to match the dataframe
        
        #adjust x and y scale so that data is always equally spaced
        ax.set_ylim([0,len(values)*multiplier])
        ax.set_xlim([0,len(columns)*multiplier])
        return ax 
    
 

        


        
      
        