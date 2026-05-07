
from ptm_pose import helpers, pose_config, project
from ptm_pose.splicing_tools.base import GenericDataset
import pandas as pd


class MAJIQ_Dataset(GenericDataset):
    """
    Given splice quantification from the MAJIQ algorithm, annotate with PTMs that are found in the differentially included regions. 

    Parameters
    ----------
    voila_tsv: str
        Path to the MAJIQ output TSV file, which contains the splice quantification data.
    samp1_name: str
        Name of the first sample in the MAJIQ output. This will be used to identify the sample in MAJIQ output and label the columns in the output dataframe.
    samp2_name: str
        Name of the second sample in the MAJIQ output. This will be used to identify the sample in MAJIQ output and label the columns in the output dataframe.
    alpha: float
        Significance threshold for filtering out non-changing lsvs. Default is 0.05.
    min_dpsi: float
        Delta PSI cutoff for filtering out non-changing lsvs. Default is 0.15.
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.


    Attributes
    ----------
    splice_data: pandas.DataFrame
        processed DataFrame containing the MAJIQ splice quantification data with each row representing a unique junction of an lsv
    """
    def __init__(self, voila_tsv_file, samp1_name, samp2_name, min_dpsi = 0, alpha = 0.05, coordinate_type = 'hg38'):
        majiq = pd.read_csv(voila_tsv_file, sep = '\t', header = 10)
        #separate lsvs with multiple junctions into unique rows
        lsv_cols = ['mean_dpsi_per_lsv_junction', 'probability_changing', 'probability_non_changing', f'{samp1_name}_mean_psi', f'{samp2_name}_mean_psi', 'de_novo_junctions', 'junctions_coords']
        for col in lsv_cols:
            majiq[col] = majiq[col].str.split(';')
        majiq = majiq.explode(column = lsv_cols)

        #convert to numeric
        for col in lsv_cols:
            if col != 'junctions_coords' and col != 'de_novo_junctions':
                majiq[col] = pd.to_numeric(majiq[col])

        #extract start and stop of junctions from junctions_coords
        majiq['junction_start'] = majiq['junctions_coords'].str.split('-').str[0].astype(int)
        majiq['junction_end'] = majiq['junctions_coords'].str.split('-').str[1].astype(int)


        super().__init__(splice_data=majiq, min_dpsi=min_dpsi, alpha=alpha, coordinate_type=coordinate_type, chromosome_col = 'seqid', strand_col = 'strand', region_start_col = 'junction_start', region_end_col = 'junction_end', dpsi_col = 'mean_dpsi_per_lsv_junction', sig_col = 'probability_non_changing', event_id_col = 'lsv_id', start_coordinate_system = '1-based')

        self.samp1_name = samp1_name
        self.samp2_name = samp2_name

    def run_pose(self, extra_cols = None, PROCESSES = 1, **kwargs):
        """
        Project PTMs onto the splice events in the dataset. This will identify PTMs that are found in the differentially included regions of the splice events. Currently cannot identify altered flanking sequences from data. Independent of NEASE analysis, which can be run separately with the run_nease function.

        Parameters
        ----------
        extra_cols: list of str
            List of additional columns from the splice_data dataframe to include in the output spliced_ptms dataframe. These columns will be included in the spliced_ptms dataframe and can be used for filtering or analysis. Default is None, which means no additional columns will be included.
        PROCESSES: int
            Number of processes to use for parallel processing. Default is 1, which means no parallel processing will be used.
        **kwargs: additional keyword arguments
            Additional keyword arguments to use for filtering the splice events/PTMs. These will be passed to the project_ptms_generic function, which will use them to filter the splice events before projecting PTMs onto them. For example, if you want to filter for only exon skipping events, you can pass event_type = 'exon_skipping' as a keyword argument. The available keyword arguments will depend on the implementation of the project_ptms_generic function, but some examples include event_type, gene_name, and dpsi_cutoff.
        """
        #check for any keyword arguments to use for filtering
        self.project_ptms_generic(extra_cols = extra_cols, PROCESSES = PROCESSES, **kwargs)

    def run_nease(self):
        """
        Run NEASE analysis for the splice events in the dataset. Independent of POSE analysis.
        """
        self.run_nease_generic()