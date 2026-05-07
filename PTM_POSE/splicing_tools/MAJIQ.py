
from ptm_pose import helpers, pose_config, project
from ptm_pose.splicing_tools.base import GenericDataset
import pandas as pd


def project_ptms_onto_MAJIQ(voila_tsv, samp1_name, samp2_name, alpha = 0.05, dpsi_cutoff = 0.1, coordinate_type = 'hg38', ptm_coordinates = None, **kwargs):
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
    dpsi_cutoff: float
        Delta PSI cutoff for filtering out non-changing lsvs. Default is 0.15.
    coordinate_type: str
        indicates the coordinate system used for the start and end positions. Either hg38 or hg19. Default is 'hg38'.
    ptm_coordinates: pandas.DataFrame
        DataFrame containing PTM coordinates and information. If not provided, will use the default PTM coordinates from the pose_config module.
    **kwargs: additional keyword arguments
        Additional keyword arguments to pass to the find_ptms_in_many_regions function, which will be fed into the `filter_ptms()` function from the helper module. These will be used to filter ptms with lower evidence. For example, if you want to filter PTMs based on the number of MS observations, you can add 'min_MS_observations = 2' to the kwargs. This will filter out any PTMs that have less than 2 MS observations. See the `filter_ptms()` function for more options.

    Returns
    -------
    spliced_ptms: pandas.DataFrame
        DataFrame containing PTMs that are found in the differentially included regions of the splice events.
    majiq: pandas.DataFrame
        DataFrame containing the MAJIQ splice quantification data, with additional columns for the PTMs that are found in the differentially included regions.
    """
    #load ptm data from config if not provided
    if ptm_coordinates is None:
        ptm_coordinates = pose_config.ptm_coordinates.copy()
    

    #check for any keyword arguments to use for filtering
    if kwargs:
        filter_arguments = helpers.extract_filter_kwargs(**kwargs)
        #check any excess unused keyword arguments, report them
        helpers.check_filter_kwargs(filter_arguments)
        #filter ptm coordinates file to include only ptms with desired evidence
        ptm_coordinates = helpers.filter_ptms(ptm_coordinates, **filter_arguments)



    majiq = pd.read_csv(voila_tsv, sep = '\t', header = 10)

    lsv_cols = ['mean_dpsi_per_lsv_junction', 'probability_changing', 'probability_non_changing', f'{samp1_name}_mean_psi', f'{samp2_name}_mean_psi', 'de_novo_junctions', 'junctions_coords']

    #separate lsvs with multiple junctions into unique rows
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

    #filter out non-changing lsvs
    majiq['probability_non_changing']
    majiq = majiq[majiq['probability_non_changing'] < alpha]
    majiq = majiq[majiq['mean_dpsi_per_lsv_junction'].abs() > dpsi_cutoff]

    majiq, spliced_ptms = project.project_ptms_onto_splice_events(majiq, ptm_coordinates = ptm_coordinates, chromosome_col = 'seqid', strand_col = 'strand', region_start_col = 'junction_start', region_end_col = 'junction_end', coordinate_type = coordinate_type, gene_col = 'gene_name', dPSI_col = 'mean_dpsi_per_lsv_junction', event_id_col = 'lsv_id', sig_col = 'probability_non_changing', extra_cols = [f'{samp1_name}_mean_psi', f'{samp2_name}_mean_psi', 'probability_changing', 'exons_coords'])
    return spliced_ptms, majiq


class MAJIQ_Dataset(GenericDataset):
    def __init__(self, voila_tsv_file, samp1_name, samp2_name, min_dpsi = 0, alpha = 0.05, dpsi_col = 'mean_dpsi_per_lsv_junction', coordinate_type = 'hg38'):
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
        #check for any keyword arguments to use for filtering
        self.project_ptms_generic(extra_cols = extra_cols, PROCESSES = PROCESSES, **kwargs)

    def run_nease(self):
        self.run_nease_generic()