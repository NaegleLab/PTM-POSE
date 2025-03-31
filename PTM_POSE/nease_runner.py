import pandas as pd
import numpy as np

from ptm_pose import pose_config, helpers
try:
    import nease
except ImportError:
    nease = None

def process_data_for_nease(splice_data, region_start_col = 'EXON START', gene_col = 'Gene', chromosome_col = 'chr', strand_col = 'strand', gene_id_type = 'name', region_end_col = 'EXON_END', dpsi_col = 'dPSI', coordinate_type = 'hg38'):
    """
    Process the data for NEASE analysis

    Parameters
    ----------
    splice_data : pd.DataFrame
        dataframe containing splice event/isoform information, including the genomic coordinates of event/isofrom
    region_start_col : str
        
        column name for the start position of the region (default is 'EXON START')
    region_end_col : str
        
        column name for the end position of the region (default is 'EXON_END')
    gene_col : str
        
        column name for the gene identifier (default is 'Gene'). This should be either the gene name, Ensembl gene ID, or UniProt ID, and this should match the ID type specified in the gene_id_type parameter.
    gene_id_type : str
        
        type of gene identifier used in the gene_col. Options are 'Ensembl', 'name', or 'uniprot'. Default is 'name'.
    chromosome_col : str
        column name containing the chromosome of event. Default is 'chr'
    strand_col : str
        column name containing the strand of the event. Default is 'strand'
    dpsi_col : str or None
        column name containing the deltaPSI of the event (optional, set to None if not provided)
    coordinate_type : str
        coordinate system of the genomic coordinates provided in the splice_data. Options are 'hg38', 'hg19', 'hg18'. Default is 'hg38'.
    
    """
    if nease is None:
        raise ImportError("NEASE module is not installed. Please install it to use this module.")


    #if gene id is not Ensembl gene id, convert it to Ensembl gene id
    if gene_id_type == 'Ensembl':
        #check if gene ids are actually Ensembl IDs
        if not splice_data[gene_col].values[0].startswith('ENSG'):
            raise ValueError("Gene IDs are not Ensembl IDs, please either convert them to Ensembl IDs or specify gene_id_type as 'name' or 'uniprot'")
        if 'Gene ensembl ID' in splice_data.columns:
            splice_data = splice_data.drop(columns = ['Gene ensembl ID'])
        splice_data = splice_data.rename(columns={gene_col: 'Gene ensembl ID'})
    elif gene_id_type == 'name':
        #convert gene names to Ensembl gene ids (first converting to uniprot ids)
        if not hasattr(pose_config, 'genename_to_uniprot'):
            pose_config.genename_to_uniprot = pose_config.flip_uniprot_dict(pose_config.uniprot_to_genename)
        #convert to uniprot
        splice_data['Gene'] = splice_data[gene_col].map(pose_config.genename_to_uniprot)
        #then convert to Ensembl gene ids
        splice_data['Gene ensembl ID'] = splice_data['Gene'].apply(lambda x: pose_config.uniprot_to_geneid[x].split(' ')[0] if x in pose_config.uniprot_to_geneid else np.nan)
        splice_data = splice_data.dropna(subset = ['Gene ensembl ID'])
    elif gene_id_type == 'uniprot':
        #convert uniprot ids to Ensembl gene ids
        splice_data['Gene ensembl ID'] = splice_data[gene_col].apply(lambda x: pose_config.uniprot_to_geneid[x].split(' ')[0])
    else:
        raise ValueError("gene_id_type must be 'Ensembl', 'name', or 'uniprot'")


    #if coord type is not hg38, convert to hg38
    if coordinate_type != 'hg38':
        print('Converting genomic coordinates to hg38')
        #convert genomic coordinates to hg38 for region start
        splice_data = helpers.convert_genomic_coordinates_df(splice_data, from_coord = coordinate_type, to_coord = 'hg38', loc_col = region_start_col, chromosome_col = chromosome_col, strand_col = strand_col, output_col=region_start_col)
        #convert genomic coordinates to hg38 for region end
        splice_data = helpers.convert_genomic_coordinates_df(splice_data, from_coord = coordinate_type, to_coord = 'hg38', loc_col = region_end_col, chromosome_col = chromosome_col, strand_col = strand_col, output_col=region_end_col)

        original_shape = splice_data.shape[0]
        #remove errors in conversion
        splice_data = splice_data[~(splice_data[region_start_col] == -1) & ~(splice_data[region_end_col] == -1)]
        splice_data = splice_data.replace([np.inf, -np.inf], np.nan) #replace inf values with nan
        #remove rows with NaN values
        splice_data = splice_data.dropna(subset = [region_start_col, region_end_col])
        if original_shape != splice_data.shape[0]:
            print(f"Removed {original_shape - splice_data.shape[0]} rows that failed genomic coordinate conversion")

    #extract only needed columns
    if dpsi_col is not None:
        splice_data = splice_data[['Gene ensembl ID', region_start_col, region_end_col, dpsi_col]].copy()
            #rename columns
        splice_data = splice_data.rename(columns={region_start_col: 'EXON START', region_end_col: 'EXON END', dpsi_col: 'dPSI'})
    else:
        splice_data = splice_data[['Gene ensembl ID', region_start_col, region_end_col]].copy()
            #rename columns
        splice_data = splice_data.rename(columns={region_start_col: 'EXON START', region_end_col: 'EXON END'})

    #force datatypes
    splice_data = splice_data.astype({'EXON START': 'str', 'EXON END': 'str','Gene ensembl ID':'str'})
    
    return splice_data

def run_nease(splice_data, region_start_col = 'EXON START', gene_col = 'Gene', chromosome_col = 'chr', strand_col = 'strand', gene_id_type = 'name', region_end_col = 'EXON_END', dpsi_col = 'dPSI', coordinate_type = 'hg38', remove_non_in_frame = False, only_divisible_by_3 = False):
    """
    Given a dataframe containing splice event/isoform information, process and run NEASE analysis

    
    Parameters
    ----------
    splice_data : pd.DataFrame
        dataframe containing splice event/isoform information, including the genomic coordinates of event/isofrom
    region_start_col : str
        
        column name for the start position of the region (default is 'EXON START')
    region_end_col : str
        
        column name for the end position of the region (default is 'EXON_END')
    gene_col : str
        
        column name for the gene identifier (default is 'Gene'). This should be either the gene name, Ensembl gene ID, or UniProt ID, and this should match the ID type specified in the gene_id_type parameter.
    gene_id_type : str
        
        type of gene identifier used in the gene_col. Options are 'Ensembl', 'name', or 'uniprot'. Default is 'name'.
    chromosome_col : str
        column name containing the chromosome of event. Default is 'chr'
    strand_col : str
        column name containing the strand of the event. Default is 'strand'
    dpsi_col : str or None
        column name containing the deltaPSI of the event (optional, set to None if not provided)
    coordinate_type : str
        coordinate system of the genomic coordinates provided in the splice_data. Options are 'hg38', 'hg19', 'hg18'. Default is 'hg38'.
    remove_non_in_frame : bool
        whether to remove events that are expected to not be in frame. Default is False. If True, only events that are in frame will be included in the output.
    only_divisible_by_3 : bool
        whether to remove events that are not multiple of 3. Default is False.
    """
    splice_data = process_data_for_nease(splice_data, region_start_col = region_start_col, gene_col = gene_col, chromosome_col = chromosome_col, strand_col = strand_col, gene_id_type = gene_id_type, region_end_col = region_end_col, dpsi_col = dpsi_col, coordinate_type = coordinate_type)
    nease_output = nease.run(splice_data, organism = 'Human', remove_non_in_frame = remove_non_in_frame, only_divisible_by_3 = only_divisible_by_3)
    return nease_output

def save_nease(nease_output, odir, file_type = 'excel'):

    """
    Write to excel file with multiple tabs

    Parameters
    ----------
    nease_ouput : nease object
        output from nease.run() / ptm_pose.nease_runner.run_nease()
    odir : str
        output directory
    file_type : str
        file type to save as, 'excel', 'csv', or 'tsv'. Default is excel, which saves each type of data in a separate sheet
    """
    #interactions
    try:
        interactions = nease_output.get_edges()
    except AttributeError:
        interactions = pd.DataFrame()
    
    try:
        slims = nease_output.get_slims()
    except AttributeError:
        slims = pd.DataFrame()

    try:
        residues = nease_output.get_residues()
    except AttributeError:
        residues = pd.DataFrame()

    try:
        domains = nease_output.get_domains()
    except AttributeError:
        domains = pd.DataFrame()

    if file_type == 'excel':
        with pd.ExcelWriter(odir + 'nease_output.xlsx') as writer:
            interactions.to_excel(writer, sheet_name='Interactions', index = False)
            slims.to_excel(writer, sheet_name='ELM', index = False)
            residues.to_excel(writer, sheet_name='PDB', index = False)
            domains.to_excel(writer, sheet_name='Domains', index = False)
    elif file_type == 'csv':
        interactions.to_csv(odir + 'nease_interactions.csv', index = False)
        slims.to_csv(odir + 'nease_slims.csv', index = False)
        residues.to_csv(odir + 'nease_residues.csv', index = False)
        domains.to_csv(odir + 'nease_domains.csv', index = False)
    elif file_type == 'tsv':
        interactions.to_csv(odir + 'nease_interactions.tsv', sep = '\t', index = False)
        slims.to_csv(odir + 'nease_slims.tsv', sep = '\t', index = False)
        residues.to_csv(odir + 'nease_residues.tsv', sep = '\t', index = False)
        domains.to_csv(odir + 'nease_domains.tsv', sep = '\t', index = False)
    else:
        raise ValueError("file_type must be 'excel', 'csv', or 'tsv'")
    
def load_nease(odir, file_type = 'excel'):
    """
    Load nease output

    Parameters
    ----------
    odir : str
        location of nease output files
    file_type : str
        type of file nease output is saved as. Options are 'excel', 'csv', or 'tsv'. Default is 'excel'.
    """
    if file_type == 'excel':
        interactions = pd.read_excel(odir + 'nease_output.xlsx', sheet_name='Interactions')
        slims = pd.read_excel(odir + 'nease_output.xlsx', sheet_name='ELM')
        residues = pd.read_excel(odir + 'nease_output.xlsx', sheet_name='PDB')
        domains = pd.read_excel(odir + 'nease_output.xlsx', sheet_name='Domains')
    elif file_type == 'csv':
        interactions = pd.read_csv(odir + 'nease_interactions.csv')
        slims = pd.read_csv(odir + 'nease_slims.csv')
        residues = pd.read_csv(odir + 'nease_residues.csv')
        domains = pd.read_csv(odir + 'nease_domains.csv')
    elif file_type == 'tsv':
        interactions = pd.read_csv(odir + 'nease_interactions.tsv', sep = '\t')
        slims = pd.read_csv(odir + 'nease_slims.tsv', sep = '\t')
        residues = pd.read_csv(odir + 'nease_residues.tsv', sep = '\t')
        domains = pd.read_csv(odir + 'nease_domains.tsv', sep = '\t')


    return interactions, slims, residues, domains





        
