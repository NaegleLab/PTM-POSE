import pandas as pd
import numpy as np
import re
import os

from ptm_pose import pose_config, helpers


#dictionaries for converting modification codes to modification names in PhosphoSitePlus data
mod_shorthand_dict = {'p': 'Phosphorylation', 'ca':'Caspase Cleavage', 'hy':'Hydroxylation', 'sn':'S-Nitrosylation', 'ng':'Glycosylation', 'ub': 'Ubiquitination', 'pa': "Palmitoylation",'ne':'Neddylation','sc':'Succinylation', 'sm': 'Sumoylation', 'ga': 'Glycosylation', 'gl': 'Glycosylation', 'ac': 'Acetylation', 'me':'Methylation', 'm1':'Methylation', 'm2': 'Dimethylation', 'm3':'Trimethylation'}
residue_dict = {'P': 'proline', 'Y':'tyrosine', 'S':'serine', 'T':'threonine', 'H':'histidine', 'D':'aspartic acid', 'I':'isoleucine', 'K':'lysine', 'R':'arginine', 'G':'glycine', 'N':'asparagine', 'M':'methionine'}
annotation_col_dict = {'PhosphoSitePlus':{'Function':'PSP:ON_FUNCTION', 'Process':'PSP:ON_PROCESS', 'Interactions':'PSP:ON_PROT_INTERACT', 'Disease':'PSP:Disease_Association', 'Kinase':'PSP:Kinase','Perturbation':'PTMsigDB:PSP-PERT'},
                        'ELM':{'Interactions':'ELM:Interactions', 'Motif Match':'ELM:Motif Matches'},
                        'PTMcode':{'Intraprotein':'PTMcode:Intraprotein_Interactions', 'Interactions':'PTMcode:Interprotein_Interactions'},
                        'PTMInt':{'Interactions':'PTMInt:Interactions'},
                        'RegPhos':{'Kinase':'RegPhos:Kinase'},
                        'DEPOD':{'Phosphatase':'DEPOD:Phosphatase'},
                        'PTMsigDB': {'WikiPathway':'PTMsigDB:PATH-WP', 'NetPath':'PTMsigDB:PATH-NP','mSigDB':'PTMsigDB:PATH-BI', 'Pertubation (DIA2)':'PTMsigDB:PERT-P100-DIA2', 'Perturbation (DIA)': 'PTMsigDB:PERT-P100-DIA', 'Perturbation (PRM)':'PTMsigDB:PERT-P100-PRM', 'Kinase':'PTMsigDB:Kinase-iKiP'}}



def add_custom_annotation(spliced_ptms, annotation_data, source_name, annotation_type, annotation_col, accession_col = 'UniProtKB Accession', residue_col = 'Residue', position_col = 'PTM Position in Canonical Isoform'):
    """
    Add custom annotation data to spliced_ptms or altered flanking sequence dataframes

    Parameters
    ----------
    annotation_data: pandas.DataFrame
        Dataframe containing the annotation data to be added to the spliced_ptms dataframe. Must contain columns for UniProtKB Accession, Residue, PTM Position in Canonical Isoform, and the annotation data to be added
    source_name: str
        Name of the source of the annotation data, will be used to label the columns in the spliced_ptms dataframe
    annotation_type: str
        Type of annotation data being added, will be used to label the columns in the spliced_ptms dataframe
    annotation_col: str
        Column name in the annotation data that contains the annotation data to be added to the spliced_ptms dataframe
    

    Returns
    -------
    spliced_ptms: pandas.DataFrame
        Contains the PTMs identified across the different splice events with an additional column for the custom annotation data
    """
    #check if annotation data contains the annotation col
    if isinstance(annotation_col, str):
        if annotation_col not in annotation_data.columns:
            raise ValueError(f'Could not find column indicated to contain {annotation_col} in annotation data. Please either change the name of your annotation data column with this information or indicate the correct column name with the annotation_col parameter')
        else:
            #make annotation col name based on source and annotation type
            annotation_col_name = source_name + ':' + annotation_type
            annotation_data = annotation_data.rename(columns = {annotation_col: annotation_col_name})
    else:
        raise ValueError('annotation_col must be a string indicating column with annotation data to be added to the spliced_ptms dataframe')

    #check to make sure annotation data has the necessary columns
    if not all([x in annotation_data.columns for x in [accession_col, residue_col, position_col]]):
        raise ValueError(f'Could not find columns containing ptm information: {accession_col}, {residue_col}, and {position_col}. Please either change the name of your annotation data columns containing this information or indicate the correct column names with the accession_col, residue_col, and position_col parameters')

    #if splice data already has the annotation columns, remove them
    if annotation_col_name in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = [annotation_col_name])

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(annotation_data, how = 'left', left_on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'], right_on = [accession_col, residue_col, position_col])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms or annotation dataframe')
    
    #report the number of PTMs identified
    num_ptms_with_custom_data = spliced_ptms.dropna(subset = annotation_col).groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"{source_name} {annotation_type} data added: {num_ptms_with_custom_data} PTMs in dataset found with {source_name} {annotation_type} information")

    return spliced_ptms

def add_PSP_regulatory_site_data(spliced_ptms, file = 'Regulatory_sites.gz', report_success = True):
    """
    Add functional information from PhosphoSitePlus (Regulatory_sites.gz) to spliced_ptms dataframe from project_ptms_onto_splice_events() function

    Parameters
    ----------
    file: str
        Path to the PhosphoSitePlus Regulatory_sites.gz file. Should be downloaded from PhosphoSitePlus in the zipped format

    Returns
    -------
    spliced_ptms: pandas.DataFrame
        Contains the PTMs identified across the different splice events with additional columns for regulatory site information, including domains, biological process, functions, and protein interactions associated with the PTMs
    """
    #check to make sure file exists
    check_file(file, expected_extension='.gz')

    #read in the kinase substrate data and add to spliced ptm info
    regulatory_site_data = pd.read_csv(file, sep = '\t', header = 2, on_bad_lines='skip',compression = 'gzip')
    regulatory_site_data = regulatory_site_data.rename(columns = {'ACC_ID':'UniProtKB Accession'})
    #drop extra modification information that is not needed
    regulatory_site_data['Residue'] = regulatory_site_data['MOD_RSD'].apply(lambda x: x.split('-')[0][0])
    regulatory_site_data['PTM Position in Canonical Isoform'] = regulatory_site_data['MOD_RSD'].apply(lambda x: int(x.split('-')[0][1:]))
    #add modification type
    regulatory_site_data['Modification Class'] = regulatory_site_data['MOD_RSD'].apply(lambda x: mod_shorthand_dict[x.split('-')[1]])

    #restrict to human data
    regulatory_site_data = regulatory_site_data[regulatory_site_data['ORGANISM'] == 'human']
    regulatory_site_data = regulatory_site_data[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'ON_PROCESS', 'ON_PROT_INTERACT', 'ON_OTHER_INTERACT', 'ON_FUNCTION']].drop_duplicates()
    
    #group like modifications into a single column
    regulatory_site_data = regulatory_site_data.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).agg(lambda x: '; '.join([y for y in x if y == y])).reset_index()
    regulatory_site_data = regulatory_site_data.replace('', np.nan)
    
    #add 'PSP:' in front of each column
    regulatory_site_data.columns = ['PSP:' + x if x not in ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'] else x for x in regulatory_site_data.columns]
    
    #if splice data already has the annotation columns, remove them
    if 'PSP:ON_FUNCTION' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['PSP:ON_FUNCTION', 'PSP:ON_PROCESS', 'PSP:ON_PROT_INTERACT', 'PSP:ON_OTHER_INTERACT'])

    #explode dataframe on modifications
    if spliced_ptms['Modification Class'].str.contains(';').any():
        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)

    #merge with spliced_ptm info
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(regulatory_site_data, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')

    
    #report the number of ptms with motif data
    if report_success:
        num_ptms_with_known_function = spliced_ptms.dropna(subset = 'PSP:ON_FUNCTION').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_known_process = spliced_ptms.dropna(subset = 'PSP:ON_PROCESS').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_known_interaction = spliced_ptms.dropna(subset = 'PSP:ON_PROT_INTERACT').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        print(f"PhosphoSitePlus regulatory_site information added:\n\t ->{num_ptms_with_known_function} PTMs in dataset found associated with a molecular function \n\t ->{num_ptms_with_known_process} PTMs in dataset found associated with a biological process\n\t ->{num_ptms_with_known_interaction} PTMs in dataset found associated with a protein interaction")
    return spliced_ptms

def add_PSP_kinase_substrate_data(spliced_ptms, file = 'Kinase_Substrate_Dataset.gz', report_success = True):
    """
    Add kinase substrate data from PhosphoSitePlus (Kinase_Substrate_Dataset.gz) to spliced_ptms dataframe from project_ptms_onto_splice_events() function

    Parameters
    ----------
    file: str
        Path to the PhosphoSitePlus Kinase_Substrate_Dataset.gz file. Should be downloaded from PhosphoSitePlus in the zipped format

    Returns
    -------
    spliced_ptms: pandas.DataFrame
        Contains the PTMs identified across the different splice events with an additional column indicating the kinases known to phosphorylate that site (not relevant to non-phosphorylation PTMs)

    """
    #check to make sure provided file exists
    check_file(file, expected_extension='.gz')

    #load data
    ks_dataset = pd.read_csv(file, sep = '\t', header = 2, on_bad_lines='skip',compression = 'gzip', encoding = "cp1252")
    #restrict to human data
    ks_dataset = ks_dataset[ks_dataset['KIN_ORGANISM'] == 'human']
    ks_dataset = ks_dataset[ks_dataset['SUB_ORGANISM'] == 'human']

    ks_dataset = ks_dataset[['GENE', 'SUB_ACC_ID', 'SUB_MOD_RSD']].groupby(['SUB_ACC_ID', 'SUB_MOD_RSD']).agg(';'.join).reset_index()
    ks_dataset.columns = ['UniProtKB Accession', 'Residue', 'PSP:Kinase']

    #separate residue and position
    ks_dataset['PTM Position in Canonical Isoform'] = ks_dataset['Residue'].apply(lambda x: int(x[1:]))
    ks_dataset['Residue'] = ks_dataset['Residue'].apply(lambda x: x[0])

    
    #if splice data already has the annotation columns, remove them
    if 'PSP:Kinase' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['PSP:Kinase'])

    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(ks_dataset, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')
    
    
        #report the number of ptms with kinase substrate information
    if report_success:
        num_ptms_with_KS = spliced_ptms.dropna(subset = 'PSP:Kinase').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        print(f"PhosphoSitePlus kinase-substrate interactions added: {num_ptms_with_KS} phosphorylation sites in dataset found associated with a kinase in PhosphoSitePlus")
    return spliced_ptms

def add_PSP_disease_association(spliced_ptms, file = 'Disease-associated_sites.gz', report_success = True):
    """
    Process disease asociation data from PhosphoSitePlus (Disease-associated_sites.gz), and add to spliced_ptms dataframe from project_ptms_onto_splice_events() function

    Parameters
    ----------
    file: str
        Path to the PhosphoSitePlus Kinase_Substrate_Dataset.gz file. Should be downloaded from PhosphoSitePlus in the zipped format

    Returns
    -------
    spliced_ptms: pandas.DataFrame
        Contains the PTMs identified across the different splice events with an additional column indicating the kinases known to phosphorylate that site (not relevant to non-phosphorylation PTMs)

    """
    #check to make sure provided file exists
    check_file(file, expected_extension='.gz')

    #load data
    disease_associated_sites = pd.read_csv(file, sep = '\t', header = 2, on_bad_lines='skip',compression = 'gzip')
    disease_associated_sites = disease_associated_sites[disease_associated_sites['ORGANISM'] == 'human']

    #removes sites without a specific disease annotation
    disease_associated_sites = disease_associated_sites.dropna(subset = ['DISEASE'])

    #drop extra modification information that is not needed
    #drop extra modification information that is not needed
    disease_associated_sites['Residue'] = disease_associated_sites['MOD_RSD'].apply(lambda x: x.split('-')[0][0])
    disease_associated_sites['PTM Position in Canonical Isoform'] = disease_associated_sites['MOD_RSD'].apply(lambda x: int(x.split('-')[0][1:]))
    #add modification type
    disease_associated_sites['Modification Class'] = disease_associated_sites['MOD_RSD'].apply(lambda x: mod_shorthand_dict[x.split('-')[1]])
    #if phosphorylation, add specific residue
    disease_associated_sites['Modification Class'] = disease_associated_sites.apply(lambda x: x['Modification Class'] + residue_dict[x['Residue'][0]] if x['Modification Class'] == 'Phospho' else x['Modification Class'], axis = 1)
    #change O-GalNac occurring on N to N-glycosylation
    disease_associated_sites['Modification Class'] = disease_associated_sites.apply(lambda x: 'N-Glycosylation' if x['Modification Class'] == 'O-Glycosylation' and x['Residue'][0] == 'N' else x['Modification Class'], axis = 1)


    #combine disease and alteration
    disease_associated_sites['ALTERATION'] = disease_associated_sites.apply(lambda x: x['DISEASE']+'->'+x['ALTERATION'] if x['ALTERATION'] == x['ALTERATION'] else x['DISEASE'], axis = 1)
    #grab only necessary columns and rename
    disease_associated_sites = disease_associated_sites[['ACC_ID', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'ALTERATION']]
    disease_associated_sites.columns = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'PSP:Disease_Association']

    #aggregate multiple disease associations
    disease_associated_sites = disease_associated_sites.groupby(['UniProtKB Accession', 'Residue','PTM Position in Canonical Isoform', 'Modification Class']).agg(';'.join).reset_index()

    #if splice data already has the annotation columns, remove them
    if 'PSP:Disease_Association' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['PSP:Disease_Association'])

    #explode dataframe on modifications
    if spliced_ptms['Modification Class'].str.contains(';').any():
        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)


    #merge with spliced_ptm info
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(disease_associated_sites, how = 'left', on = ['UniProtKB Accession', 'Residue','PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')
    
    #
    #report the number of ptms with motif data
    if report_success:
        num_ptms_with_disease = spliced_ptms.dropna(subset = 'PSP:Disease_Association').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        print(f"PhosphoSitePlus disease associations added: {num_ptms_with_disease} PTM sites in dataset found associated with a disease in PhosphoSitePlus")
    
    
    return spliced_ptms


def add_ELM_interactions(spliced_ptms, file = None, report_success =True):
    """
    Given a spliced ptms dataframe from the project module, add ELM interaction data to the dataframe
    """
    #load data
    if file is None:
        elm_interactions = pd.read_csv('http://elm.eu.org/interactions/as_tsv', sep = '\t', header = 0)
    else:
        check_file(file, expected_extension='.tsv')
        elm_interactions = pd.read_csv(file, sep = '\t', header = 0)

    elm_interactions = elm_interactions[(elm_interactions['taxonomyElm'] == '9606(Homo sapiens)') & (elm_interactions['taxonomyDomain'] == '9606(Homo sapiens)')]

    elm_list = []
    elm_type = []
    elm_interactor = []
    for i, row in spliced_ptms.iterrows():
        #grab ptm location from residue column (gives residue and position (S981), so need to remove residue and convert to int)
        ptm_loc = int(row['PTM Position in Canonical Isoform']) if row['PTM Position in Canonical Isoform'] == row['PTM Position in Canonical Isoform'] and row['PTM Position in Canonical Isoform'] != 'none' else None

        #if data does not have position information, move to the next
        if ptm_loc is None:
            elm_list.append(np.nan)
            elm_type.append(np.nan)
            elm_interactor.append(np.nan)
            continue

        #find if any of the linear motifs match ptm loc
        protein_match = row['UniProtKB Accession'] == elm_interactions['interactorElm']
        region_match = (ptm_loc >= elm_interactions['StartElm'])  & (ptm_loc <=elm_interactions['StopElm'])
        elm_subset_motif = elm_interactions[protein_match & region_match]
        #if any interactions were found, record and continue to the next (assumes a single ptm won't be found as both a SLiM and domain)
        if elm_subset_motif.shape[0] > 0:
            elm_list.append(';'.join(elm_subset_motif['Elm'].values))
            elm_type.append('SLiM')
            elm_interactor.append(';'.join(elm_subset_motif['interactorDomain'].values))
            continue


        #domain
        protein_match = row['UniProtKB Accession'] == elm_interactions['interactorDomain']
        region_match = (ptm_loc >= elm_interactions['StartDomain'])  & (ptm_loc <=elm_interactions['StopDomain'])
        elm_subset_domain = elm_interactions[protein_match & region_match]
        #if any interactions were found, record and continue to the next (assumes a single ptm won't be found as both a SLiM and domain)
        if elm_subset_domain.shape[0] > 0:
            elm_list.append(';'.join(elm_subset_domain['Elm'].values))
            elm_type.append('Domain')
            elm_interactor.append(';'.join(elm_subset_domain['interactorElm'].values))
            continue

        #if no interactions wer found, record as np.nan
        elm_list.append(np.nan)
        elm_type.append(np.nan)
        elm_interactor.append(np.nan)

    spliced_ptms['ELM:Interactions'] = elm_interactor
    spliced_ptms['ELM:Location of PTM for Interaction'] = elm_type
    spliced_ptms['ELM:Motifs Associated with Interactions'] = elm_list
    
    #report the number of ptms with motif data
    if report_success:
        num_ptms_with_ELM_instance = spliced_ptms.dropna(subset = 'ELM:Interactions').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).size().shape[0]
        print(f"ELM interaction instances added: {num_ptms_with_ELM_instance} PTMs in dataset found associated with at least one known ELM instance")
    return spliced_ptms


def add_ELM_matched_motifs(spliced_ptms, flank_size = 7, file = None, report_success = True):
    if file is None:
        elm_classes = pd.read_csv('http://elm.eu.org/elms/elms_index.tsv', sep = '\t', header = 5)
    else:
        check_file(file, expected_extension='.tsv')
        elm_classes = pd.read_csv(file, sep = '\t', header = 5)

    ptm_coordinates = pose_config.ptm_coordinates.copy()
    #create corresponding label for ptm_coordinate data
    ptm_coordinates['PTM Label'] = ptm_coordinates['UniProtKB Accession'] + '_' + ptm_coordinates['Residue'] + ptm_coordinates['PTM Position in Canonical Isoform'].apply(lambda x: int(float(x)) if x == x else np.nan).astype(str)
    
    match_list = []
    for i, row in spliced_ptms.iterrows():
        matches = []
        #grab ptm information
        #grab flanking sequence for the ptm
        loc = int(row["PTM Position in Canonical Isoform"]) if row['PTM Position in Canonical Isoform'] == row['PTM Position in Canonical Isoform'] else np.nan
        ptm = row['UniProtKB Accession'] + '_' + row['Residue'] + str(loc)

        
        if ptm in ptm_coordinates['PTM Label'].values:
            ptm_flanking_seq = ptm_coordinates.loc[ptm_coordinates['PTM Label'] == ptm, 'Expected Flanking Sequence'].values[0]
            #make sure flanking sequence is present
            if isinstance(ptm_flanking_seq, str):

                #default flanking sequence is 10, if requested flanking sequence is different, then adjust
                if flank_size > 10:
                    raise ValueError('Flanking size must be equal to or less than 10')
                elif flank_size < 10:
                    ptm_flanking_seq = ptm_flanking_seq[10-flank_size:10+flank_size]

                for j, elm_row in elm_classes.iterrows():
                    reg_ex = elm_row['Regex']
                    if re.search(reg_ex, ptm_flanking_seq) is not None:
                        matches.append(elm_row['ELMIdentifier'])

                match_list.append(';'.join(matches))
            else:
                match_list.append(np.nan)
        else:
            #print(f'PTM {ptm} not found in PTM info file')
            match_list.append(np.nan)
    
    spliced_ptms['ELM:Motif Matches'] = match_list

    #report the number of ptms with motif data
    if report_success:
        num_ptms_with_matched_motif = spliced_ptms.dropna(subset = 'ELM:Motif Matches').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).size().shape[0]
        print(f"ELM Class motif matches found: {num_ptms_with_matched_motif} PTMs in dataset found with at least one matched motif")
    return spliced_ptms

def add_PTMInt_data(spliced_ptms, file = None, report_success = True):
    """
    Given spliced_ptms data from project module, add PTMInt interaction data, which will include the protein that is being interacted with, whether it enchances or inhibits binding, and the localization of the interaction. This will be added as a new column labeled PTMInt:Interactions and each entry will be formatted like 'Protein->Effect|Localization'. If multiple interactions, they will be separated by a semicolon
    """
    #load file
    if file is None:
        PTMint = pd.read_csv('https://ptmint.sjtu.edu.cn/data/PTM%20experimental%20evidence.csv')
    else:
        check_file(file, expected_extension='.csv')
        PTMint = pd.read_csv(file)

    PTMint = PTMint.rename(columns={'Uniprot':'UniProtKB Accession', 'AA':'Residue', 'Site':'PTM Position in Canonical Isoform'})
    #PTMint['Site'] = PTMint['AA'] + PTMint['Site'].astype(str)
    PTMint['PTMInt:Interaction'] = PTMint['Int_gene']+'->'+PTMint['Effect']
    PTMint = PTMint[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'PTMInt:Interaction']]
    #PTMint['PTM Position in Canonical Isoform'] = PTMint['PTM Position in Canonical Isoform'].astype(str)

    #aggregate PTMint data on the same PTMs
    PTMint = PTMint.groupby(['UniProtKB Accession','Residue','PTM Position in Canonical Isoform'], as_index = False).agg(';'.join)

    #if splice data already has the annotation columns, remove them
    if 'PTMInt:Interaction' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['PTMInt:Interaction'])

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(PTMint[['UniProtKB Accession','Residue','PTM Position in Canonical Isoform', 'PTMInt:Interaction']], on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'], how = 'left')
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')

    #report the number of PTMs identified
    if report_success:
        num_ptms_with_PTMInt_data = spliced_ptms.dropna(subset = 'PTMInt:Interaction').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
        print(f"PTMInt data added: {num_ptms_with_PTMInt_data} PTMs in dataset found with PTMInt interaction information")

    return spliced_ptms
    #delete source PTMint data
    #os.remove(pdir + './Data/PTM_experimental_evidence.csv')

#def add_PTMcode_intraprotein(spliced_ptms, fname = None, report_success = True):
#    #load ptmcode info
#    if fname is None:
#        ptmcode = pd.read_csv('https://ptmcode.embl.de/data/PTMcode2_associations_within_proteins.txt.gz', sep = '\t', header = 2, compression='gzip')
#    else:
#        check_file(fname, expected_extension = '.gz')
#        ptmcode = pd.read_csv(fname, sep = '\t', header = 2, compression = 'gzip')
#    
#    #grab humn data
#    ptmcode = ptmcode[ptmcode['Species'] == 'Homo sapiens']
#
#    #add gene name to data
#    translator = pd.DataFrame(pose_config.uniprot_to_genename, index = ['Gene']).T
#    translator['Gene'] = translator['Gene'].apply(lambda x: x.split(' '))
#    translator = translator.explode('Gene')
#    translator = translator.reset_index()
#    translator.columns = ['UniProtKB/Swiss-Prot ID', 'Gene name']
#
#    #add uniprot ID information
#    ptmcode = ptmcode.merge(translator.dropna().drop_duplicates(), left_on = '## Protein', right_on = 'Gene name', how = 'left')
#
#    #convert modification names to match annotation data
#    convert_dict = {'Adp ribosylation': 'ADP Ribosylation', 'Glutamine deamidation':'Deamidation'}
#    new_mod_names = []
#    failed_mod = []
#    mod_list = ptmcode['PTM1'].unique()
#    for mod in mod_list:
#        mod = mod.capitalize()
#        if 'glycosylation' in mod: #if glycosylation, group into one gorup
#            new_mod_names.append('Glycosylation')
#        elif mod in pose_config.modification_conversion['Modification Class'].values: #if already in modification class data, keep
#            new_mod_names.append(mod)
#        elif mod in convert_dict.keys():
#            new_mod_names.append(convert_dict[mod])
#        else:
#            try:
#                new_mod = pose_config.modification_conversion[pose_config.modification_conversion['Modification'] == mod].values[0][0]
#                new_mod_names.append(new_mod)
#            except:
#                failed_mod.append(mod)
#                new_mod_names.append(mod)
#    conversion_df = pd.DataFrame({'PTM1':mod_list, 'Modification Class':new_mod_names})
#
#    #add new modification labels to data
#    ptmcode = ptmcode.merge(conversion_df, on = 'PTM1', how = 'left')
#    
#    #groupby by PTM1 and rename to match column names in annotation data
#    ptmcode = ptmcode[['UniProtKB/Swiss-Prot ID', 'Modification Class', 'Residue1', 'Residue2']].dropna(subset = 'UniProtKB/Swiss-Prot ID')
#    ptmcode = ptmcode.groupby(['UniProtKB/Swiss-Prot ID', 'Modification Class', 'Residue1'])['Residue2'].agg(';'.join).reset_index()
#    ptmcode = ptmcode.rename(columns = {'UniProtKB/Swiss-Prot ID':'UniProtKB Accession', 'Residue1':'Residue', 'Residue2':'PTMcode:Intraprotein_Interactions'})
#    
#    #separate residue information into separate columns, one for amino acid and one for position
#    ptmcode['PTM Position in Canonical Isoform'] = ptmcode['Residue'].apply(lambda x: int(x[1:]))
#    ptmcode['Residue'] = ptmcode['Residue'].apply(lambda x: x[0])
#
#        #if splice data already has the annotation columns, remove them
#    if 'PTMcode:Intraprotein_Interactions' in spliced_ptms.columns:
#        spliced_ptms = spliced_ptms.drop(columns = ['PTMcode:Intraprotein_Interactions'])
#
#    #explode dataframe on modifications
#    if spliced_ptms['Modification Class'].str.contains(';').any():
#        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
#        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)
#
#    #add to splice data
#    original_data_size = spliced_ptms.shape[0]
#    spliced_ptms = spliced_ptms.merge(ptmcode, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
#    if spliced_ptms.shape[0] != original_data_size:
#        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
#    
#    #report the number of PTMs identified
#    if report_success:
#        num_ptms_with_PTMcode_data = spliced_ptms.dropna(subset = 'PTMcode:Intraprotein_Interactions').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
#        print(f"PTMcode intraprotein interactions added: {num_ptms_with_PTMcode_data} PTMs in dataset found with PTMcode intraprotein interaction information")
#
#    return spliced_ptms

def extract_ids_PTMcode(df, col = '## Protein1'):

    #add gene name to data
    name_to_uniprot = pd.DataFrame(pose_config.uniprot_to_genename, index = ['Gene']).T
    name_to_uniprot['Gene'] = name_to_uniprot['Gene'].apply(lambda x: x.split(' ') if x == x else np.nan)
    name_to_uniprot = name_to_uniprot.explode('Gene')
    name_to_uniprot = name_to_uniprot.reset_index()
    name_to_uniprot.columns = ['UniProtKB/Swiss-Prot ID', 'Gene name']
    name_to_uniprot = name_to_uniprot.drop_duplicates(subset = 'Gene name', keep = False)

    #protein name is provided as either ensemble gene id or gene name check for both
    df = df.merge(pose_config.translator[['Gene stable ID']].reset_index().dropna().drop_duplicates(), left_on = col, right_on = 'Gene stable ID', how = 'left')
    df = df.rename(columns = {'index': 'From_ID'})
    df = df.merge(name_to_uniprot, left_on = col, right_on = 'Gene name', how = 'left')
    df = df.rename(columns = {'UniProtKB/Swiss-Prot ID': 'From_Name'})

    #grab unique id from 'From_ID' and 'From_Name' column, if available
    uniprot_ids = df['From_Name'].combine_first(df['From_ID'])
    return uniprot_ids.values

def add_PTMcode_interprotein(spliced_ptms, fname = None, report_success = True):
    if fname is None:
        ptmcode = pd.read_csv('https://ptmcode.embl.de/data/PTMcode2_associations_between_proteins.txt.gz', sep = '\t', header = 2, compression = 'gzip')
    else:
        check_file(fname, expected_extension = '.gz')
        ptmcode = pd.read_csv(fname, sep = '\t', header = 2, compression='gzip')

    #grab human interactions
    ptmcode = ptmcode[ptmcode['Species'] == 'Homo sapiens']
    #ignore intraprotein interactions
    ptmcode = ptmcode[ptmcode['## Protein1'] != ptmcode['Protein2']]

    #get uniprot id for primary protein and interacting protein
    ptmcode['UniProtKB Accession'] = extract_ids_PTMcode(ptmcode, '## Protein1')
    ptmcode['Interacting Protein'] = extract_ids_PTMcode(ptmcode, 'Protein2')

    ptmcode = ptmcode.dropna(subset = ['UniProtKB Accession', 'Interacting Protein'])
    #remove duplicate proteins (some entries have different ids but are actually the same protein)
    ptmcode = ptmcode[ptmcode['UniProtKB Accession'] != ptmcode['Interacting Protein']]

    #aggregate interactions
    ptmcode['Interacting Residue'] = ptmcode['Interacting Protein'] + '_' + ptmcode['Residue2']


    #convert modification names
    convert_dict = {'Adp ribosylation': 'ADP Ribosylation', 'Glutamine deamidation':'Deamidation'}
    new_mod_names = []
    failed_mod = []
    mod_list = ptmcode['PTM1'].unique()
    for mod in mod_list:
        mod = mod.capitalize()
        if 'glycosylation' in mod:
            new_mod_names.append('Glycosylation')
        elif mod in pose_config.modification_conversion['Modification Class'].values:
            new_mod_names.append(mod)
        elif mod in convert_dict.keys():
            new_mod_names.append(convert_dict[mod])
        else:
            try:
                new_mod = pose_config.modification_conversion[pose_config.modification_conversion['Modification'] == mod].values[0][0]
                new_mod_names.append(new_mod)
            except:
                failed_mod.append(mod)
                new_mod_names.append(mod)
    conversion_df = pd.DataFrame({'PTM1':mod_list, 'Modification Class':new_mod_names})

    ptmcode = ptmcode.merge(conversion_df, on = 'PTM1', how = 'left')


    ptmcode = ptmcode.rename(columns = {'Residue1':'Residue'})
    ptmcode = ptmcode.groupby(['UniProtKB Accession', 'Residue', 'Modification Class'])['Interacting Residue'].agg(';'.join).reset_index()
    ptmcode = ptmcode.rename(columns = {'UniProtKB/Swiss-Prot ID':'UniProtKB Accession', 'Residue1':'Residue', 'Interacting Residue':'PTMcode:Interprotein_Interactions'})

    #separate residue information into separate columns, one for amino acid and one for position
    ptmcode['PTM Position in Canonical Isoform'] = ptmcode['Residue'].apply(lambda x: float(x[1:]))
    ptmcode['Residue'] = ptmcode['Residue'].apply(lambda x: x[0])

            #if splice data already has the annotation columns, remove them
    if 'PTMcode:Interprotein_Interactions' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['PTMcode:Interprotein_Interactions'])

        #explode dataframe on modifications
    if spliced_ptms['Modification Class'].str.contains(';').any():
        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(ptmcode, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    if report_success:
        num_ptms_with_PTMcode_data = spliced_ptms.dropna(subset = 'PTMcode:Interprotein_Interactions').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
        print(f"PTMcode interprotein interactions added: {num_ptms_with_PTMcode_data} PTMs in dataset found with PTMcode interprotein interaction information")

    return spliced_ptms

def extract_positions_from_DEPOD(x):
    """
    Given string object consisting of multiple modifications in the form of 'Residue-Position' separated by ', ', extract the residue and position. Ignore any excess details in the string.
    """
    x = x.split('[')[0].split(', ')
    #for each residue in list, find location of 'Ser', 'Thr' and 'Tyr' in the string (should either have '-' or a number immediately after it)
    new_x = []
    for item in x:
        #determine type of modification
        if 'Ser' in item:
            loc = [match.start() for match in re.finditer('Ser', item)]
            res = 'S'
        elif 'Thr' in item:
            loc = [match.start() for match in re.finditer('Thr', item)]
            res = 'T'
        elif 'Tyr' in item:
            loc = [match.start() for match in re.finditer('Tyr', item)]
            res = 'Y'
        elif 'His' in item:
            loc = [match.start() for match in re.finditer('His', item)]
            res = 'H'
        else:
            loc = -1

        #check if multiple locations were found, if so grab last entry
        if loc == -1:
            item = np.nan
            make_string = False
        elif len(loc) > 1:
            make_string = True
            loc = loc[-1]
        else:
            loc = loc[0]
            make_string = True
        
        #find integer
        if make_string:
            if '-' in item[loc:]:
                item = item.split('-')
                item = res + item[1].strip()
            else:
                item = item[loc+3:]
                item = res + item

        new_x.append(item)
    
    return new_x

def add_DEPOD_phosphatase_data(spliced_ptms, report_success = True):

    #download data
    depod1 = pd.read_excel('https://depod.bioss.uni-freiburg.de/download/PPase_protSubtrates_201903.xls', sheet_name='PSprots')
    depod2 = pd.read_excel('https://depod.bioss.uni-freiburg.de/download/PPase_protSubtrates_newPairs_201903.xls', sheet_name = 'newPSprots')
    depod = pd.concat([depod1, depod2])

    #remove any rows with missing sit information
    depod = depod.dropna(subset = 'Dephosphosites')

    #remove excess annotations that make parsing difficult
    depod['Dephosphosites'] = depod['Dephosphosites'].apply(lambda x: x.split('[')[0])
    depod['Dephosphosites'] = depod['Dephosphosites'].apply(lambda x: x.split('(')[0])
    depod['Dephosphosites'] = depod['Dephosphosites'].apply(lambda x: x.split(';')[0])
    depod['Dephosphosites'] = depod['Dephosphosites'].apply(lambda x: x.split('in')[0])
    depod['Dephosphosites'] = depod['Dephosphosites'].str.replace('in ref.', '')

    #separate individual sites
    depod['Dephosphosites'] = depod['Dephosphosites'].str.split(',')
    depod = depod.explode('Dephosphosites')
    depod = depod[(~depod['Dephosphosites'].str.contains('Isoform')) & (~depod['Dephosphosites'].str.contains('isoform'))]

    #process dephosphosite strings to extract residue and position and explode so that each phosphosite is its own row
    depod['Dephosphosites'] = depod['Dephosphosites'].apply(extract_positions_from_DEPOD)
    depod = depod.explode('Dephosphosites')

    #separate multiple substrate accessions into their own rows (many of these link back to the same ID, but will keep just in case)
    depod['Substrate accession numbers'] = depod['Substrate accession numbers'].str.split(' ')
    depod = depod.explode('Substrate accession numbers')
    depod = depod.dropna(subset = ['Substrate accession numbers'])

    #extract only needed information and add phosphorylation as modification type
    #extract only needed information and add phosphorylation as modification type
    depod['Residue'] = depod['Dephosphosites'].apply(lambda x: x[0] if x == x else np.nan)
    depod['PTM Position in Canonical Isoform'] = depod['Dephosphosites'].apply(lambda x: int(x[1:]) if x == x else np.nan)
    depod = depod.rename({'Substrate accession numbers': 'UniProtKB Accession', 'Phosphatase entry names':'DEPOD:Phosphatase'}, axis = 1)
    depod = depod[['DEPOD:Phosphatase', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']]
    depod['Modification Class'] = 'Phosphorylation'

    #combine on the same PTM
    depod = depod.drop_duplicates()
    depod = depod.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'], as_index = False)['DEPOD:Phosphatase'].agg(';'.join)

        #if splice data already has the annotation columns, remove them
    if 'DEPOD:Phosphatase' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['DEPOD:Phosphatase'])

        #explode dataframe on modifications
    if spliced_ptms['Modification Class'].str.contains(';').any():
        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(depod, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    if report_success:
        num_ptms_with_PTMcode_data = spliced_ptms.dropna(subset = 'DEPOD:Phosphatase').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
        print(f"DEPOD Phosphatase substrates added: {num_ptms_with_PTMcode_data} PTMs in dataset found with Phosphatase substrate information")

    return spliced_ptms

def add_RegPhos_data(spliced_ptms, file = None, report_success = True):
    if file is None:
        regphos = pd.read_csv('http://140.138.144.141/~RegPhos/download/RegPhos_Phos_human.txt', sep = '\t', dtype = {'position':int, 'description':str,'catalytic kinase':str, 'reference':'str'})
    else:
        check_file(file, expected_extension = '.txt')
        regphos = pd.read_csv(file, sep = '\t')

    regphos = regphos.dropna(subset = 'catalytic kinase')
    #regphos['Residue'] = regphos['code'] + regphos['position'].astype(str)
    regphos = regphos.rename(columns = {'code': 'Residue', 'position':'PTM Position in Canonical Isoform', 'AC': 'UniProtKB Accession', 'catalytic kinase': 'RegPhos:Kinase'})
    regphos['Modification Class'] = 'Phosphorylation'
    regphos = regphos[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'RegPhos:Kinase']].dropna()
    regphos = regphos.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).agg(';'.join).reset_index()

    #if splice data already has the annotation columns, remove them
    if 'RegPhos:Kinase' in spliced_ptms.columns:
        spliced_ptms = spliced_ptms.drop(columns = ['RegPhos:Kinase'])

    #explode dataframe on modifications
    if spliced_ptms['Modification Class'].str.contains(';').any():
        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(regphos, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    if report_success:
        num_ptms_with_regphos_data = spliced_ptms.dropna(subset = 'RegPhos:Kinase').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).size().shape[0]
        print(f"RegPhos kinase-substrate data added: {num_ptms_with_regphos_data} PTMs in dataset found with kinase-substrate information")

    return spliced_ptms


def add_PTMsigDB_data(spliced_ptms, file = None, report_success = True):
    #if file is None:
    #    ptmsigdb = pd.read_excel('https://proteomics.broadapps.org/ptmsigdb/_w_8b062d9e/appff37efd164a676afcc8e6e42e6058e01/session/a2b28c4ed29deadd6779fdd26aec33c1/download/download.xlsx?w=8b062d9e', sheet_name = 'human')
    #else:
    check_file(file, expected_extension = '.xlsx')
    ptmsigdb = pd.read_excel(file, sheet_name = 'human')


    ptmsigdb['UniProtKB Accession'] = ptmsigdb['site.uniprot'].str.split(';').str[0]
    ptmsigdb['Residue'] = ptmsigdb['site.uniprot'].str.split(';').str[1].str[0]
    ptmsigdb['PTM Position in Canonical Isoform'] = ptmsigdb['site.uniprot'].apply(lambda x: int(x.split(';')[1].split('-')[0][1:]))

    #filter out excess information in some of the site.ptm column, then convert to modification class details
    ptmsigdb['site.ptm'] = ptmsigdb['site.ptm'].apply(lambda x: x.split(';')[1].split('-')[1] if ';' in x else x)
    ptmsigdb['Modification Class'] = ptmsigdb['site.ptm'].map(mod_shorthand_dict)

    #combine signature and direction for annotation column
    ptmsigdb['Signature'] = ptmsigdb['signature'] +'->'+ ptmsigdb['site.direction']

    #drop unneeded columns
    ptmsigdb = ptmsigdb[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'Signature', 'category']]
    ptmsigdb['Signature'] = ptmsigdb.apply(lambda x: x['Signature'].replace(x['category'] + '_', ''), axis = 1)
    ptmsigdb['category'] = 'PTMsigDB:' + ptmsigdb['category'] 
    ptmsigdb = ptmsigdb.drop_duplicates()

        #convert to pivot table with each category being a separate column
    ptmsigdb = ptmsigdb.pivot_table(index = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'], columns = 'category', values = 'Signature', aggfunc=';'.join).reset_index()

    #remove psp data if it is already in spliced ptms
    if 'PSP:Kinase' in spliced_ptms.columns:
        ptmsigdb = ptmsigdb.drop(columns = 'PTMsigDB:KINASE-PSP')

    if 'PSP:Disease_Association' in spliced_ptms.columns:
        ptmsigdb = ptmsigdb.drop(columns = 'PTMsigDB:DISEASE-PSP')


    #if splice data already has the annotation columns, remove them
    if 'PTMsigDB:PATH-BI' in spliced_ptms.columns:
        cols_in_data = [col for col in spliced_ptms.columns if 'PTMsigDB' in col]
        spliced_ptms = spliced_ptms.drop(columns = cols_in_data)


    #explode dataframe on modifications
    if spliced_ptms['Modification Class'].str.contains(';').any():
        spliced_ptms['Modification Class'] = spliced_ptms['Modification Class'].str.split(';')
        spliced_ptms = spliced_ptms.explode('Modification Class').reset_index(drop = True)

    #merge with spliced_ptm info
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(ptmsigdb, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')


    #report the number of ptms with motif data
    if report_success:
        num_ptms_with_ikip = spliced_ptms.dropna(subset = 'PTMsigDB:KINASE-iKiP').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_path_bi = spliced_ptms.dropna(subset = 'PTMsigDB:PATH-BI').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_path_np= spliced_ptms.dropna(subset = 'PTMsigDB:PATH-NP').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_path_wp = spliced_ptms.dropna(subset = 'PTMsigDB:PATH-WP').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_dia_pert = spliced_ptms.dropna(subset = 'PTMsigDB:PERT-P100-DIA').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_dia2_pert = spliced_ptms.dropna(subset = 'PTMsigDB:PERT-P100-DIA2').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_prm_pert = spliced_ptms.dropna(subset = 'PTMsigDB:PERT-P100-PRM').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        num_ptms_with_psp_pert = spliced_ptms.dropna(subset = 'PTMsigDB:PERT-PSP').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).size().shape[0]
        print(f"PTMsigDB added:\n\t ->{num_ptms_with_ikip} PTMs associated with kinases in iKiP\n\t ->{num_ptms_with_path_wp} PTMs associated with molecular pathway signatures from WikiPathways\n\t ->{num_ptms_with_path_np} PTMs associated with molecular pathway signatures from NetPath\n\t ->{num_ptms_with_psp_pert} PTMs with PhosphoSitePlus perturbations\n\t ->{num_ptms_with_dia_pert} with perturbations in LINCS P1000 DIA dataset \n\t ->{num_ptms_with_dia2_pert} with perturbations in LINCS P1000 DIA2 dataset\n\t ->{num_ptms_with_prm_pert} with perturbations in LINCS P1000 PRM dataset")
    return spliced_ptms



######### Functions for combining annotations from multiple sources ########

def convert_PSP_label_to_UniProt(label):
    """
    Given a label for an interacting protein from PhosphoSitePlus, convert to UniProtKB accession. Required as PhosphoSitePlus interactions are recorded in various ways that aren't necessarily consistent with other databases (i.e. not always gene name)

    Parameters
    ----------
    label: str
        Label for interacting protein from PhosphoSitePlus
    """
    if not hasattr(pose_config, 'genename_to_uniprot'):
        #using uniprot to gene name dict, construct dict to go the other direction (gene name to uniprot id)
        pose_config.genename_to_uniprot = pose_config.flip_uniprot_dict(pose_config.uniprot_to_genename)


    #remove isoform label if present
    if label in pose_config.genename_to_uniprot: #if PSP name is gene name found in uniprot
        return pose_config.genename_to_uniprot[label]
    elif label.upper() in pose_config.genename_to_uniprot:
        return pose_config.genename_to_uniprot[label.upper()]
    elif label.split(' ')[0].upper() in pose_config.genename_to_uniprot:
        return pose_config.genename_to_uniprot[label.split(' ')[0].upper()]
    elif label.replace('-', '').upper() in pose_config.genename_to_uniprot:
        return pose_config.genename_to_uniprot[label.replace('-', '').upper()]
    elif label in pose_config.psp_name_dict: # if PSP name is not gene name, but is in conversion dictionary
        return pose_config.psp_name_dict[label]
    else: #otherwise note that gene was missed
        return np.nan
        #missed_genes.append(gene)

def extract_interaction_details(interaction, column = "PSP:ON_PROT_INTERACT"):

    interaction_types = {'PTMcode:Interprotein_Interactions':'INDUCES', 'PSP:Kinase':'REGULATES', 'DEPOD:Phosphatase':'REGULATES', 'RegPhos:Kinase':'REGULATES', 'Combined:Kinase':'REGULATES', 'ELM:Interactions':'UNCLEAR'}
    if column == 'PSP:ON_PROT_INTERACT':
        type = interaction.split('(')[1].split(')')[0]
        protein = interaction.split('(')[0].strip(' ')
    elif column == 'PTMInt:Interaction':
        ptmint_type_conversion = {'Inhibit':'DISRUPTS', 'Enhance':"INDUCES"}
        type = ptmint_type_conversion[interaction.split('->')[1]]
        protein = interaction.split('->')[0]
    elif column == 'PTMcode:Interprotein_Interactions':
        type = 'INDUCES'
        protein = interaction.split('_')[0]
    else:
        type = interaction_types[column]
        protein = interaction

    return type, protein

def unify_interaction_data(spliced_ptms, interaction_col, name_dict = {}):
    """
    Given spliced ptm data and a column containing interaction data, extract the interacting protein, type of interaction, and convert to UniProtKB accession. This will be added as a new column labeled 'Interacting ID'

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe containing PTM data
    interaction_col: str
        column containing interaction information from a specific database
    name_dict: dict
        dictionary to convert names within given database to UniProt IDs. For cases when name is not necessarily one of the gene names listed in UniProt

    Returns
    -------
    interact: pd.DataFrame
        Contains PTMs and their interacting proteins, the type of influence the PTM has on the interaction (DISRUPTS, INDUCES, or REGULATES)
    """
    if not hasattr(pose_config, 'genename_to_uniprot'):
        #using uniprot to gene name dict, construct dict to go the other direction (gene name to uniprot id)
        pose_config.genename_to_uniprot = pose_config.flip_uniprot_dict(pose_config.uniprot_to_genename)

    #extract PSP data from annotated PTMs, separate cases in which single PTM has multipe interactions
    data_cols = [col for col in spliced_ptms.columns if col in ['Significance', 'dPSI']]
    interact = spliced_ptms.dropna(subset = interaction_col)[['Gene', 'UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class',interaction_col] + data_cols]
    if interact.empty:
        print(f"No PTMs associated with {interaction_col}")
        return interact
    
    interact[interaction_col] = interact[interaction_col].apply(lambda x: x.split(';'))
    interact = interact.explode(interaction_col)

    #extract protein and type of interaction (currently for phosphosite plus)
    type = []
    protein = []
    for i, row in interact.iterrows():
        processed = extract_interaction_details(row[interaction_col], interaction_col)
        type.append(processed[0])
        protein.append(processed[1])
    interact['Type']  = type
    interact['Interacting Protein'] = protein
        

    #convert interacting protein to uniprot id for databases that are not reported in uniprot ids
    if interaction_col not in ['PTMcode:Interprotein_Interactions', 'ELM:Interactions']:
        interacting_id = []
        missed_genes = []
        for gene in interact['Interacting Protein']:
            #remove isoform label if present
            if gene in pose_config.genename_to_uniprot: #if PSP name is gene name found in uniprot
                interacting_id.append(pose_config.genename_to_uniprot[gene])
            elif gene.upper() in pose_config.genename_to_uniprot:
                interacting_id.append(pose_config.genename_to_uniprot[gene.upper()])
            elif gene.split(' ')[0].upper() in pose_config.genename_to_uniprot:
                interacting_id.append(pose_config.genename_to_uniprot[gene.split(' ')[0].upper()])
            elif gene.replace('-', '').upper() in pose_config.genename_to_uniprot:
                interacting_id.append(pose_config.genename_to_uniprot[gene.replace('-', '').upper()])
            elif gene in name_dict: # if PSP name is not gene name, but is in conversion dictionary
                interacting_id.append(name_dict[gene])
            else: #otherwise note that gene was missed
                interacting_id.append(np.nan)
                missed_genes.append(gene)

        #save information
        interact['Interacting ID'] = interacting_id
        interact = interact.dropna(subset = 'Interacting ID')


        #check if there multiple in one row
        if interact['Interacting ID'].str.contains(';').any():
            interact['Interacting ID'] = interact['Interacting ID'].apply(lambda x: x.split(';'))
            interact = interact.explode('Interacting ID')
    else:
        interact['Interacting ID'] = interact['Interacting Protein']
    

    interact['Interacting ID'] = interact['Interacting ID'].apply(lambda x: x.split(' ')[0] if x == x else np.nan)
    interact = interact.explode('Interacting ID')
    interact = interact.dropna(subset = 'Interacting ID')
    interact = interact.drop(columns = interaction_col).drop_duplicates()

    return interact

def add_annotation(spliced_ptms, database = 'PhosphoSitePlus', annotation_type = 'Function', file = None, check_existing = False):
    """
    Given a desired database and annotation type, add the corresponding annotation data to the spliced ptm dataframe

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe containing PTM data
    database: str
        Database to extract annotation data from. Options include 'PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD'
    annotation_type: str
        Type of annotation to extract. Options include 'Function', 'Process', 'Interactions', 'Disease', 'Kinase', 'Phosphatase', but depend on the specific database (run analyze.get_annotation_categories())
    file: str
        File path to annotation data. If None, will download from online source, except for PhosphoSitePlus (due to licensing restrictions)
    """
    if check_existing:
        annot_col = annotation_col_dict[database][annotation_type]
        if annot_col in spliced_ptms.columns:
            print(f"Annotation data for {database} {annotation_type} already present in provided dataframe, skipping. If you would like to update annotation data, set check_existing = False")
            return spliced_ptms

    if database == "PhosphoSitePlus":
        if annotation_type in ['Function', 'Process', 'Interactions']:
            check_file(file, expected_extension='.gz')
            spliced_ptms = add_PSP_regulatory_site_data(spliced_ptms, file = file)
        elif annotation_type == 'Disease':
            check_file(file, expected_extension='.gz')
            spliced_ptms = add_PSP_disease_association(spliced_ptms, file = file)
        elif annotation_type == 'Kinase':
            check_file(file, expected_extension='.gz')
            spliced_ptms = add_PSP_kinase_substrate_data(spliced_ptms, file = file)
        else:
            raise ValueError(f"Annotation type {annotation_type} not recognized for PhosphoSitePlus")
    elif database == 'PTMcode':
        #if annotation_type == 'Intraprotein':
        #    if file is not None:
        #        check_file(file, expected_extension='.gz')
        #        spliced_ptms = add_PTMcode_intraprotein(spliced_ptms, file = file)
        #    else:
        #        spliced_ptms = add_PTMcode_intraprotein(spliced_ptms)
        if annotation_type == 'Interactions':
            if file is not None:
                check_file(file, expected_extension='.gz')
                spliced_ptms = add_PTMcode_interprotein(spliced_ptms, file = file)
            else:
                spliced_ptms = add_PTMcode_interprotein(spliced_ptms)
        else:
            raise ValueError(f"Annotation type {annotation_type} not recognized for PTMcode")
    elif database == 'PTMInt':
        if annotation_type == 'Interactions':
            if file is not None:
                check_file(file, expected_extension='.csv')
                spliced_ptms = add_PTMInt_data(spliced_ptms, file = file)
            else:
                spliced_ptms = add_PTMInt_data(spliced_ptms)
        else:
            raise ValueError(f"Annotation type {annotation_type} not recognized for PTMInt")
    elif database == 'RegPhos':
        if annotation_type == 'Kinase':
            if file is not None:
                check_file(file, expected_extension='.txt')
                spliced_ptms = add_RegPhos_data(spliced_ptms, file = file)
            else:
                spliced_ptms = add_RegPhos_data(spliced_ptms)
        else:
            raise ValueError(f"Annotation type {annotation_type} not recognized for RegPhos")
    elif database == 'DEPOD':
        if annotation_type == 'Phosphatase':
            spliced_ptms = add_DEPOD_phosphatase_data(spliced_ptms, file = file)
        else:
            raise ValueError(f"Annotation type {annotation_type} not recognized for RegPhos")
    elif database == 'Combined':
        if annotation_type == 'Kinase':
            if 'PSP:Kinase' not in spliced_ptms.columns:
                raise ValueError("PhosphoSitePlus kinase data not found in spliced PTM dataframe, please annotate with this first")
            if 'RegPhos:Kinase' not in spliced_ptms.columns:
                spliced_ptms = add_RegPhos_data(spliced_ptms)
            spliced_ptms = combine_KS_data(spliced_ptms)
        elif annotation_type == 'Interactions':
            spliced_ptms = combine_interaction_data(spliced_ptms)
    else:
        raise ValueError(f"Database {database} not recognized")
    
    return spliced_ptms
    

def combine_interaction_data(spliced_ptms, interaction_databases = ['PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD', 'ELM'], include_enzyme_interactions = True):
    """
    Given annotated spliced ptm data, extract interaction data from various databases and combine into a single dataframe. This will include the interacting protein, the type of interaction, and the source of the interaction data

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Dataframe containing PTM data and associated interaction annotations from various databases
    interaction_databases: list
        List of databases to extract interaction data from. Options include 'PhosphoSitePlus', 'PTMcode', 'PTMInt', 'RegPhos', 'DEPOD'. These should already have annotation columns in the spliced_ptms dataframe, otherwise they will be ignored. For kinase-substrate interactions, if combined column is present, will use that instead of individual databases
    include_enzyme_interactions: bool
        If True, will include kinase-substrate and phosphatase interactions in the output dataframe

    Returns
    -------
    interact_data: list
        List of dataframes containing PTMs and their interacting proteins, the type of influence the PTM has on the interaction (DISRUPTS, INDUCES, or REGULATES), and the source of the interaction data

    """
    interact_data = []
    combined_added = False
    for database in interaction_databases:
        if database == 'PhosphoSitePlus' and 'PSP:ON_PROT_INTERACT' in spliced_ptms.columns:
            if not spliced_ptms['PSP:ON_PROT_INTERACT'].isna().all():
                print('PhosphoSitePlus regulatory site data found and added')
                interact = unify_interaction_data(spliced_ptms, 'PSP:ON_PROT_INTERACT', pose_config.psp_name_dict)
                interact['Source'] = database
                interact_data.append(interact)

  
        if database == 'PTMcode' and 'PTMcode:Interprotein_Interactions' in spliced_ptms.columns:
            if not spliced_ptms['PTMcode:Interprotein_Interactions'].isna().all():
                print('PTMcode data found and added')
                interact = unify_interaction_data(spliced_ptms, 'PTMcode:Interprotein_Interactions')
                interact['Source'] = database
                interact_data.append(interact)
        if database == 'PTMInt' and 'PTMInt:Interaction' in spliced_ptms.columns:
            if not spliced_ptms['PTMInt:Interaction'].isna().all():
                print('PTMInt data found and added')
                interact = unify_interaction_data(spliced_ptms, 'PTMInt:Interaction')
                interact['Source'] = database
                interact_data.append(interact)
        if database == 'ELM' and 'ELM:Interactions' in spliced_ptms.columns:
            if not spliced_ptms['ELM:Interactions'].isna().all():
                print('ELM data found and added')
                interact = unify_interaction_data(spliced_ptms, 'ELM:Interactions')
                interact['Source'] = database
                interact_data.append(interact)
        
        if include_enzyme_interactions:
            #dictionary to convert kinase names to gene names
            ks_genes_to_uniprot = {'ABL1(ABL)':'P00519', 'ACK':'Q07912', 'AURC':'Q9UQB9', 'ERK1(MAPK3)':'P27361','ERK2(MAPK1)':'P28482',  'ERK5(MAPK7)':'Q13164','JNK1(MAPK8)':'P45983', 'CK1A':'P48729', 'JNK2(MAPK9)':'P45984', 'JNK3(MAPK10)':'P53779', 'P38A(MAPK14)':'Q16539','P38B(MAPK11)':'Q15759', 'P38G(MAPK12)':'P53778','P70S6K' :'Q9UBS0', 'PAK':'Q13153', 'PKCZ':'Q05513', 'CK2A':'P19784', 'ABL2':'P42684', 'AMPKA1':'Q13131', 'AMPKA2':'Q13131', 'AURB':'Q96GD4', 'CAMK1A':'Q14012', 'CDC42BP':'Q9Y5S2','CK1D':'P48730','CK1E':'P49674','CK2B':'P67870','DMPK1':'Q09013', 'DNAPK':'P78527','DSDNA KINASE':'P78527', 'EG3 KINASE':'P49840','ERK3(MAPK6)':'Q16659','GSK3':'P49840', 'MRCKA':'Q5VT25', 'P38D(MAPK13)':'O15264','P70S6KB':'Q9UBS0','PDKC':'P78527','PKCH':'P24723','PKCI':'P41743','PKCT':'Q04759','PKD3':'O94806','PKG1':'Q13976','PKG2':'Q13237','SMMLCK':'Q15746'}
            if 'Combined:Kinase' in spliced_ptms.columns and not combined_added:
                if not spliced_ptms['Combined:Kinase'].isna().all():
                    print('Combined kinase-substrate data found and added')
                    interact = unify_interaction_data(spliced_ptms, 'Combined:Kinase', ks_genes_to_uniprot)
                    interact['Source'] = 'PSP/RegPhos'
                    interact_data.append(interact)
                    combined_added = True
            elif 'Combined:Kinase' not in spliced_ptms.columns:
                if 'RegPhos:Kinase' in spliced_ptms.columns and database == 'RegPhos':
                    if not spliced_ptms['RegPhos:Kinase'].isna().all():
                        print('RegPhos kinase-substrate data found and added')
                        interact = unify_interaction_data(spliced_ptms, 'RegPhos:Kinase', ks_genes_to_uniprot)
                        interact['Source'] = database
                        interact_data.append(interact)
                if 'PSP:Kinase' in spliced_ptms.columns and database == 'PhosphoSitePlus':
                    if not spliced_ptms['PSP:Kinase'].isna().all():
                        print('PhosphoSitePlus kinase-substrate data found and added')
                        interact = unify_interaction_data(spliced_ptms, 'PSP:Kinase', ks_genes_to_uniprot)
                        interact['Source'] = database
                        interact_data.append(interact)

            if database == 'DEPOD' and 'DEPOD:Phosphatase' in spliced_ptms.columns:
                if not spliced_ptms['DEPOD:Phosphatase'].isna().all():
                    print('DEPOD phosphatase-substrate data found and added')
                    interact = unify_interaction_data(spliced_ptms, 'DEPOD:Phosphatase')
                    interact['Source'] = database
                    interact_data.append(interact)
    
    if len(interact_data) > 0:
        interact_data = pd.concat(interact_data)
        extra_cols = [col for col in interact_data.columns if col in ['dPSI', 'Significance']]
        interact_data = interact_data.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Interacting ID', 'Type']+extra_cols, dropna = False, as_index = False)['Source'].apply(helpers.join_unique_entries)
    
        #convert uniprot ids back to gene names for interpretability
        ptm_gene = []
        interacting_gene = []
        for i, row in interact_data.iterrows():
            ptm_gene.append(pose_config.uniprot_to_genename[row['UniProtKB Accession'].split('-')[0]].split(' ')[0]) if row['UniProtKB Accession'].split('-')[0] in pose_config.uniprot_to_genename else ptm_gene.append(row['UniProtKB Accession'])
            interacting_gene.append(pose_config.uniprot_to_genename[row['Interacting ID'].split('-')[0]].split(' ')[0]) if row['Interacting ID'].split('-')[0] in pose_config.uniprot_to_genename else interacting_gene.append(row['Interacting ID'])
        interact_data['Modified Gene'] = ptm_gene
        interact_data["Interacting Gene"] = interacting_gene
  
  
        return interact_data.drop_duplicates()
    else:
        return pd.DataFrame()



def combine_KS_data(spliced_ptms, ks_databases = ['PhosphoSitePlus', 'RegPhos'], regphos_conversion = {'ERK1(MAPK3)':'MAPK3', 'ERK2(MAPK1)':'MAPK1', 'JNK2(MAPK9)':'MAPK9','CDC2':'CDK1', 'CK2A1':'CSNK2A1', 'PKACA':'PRKACA', 'ABL1(ABL)':'ABL1'}):
    """
    Given spliced ptm information, combine kinase-substrate data from multiple databases (currently support PhosphoSitePlus and RegPhos), assuming that the kinase data from these resources has already been added to the spliced ptm data. The combined kinase data will be added as a new column labeled 'Combined:Kinase'

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Spliced PTM data from project module
    ks_databases: list
        List of databases to combine kinase data from. Currently support PhosphoSitePlus and RegPhos
    regphos_conversion: dict
        Allows conversion of RegPhos names to matching names in PhosphoSitePlus.

    Returns
    -------
    splicde_ptms: pd.DataFrame
        Spliced PTM data with combined kinase data added
    
    """
    if not hasattr(pose_config, 'genename_to_uniprot'):
        pose_config.genename_to_uniprot = pose_config.flip_uniprot_dict(pose_config.uniprot_to_genename)

    ks_data = []
    for i, row in spliced_ptms.iterrows():
        combined = []
        for db in ks_databases:
            if db == 'PhosphoSitePlus':
                psp = row['PSP:Kinase'].split(';') if row['PSP:Kinase'] == row['PSP:Kinase'] else []
                #convert PSP names to a common name (first gene name provided by uniprot)
                psp = [pose_config.uniprot_to_genename[pose_config.genename_to_uniprot[kin]].split(' ')[0]  if kin in pose_config.genename_to_uniprot else kin for kin in psp]
                combined += psp
            elif db == 'RegPhos':
                regphos = row['RegPhos:Kinase'].split(';') if row['RegPhos:Kinase'] == row['RegPhos:Kinase'] else []
                for i, rp in enumerate(regphos):
                    if rp in pose_config.genename_to_uniprot:
                        regphos[i] = pose_config.uniprot_to_genename[pose_config.genename_to_uniprot[rp]].split(' ')[0]
                    elif rp.split('(')[0] in pose_config.genename_to_uniprot:
                        regphos[i] = pose_config.uniprot_to_genename[pose_config.genename_to_uniprot[rp.split('(')[0]]].split(' ')[0]
                    elif rp.upper() in regphos_conversion:
                        regphos[i] = regphos_conversion[rp.upper()]
                    else:
                        regphos[i] = rp.upper()
                combined += regphos


        if len(combined) > 0:
            ks_data.append(';'.join(set(combined)))
        else:
            ks_data.append(np.nan)

    spliced_ptms['Combined:Kinase'] = ks_data
    return spliced_ptms


def check_file(fname, expected_extension = '.tsv'):
    """
    Given a file name, check if the file exists and has the expected extension. If the file does not exist or has the wrong extension, raise an error.

    Parameters
    ----------
    fname: str
        File name to check
    expected_extension: str
        Expected file extension. Default is '.tsv'
    """
    if fname is None:
        raise ValueError('Annotation file path must be provided')
    if not os.path.exists(fname):
        raise ValueError(f'File {fname} not found')
    
    if not fname.endswith(expected_extension):
        raise ValueError(f'File {fname} does not have the expected extension ({expected_extension})')
    




def annotate_ptms(spliced_ptms, psp_regulatory_site_file = None, psp_ks_file = None, psp_disease_file = None, elm_interactions = False, elm_motifs = False, PTMint = False, PTMcode_interprotein = False, DEPOD = False, RegPhos = False, ptmsigdb_file = None, interactions_to_combine = ['PTMcode', 'PhosphoSitePlus', 'RegPhos', 'PTMInt'], kinases_to_combine = ['PhosphoSitePlus', 'RegPhos'], combine_similar = True):
    """
    Given spliced ptm data, add annotations from various databases. The annotations that can be added are the following:
    - PhosphoSitePlus 
        - regulatory site data (file must be provided)
        - kinase-substrate data (file must be provided)
        - disease association data (file must be provided)
    - ELM 
        - interaction data (can be downloaded automatically or provided as a file)
        - motif matches (elm class data can be downloaded automatically or provided as a file)
    - PTMInt
        - interaction data (will be downloaded automatically)
    - PTMcode
        - intraprotein interactions (can be downloaded automatically or provided as a file)
        - interprotein interactions (can be downloaded automatically or provided as a file)
    - DEPOD
        - phosphatase-substrate data (will be downloaded automatically)
    - RegPhos
        - kinase-substrate data (will be downloaded automatically)

    Parameters
    ----------
    spliced_ptms: pd.DataFrame
        Spliced PTM data from project module
    psp_regulatory_site_file: str
        File path to PhosphoSitePlus regulatory site data
    psp_ks_file: str
        File path to PhosphoSitePlus kinase-substrate data
    psp_disease_file: str
        File path to PhosphoSitePlus disease association data
    elm_interactions: bool or str
        If True, download ELM interaction data automatically. If str, provide file path to ELM interaction data
    elm_motifs: bool or str
        If True, download ELM motif data automatically. If str, provide file path to ELM motif data
    PTMint: bool
        If True, download PTMInt data automatically
    PTMcode_intraprotein: bool or str
        If True, download PTMcode intraprotein data automatically. If str, provide file path to PTMcode intraprotein data
    PTMcode_interprotein: bool or str
        If True, download PTMcode interprotein data automatically. If str, provide file path to PTMcode interprotein data
    DEPOD: bool
        If True, download DEPOD data automatically
    RegPhos: bool
        If True, download RegPhos data automatically
    ptmsigdb_file: str
        File path to PTMsigDB data
    interactions_to_combine: list
        List of databases to combine interaction data from. Default is ['PTMcode', 'PhosphoSitePlus', 'RegPhos', 'PTMInt']
    kinases_to_combine: list
        List of databases to combine kinase-substrate data from. Default is ['PhosphoSitePlus', 'RegPhos']
    combine_similar: bool
        Whether to combine annotations of similar information (kinase, interactions, etc) from multiple databases into another column labeled as 'Combined'. Default is True
    """
    if psp_regulatory_site_file is not None:
        try:
            check_file(psp_regulatory_site_file, expected_extension='.gz')
            spliced_ptms = add_PSP_regulatory_site_data(spliced_ptms, file = psp_regulatory_site_file)
        except Exception as e:
            raise RuntimeError(f'Error adding PhosphoSitePlus regulatory site data. Error message: {e}')
    if psp_ks_file is not None:
        try:    
            check_file(psp_ks_file, expected_extension='.gz')
            spliced_ptms = add_PSP_kinase_substrate_data(spliced_ptms, file = psp_ks_file)
        except Exception as e:
            raise RuntimeError(f'Error adding PhosphoSitePlus kinase-substrate data. Error message: {e}')
    if psp_disease_file is not None:
        try:
            check_file(psp_disease_file, expected_extension='.gz')
            spliced_ptms = add_PSP_disease_association(spliced_ptms, file = psp_disease_file)
        except Exception as e:
            raise RuntimeError(f'Error adding PhosphoSitePlus disease association data. Error message: {e}')
    if elm_interactions:
        try:
            if isinstance(elm_interactions, bool):
                spliced_ptms = add_ELM_interactions(spliced_ptms)
            elif isinstance(elm_interactions, str):
                check_file(elm_interactions, expected_extension='.tsv')
                spliced_ptms = add_ELM_interactions(spliced_ptms, file = elm_interactions)
            else:
                raise ValueError('elm_interactions must be either a boolean (download elm data automatically, slower) or a string (path to elm data tsv file, faster)')
        except Exception as e:
            raise RuntimeError(f'Error adding ELM interaction data. Error message: {e}')
    if elm_motifs:
        try:
            if isinstance(elm_motifs, bool):
                spliced_ptms = add_ELM_matched_motifs(spliced_ptms)
            elif isinstance(elm_motifs, str):
                check_file(elm_motifs, expected_extension='.tsv')
                spliced_ptms = add_ELM_matched_motifs(spliced_ptms, file = elm_motifs)
            else:
                raise ValueError('elm_interactions must be either a boolean (download elm data automatically, slower) or a string (path to elm data tsv file, faster)')
        except Exception as e:
            raise RuntimeError(f'Error adding ELM motif matches. Error message: {e}')
    if PTMint:
        try:
            if isinstance(PTMint, bool):
                spliced_ptms = add_PTMInt_data(spliced_ptms)
            elif isinstance(PTMint, str):
                check_file(PTMint, expected_extension='.csv')
                spliced_ptms = add_PTMInt_data(spliced_ptms, file = PTMint)
            else:
                raise ValueError('PTMint must be either a boolean (download PTMInt data automatically, slower) or a string (path to PTMInt data csv file, faster)')
        except Exception as e:
            raise RuntimeError(f'Error adding PTMInt interaction data. Error message: {e}')
    #if PTMcode_intraprotein:
    #    try:
    #        if isinstance(PTMcode_intraprotein, bool):
    #            spliced_ptms = add_PTMcode_intraprotein(spliced_ptms)
    #        elif isinstance(PTMcode_intraprotein, str):
    #            check_file(PTMcode_intraprotein, expected_extension='.gz')
    #            spliced_ptms = add_PTMcode_intraprotein(spliced_ptms, fname = PTMcode_intraprotein)
    #        else:
    #            raise ValueError('PTMcode_intraprotein must be either a boolean (download PTMcode data automatically, slower) or a string (path to PTMcode data file, faster)')
    #    except Exception as e:
    #        print(f'Error adding PTMcode intraprotein interaction data. Error message: {e}')
    if PTMcode_interprotein:
        try:
            if isinstance(PTMcode_interprotein, bool):
                spliced_ptms = add_PTMcode_interprotein(spliced_ptms)
            elif isinstance(PTMcode_interprotein, str):
                check_file(PTMcode_interprotein, expected_extension='.gz')
                spliced_ptms = add_PTMcode_interprotein(spliced_ptms, fname = PTMcode_interprotein)
            else:
                raise ValueError('PTMcode_interprotein must be either a boolean (download PTMcode data automatically, slower) or a string (path to PTMcode data file, faster)')
        except Exception as e:
            raise RuntimeError(f'Error adding PTMcode interprotein interaction data. Error message: {e}')
    if DEPOD:
        try:
            spliced_ptms = add_DEPOD_phosphatase_data(spliced_ptms)
        except Exception as e:
            raise RuntimeError(f'Error adding DEPOD phosphatase data. Error message: {e}')
    if RegPhos:
        try:
            if isinstance(RegPhos, str):
                check_file(RegPhos, expected_extension='.txt')
                spliced_ptms = add_RegPhos_data(spliced_ptms, file = RegPhos)
            else:
                spliced_ptms = add_RegPhos_data(spliced_ptms)
        except Exception as e:
            raise RuntimeError(f'Error adding RegPhos kinase substrate data data. Error message: {e}')
    if ptmsigdb_file is not None:
        try:
            spliced_ptms = add_PTMsigDB_data(spliced_ptms, file = ptmsigdb_file)
        except Exception as e:
            raise RuntimeError(f'Error adding PTMsigDB data. Error message: {e}')

    if combine_similar:
        interaction_cols = ['PTMcode:Interprotein_Interactions', 'PSP:ON_PROT_INTERACT', 'PSP:Kinase', 'PTMInt:Interaction', 'RegPhos:Kinase', 'DEPOD:Phosphatase']
        if set(interaction_cols).intersection(spliced_ptms.columns) != 0:
            print('\nCombining interaction data from multiple databases')
            interact = combine_interaction_data(spliced_ptms, interaction_databases = interactions_to_combine)
            if not interact.empty:
                interact['Combined:Interactions'] = interact['Interacting Gene']+'->'+interact['Type']
                interact = interact.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'], dropna = False, as_index = False)['Combined:Interactions'].apply(lambda x: ';'.join(np.unique(x)))
                if 'Combined:Interactions' in spliced_ptms.columns:
                    spliced_ptms = spliced_ptms.drop(columns = ['Combined:Interactions'])
    
                spliced_ptms = spliced_ptms.merge(interact, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'])
            else:
                spliced_ptms['Combined:Interactions'] = np.nan

        #check for what kinase data is available
        spliced_ptms = combine_KS_data(spliced_ptms, ks_databases=kinases_to_combine) #add combined kinase column


    return spliced_ptms


