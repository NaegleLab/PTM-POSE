import pandas as pd
import numpy as np
import re


def add_PSP_regulatory_site_data(spliced_ptms, file = 'Regulatory_sites.gz'):
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
    #read in the kinase substrate data and add to spliced ptm info
    regulatory_site_data = pd.read_csv(file, sep = '\t', header = 2, on_bad_lines='skip',compression = 'gzip')
    #drop extra modification information that is not needed
    regulatory_site_data['Residue'] = regulatory_site_data['MOD_RSD'].apply(lambda x: x.split('-')[0])
    #restrict to human data
    regulatory_site_data = regulatory_site_data[regulatory_site_data['ORGANISM'] == 'human']
    regulatory_site_data = regulatory_site_data[['ACC_ID', 'Residue', 'DOMAIN', 'ON_PROCESS', 'ON_PROT_INTERACT', 'ON_OTHER_INTERACT', 'ON_FUNCTION']]
    #add 'PSP:' in front of each column
    regulatory_site_data.columns = ['PSP:' + x for x in regulatory_site_data.columns]

    #merge with spliced_ptm info
    spliced_ptms = spliced_ptms.merge(regulatory_site_data, how = 'left', left_on = ['UniProtKB Accession', 'Residue'], right_on = ['PSP:ACC_ID', 'PSP:Residue'])

    spliced_ptms = spliced_ptms.drop(columns = ['PSP:ACC_ID', 'PSP:Residue'], axis = 1)
    return spliced_ptms

def add_PSP_kinase_substrate_data(spliced_ptms, file = 'Kinase_Substrate_Dataset.gz'):
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
    ks_dataset = pd.read_csv(file, sep = '\t', header = 2, on_bad_lines='skip',compression = 'gzip', encoding = "cp1252")
    #restrict to human data
    ks_dataset = ks_dataset[ks_dataset['KIN_ORGANISM'] == 'human']
    ks_dataset = ks_dataset[ks_dataset['SUB_ORGANISM'] == 'human']

    ks_dataset = ks_dataset[['GENE', 'SUB_ACC_ID', 'SUB_MOD_RSD']].groupby(['SUB_ACC_ID', 'SUB_MOD_RSD']).agg(';'.join).reset_index()
    ks_dataset.columns = ['UniProtKB Accession', 'Residue', 'PSP:Kinase']

    spliced_ptms = spliced_ptms.merge(ks_dataset, how = 'left', on = ['UniProtKB Accession', 'Residue'])
    return spliced_ptms

def add_PSP_disease_association(spliced_ptms, file = 'Disease-associated_sites.gz'):
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
    disease_associated_sites = pd.read_csv(file, sep = '\t', header = 2, on_bad_lines='skip',compression = 'gzip')
    disease_associated_sites = disease_associated_sites[disease_associated_sites['ORGANISM'] == 'human']

    #removes sites without a specific disease annotation
    disease_associated_sites = disease_associated_sites.dropna(subset = ['DISEASE'])

    #combine disease and alteration
    disease_associated_sites['ALTERATION'] = disease_associated_sites.apply(lambda x: x['DISEASE']+'->'+x['ALTERATION'] if x['ALTERATION'] == x['ALTERATION'] else x['DISEASE'], axis = 1)
    #remove extra information from residue column
    disease_associated_sites['MOD_RSD'] = disease_associated_sites['MOD_RSD'].apply(lambda x: x.split('-')[0])
    #grab only necessary columns and rename
    disease_associated_sites = disease_associated_sites[['ACC_ID', 'MOD_RSD', 'ALTERATION']]
    disease_associated_sites.columns = ['UniProtKB Accession', 'Residue', 'PSP:Disease_Association']

    #aggregate multiple disease associations
    disease_associated_sites = disease_associated_sites.groupby(['UniProtKB Accession', 'Residue']).agg(';'.join).reset_index()

    #merge with spliced_ptm info
    spliced_ptms = spliced_ptms.merge(disease_associated_sites, how = 'left', on = ['UniProtKB Accession', 'Residue'])
    return spliced_ptms


def add_ELM_interactions(spliced_ptms, elm_interactions_fname):
    elm_interactions = pd.read_csv(elm_interactions_fname, sep = '\t')
    elm_interactions = elm_interactions[(elm_interactions['taxonomyElm'] == '9606(Homo sapiens)') & (elm_interactions['taxonomyDomain'] == '9606(Homo sapiens)')]

    elm_list = []
    elm_type = []
    elm_interactor = []
    for i, row in spliced_ptms.iterrows():
        #grab ptm location from residue column (gives residue and position (S981), so need to remove residue and convert to int)
        ptm_loc = int(row['Residue'][1:])

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

    spliced_ptms['ELM Interactions'] = elm_list
    spliced_ptms['Location of PTM for ELM Interaction'] = elm_type
    spliced_ptms['Interacting Protein for ELM Interaction'] = elm_interactor
    return spliced_ptms

def add_ELM_matched_motifs(spliced_ptms, ptm_info, elm_classes_fname, flank_size = 7):
    elm_classes = pd.read_csv(elm_classes_fname, sep = '\t', header = 5, index_col = 0)

    match_list = []
    for i, row in spliced_ptms.iterrows():
        matches = []
        #grab flanking sequence for the ptm
        ptm = row['UniProtKB Accession'] + '_' + row['Residue']
        ptm_flanking_seq = ptm_info.loc[ptm, 'Flanking Sequence']

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
    
    spliced_ptms['ELM Motif Matches'] = match_list
    return spliced_ptms




