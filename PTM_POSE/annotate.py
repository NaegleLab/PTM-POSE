import pandas as pd
import numpy as np
import re

from PTM_POSE import POSE_config


#dictionaries for converting modification codes to modification names in PhosphoSitePlus data
psp_dict = {'p': 'Phosphorylation', 'ca':'Caspase Cleavage', 'hy':'Hydroxylation', 'sn':'S-Nitrosylation', 'ng':'Glycosylation', 'ub': 'Ubiquitination', 'pa': "Palmitoylation",'ne':'Neddylation','sc':'Succinylation', 'sm': 'Sumoylation', 'ga': 'Glycosylation', 'gl': 'Glycosylation', 'ac': 'Acetylation', 'me':'Methylation', 'm1':'Methylation', 'm2': 'Dimethylation', 'm3':'Trimethylation'}
residue_dict = {'P': 'proline', 'Y':'tyrosine', 'S':'serine', 'T':'threonine', 'H':'histidine', 'D':'aspartic acid', 'I':'isoleucine', 'K':'lysine', 'R':'arginine', 'G':'glycine', 'N':'asparagine', 'M':'methionine'}



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
    regulatory_site_data = regulatory_site_data.rename(columns = {'ACC_ID':'UniProtKB Accession'})
    #drop extra modification information that is not needed
    regulatory_site_data['Residue'] = regulatory_site_data['MOD_RSD'].apply(lambda x: x.split('-')[0][0])
    regulatory_site_data['PTM Position in Canonical Isoform'] = regulatory_site_data['MOD_RSD'].apply(lambda x: int(x.split('-')[0][1:]))
    #add modification type
    regulatory_site_data['Modification Class'] = regulatory_site_data['MOD_RSD'].apply(lambda x: psp_dict[x.split('-')[1]])

    #restrict to human data
    regulatory_site_data = regulatory_site_data[regulatory_site_data['ORGANISM'] == 'human']
    regulatory_site_data = regulatory_site_data[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'ON_PROCESS', 'ON_PROT_INTERACT', 'ON_OTHER_INTERACT', 'ON_FUNCTION']].drop_duplicates()
    
    #group like modifications into a single column
    regulatory_site_data = regulatory_site_data.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).agg(lambda x: '; '.join([y for y in x if y == y])).reset_index()
    regulatory_site_data = regulatory_site_data.replace('', np.nan)
    
    #add 'PSP:' in front of each column
    regulatory_site_data.columns = ['PSP:' + x if x not in ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'] else x for x in regulatory_site_data.columns]

    #merge with spliced_ptm info
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(regulatory_site_data, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')

    
    #report the number of ptms with motif data
    num_ptms_with_known_function = spliced_ptms.dropna(subset = 'PSP:ON_FUNCTION').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    num_ptms_with_known_process = spliced_ptms.dropna(subset = 'PSP:ON_PROCESS').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    num_ptms_with_known_interaction = spliced_ptms.dropna(subset = 'PSP:ON_PROT_INTERACT').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    #num_ptms_in_domain = spliced_ptms.dropna(subset = 'PSP:DOMAIN').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"PhosphoSitePlus regulatory_site information added:\n\t ->{num_ptms_with_known_function} PTMs in dataset found associated with a molecular function \n\t ->{num_ptms_with_known_process} PTMs in dataset found associated with a biological process\n\t ->{num_ptms_with_known_interaction} PTMs in dataset found associated with a protein interaction")
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

    #separate residue and position
    ks_dataset['PTM Position in Canonical Isoform'] = ks_dataset['Residue'].apply(lambda x: int(x[1:]))
    ks_dataset['Residue'] = ks_dataset['Residue'].apply(lambda x: x[0])

    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(ks_dataset, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')
    
    
        #report the number of ptms with kinase substrate information
    num_ptms_with_KS = spliced_ptms.dropna(subset = 'PSP:Kinase').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"PhosphoSitePlus kinase-substrate interactions added: {num_ptms_with_KS} phosphorylation sites in dataset found associated with a kinase in PhosphoSitePlus")
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

    #drop extra modification information that is not needed
    #drop extra modification information that is not needed
    disease_associated_sites['Residue'] = disease_associated_sites['MOD_RSD'].apply(lambda x: x.split('-')[0][0])
    disease_associated_sites['PTM Position in Canonical Isoform'] = disease_associated_sites['MOD_RSD'].apply(lambda x: int(x.split('-')[0][1:]))
    #add modification type
    disease_associated_sites['Modification Class'] = disease_associated_sites['MOD_RSD'].apply(lambda x: psp_dict[x.split('-')[1]])
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

    #merge with spliced_ptm info
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(disease_associated_sites, how = 'left', on = ['UniProtKB Accession', 'Residue','PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataset size changed upon merge, please make sure there are no duplicates in spliced ptms data')
    
    #
    #report the number of ptms with motif data
    num_ptms_with_disease = spliced_ptms.dropna(subset = 'PSP:Disease_Association').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"PhosphoSitePlus disease associations added: {num_ptms_with_disease} PTM sites in dataset found associated with a disease in PhosphoSitePlus")
    
    
    return spliced_ptms


def add_ELM_interactions(spliced_ptms, fname = None):
    """
    Given a spliced ptms dataframe from the project module, add ELM interaction data to the dataframe
    """
    if fname is None:
        elm_interactions = pd.read_csv('http://elm.eu.org/interactions/as_tsv', sep = '\t', header = 0)
    else:
        elm_interactions = pd.read_csv(fname, sep = '\t', header = 0)
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

    spliced_ptms['ELM:Interactions'] = elm_list
    spliced_ptms['Location of PTM for ELM Interaction'] = elm_type
    spliced_ptms['Interacting Protein for ELM Interaction'] = elm_interactor
    
    #report the number of ptms with motif data
    num_ptms_with_ELM_instance = spliced_ptms.dropna(subset = 'ELM:Interactions').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"ELM interaction instances added: {num_ptms_with_ELM_instance} PTMs in dataset found associated with at least one known ELM instance")
    return spliced_ptms


def add_ELM_matched_motifs(spliced_ptms, ptm_coordinates, flank_size = 7):
    elm_classes = pd.read_csv('http://elm.eu.org/elms/elms_index.tsv', sep = '\t', header = 5)
    #create corresponding label for ptm_coordinate data
    ptm_coordinates['PTM Label'] = ptm_coordinates['UniProtKB Accession'] + '_' + ptm_coordinates['Residue'] + ptm_coordinates['PTM Position in Canonical Isoform'].apply(lambda x: int(x) if x == x else np.nan).astype(str)
    
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
            print(f'PTM {ptm} not found in PTM info file')
            match_list.append(np.nan)
    
    spliced_ptms['ELM:Motif Matches'] = match_list

    #report the number of ptms with motif data
    num_ptms_with_matched_motif = spliced_ptms.dropna(subset = 'ELM:Motif Matches').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"ELM Class motif matches found: {num_ptms_with_matched_motif} PTMs in dataset found with at least one matched motif")
    return spliced_ptms

def add_PTMint_data(spliced_ptms):
    """
    Given spliced_ptms data from project module, add PTMInt interaction data, which will include the protein that is being interacted with, whether it enchances or inhibits binding, and the localization of the interaction. This will be added as a new column labeled PTMInt:Interactions and each entry will be formatted like 'Protein->Effect|Localization'. If multiple interactions, they will be separated by a semicolon
    """
    #identify package directory
    #pdir = os.path.dirname(__file__)
    #if no directory has been provided, download the data to directory
    #if not os.path.exists(pdir+ "/Data/PTM_experimental_evidence.csv"):
    #    download_PTMint_data(pdir+'/Data/')

    PTMint = pd.read_csv('https://ptmint.sjtu.edu.cn/data/PTM%20experimental%20evidence.csv')
    PTMint = PTMint.rename(columns={'Uniprot':'UniProtKB Accession', 'AA':'Residue', 'Site':'PTM Position in Canonical Isoform'})
    #PTMint['Site'] = PTMint['AA'] + PTMint['Site'].astype(str)
    PTMint['PTMInt:Interaction'] = PTMint['Int_gene']+'->'+PTMint['Effect']
    PTMint = PTMint[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'PTMInt:Interaction']]
    #PTMint['PTM Position in Canonical Isoform'] = PTMint['PTM Position in Canonical Isoform'].astype(str)

    #aggregate PTMint data on the same PTMs
    PTMint = PTMint.groupby(['UniProtKB Accession','Residue','PTM Position in Canonical Isoform'], as_index = False).agg(';'.join)

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(PTMint[['UniProtKB Accession','Residue','PTM Position in Canonical Isoform', 'PTMInt:Interaction']], on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform'], how = 'left')
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')

    #report the number of PTMs identified
    num_ptms_with_PTMInt_data = spliced_ptms.dropna(subset = 'PTMInt:Interaction').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"PTMInt data added: {num_ptms_with_PTMInt_data} PTMs in dataset found with PTMInt interaction information")

    return spliced_ptms
    #delete source PTMint data
    #os.remove(pdir + './Data/PTM_experimental_evidence.csv')

def add_PTMcode_intraprotein(spliced_ptms, fname = None):
    #load ptmcode info
    if fname is None:
        ptmcode = pd.read_csv('https://ptmcode.embl.de/data/PTMcode2_associations_within_proteins.txt.gz', sep = '\t', header = 2, compression='gzip')
    else:
        ptmcode = pd.read_csv(fname, sep = '\t', header = 2, compression = 'gzip')
    
    #grab humn data
    ptmcode = ptmcode[ptmcode['Species'] == 'Homo sapiens']

    #add gene name to data
    translator = pd.DataFrame(POSE_config.uniprot_to_gene_name, index = ['Gene']).T
    translator['Gene'] = translator['Gene'].apply(lambda x: x.split(' '))
    translator = translator.explode('Gene')
    translator = translator.reset_index()
    translator.columns = ['UniProtKB/Swiss-Prot ID', 'Gene name']

    #add uniprot ID information
    ptmcode = ptmcode.merge(translator.dropna().drop_duplicates(), left_on = '## Protein', right_on = 'Gene name', how = 'left')

    #convert modification names to match annotation data
    convert_dict = {'Adp ribosylation': 'ADP Ribosylation', 'Glutamine deamidation':'Deamidation'}
    new_mod_names = []
    failed_mod = []
    mod_list = ptmcode['PTM1'].unique()
    for mod in mod_list:
        mod = mod.capitalize()
        if 'glycosylation' in mod: #if glycosylation, group into one gorup
            new_mod_names.append('Glycosylation')
        elif mod in POSE_config.modification_conversion['Modification Class'].values: #if already in modification class data, keep
            new_mod_names.append(mod)
        elif mod in convert_dict.keys():
            new_mod_names.append(convert_dict[mod])
        else:
            try:
                new_mod = POSE_config.modification_conversion[POSE_config.modification_conversion['Modification'] == mod].values[0][0]
                new_mod_names.append(new_mod)
            except:
                failed_mod.append(mod)
                new_mod_names.append(mod)
    conversion_df = pd.DataFrame({'PTM1':mod_list, 'Modification Class':new_mod_names})

    #add new modification labels to data
    ptmcode = ptmcode.merge(conversion_df, on = 'PTM1', how = 'left')
    
    #groupby by PTM1 and rename to match column names in annotation data
    ptmcode = ptmcode[['UniProtKB/Swiss-Prot ID', 'Modification Class', 'Residue1', 'Residue2']].dropna(subset = 'UniProtKB/Swiss-Prot ID')
    ptmcode = ptmcode.groupby(['UniProtKB/Swiss-Prot ID', 'Modification Class', 'Residue1'])['Residue2'].agg(';'.join).reset_index()
    ptmcode = ptmcode.rename(columns = {'UniProtKB/Swiss-Prot ID':'UniProtKB Accession', 'Residue1':'Residue', 'Residue2':'PTMcode:Intraprotein_Interactions'})
    
    #separate residue information into separate columns, one for amino acid and one for position
    ptmcode['PTM Position in Canonical Isoform'] = ptmcode['Residue'].apply(lambda x: int(x[1:]))
    ptmcode['Residue'] = ptmcode['Residue'].apply(lambda x: x[0])

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(ptmcode, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    num_ptms_with_PTMcode_data = spliced_ptms.dropna(subset = 'PTMcode:Intraprotein_Interactions').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"PTMcode intraprotein interactions added: {num_ptms_with_PTMcode_data} PTMs in dataset found with PTMcode intraprotein interaction information")

    return spliced_ptms

def extract_ids_PTMcode(df, col = '## Protein1'):

    #add gene name to data
    name_to_uniprot = pd.DataFrame(POSE_config.uniprot_to_genename, index = ['Gene']).T
    name_to_uniprot['Gene'] = name_to_uniprot['Gene'].apply(lambda x: x.split(' '))
    name_to_uniprot = name_to_uniprot.explode('Gene')
    name_to_uniprot = name_to_uniprot.reset_index()
    name_to_uniprot.columns = ['UniProtKB/Swiss-Prot ID', 'Gene name']
    name_to_uniprot = name_to_uniprot.drop_duplicates(subset = 'Gene name', keep = False)

    #protein name is provided as either ensemble gene id or gene name check for both
    df = df.merge(POSE_config.translator[['Gene stable ID']].reset_index().dropna().drop_duplicates(), left_on = col, right_on = 'Gene stable ID', how = 'left')
    df = df.rename(columns = {'index': 'From_ID'})
    df = df.merge(name_to_uniprot, left_on = col, right_on = 'Gene name', how = 'left')
    df = df.rename(columns = {'UniProtKB/Swiss-Prot ID': 'From_Name'})

    #grab unique id from 'From_ID' and 'From_Name' column, if available
    uniprot_ids = df['From_Name'].combine_first(df['From_ID'])
    return uniprot_ids.values

def add_PTMcode_interprotein(spliced_ptms, fname = None):
    if fname is None:
        ptmcode = pd.read_csv('https://ptmcode.embl.de/data/PTMcode2_associations_between_proteins.txt.gz', sep = '\t', header = 2, compression = 'gzip')
    else:
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
        elif mod in POSE_config.modification_conversion['Modification Class'].values:
            new_mod_names.append(mod)
        elif mod in convert_dict.keys():
            new_mod_names.append(convert_dict[mod])
        else:
            try:
                new_mod = POSE_config.modification_conversion[POSE_config.modification_conversion['Modification'] == mod].values[0][0]
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

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(ptmcode, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    num_ptms_with_PTMcode_data = spliced_ptms.dropna(subset = 'PTMcode:Interprotein_Interactions').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"PTMcode interprotein interactions added: {num_ptms_with_PTMcode_data} PTMs in dataset found with PTMcode interprotein interaction information")

    return spliced_ptms

def extract_positions_from_DEPOD(x):
    """
    Given string object consisting of multiple modifications in the form of 'Residue-Position' separated by ', ', extract the residue and position. Ignore any excess details in the string.
    """
    x = x.split(', ')
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

def add_DEPOD_phosphatase_data(spliced_ptms):
    #download data
    depod1 = pd.read_excel('https://depod.bioss.uni-freiburg.de/download/PPase_protSubtrates_201903.xls', sheet_name='PSprots')
    depod2 = pd.read_excel('https://depod.bioss.uni-freiburg.de/download/PPase_protSubtrates_newPairs_201903.xls', sheet_name = 'newPSprots')
    depod = pd.concat([depod1, depod2])

    #remove any rows with missing sit information
    depod = depod.dropna(subset = ['Dephosphosites'])

    #process dephosphosite strings to extract residue and position and explode so that each phosphosite is its own row
    depod['Dephosphosites'] = depod['Dephosphosites'].apply(extract_positions_from_DEPOD)
    depod = depod.explode('Dephosphosites')

    #separate multiple substrate accessions into their own rows (many of these link back to the same ID, but will keep just in case)
    depod['Substrate accession numbers'] = depod['Substrate accession numbers'].str.split(' ')
    depod = depod.explode('Substrate accession numbers')
    depod = depod.dropna(subset = ['Substrate accession numbers'])

    #extract only needed information and add phosphorylation as modification type
    depod = depod.rename({'Substrate accession numbers': 'UniProtKB Accession', 'Dephosphosites': 'Residue', 'Phosphatase entry names':'DEPOD:Phosphatase'}, axis = 1)
    depod = depod[['DEPOD:Phosphatase', 'UniProtKB Accession', 'Residue']]
    depod['Modification Class'] = 'Phosphorylation'

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(depod, how = 'left', on = ['UniProtKB Accession', 'Residue', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    num_ptms_with_PTMcode_data = spliced_ptms.dropna(subset = 'DEPOD:Phosphatase').groupby(['UniProtKB Accession', 'Residue']).size().shape[0]
    print(f"DEPOD Phosphatase substrates added: {num_ptms_with_PTMcode_data} PTMs in dataset found with Phosphatase substrate information")

    return spliced_ptms

def add_RegPhos_data(spliced_ptms, fname = None):
    if fname is None:
        regphos = pd.read_csv('http://140.138.144.141/~RegPhos/download/RegPhos_Phos_human.txt', sep = '\t')
    else:
        regphos = pd.read_csv(fname, sep = '\t')

    regphos = regphos.dropna(subset = 'catalytic kinase')
    #regphos['Residue'] = regphos['code'] + regphos['position'].astype(str)
    regphos = regphos.rename(columns = {'code': 'Residue', 'position':'PTM Position in Canonical Isoform', 'AC': 'UniProtKB Accession', 'catalytic kinase': 'RegPhos:Kinase'})
    regphos['Modification Class'] = 'Phosphorylation'
    regphos = regphos[['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class', 'RegPhos:Kinase']].dropna()
    regphos = regphos.groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class']).agg(';'.join).reset_index()

    #add to splice data
    original_data_size = spliced_ptms.shape[0]
    spliced_ptms = spliced_ptms.merge(regphos, how = 'left', on = ['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform', 'Modification Class'])
    if spliced_ptms.shape[0] != original_data_size:
        raise RuntimeError('Dataframe size has changed, check for duplicates in spliced ptms dataframe')
    
    #report the number of PTMs identified
    num_ptms_with_regphos_data = spliced_ptms.dropna(subset = 'RegPhos:Kinase').groupby(['UniProtKB Accession', 'Residue', 'PTM Position in Canonical Isoform']).size().shape[0]
    print(f"RegPhos kinase-substrate data added: {num_ptms_with_regphos_data} PTMs in dataset found with kinase-substrate information")

    return spliced_ptms
