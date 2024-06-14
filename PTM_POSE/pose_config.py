import pandas as pd
import numpy as np

#base python packages
import os
import time


from ptm_pose import database_interfacing as di

#identify package directory
package_dir = os.path.dirname(os.path.abspath(__file__))
resource_dir = package_dir + '/Resource_Files/'

#download modification conversion file (allows for conversion between modificaiton subtypes and clases)
modification_conversion = pd.read_csv(resource_dir + 'modification_conversion.csv')

#load ptm_coordinates dataframe, if present
if os.path.isfile(resource_dir + 'ptm_coordinates.csv'):
    ptm_coordinates = pd.read_csv(resource_dir + 'ptm_coordinates.csv',index_col = 0, dtype = {'Chromosome/scaffold name': str, 'PTM Position in Canonical Isoform': str})
else:
    print('ptm_coordinates file not found. Please run download_ptm_coordinates() to download the file from GitHub LFS. Set save = True to save the file locally and avoid downloading in the future.')
    ptm_coordinates = None

def download_ptm_coordinates(save = False, max_retries = 5, delay = 10):
    """
    Download ptm_coordinates dataframe from GitHub Large File Storage (LFS). By default, this will not save the file locally due the larger size (do not want to force users to download but highly encourage), but an option to save the file is provided if desired

    Parameters
    ----------
    save : bool, optional
        Whether to save the file locally into Resource Files directory. The default is False.
    max_retries : int, optional
        Number of times to attempt to download the file. The default is 5.
    delay : int, optional
        Time to wait between download attempts. The default is 10.

    """
    for i in range(max_retries):
        try:
            ptm_coordinates = pd.read_csv('https://github.com/NaegleLab/PTM-POSE/raw/main/PTM_POSE/Resource_Files/ptm_coordinates.csv?download=', index_col = 0, dtype = {'Chromosome/scaffold name': str, 'PTM Position in Canonical Isoform': str})
            break
        except: 
            time.sleep(delay)
    else:
        raise Exception('Failed to download ptm_coordinates file after ' + str(max_retries) + ' attempts. Please try again.')

    

    if save:
        ptm_coordinates.to_csv(resource_dir + 'ptm_coordinates.csv')
    
    return ptm_coordinates


#load uniprot translator dataframe, process if need be
if os.path.isfile(resource_dir + 'translator.csv'):
    translator = pd.read_csv(resource_dir + 'translator.csv', index_col=0)
    uniprot_to_genename = translator['Gene name'].to_dict()
    uniprot_to_geneid = translator['Gene stable ID'].to_dict()

    #replace empty strings with np.nan
    translator = translator.replace('', np.nan)
else:
    print('Translator file not found. Downloading mapping information between UniProt and Gene Names from UniProt')
    uniprot_to_genename, uniprot_to_geneid = di.get_uniprot_to_gene()
    translator = pd.DataFrame({'Gene stable ID': uniprot_to_geneid, 'Gene name':uniprot_to_genename})
    #convert to dataframe, save to file
    translator.to_csv(resource_dir + 'translator.csv')


#additional information
#manually curated dictionary to convert phosphositeplus names that are not standard gene names to UniProt IDs
psp_name_dict = {'Actinfilin':'Q6TDP4','14-3-3 zeta':'P63104','14-3-3 epsilon':'P62258','14-3-3 sigma':'P31947','P130Cas':'P56945','ENaC-beta':'P51168','ENaC-alpha':'P37088','14-3-3 eta':'Q04917','14-3-3 beta':'P31946', '14-3-3 gamma':'P61981', '14-3-3 theta':'P27348','Securin':'O95997','GPIbA':'P07359','occludin':'Q16625','ER-beta':'Q92731','53BP1': 'Q12888','4E-T':'Q9NRA8','53BP2':'Q13625','AP-2 beta':'Q92481','APAF':'O14727','Bcl-xL':'Q07817','C/EBP-epsilon':'Q15744','CREB':'P16220','Calmodulin':'P0DP23','Cortactin':'Q14247','DNAPK':'P78527', 'Diaphanous-1':'O60610', 'ER-alpha':'P03372', 'Exportin-1':'O14980', 'Ezrin':'P15311', 'H3':'Q6NXT2','HSP70':'P0DMV8;P0DMV9','IKKG':'Q9Y6K9', 'Ig-beta':'P40259','Ku80':'P13010','LC8':'Q96FJ2', 'MRLC2V':'P10916', 'Merlin':'P35240','NFkB-p105':'P19838', 'Rb':'P06400', 'RhoGDI alpha':'P52565', 'Rhodopsin':'P08100', 'SHP-1':'P29350', 'SHP-2':'Q06124','SLP76':'Q13094','SMRT':'Q9Y618','SRC-3':'Q9Y6Q9','STI1':'Q9BPY8','Vinculin':'P18206','beclin 1':'Q14457','claspin':'Q9HAW4', 'gp130':'P40189','leupaxin':'O60711','p14ARF':'Q8N726','rubicon':'Q92622','snRNP A':'P09661','snRNP B1':'P08579','snRNP C':'P09234','syntenin':'O00560;Q9H190','talin 1':'Q9Y490', 'ubiquitin':'P0CG47', '4E-BP1':'Q13541', 'ALK2':'Q04771', 'AMPKA1':'Q13131','AurA':'O14965','AurB':'Q96GD4', 'AurC':'Q9UQB9', 'C/EBP-beta':'P17676', 'CAMK1A':'Q14012', 'CHD-3 iso3':'Q12873', 'CK1A':'P48729', 'CK2B':'P67870', 'DAT':'Q01959', 'DJ-1':'Q99497', 'DOR-1':'P41143', 'DYN1':'Q05193','Desmoplakin':'P15924', 'Exportin-4':'Q9C0E2', 'FBPase':'P09467', 'FBPase 2':'O60825', 'G-alpha':'P63096', 'G-alpha 13':'Q14344', 'G-alpha i1':'P63096', 'G-beta 1':'P62873', 'G-beta 2':'P62879', 'G6PI':'P06744', 'GM130':'Q08379', 'GR':'P04150', 'H4':'P62805', 'HP1 alpha':'P45973', 'IkB-alpha':'P25963', 'IkB-beta':'Q15653', 'PPAR-gamma':'P37231', 'Claudin-1':'O95832', 'Claudin-2':'P57739', 'Cofilin-1':'P23528', 'K14':'P02533', 'K18':'P05783', 'K5':'P13647','K8':'P05787','Ku70':'P12956', 'Moesin':'P26038','N-WASP':'O00401','Nur77':'P22736','P38A':'Q16539','P38B':'Q15759', 'P70S6KB':'P23443','PGC-1 alpha':'Q9UBK2','PKHF1':'Q96S99','P38G':'P53778','PKCI':'P41743','PKCZ':'Q05513', 'PKG1':'Q13976', 'PTP-PEST':'Q05209','Plectin-1':'Q15149','RFA2':'P15927','SERCA2':'P16615','SH2-B-beta':'Q9NRF2', 'SNAP-alpha':'P54920', 'SPT16':'Q9BXB7', 'SPT6':'Q7KZ85','STEP':'P54829','STLK3':'Q9UEW8', 'Snail1':'O95863', 'Snail2':'O43623', 'Stargazin':'P62955','Survivin':'O15392','TARP':'P09693','TK':'P04183','TOM20':'Q15388','TR-alpha':'P10827','Titin':'Q8WZ42','Vimentin':'P08670','WASP':'P42768','ZAP':'Q7Z2W4',  'Zyxin':'Q15942', 'cIAP1':'Q13490','caveolin-1':'Q03135', 'coronin 2A':'Q92828', 'desmin':'P17661','eIF2-alpha':'Q9BY44', 'eIF2-beta':'P20042', 'eIF3-alpha':'O75822', 'eIF3-eta':'P55884', 'eIF3-zeta':'O15371', 'eNOS':'P29474', 'emerin':'P50402', 'epsin 1':'Q9Y6I3', 'glutaminase':'O94925','hnRNP A1':'P09651', 'hnRNP A2/B1':'P22626', 'hnRNP A3':'P51991','hnRNP D0':'Q14103', 'hnRNP E2':'Q15366','hnRNP P2':'P35637','hnRNP U':'Q00839', 'kindlin-2':'Q96AC1', 'kindlin-3':'Q86UX7','lamin A/C':'P02545', 'mucolipin 1':'Q9GZU1','nNOS':'Q8WY41','p21Cip1':'P38936', 'p27Kip1':'P46527','p47phox':'P14598','p90RSK':'Q15418','palladin':'Q8WX93','polybromo 1':'Q86U86', 'syndecan-4':'P31431', 'tensin 1 iso1':'Q9HBL0', 'utrophin':'P46939','DKFZp686L1814':'Q6MZP7', 'EB1':'Q15691', 'EB2':'Q15555', 'G-alpha i3':'P08754','HSP20':'O14558','HSP40':'P25685', 'Hic-5':'O43294', 'Ig-alpha':'P11912', 'LC3A':'Q9H492', 'LC3B':'Q9GZQ8', 'LC3C':'Q9BXW4','NFkB-p100':'Q00653','NFkB-p65':'Q04206','Pnk1':'Q96T60', 'RPT2':'P62191','EB3':'Q9UPY8'}



def flip_uniprot_dict(uniprot_dict):
    """
    Given one of the uniprot id to gene name or gene id dictionaries, flip the dictionary so that the gene name or id is the key and the uniprot id is the value
    """
    uniprot_dict = pd.DataFrame(uniprot_dict, index = ['Gene']).T.reset_index()
    uniprot_dict['Gene'] = uniprot_dict['Gene'].str.split(' ')
    uniprot_dict = uniprot_dict.explode('Gene')
    uniprot_dict = uniprot_dict.set_index('Gene')['index'].to_dict()
    return uniprot_dict