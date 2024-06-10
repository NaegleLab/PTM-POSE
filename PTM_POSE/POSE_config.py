import pandas as pd
import numpy as np

#file processing packages
import os
import json


from PTM_POSE import database_interfacing as di

package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
resource_dir = package_dir + '/Resource_Files/'
modification_conversion = pd.read_csv(resource_dir + 'modification_conversion.csv')

if os.path.isfile(resource_dir + 'ptm_coordinates.csv'):
    ptm_coordinates = pd.read_csv(resource_dir + 'ptm_coordinates.csv',index_col = 0, dtype = {'Chromosome/scaffold name': str, 'PTM Position in Canonical Isoform': str})
else:
    print('ptm_coordinates file not found. Please run download_ptm_coordinates() to download the file from GitHub LFS. Set save = True to save the file locally and avoid downloading in the future.')
    ptm_coordinates = None

def download_ptm_coordinates(save = False):
    """
    Download ptm_coordinates dataframe from GitHub Large File Storage (LFS). By default, this will not save the file locally due the larger size (do not want to force users to download but highly encourage), but an option to save the file is provided if desired
    """
    ptm_coordinates = pd.read_csv('https://github.com/NaegleLab/PTM-POSE/raw/main/Resource_Files/ptm_coordinates.csv?download=', index_col = 0, dtype = {'Chromosome/scaffold name': str, 'PTM Position in Canonical Isoform': str})
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
    print('Translator file not found. Downloading mapping information between UniProt and Gene Names from pybiomart')
    uniprot_to_genename, uniprot_to_geneid = di.get_uniprot_to_gene()
    translator = pd.DataFrame({'Gene stable ID': uniprot_to_geneid, 'Gene name':uniprot_to_genename})
    #convert to dataframe, save to file
    translator.to_csv(resource_dir + 'translator.csv')

