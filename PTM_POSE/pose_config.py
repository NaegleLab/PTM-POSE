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
    print('Translator file not found. Downloading mapping information between UniProt and Gene Names from pybiomart')
    uniprot_to_genename, uniprot_to_geneid = di.get_uniprot_to_gene()
    translator = pd.DataFrame({'Gene stable ID': uniprot_to_geneid, 'Gene name':uniprot_to_genename})
    #convert to dataframe, save to file
    translator.to_csv(resource_dir + 'translator.csv')

