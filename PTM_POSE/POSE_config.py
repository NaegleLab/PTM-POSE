import pandas as pd
import numpy as np

#file processing packages
import os
import json


from PTM_POSE import database_interfacing as di

package_dir = os.path.dirname(os.path.abspath(__file__))
modification_conversion = pd.read_csv(package_dir + '/../Resource_Files/modification_conversion.csv')



def download_translator():
    """
    Download information from biomart to convert between database IDs
    """
    chromosomes = ['X', '20', '1', '6', '3', '7', '12', '11', '4', '17', '2', '16',
       '8', '19', '9', '13', '14', '5', '22', '10', 'Y', '18', '15', '21',
       'MT']

    dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl',
                   host='http://www.ensembl.org')
    
    #load ID data that relates Ensembl to UniProt
    translator = dataset.query(attributes=['ensembl_gene_id','external_gene_name',
                                       'uniprotswissprot'],
             filters = {'transcript_biotype':'protein_coding',
             'chromosome_name':chromosomes})

    translator = translator.drop_duplicates()
    return translator


#load uniprot translator dataframe, process if need be
if os.path.isfile(package_dir + '/../Resource_Files/translator.csv'):
    translator = pd.read_csv(package_dir + '/../Resource_Files/translator.csv')
else:
    print('Translator file not found. Downloading mapping information between UniProt and Gene Names from pybiomart')
    translator = download_translator()

    #save to processed data directory
    translator.to_csv(package_dir + '/../Resource_Files/translator.csv')


#extract dictionary for converting from uniprot ids to gene names
if os.path.isfile(package_dir + '/../Resource_Files/uniprot_to_gene_name.json'):
    with open(package_dir + '/../Resource_Files/uniprot_to_gene_name.json', 'r') as f:
        uniprot_to_gene = json.load(f)
else:
    uniprot_to_gene = di.get_uniprot_gene_names()
    
    #save dictionary as json file
    with open(package_dir + '/../Resource_Files/uniprot_to_gene_name.json', 'w') as f:
        json.dump(uniprot_to_gene, f)


