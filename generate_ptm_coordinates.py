
import numpy as np
import pandas as pd
from ExonPTMapper import mapping

import datetime
import pyliftover
from tqdm import tqdm


phosphosite_data = pd.read_csv('../../Database_Information/PhosphoSitePlus/reformatted_sites/phosphositeplus_data.csv', index_col = 0)
mapper = mapping.PTM_mapper()

mapper.find_ptms_all()
mapper.mapPTMs_all()


ptm_coordinates = mapper.ptm_coordinates.copy()
new_coords = []
liftover_object = pyliftover.LiftOver('hg38','hg19')
for i, row in tqdm(ptm_coordinates.iterrows(), total = ptm_coordinates.shape[0], desc = 'Converting from hg38 to hg19 coordinates'):
    new_coords.append(mapping.convert_to_new_coordinates(row['Gene Location (hg38)'], row['Chromosome/scaffold name'], row['Strand'], liftover_object = liftover_object))
ptm_coordinates['Gene Location (hg19)'] = new_coords

new_coords = []
liftover_object = pyliftover.LiftOver('hg19','hg18')
for i, row in tqdm(ptm_coordinates.iterrows(), total = ptm_coordinates.shape[0], desc = 'Converting from hg19 to hg18 coordinates'):
    new_coords.append(mapping.convert_to_new_coordinates(row['Gene Location (hg19)'], row['Chromosome/scaffold name'], row['Strand'], liftover_object = liftover_object))
ptm_coordinates['Gene Location (hg18)'] = new_coords

ptm_coordinates['PTM Position in Canonical Isoform'] = ptm_coordinates['PTM Position in Canonical Isoform'].str.split(';')
ptm_coordinates['UniProtKB Accession'] = ptm_coordinates['UniProtKB Accession'].str.split(';')
ptm_coordinates = ptm_coordinates.explode(['UniProtKB Accession', 'PTM Position in Canonical Isoform']).reset_index()

ptm_coordinates.reset_index().to_csv('../PTM_POSE/ptm_pose/Resource_Files/ptm_coordinates.csv', index = False)

#write to text file indicating when the data was last updated
with open('./ptm_pose/Resource_Files/last_updated.txt', 'w') as f:
    f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))