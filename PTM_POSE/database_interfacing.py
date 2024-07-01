import re

#packages for web interfacing
import requests
from requests.adapters import HTTPAdapter, Retry
import re

import numpy as np

from ptm_pose import project


#UniProt accession services adapted from suggested python code on UniProt website
def establish_session():
    """
    Establish a session for interfacing with the UniProt REST API. Taken from example REST API code provided by UniProt
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session, re_next_link

def get_next_link(headers, re_next_link):
    """
    Given a header, return the next link associated with the header. Taken from example REST API code provided by UniProt
    """
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url, session, re_next_link):
    """
    Given a URL, get the batch of data associated with the URL. Taken from example REST API code provided by UniProt
    """
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers, re_next_link)

def get_uniprot_to_gene(genename = True, geneid = True):
    """
    Construct a dictionary for converting from UniProt IDs to any gene names associated with that ID. Do this for all human, reviewed uniprot ids

    Returns
    -------
    dict
        Dictionary where keys are UniProt IDs and values are gene names associated with that ID separated by a space
    """
    #start up session for interfacting with rest api
    session, re_next_link = establish_session()
    if not genename and not geneid:
        raise ValueError('Must request at least one of genename or geneid')
    
    #establish query url
    fields = ['accession'] + ['gene_names'] * genename + ['xref_ensembl'] * geneid 
    fields = ','.join(fields)
    url =  f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606&format=tsv&fields={fields}_full&size=500"

    #run qeury in batch sizes of 500 and extract info
    uni_to_genename = {}
    uni_to_geneid = {}
    for batch, total in get_batch(url, session, re_next_link):
        for line in batch.text.splitlines()[1:]:
            results = line.split('\t')
            if genename and geneid:
                primaryAccession, gene_names, gene_id = results
            elif genename:
                primaryAccession, gene_names = results
            elif geneid:
                primaryAccession, gene_id = results
            

            if genename:
                uni_to_genename[primaryAccession] = gene_names
            if geneid:
                #ensembl reference exists, extract gene id from reference list and combine unique entries
                if len(gene_id) > 0:
                    gene_id = ' '.join(np.unique([id.split('.')[0].strip() for id in gene_id.split(';') if 'ENSG' in id]))
                uni_to_geneid[primaryAccession] = gene_id


    #return information
    if genename and geneid:
        return uni_to_genename, uni_to_geneid
    elif genename:
        return uni_to_genename
    elif geneid:
        return uni_to_geneid


def get_region_sequence(chromosome, strand, region_start, region_end, coordinate_type = 'hg38'):
    """
    Given a genomic region, return the DNA sequence associated with the region. Adapted from example REST API code provided by Ensembl

    Parameters
    ----------
    chromosome: str
        The chromosome of the region
    strand: str
        The strand of the region. Can be 1 (forward strand) or -1 (reverse strand)
    region_start: int
        The start of the region
    region_end: int
        The end of the region
    coordinate_type: str
        The coordinate system version to use. Can be 'hg38' or 'hg19'

    Returns
    -------
    str
        DNA sequence associated with the region
    """
    if coordinate_type == 'hg38':
        coord_system_version = 'GRCh38'
    elif coordinate_type == 'hg19':
        coord_system_version = 'GRCh37'


    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{region_start}..{region_end}:{strand}?coord_system_version={coord_system_version}"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    
    if not r.ok:
        r.raise_for_status()
        return None
    else:
        return r.text

def get_region_sequences_from_list(regions_list, coordinate_type = 'hg38'):
    """
    Given a list of genomic regions, return the DNA sequences associated with the regions. Adapted from example REST API code provided by Ensembl

    Parameters
    ----------
    regions_list: list of tuples
        Contain information about each different region. Each tuple in list should be in the format (chromosome, strand, region_start, region_end)
    coordinate_type: str
        The coordinate system version to use. Can be 'hg38' or 'hg19'
    
    Returns
    -------
    list of str
        List of DNA sequences associated with the regions in the same order as the inputted list
    """
    if coordinate_type == 'hg38':
        coord_system_version = 'GRCh38'
    elif coordinate_type == 'hg19':
        coord_system_version = 'GRCh37'
 
    region_list_str = '['
    region_coords = []
    for i,region_info in enumerate(regions_list):
        strand = project.convert_strand_symbol(region_info[1])
        coord = f'{region_info[0]}:{region_info[2]}..{region_info[3]}:{strand}'
        region_coords.append(coord)
        if i == len(regions_list) - 1:
            region_list_str += f'"{coord}"'
        else:
            region_list_str += f'"{coord}",'
    region_list_str = region_list_str + ']'


    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    r = requests.post(server+ext, headers=headers, data='{ "regions" : %s}' % region_list_str, params = {'coord_system_version':coord_system_version})
    
    if not r.ok:
        r.raise_for_status()
        return None
    else:
        decoded = r.json()

        #extract sequences, making sure they are in the same order as the inputted list
        seq_list = []
        for region in region_coords:
            #find seq info associated with query
            for result in decoded:
                if result['query'] == region:
                    seq_list.append(result['seq'])
                    break

        #return sequences
        return seq_list

