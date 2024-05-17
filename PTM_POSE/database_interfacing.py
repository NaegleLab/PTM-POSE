import re

#packages for web interfacing
import requests
from requests.adapters import HTTPAdapter, Retry
import re

import project



#UniProt accession services adapted from suggested python code on UniProt website
def establish_session():
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session, re_next_link

def get_next_link(headers, re_next_link):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url, session, re_next_link):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers, re_next_link)

def get_uniprot_gene_names():
    """
    Construct a dictionary for converting from UniProt IDs to any gene names associated with that ID. Do this for all human, reviewed uniprot ids
    """
    #start up session for interfacting with rest api
    session, re_next_link = establish_session()

    url =  "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606&format=tsv&fields=accession,gene_names&size=500"
    id_to_gene = {}
    for batch, total in get_batch(url, session, re_next_link):
        for line in batch.text.splitlines()[1:]:
            primaryAccession, gene_names = line.split('\t')
            id_to_gene[primaryAccession] = gene_names
    return id_to_gene


def get_region_sequence(chromosome, strand, region_start, region_end, coordinate_type = 'hg38'):
    """
    Given a genomic region, return the DNA sequence associated with the region. Adapted from example REST API code provided by Ensembl
    """
    if coordinate_type == 'hg38':
        coord_system_version = 'GRCh38'
    elif coordinate_type == 'hg19':
        coord_system_version = 'GRCh37'

    #check if information is provided for a single gene or multiple genes
    #if isinstance(chromosome, list) and isinstance(strand, list) and isinstance(region_start, list) and isinstance(region_end, list):
    #    if len(chromosome) == len(strand) == len(region_start) == len(region_end):
    #        for i in range(len(chromosome)):

    #    else:
    #        raise ValueError('Length of region info lists is different, please make sure chromosome, strand, region_start, and region_end are all the same length')
        

    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{region_start}..{region_end}:{strand}?coord_system_version={coord_system_version}"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    
    if not r.ok:
        r.raise_for_status()
        return None
    else:
        return r.text

def get_region_sequences_from_list(regions_list, coordinate_type = 'hg38'):
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

