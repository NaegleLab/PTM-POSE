import re

#packages for web interfacing
import requests
from requests.adapters import HTTPAdapter, Retry
import re



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

