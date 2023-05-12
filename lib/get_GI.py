from Bio import Entrez


def get_GI(tmp, gene_name_list_file, gene_id_name_file):
    # Get gene id and other information from ncbi dataset(api)
    ## Generate gene_search_list from gene_name_list
    organism_of_interest = 'Mus musculus'
    n_type_of_interest = 'mRNA'
    with open(tmp+gene_name_list_file) as f:
        gene_name_list = f.read().splitlines()
    gene_search_list = [', '.join([name, organism_of_interest, n_type_of_interest]) for name in gene_name_list]


    ## Get gene id list using Entrez.esearch
    id_list = []
    for gene_search in gene_search_list:
        Entrez.email = "1418767067@qq.com"
        handle = Entrez.esearch(db="nuccore", term=gene_search)
        record = Entrez.read(handle)
        handle.close()
        id_list += record["IdList"][:1] # set number of search results to read
    with open(tmp+gene_id_name_file, 'w') as f:
        f.write('\n'.join(id_list))