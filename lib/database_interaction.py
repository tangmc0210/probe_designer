from Bio import Entrez


# NCBI_database
def get_GI(tmp, gene_name_list_file, gene_id_name_file):
    # Get gene id and other information from ncbi dataset(api)
    ## Generate gene_search_list from gene_name_list
    organism_of_interest = "Mus musculus"
    n_type_of_interest = "mRNA"
    with open(tmp + gene_name_list_file) as f:
        gene_name_list = f.read().splitlines()
    gene_search_list = [
        ", ".join([name, organism_of_interest, n_type_of_interest])
        for name in gene_name_list
    ]

    ## Get gene id list using Entrez.esearch
    id_list = []
    for gene_search in gene_search_list:
        Entrez.email = "1418767067@qq.com"
        handle = Entrez.esearch(db="nuccore", term=gene_search)
        record = Entrez.read(handle)
        handle.close()
        id_list += record["IdList"][:1]  # set number of search results to read
    with open(tmp + gene_id_name_file, "w") as f:
        f.write("\n".join(id_list))


def get_genbank_from_GI(tmp, gene_id_name_file, genebank_file="3_gene_seq_in_file.gb"):
    ## Read id_list from existing file
    with open(tmp + gene_id_name_file, "r") as f:
        id_list = f.read().split("\n")

    # Get the genbank file of each gene by search for id list
    fetch_per_round = 3
    round = -(-len(id_list) // fetch_per_round)
    for i in range(round):
        id_list_per_round = id_list[i * fetch_per_round : (i + 1) * fetch_per_round]
        Entrez.email = "1418767067@qq.com"
        handle = Entrez.efetch(
            db="nuccore",
            strand=1,  # plus if strand=1
            id=id_list_per_round,
            rettype="gbwithparts",
            retmode="text",
        )
        seq_record = handle.read()
        handle.close()
        print(i + 1, "{:.2f} %".format((i + 1) / round * 100))
        with open(tmp + genebank_file, "a") as f:
            f.write(seq_record)


# Ensembl_database
import requests


def enn_lookup(gene="BRCA1", species="homo_sapiens"):
    server = "http://rest.ensembl.org"
    ext = f"/lookup/symbol/{species}/{gene}?"
    options = ";".join(
        ["type=cds", "multiple_sequences=true", "content-type=application/json"]
    )
    response = requests.get(url=server + ext + options)
    lookup = response.json()

    id = lookup["id"]
    ext = f"/sequence/id/{id}?"
    options = ";".join(
        [
            f"species={species}",
            "type=cds",
            "multiple_sequences=true",
            "content-type=application/json",
        ]
    )
    response = requests.get(url=server + ext + options)
    sequence = response.json()

    return lookup, sequence


def fetch_seq(gene, species, method="ensembl"):
    if method == "ensembl":
        lookup, sequences = enn_lookup(gene=gene, species=species)
    return sequences
