from Bio import Entrez


def get_genbank_from_GI(tmp, gene_id_name_file, genebank_file='3_gene_seq_in_file.gb'):
    ## Read id_list from existing file
    with open(tmp+gene_id_name_file, 'r') as f:
        id_list = f.read().split('\n')

    # Get the genbank file of each gene by search for id list
    fetch_per_round = 3
    round = -(-len(id_list) // fetch_per_round)
    for i in range(round):
        id_list_per_round = id_list[i*fetch_per_round: (i+1)*fetch_per_round]
        Entrez.email = "1418767067@qq.com"
        handle = Entrez.efetch(db="nuccore", strand=1,  # plus if strand=1
                            id=id_list_per_round, rettype="gbwithparts", retmode="text")
        seq_record = handle.read()
        handle.close()
        print(i+1, '{:.2f} %'.format((i+1)/round*100))
        with open(tmp+genebank_file, 'a') as f:
            f.write(seq_record)