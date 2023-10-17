import os
import time
import pandas as pd
from lib.Entrez_GI_and_genbank import get_genbank_from_GI
from lib.search_method import site_searcher
from lib.blast import blast_local, extract_blast_inf
from lib.seq_selection import select_wanted
from lib.results_merge import merge


workdir = "~/Github/Probe_designer/dataset"
gene_name_file = "immune.xlsx"

_ = pd.read_excel(gene_name_file)
total_gene_name_list = list(set(_["gene_name"]))
os.chdir(workdir)
gene_per_round = 30

for round in range(0, -(-len(total_gene_name_list) // gene_per_round)):
    current_time = time.localtime()
    formatted_time = time.strftime("%Y%m%d_%H%M%S", current_time)
    tmp = "./results/" + formatted_time + "/tmp/"
    output = "./results/" + formatted_time + "/"

    # Initiation of array
    binding_site_FOIs = [
        "accession",
        "gene_name",
        "mol_type",
        "organism",
        "binding",
        "Tm_l",
        "Tm_r",
        "wanted",
    ]
    align_FOIs = ["align_num", "align_accession", "align_descrip", "plus/minus"]
    tmp_output_pd = pd.DataFrame(columns=binding_site_FOIs + align_FOIs)

    gene_name_list = total_gene_name_list[
        round * gene_per_round, (round + 1) * gene_per_round
    ]

    # read gene name list
    with open(tmp + "1_gene_name_list.txt", "w") as f:
        f.write("\n".join(gene_name_list))

    # create gene id list file
    with open(tmp + "2_id_list.txt", "w") as f:
        f.write("")

    while True:
        ready = input('ready for "./results/time/tmp/2_id_list.txt" file? [True/False]')
        if ready:
            try:
                # read gene ID list
                with open(tmp + "2_id_list.txt", "r") as f:
                    id_list = f.read().split("\n")
                    break
            except:
                print('"./results/time/tmp/2_id_list.txt" file not found')
        else:
            print('please prepare for "./results/time/tmp/2_id_list.txt" file')

    # Get the genbank file of each gene by search for id list
    get_genbank_from_GI(tmp, "2_id_list.txt", genebank_file="3_gene_seq_in_file.gb")

    # Search binding sites on mRNA sequence
    tmp_output_pd = site_searcher(
        tmp_output_pd,
        BDS_len=40,
        max_num=30,
        gene_name_list=gene_name_list,
        tmp=tmp,
        genbank_file="3_gene_seq_in_file.gb",
        BDS_file_out_dir="./pre_binding_dir",
        pre_binding_file_suffix="_pre_binding.fasta",
        total_pre_binding_file_name="_total_pre_binding.fasta",
        pre_binding_num_file="4_pre_binding_num.json",
    )

    # sblast
    # # perform local blast
    # blast_local(
    #     tmp,
    #     BDS_file_out_dir="./pre_binding_dir",
    #     total_pre_binding_file_name="_total_pre_binding.fasta",
    #     blast_results_file="5_blast_results.xml",
    #     db="refseq_rna",
    #     task="megablast",
    # )

    while True:
        ready = input(
            'ready for "./results/time/tmp/5_blast_results.xml" file? [True/False]'
        )
        if ready:
            try:
                ## Extract interested information from blast_results
                tmp_output_pd = extract_blast_inf(
                    tmp=tmp,
                    blast_results_file="5_blast_results.xml",
                    tmp_output_pd=tmp_output_pd,
                )
                break
            except:
                print(
                    '"./results/time/tmp/5_blast_results.xml" file not found or file break, please download again'
                )
        else:
            print('please prepare for "./results/time/tmp/5_blast_results.xml" file')

    # ## Extract interested information from blast_results
    # tmp_output_pd = extract_blast_inf(
    #     tmp=tmp, blast_results_file="5_blast_results.xml", tmp_output_pd=tmp_output_pd
    # )

    # select wanted BDS

    gene_name_list_tosearch = select_wanted(
        tmp_output_pd,
        tmp,
        output,
        gene_name_list,
        gene_name_list_tosearch="gene_name_list_tosearch.txt",
    )

# merge
merge(workdir=workdir, resultsdir="./results/", gene_name_file=gene_name_file)
