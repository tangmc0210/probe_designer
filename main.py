import os
import time
import pandas as pd
from lib.get_genbank_from_GI import get_genbank_from_GI
from lib.binding_site_searcher import site_searcher
from lib.blast import blast_local, extract_blast_inf
from lib.select_wanted_BS import select_wanted
from lib.merge import merge


workdir = 'C:/Users/14187/Documents/Github/Probe_designer/dataset'

current_time = time.localtime()
formatted_time = time.strftime("%Y%m%d_%H%M%S", current_time)
tmp = './results/' + formatted_time + '/tmp/'
output = './results/' + formatted_time + '/'

os.chdir(workdir)

# Initiation of array
binding_site_FOIs = ['accession', 'gene_name', 'mol_type', 'organism', 'binding', 'wanted']
align_FOIs = ['align_num', 'align_accession', 'align_descrip', 'plus/minus']
tmp_output_pd = pd.DataFrame(columns=binding_site_FOIs+align_FOIs)


# read gene name list
with open(tmp+'1_gene_name_list.txt') as f:
    gene_name_list = f.read().splitlines()

# read gene ID list
with open(tmp+'2_id_list.txt', 'r') as f:
    id_list = f.read().split('\n')

# Get the genbank file of each gene by search for id list
get_genbank_from_GI(tmp, '2_id_list.txt', genebank_file='3_gene_seq_in_file.gb')

# Search binding sites on mRNA sequence
tmp_output_pd = site_searcher(tmp_output_pd, BDS_len=40, max_num=30, 
                            tmp=tmp ,genbank_file='3_gene_seq_in_file.gb', BDS_file_out_dir='./pre_binding_dir', 
                            pre_binding_file_suffix='_pre_binding.fasta', total_pre_binding_file_name='_total_pre_binding.fasta', 
                            pre_binding_num_file='4_pre_binding_num.json')


# blast
## perform local blast
blast_local(tmp, BDS_file_out_dir='./pre_binding_dir', total_pre_binding_file_name='_total_pre_binding.fasta', blast_results_file='5_blast_results.xml')

## Extract interested information from blast_results
tmp_output_pd = extract_blast_inf(s, tmp, blast_results_file='5_blast_results.xml')


# select wanted BDS
tmp_output_pd = select_wanted(tmp_output_pd, tmp, output, 
                  gene_name_list, gene_name_list_tosearch='gene_name_list_tosearch.txt')

# merge
merge(workdir=str ,resultsdir='./results/', gene_name_file='immune.xlsx')