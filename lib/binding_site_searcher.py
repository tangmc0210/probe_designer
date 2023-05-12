import os
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def site_searcher(tmp_output_pd, genbank_file=str, tmp='./tmp/' ,BDS_file_out_dir='./pre_binding_dir', 
                  pre_binding_file_suffix='_pre_binding.fasta', total_pre_binding_file_name='_total_pre_binding.fasta', 
                  pre_binding_num_file='4_pre_binding_num.json',
                  BDS_len=40, max_num=30):
    
    try: os.mkdir(BDS_file_out_dir)
    except: pass

    pre_binding_num = {}
    pos, length = 0, 0
    for gene_seq_in in SeqIO.parse(tmp+genbank_file, 'genbank'):
        # get information of gene
        id = gene_seq_in.id # get seq id
        # get gene_name
        if gene_seq_in.features:
            for feature in gene_seq_in.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ['NAN'])[0]

        # get molecule_type
        mol_type = gene_seq_in.annotations['molecule_type']
        # get organism
        organism = gene_seq_in.annotations['organism'] 
        # get minus seq
        translib = {'A':'T','T':'A','C':'G','G':'C'}
        try: 
            seq_minus = [translib[i] for i in str(gene_seq_in.seq)]
            seq = ''.join(list(reversed(seq_minus)))
        except: seq = str(gene_seq_in.seq)
        

        # set start point and pre_binding_num
        length = len(seq)
        pre_binding_num_tmp = min(length//BDS_len, max_num)
        st = length // 2 - pre_binding_num_tmp * BDS_len // 2
        # generate and write pre_binding for each gene in a fasta file
        record_list = []
        file_out = f'{BDS_file_out_dir}{id}{pre_binding_file_suffix}'

        for i in range(pre_binding_num_tmp):
            pre_binding_tmp = seq[st+i*BDS_len:st+i*BDS_len+BDS_len]
            record_list.append(SeqRecord(
                Seq(pre_binding_tmp),
                id='pre_binding'+str(i), 
                description='|'.join([id,gene_name,organism,mol_type])
            ))
            # add information about binding sites to tmp_output_pd
            tmp_output_pd.loc[i+pos, ['accession', 'gene_name', 'mol_type', 'organism', 'binding']] \
                = [id, gene_name, mol_type, organism, pre_binding_tmp]
            
        pos += pre_binding_num_tmp
        
        # write pre_binding to files
        with open(file_out, "w") as f:
            for new_record in record_list:
                SeqIO.write(new_record, f, "fasta")
        with open(BDS_file_out_dir+total_pre_binding_file_name, "a") as handle:
            for new_record in record_list:
                SeqIO.write(new_record, handle, "fasta")

        # record the num of pre_binding for each gene
        pre_binding_num[id] = [gene_name, pre_binding_num_tmp]

    with open(tmp+pre_binding_num_file, "w") as f:
        json.dump(pre_binding_num, f)

    return tmp_output_pd