import pandas as pd


def select_wanted(tmp_output_pd, tmp, output, 
                  gene_name_list, gene_name_list_tosearch):
    
    tmp_output_pd['wanted'] = [True]*len(tmp_output_pd)
    
    # sieve for the suitable binding site
    gene_name_list_out = [i for i in gene_name_list]
    
    for i in range(len(tmp_output_pd)):
        # check gene_name
        gene_name = tmp_output_pd.loc[i, 'gene_name']
        if gene_name not in gene_name_list:
            tmp_output_pd.loc[i, 'wanted'] = False
        else:
            try: gene_name_list_out.remove(gene_name)
            except: pass

        # check DNA or mRNA type
        if tmp_output_pd.loc[i, 'wanted'] == True:
            if tmp_output_pd.loc[i, 'mol_type'] != 'mRNA':
                tmp_output_pd.loc[i, 'wanted'] = False
            
        # check gene_organism name
        if tmp_output_pd.loc[i, 'wanted'] == True:
            spe_ori, gene_ori = tmp_output_pd.loc[i, 'organism'], tmp_output_pd.loc[i, 'gene_name']        
            descrip = tmp_output_pd.loc[i, 'align_descrip'].split('|')
            for des in descrip:
                if gene_ori not in des and spe_ori in des:
                    tmp_output_pd.loc[i, 'wanted'] = False
                    break

        # check plus/minus
        if tmp_output_pd.loc[i, 'wanted'] == True:
            if pd.isnull(tmp_output_pd.loc[i, 'plus/minus']):
                tmp_output_pd.loc[i, 'wanted'] = False
            else:    
                pm_list = tmp_output_pd.loc[i, 'plus/minus'].split(',')
                if '-1' not in pm_list:
                    tmp_output_pd.loc[i, 'wanted'] = False
        
        # check G 40%-70%, non consective 5 base
        if tmp_output_pd.loc[i, 'wanted'] == True:
            seq = tmp_output_pd.loc[i, 'binding']
            if 'GGGGG' in seq:
                tmp_output_pd.loc[i, 'wanted'] = False
            else:
                G_per = seq.count('G')/len(seq)
                if G_per < 0.3 or G_per > 0.7:
                    tmp_output_pd.loc[i, 'wanted'] = False

    # write the gene name to search next round
    with open(output+gene_name_list_tosearch, 'w') as f:
        f.write('\n'.join(gene_name_list_out))

    # write the whole information of interest to a excel file in tmp dir
    tmp_output_pd.to_excel(tmp+'6_probes_sieve.xlsx')

    # get the sub dataframe of wanted probes
    output_df = tmp_output_pd[tmp_output_pd['wanted'] == True]
    # write the output to a xlsx file
    output_df.to_excel(output+'probes_wanted.xlsx')