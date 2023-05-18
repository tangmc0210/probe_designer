# probe_designer

## Brief introduction

This program is intended to design padlock probes for spatial transcriptome. When input with gene list you are interested in, the program will search for the corresponding mRNA sequence for the gene and, with some strategy, search for binding sites of specific length you assign.

## Quickstart

### Create environment using conda or pip

```powershell
conda create --name <env_name> --file requirements.txt # where <env_name> could be probe_designer for example
pip install -r requirements.txt
```

### Prepare your dataset

Create a directory for you task, like *exmple_dataset*.

Prepare a .xlsx file including column "gene_name" which include target genes as *exmple_dataset/gene_name.xlsx*.

*Example dataset:*

| gene_name |... |
| --------  |--- |
| CD3D      |... |
| CD4      |... |
| CD8A      |... |
| ...       |... |

### Change work dir in main.py

`workdir = '~/Github/Probe_designer/dataset'`

## Binding site search strategy

Current strategy for binding site search is bruteforce. The program will begin in the middle of a mRNA sequence and extend with a step length(length of mRNA sequence divided by binding site length) to both side. You will get n binding site sequences, where n is the maxium of your intended number and the max number you can get from a mRNA sequence.

After that, the binding site sequences will be selected by gene and species specificity(by blast results), G content and consecutive G, which will make a lot of misses in the seq process because we use bicolor coding and G is double off.

## The process and tmp file in the process

### Get gene_id from ncbi dataset

Given the gene name list, the program search for target sequence using entrez provided by ncbi. However, the search results are usually mismatched so I recommand to search on website and choose target GI for next step analysis.

This step returns file `results/time/tmp/1_id_list.txt`.

### Get genbank file of each mRNA from ncbi dataset

The program downloads rna seq from ncbi according to GI number given in tmp file "1_id_list.txt". So prepare the id_list file yourself by searching on web if you find it awful using auto search in previous step.

This step returns file `results/time/tmp/2_gene_seq_in_file.gb`

### Binding site searcher

As mentioned in **Binding site search strategy**, current search stratagy is bruteforce. The program will begin in the middle of a mRNA sequence and extend with a step length(length of mRNA sequence divided by binding site length) to both side.

This part is most likely to be improved if you want to get better hybriding performance.

This step returns a series of .fasta file as `results/time/tmp/pre_biding/.fasta` where `_total_pre_binding.fasta` is to blast.

### Blast and decoding of blast_results

You can perform local blast or web blast on `_total_pre_binding.fasta` to get `results/time/tmp/4_blast_results.xml` file for decoding. You will need to build your local blast if you need to perform blast in this code. Otherwise, you can get blast results by uploading the file to blast web and download the blast results as `results/time/tmp/4_blast_results.xml`.

*Tip1: web blast length is limited to about 1200 sequences once, so split the gene name list if target gene number is too large.*

*Tip2: How to build local blast: [Download BLAST Software and Databases](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)*

### Select wanted binding site by specific rules

Current rules are divided into two kinds: pre_blast rules and post_blast rules. Pre_blast rules work before blast the potential binding site sequences and can often exempt most of sequences while post_blast rules work after blast and are intended to get sequence specifity information.

*Present pre_blast rules:*

1. G percentage is between 40% and 70%;
2. G is not consecutive up to 5 or more;
3. Sequence type is mRNA (judge by genbank file);

*Present post_blast rules:*

1. Binding site sequence is specific to this gene in interested organisms (judge by blast results);
2. Binding site sequence is plus/minus to target gene;

### Export *not-found* gene file for next round search

Considering the error occurred in the id or sequencing searching, sift out in selection, the files of each step are saved as tmp file in `results/time/gene_name_list_tosearch`.

### Merge the results got from above

run *merge.ipynb* to merge results got above and return a xlsx file of wanted binding sites for each gene.

## Perspective and plans

1. Improve the searching strategy to get better binding performance.