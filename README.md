# Probe_designer

## Brief introduction
This program is intended to design padlock probe for spatial transcriptome. When input with gene list you are interested in, the program will search for the corresponding mRNA sequence for the gene and, with some strategy, search for binding sites of specific length you assign.

## Binding site search strategy
Current strategy for binding site search is bruteforce. The program will begin in the middle of a mRNA sequence and extend with a step length(length of mRNA sequence divided by binding site length) to both side. You will get n binding site sequences, where n is the maxium of your intended number and the max number you can get from a mRNA sequence.

After that, the binding site sequences will be selected by blast(specificity), G content and consecutive G, which will make a lot of misses in the seq process because we use bicolor coding and G is double off.

## The process and tmp file in the process
### Get gene_id from ncbi dataset
By the gene name list given, from 

### Get genbank file of each mRNA from ncbi dataset

### Binding site searcher

### Blast and decoding of blast_results
*tip: blast length is limited to about 1200 sequences once, so split the gene name list if target gene number is too large.*

### Select wanted binding site by specific rules
*Present rules*:
1. sequence type is mRNA (judge by genbank file);
2. binding site sequence is specific to this gene in interested organisms (judge by blast results);
3. G percentage is between 40% and 70%;
4. G is not consecutive up to 5 or more;
5. binding site sequence is plus/minus to target gene; 

### Export *not-found* gene file for next round search
Considering the error occurred in the id or sequencing searching, sift out in selection, the files of each step are saved as tmp file in "results/time" folder.