# TLR_Analysis
Scripts for Alignment and Analysis of TLR sequences related to:
"A novel bioinformatic approach to characterise toll-like receptor gene diversity in threatened birds,"
a master's thesis characterising TLR gene diversity in tūturuatu/shore plover (Thinornis novaeseelandiae), 
kākāriki karaka/orange-fronted parakeet (Cyanoramphus malherbi), and kakī/black stilt (Himantopus novaezelandiae).
Thesis Permanent Link: https://ir.canterbury.ac.nz/handle/10092/101832.

Associated Peer-reviewed Publications:

All scripts written by @molmagid

# Files
[0_blast_tlr_search.sh](https://github.com/molmagid/TLR_Analysis/blob/main/0_blast_tlr_search.sh)
A script to find TLR sequences within a species' genome, using reference TLR sequences from related species
[1_trimming.sh](https://github.com/molmagid/TLR_Analysis/blob/main/1_trimming.sh)
Script to trim raw fastq files and analyse them with fastqc prior to aligment /t
[2_alignment_pipeline.sh](https://github.com/molmagid/TLR_Analysis/blob/15b320d9862f74d9636d6c373cdf5c702a46857a/2_alignment_pipeline.sh)
Script to align trimmed fastq files to species reference genome /t
[3_filtering.sh](https://github.com/molmagid/TLR_Analysis/blob/main/3_filtering.sh)
Script to filter population vcf file /t 
[4_haplotype_creation.sh](https://github.com/molmagid/TLR_Analysis/blob/main/4_haplotype_creation.sh)
Script to phase and create haplotypes for all individuals in the population vcf file
