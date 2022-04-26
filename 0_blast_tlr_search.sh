#!/bin/sh -e
#script for align reference species TLR sequences sourced from NCBI 
#to search within the reference genome of a species of interest
#uses BLAST+ version 2.2.18


#make blast database for reference genome
makeblastdb -dbtype nucl species_genome.fasta

#blast each reference TLR fasta against genome database
for file in *TLR*.fasta
do
base=$(basename file .fasta)
blastn  -query $file -db species_genome.fasta -out ${base}.txt
done

#blast each reference TLR protein against genome database
for file in *TLR_protein*.fasta
do
base=$(basename file .fasta)
tblastn -query $file -db species_genome.fasta -out ${base}.txt
done
