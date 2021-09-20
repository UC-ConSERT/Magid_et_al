#!/bin/sh -e

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
