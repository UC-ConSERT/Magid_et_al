#!/bin/bash -e

#phase each haplotype into a filtered vcf file with Beagle v. 5.2
java -jar /usr/bin/beagle.28Jun21.220.jar gt=#filtered_vcf_file out=#out_vcf_file_name

#create consensus sequences with all ALT variants in a population using samtools faidx
#this will create a file for comparison to ref genome sequence
samtools faidx ref_genome.fa '#tlr_region chr:start-end'| bcftools consensus phased_tlr.vcf.gz -M N >> TLR_variants.fasta

#creates TLR haplotypes for all individuals and output sequences into a MEGA file format to use in DNAsp
counter=1 #counter adds unique number labels to each individual for MEGA format
echo "#MEGA" > tlr_haplotypes.meg
#for loop to output haplotypes for each individual bam file within the final vcf for one TLR region
for file in final_bams/*_merged.aligned.sorted.bam #loop through the bam files
do
base=$(basename $file _merged.aligned.sorted.bam)
echo "#"${counter}"-1" >> tlr_haplotypes.meg
#outputs sequence for phased haplotype 1 for one individual
samtools faidx ref_genome.fa #tlr_region | bcftools consensus phased_tlr.vcf.gz -s $base -M N -H 1pIu >> tlr_haplotypes.meg
echo "#"${counter}"-2" >> tlr_haplotypes.meg
#outputs sequence for phased haplotype 2 for one individual
samtools faidx ref_genome.fa #tlr_region | bcftools consensus phased_tlr.vcf.gz -s $base -M N -H 2pIu >> tlr_haplotypes.meg
counter=$(($counter+1))
done
