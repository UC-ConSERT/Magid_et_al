#!/bin/bash -e

#phase each haplotype into a filtered vcf file
java -jar /usr/bin/beagle.28Jun21.220.jar gt=#filtered_vcf_file out=#out_vcf_file

#create consensus sequences with all ALT variants for comparison to ref sequence
samtools faidx ref_genome.fa #tlr_region | bcftools consensus phased_tlr.vcf.gz -M N >> TLR_variants.fasta

#creates TLR haplotypes for all individuals and output into a MEGA file format
counter=1
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
