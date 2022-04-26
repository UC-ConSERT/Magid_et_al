#!/bin/bash -e
#script for aligning individual samples to species reference genome to create a population dataset

ref= #reference genome for alignment
datadir= #directory with trimmed fastq data
samdir= #sam file directory
bamdir= #bam file directory
bcf_file= #bcf file directory
fq1=_val_1.fq.gz #Read 1 suffix
fq2=_val_2.fq.gz #Read 2 suffix
platform="Illumina"

#now, retrieving read group and instrument information.
for samp in ${datadir}*${fq1} #remember to be explicit with file location
do
        #extract sample information from fastq name
	base=$(basename $samp _val_1.fq.gz)
	infoline=$(zcat ${samp} | head -n 1)
	instrument=`echo ${infoline} | cut -d ':' -f1`
	instrumentrun=`echo $infoline | cut -d ':' -f2`
	flowcell=`echo $infoline | cut -d ':' -f3`
	lane=`echo $infoline | cut -d ':' -f4`
	index=`echo $infoline | cut -d ':' -f10`

	#now to incorporate sample information into the alignment
	rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
	rgpl="PL:${platform}"
	rgpu="PU:${flowcell}.${lane}"
	rglb="LB:${base}_library1"
	rgsm="SM:${base}"

        echo "Aligning reads for $base" #align paired sample reads with ref genome using bwa v 0.7.17
        time bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 64 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam
	time bwa mem -M -t 64 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam
	echo "Converting sam file to bam file for $base" #convert sam file to bam file with SAMtools v 1.10
	time samtools view -@ 16 -T $ref -b ${samdir}${base}.sam > ${bamdir}${base}.bam

	echo "Sorting and indexing file"
	samtools sort -@ 16 -o ${bamdir}${base}.aligned.sorted.bam ${bamdir}${base}.bam #sorting file with samtools
	samtools index -b ${bamdir}${base}.aligned.sorted.bam #indexing file with samtools	
	rm ${samdir}${base}.sam #remove intermediate samtools file
done

#chunk bam files into 12 pieces using custom perl script 
#@Lanilen/SubSampler_SNPcaller/split_bamfiles_tasks.pl
#chunked files will help mpileup run faster
ls ${bamdir}*aligned.sorted.bam > ${bamdir}bam_list.txt
perl /data/SubSampler_SNPcaller/split_bamfiles_tasks.pl -b ${bamdir}bam_list.txt -g $ref -n 12 -o ${bamdir}/chunks | parallel -j 12 {}


#run bcftools mpileup in parallel on chunks of bam files with BCFtools v 1.11
for (( i=1; i<=12; i++ )); do
        bcftools mpileup -E -O b -f $ref -a AD,ADF,DP,ADR,SP -o ${bcf_file}OFK_${i}_raw.bcf ${bamdir}/chunks/${i}/* &
done
wait
echo “mpileup is done running”

#SNP calling on bcf files with bcftools call
for file in ${bcf_file}*.bcf
do
base=$(basename $file .bcf)
bcftools call $file -mv -O v -o ${bcf_file}${base}_VariantCalls.vcf &
done
wait
echo “variant calling is complete”


#prepare files for filtering
for file in ${bcf_file}*.vcf
do
base=$(basename $file .vcf)
#reheader each chunked bcf so it has the same sample names
bcftools reheader -s ${bamdir}OFK_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf
wait
#put bcf files names into a list for concatenation
ls ${bcf_file}${base}_reheader.bcf >> ${bcf_file}list_of_bcf.txt 
done

#concatenate the chunked bcf files into a whole population bcf
bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O b -o ${bcf_file}OFK_VariantCalls_concat.bcf --threads 16 
echo “population bcf file is ready for filtering!”


