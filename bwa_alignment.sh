#!/bin/bash -e
ref=/data/Tuturuatu/ref_genomes/Maui_TLR.fasta #Reference genome for alignment
rawdata=/data/Tuturuatu/data/raw_data/
datadir=/data/Tuturuatu/data/trimmed_galore/ #Directory with fastq data
samdir=/data/Tuturuatu/tlr_results/sam/ #Sam file output
bamdir=/data/Tuturuatu/tlr_results/bam/ #Bam file output
bcf_file=/data/Tuturuatu/tlr_results/bcf/ #bcf file output
fq1=_val_1.fq.gz #Read 1 suffix
fq2=_val_2.fq.gz #Read 2 suffix
platform="Illumina"

#Trim files before alignment
for file in ${rawdata}*_R1_001.fastq.gz
do
	base=$(basename ${file} _R1_001.fastq.gz)
	/usr/bin/TrimGalore-0.6.6/trim_galore --paired --2colour 20 --basename $base -o ${datadir} $file ${rawdata}${base}_R2_001.fastq.gz
done
wait
echo "done trimming"

for file in ${datadir}*.txt
do
ls $file >> trimming_reports.txt
done

#First index the reference genome
time bwa index $ref

#Now, retrieving read group and instrument information.
for samp in ${datadir}*${fq1} #Remember to be explicit with file location
do
        base=$(basename $samp _val_1.fq.gz)
	infoline=$(zcat ${samp} | head -n 1)
	instrument=`echo ${infoline} | cut -d ':' -f1`
	instrumentrun=`echo $infoline | cut -d ':' -f2`
	flowcell=`echo $infoline | cut -d ':' -f3`
	lane=`echo $infoline | cut -d ':' -f4`
	index=`echo $infoline | cut -d ':' -f10`

	#Now to incorporate this information into the alignment
	rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
	rgpl="PL:${platform}"
	rgpu="PU:${flowcell}.${lane}"
	rglb="LB:${base}_library1"
	rgsm="SM:${base}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 64 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam
	time bwa mem -M -t 64 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam
	echo "Converting sam file to bam file for $base" 
	time samtools view -@ 16 -T $ref -b ${samdir}${base}.sam > ${bamdir}${base}.bam

	echo "Aligning and indexing file"
	samtools sort -@ 16 -o ${bamdir}${base}.aligned.sorted.bam ${bamdir}${base}.bam
	samtools index -@ 16 -b ${bamdir}${base}.aligned.sorted.bam	
	rm ${samdir}${base}.sam
done

#chunk bam files for mpileup
ls ${bamdir}*aligned.sorted.bam > ${bamdir}OFK_bam_list.txt
perl /data/SubSampler_SNPcaller/split_bamfiles_tasks.pl -b ${bamdir}OFK_bam_list.txt -g $ref -n 12 -o /data/OFK/tlr_results/chunks | parallel -j 12 {}

#run mpileup on chunks of bam files
for (( i=1; i<=12; i++ )); do
        bcftools mpileup -E -O b -f $ref -a AD,ADF,DP,ADR,SP -o ${bcf_file}OFK_${i}_raw.bcf /data/OFK/tlr_results/chunks/${i}/* &
done
wait
echo “mpileup is done running”

#variant calling on bcf files
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
bcftools reheader -s ${bamdir}OFK_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.vcf
bgzip ${bcf_file}${base}_reheader.vcf
bcftools index ${bcf_file}${base}_reheader.vcf.gz
ls ${bcf_file}${base}_reheader.vcf.gz >> ${bcf_file}list_of_vcf.txt
done

#concatenate the chunked vcf files
bcftools concat --file-list ${bcf_file}list_of_vcf.txt -O v -o ${bcf_file}OFK_VariantCalls_concat.vcf --threads 16 
echo “vcf file is ready for filtering!”


