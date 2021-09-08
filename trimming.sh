#!/bin/bash -e
rawdata= #directory with raw fastq data
datadir= #directory with fastq data

#run Fastqc on raw data files
for file in ${rawdata}*.gz
do
fastqc $file
done

#trim files before alignment- change parameters based raw file results
for file in ${rawdata}*_R1_001.fastq.gz
do
       base=$(basename ${file} _R1_001.fastq.gz)
       /usr/bin/TrimGalore-0.6.6/trim_galore /
       --paired --nextera --cores 6 --nextseq 28 /
       --length 50 --three_prime_clip_R1 5 --three_prime_clip_R2 5 /
       --fastqc --2colour 20 --clip_R1 20 --clip_R2 20 ${file} /
       ${rawdata}${base}_R2_001.fastq.gz -o ${datadir} --basename ${base}
done
wait
echo "done trimming"

#collect trimming reports into one file
for file in ${datadir}*.txt
do
cat $file >> trimming_reports.txt 
done


#run Fastqc on trimmed data files
for file in ${datadir}*.gz
do
fastqc $file
done
#examine fastqc output to determine whether additional trimming is needed
