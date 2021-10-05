#!/bin/sh

work= #directory where files to filter are
vcf_out=${work}filter_trial/
noLD=${vcf_out}/noLD
LD=${vcf_out}/LD_filter

#for loop to filter file with different values for parameters including
#missingness, depth, and GQ
for bcf in ${work}*_concat.bcf
do
    base=$(basename ${bcf} _concat.bcf)
    for i in {4..5} #filtering files for 4x and 5x depth
    do
        echo "Filtering SNPs for ${base}...." 
        stats --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_noMinGQ.bcf \
            --minDP ${i} \
            --maxDP 200  \ 
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        stats --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_noMinGQ.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        stats --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ10.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 10 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        stats --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ10.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 10 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        stats --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ20.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        stats --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ20.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all
    done;
done

echo "Filtering for Linkage parameters..."
#for loop to filter previous filtered files for linkage
for bcf in ${noLD}*.bcf
do
  base=$(basename ${bcf} .bcf)
    echo "Running light LD pruning at 0.8 for ${base}...."
   bcftools +prune \
        -l 0.8 \
        -w 1000 \
        -O b \
        -o ${LD}${base}_0.8LD_VariantCalls.bcf \
        ${bcf} &
    echo "Running moderate LD pruning at 0.6 for ${base}...."
    bcftools +prune \
        -l 0.6 \
        -w 1000 \
        -O b \
        -o ${LD}${base}_0.6LD_VariantCalls.bcf \
        ${bcf} &
    echo "Running strong LD pruning at 0.4 for ${base}...."
    bcftools +prune \
        -l 0.4 \
        -w 1000 \
        -O b \
        -o ${LD}${base}_0.4LD_VariantCalls.bcf \
        ${bcf}
done


mkdir /data/Tuturuatu/wild_data/results/bcf/stats/

#calculating statistics for no linkage filtered files
for file in ${noLD}*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --site-depth &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --depth &
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-site &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-indv &
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --het
done

#calculating statistics for linkage filtered files
for file in ${LD}*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --site-depth &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --depth &
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \        --missing-site &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-indv &
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --het
done
