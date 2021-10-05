#!/bin/bash -e

#phase haplotypes in filtered vcf file
java -jar /usr/bin/beagle.28Jun21.220.jat gt=#filtered_vcf_file out=#out_vcf_file

#
