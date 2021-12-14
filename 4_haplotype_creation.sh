#!/bin/bash -e

#phase haplotypes in filtered vcf file
java -jar /usr/bin/beagle.28Jun21.220.jar gt=#filtered_vcf_file out=#out_vcf_file

#
