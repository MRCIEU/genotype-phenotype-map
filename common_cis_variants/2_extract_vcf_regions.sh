#!/bin/bash

vcftools --vcf input.vcf --bed bed_file_describing_the_range.bed --out output_prefix --recode --keep-INFO-all

vcftools --regions chr:beg-end
