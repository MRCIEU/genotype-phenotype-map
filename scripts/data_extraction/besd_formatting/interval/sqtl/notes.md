converting interval sqtl data


Data downloaded from here: https://zenodo.org/records/10354433 

Do this per chromosome. 

colnames in eqlt summ stats: 

- phenotype_id - splice id - keep 
- variant_id - rsid
- splicemid_distance - same as tss - keep 
- af - eaf - keep 
- ma_samples
- ma_count
- pval_nominal - keep 
- slope - keep 
- slope_se - keep 
- chr - keep 
- pos_b38 - keep 
- effect_allele - keep
- other_allele - keep 
- pos_b7

summary file columns: - using this as a reference for tagged gene
- phenotype_id - keep 
- gene_id - keep 
- gene_name - keep 
- os
- variant_rsid
- splicemid_distance
- ma_samples
- ma_count
- af
- slope
- slope_se
- pval_nominal
- pval_beta
- qval
- qval_sig
- pval_nominal_threshold


columns needed for esd: 
- Chr
- SNP
- Bp
- A1
- A2
- Freq
- Beta
- se
- p


columns needed for flist: 
- Chr
- ProbeID
- GeneticDistance
- ProbeBp
- Gene
- Orientation
- PathOfEsd


#### continuing 
I think just need to upate the flist drafts then will be ready to convert to besd. 

renaming the jsons:

 for file in INTERVAL_eQTL_nominal_chr*.json; do
  mv "$file" "${file/eQTL/sQTL}"
done