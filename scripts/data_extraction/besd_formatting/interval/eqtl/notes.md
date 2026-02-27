Formatting INTERVAL eQTL data for besd. 

Data downloaded from here: https://zenodo.org/records/10354433 

Do this per chromosome. 

colnames in eqlt summ stats: 

- phenotype_id (ENSG)	- keep 
- variant_id - risd
- tss_distance - keep 
- af (I believe this is eaf, after checking a few on their portal online and ferqs >0.5) - keep 
- ma_samples
- ma_count
- pval_nominal - keep 
- slope - keep 
- slope_se - keep 
- chr - keep 
- pos_b38 - keep 
- effect_allele - keep 
- other_allele - keep
- pos_b37

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



#### format eqtl dat
- run script on raw dat
- update flists 
- convert to besd 

