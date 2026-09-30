Formatting DICE eQTL data for besd.

Data downloaded from: https://dice-database.org/downloads (DICE_DB_1, blood, unfiltered eQTLs).

Do this per cell type. The unfiltered VCFs contain every cis SNP-gene pair (all data
in the cis region), matching the onek1k besd.

#### Raw VCF format

- `CHROM` (`chrN`), `POS` (GRCh37), `REF`, `ALT`
- `INFO`: `Gene` (ENSG with version), `GeneSymbol`, `Pvalue`, `Beta`, `Statistic`, `FDR`

There is no `SE` or allele frequency in the VCF, so:

- `se = abs(Beta / Statistic)` (the `Statistic` is the MatrixEQTL ANOVA t-statistic)
- `Freq` is written as `NA`; the pipeline fills `EAF` from the LD reference panel

#### Steps

1. Run `format_dice_eqtl.R` (inside the `mrcieu/genotype-phenotype-map` image). For
   each VCF in `/local-scratch/data/raw/dice` it:
   - adds the missing `##contig` header lines (from `hg19.fa.fai`),
   - lifts the VCF from GRCh37 to GRCh38 with `bcftools +liftover`, using the same
     `hg19.fa`/`hg38.fa`/`hg19ToHg38` resources as the pipeline,
   - streams `bcftools query` output into R and reuses
     `standardise_alleles()` (and the SNP naming) from
     `pipeline_steps/common_extraction_functions.R`,
   - attaches GRCh38 gene coordinates from `gene_name_map_full.tsv`,
   - writes one `.esd` per gene plus a `.flist` and `.json` per cell type.

   Pass a cell type prefix (e.g. `B_CELL_NAIVE`) as an argument to process a single VCF.

2. Run `make_besd.sh` to build the `.besd`/`.epi`/`.esi` files with
   `smr --eqtl-flist <flist> --make-besd` and copy the json across.

Output: `/local-scratch/data/hg38/dice/<cell_type_prefix>.*`

#### Allele handling

`bcftools +liftover` records when it swaps REF/ALT in `INFO/SWAP` (this happens
~0.05% of the time). Because the DICE `Beta` is relative to the original ALT, the
beta is negated for `SWAP=1`. Multi-allelic results (`SWAP=-1`, a new reference
allele) are dropped since they no longer map onto a single effect allele.
`INFO/FLIP` (strand flip) does not change the effect direction, so no action.

The hg19->hg38 chain maps some terminal repeat regions across chromosomes (e.g. the
end of chr2 onto chr1/chr8). Because DICE is cis, variants are kept only when their
chromosome matches the gene's chromosome from `gene_name_map_full.tsv` (or the modal
chromosome for genes absent from the map). Without this, a gene ends up with variants
on several chromosomes and `smr --make-besd` fails with "Duplicated ESD file name".

#### Cell type mapping

| DICE file | cell_type | sample size |
|---|---|---|
| B_CELL_NAIVE | B IN | 91 |
| MONOCYTES | Mono C | 91 |
| M2 | Mono NC | 90 |
| NK | NK | 90 |
| TREG_MEM | Treg Mem | 89 |
| CD4_NAIVE | CD4 NC | 88 |
| CD4_STIM | CD4 STIM | 88 |
| TREG_NAIVE | Treg Naive | 89 |
| TFH | TFH | 89 |
| TH1 | TH1 | 81 |
| THSTAR | THSTAR | 88 |
| TH17 | TH17 | 89 |
| TH2 | TH2 | 89 |
| CD8_NAIVE | CD8 NC | 89 |
| CD8_STIM | CD8 STIM | 88 |

New labels (`CD4 STIM`, `CD8 STIM`, `TFH`, `TH1`, `TH2`, `TH17`, `THSTAR`,
`Treg Mem`, `Treg Naive`) are added to `cell_types` in `pipeline_steps/constants.R`.

#### columns needed for esd

Chr, SNP, Bp, A1, A2, Freq, Beta, se, p

#### columns needed for flist

Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd
