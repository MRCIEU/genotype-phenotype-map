library(here)
library(dplyr)
library(GenomeInfoDb)

# hg38 ld blocks
blocks <- read.table(here("pipeline_steps/data/ld_regions_hg38_updated.tsv"), header = T)

# chromosome lengths
hg38_chrs <- getChromInfoFromEnsembl("GRCh38", assembled.molecules.only = T) |> 
    filter(name %in% seq(1:22)) |> 
    arrange(as.numeric(name))

# start block
blocks_start <- blocks |> group_by(chr,ancestry) |> 
    filter(start == min(start)) |>
    mutate(end = start, start = 0)

# end block
blocks_end <- blocks |> group_by(chr,ancestry) |> 
    filter(end == max(end)) |>
    mutate(start = end, end = hg38_chrs[match(chr,hg38_chrs$name),"length"])

# check chr length > previous block end
# table(blocks$end - blocks$start > 0) 

blocks_updated <- rbind(blocks, blocks_start, blocks_end) |>
    arrange(ancestry,chr,start)

write.table(blocks_updated, here("pipeline_steps/data/ld_regions_hg38_fulllength.tsv"), 
    row.names = F, col.names = T, quote = F, sep = "\t")