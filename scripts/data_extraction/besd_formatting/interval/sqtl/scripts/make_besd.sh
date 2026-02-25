#!/bin/bash

in_dir="/local-scratch/data/tmp_processes/interval/sqtl/updated_flists"
out_dir="/local-scratch/data/hg38/interval/sqtl"


mkdir -p "$out_dir"

# Function to convert to besd
process_file() {
    local flist="$1"
    local out_dir="$2"


    local file_name=$(basename "$flist")
    local base_name="${file_name%.*}"

    # output path
    local besd_out="$out_dir/${base_name}"

    # Convert to besd
    echo "Converting $flist to besd..."
    smr --eqtl-flist "$flist" --make-besd --out "$besd_out"
    echo "$flist conversion complete, output written to $besd_out"
}

export -f process_file  # Export for GNU parallel
export out_dir          # Export variable


find "$in_dir" -type f | parallel -j 10 process_file {} "$out_dir"