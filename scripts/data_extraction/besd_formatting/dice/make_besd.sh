#!/bin/bash

in_dir="/local-scratch/data/tmp_processes/dice/flists"
out_dir="/local-scratch/data/hg38/dice"
jobs="${JOBS:-4}"

mkdir -p "$out_dir/logs"

# Function to convert to besd
process_file() {
    local flist="$1"
    local out_dir="$2"

    local base_name="$(basename "${flist%.*}")"
    local besd_out="$out_dir/${base_name}"
    local log="$out_dir/logs/${base_name}.log"

    # Skip cell types that have already been converted
    if [[ -f "$besd_out.besd" && -f "$besd_out.epi" && -f "$besd_out.esi" && -f "$out_dir/${base_name}.json" ]]; then
        echo "SKIP: ${base_name} (already exists)"
        return 0
    fi

    # Convert to besd (quiet: only report pass/fail, full output in the log)
    if smr --eqtl-flist "$flist" --make-besd --out "$besd_out" > "$log" 2>&1 &&
        [[ -f "$besd_out.besd" ]]; then
        cp "${flist%.*}.json" "$out_dir/${base_name}.json"
        echo "OK: ${base_name}"
    else
        echo "FAILED: ${base_name} (see $log)"
        tail -n 20 "$log"
        return 1
    fi
}

export -f process_file  # Export for xargs bash subshells
export out_dir          # Export variable

find "$in_dir" -type f -name "*.flist" -print0 |
    xargs -0 -P "$jobs" -I {} bash -c 'process_file "$1" "$2"' _ {} "$out_dir"
