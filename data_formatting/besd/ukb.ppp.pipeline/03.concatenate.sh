#!/bin/bash

#GZ_FILES="/home/rj18633/scratch/gp.map/data/besd.formatting/test.out/v1/v3/snps.updated"
#ESD_DIR="/home/rj18633/scratch/gp.map/data/besd.formatting/test.out/v1/v3/esd"
gz_files="$SNPS_UPDATED_DIR"
esd_dir="$ESD_DIR"


HEADER="Chr\tSNP\tBp\tA1\tA2\tFreq\tBeta\tse\tp"

# Ensure esd_dir exists
mkdir -p "$esd_dir"

# Process each subdirectory in gzfiles
for SUBDIR in "$gz_files"/*; do
    if [ -d "$SUBDIR" ]; then
        SUBDIR_NAME=$(basename "$SUBDIR")
        OUTFILE="$esd_dir/$SUBDIR_NAME.esd"

        # Start the output file with the header
        echo -e "$HEADER" > "$OUTFILE"

        # Concatenate each .gz file in the subdirectory
        for GZFILE in "$SUBDIR"/*.gz; do
            if [ -f "$GZFILE" ]; then
                # Remove the header from the .gz file and append to OUTFILE
                gunzip -c "$GZFILE" | tail -n +2 >> "$OUTFILE"
                # Remove the original .gz file
                rm "$GZFILE"
            fi
        done
    fi
done