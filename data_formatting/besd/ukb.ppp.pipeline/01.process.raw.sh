#!/bin/bash


# Use environment variables for input and output directories instead
input_dir="$INPUT_DIR"
output_dir="$PROCESSED_DIR"



# Function to process a single .gz file
process_gz_file() {
    local file=$1
    local output=$2
    base_name=$(basename "$file" .gz)
    output_file="$output/${base_name}_processed.gz"
    
    zcat "$file" | awk 'NR > 1 {
        # Split column 3 by ":" and keep the first two parts
        split($3, parts, ":")
        $3 = parts[1] ":" parts[2]
        base_position = parts[2]
        
        # get p-value from log10p 
        log10p = $13
        p_value = 10^(-log10p)

        # Print fields separated by tabs
        print $1 "\t" $3 "\t" base_position "\t" $5 "\t" $4 "\t" $6 "\t" $10 "\t" $11 "\t" p_value
    }' | gzip > "$output_file"
    
    echo "Processed $file and saved to $output_file"
}

# Function to concatenate processed .gz files and add header
concatenate_files() {
    local dir=$1
    local output_file_esd="$dir/$(basename "$dir").esd"
    local temp_file=$(mktemp)

    header="Chr\tSNP\tBp\tA1\tA2\tFreq\tBeta\tse\tp"

    # write header to the temp file
    echo -e "$header" > "$temp_file"
    
    # loop through all processed .gz files and concatenate them
    for file in "$dir"/*.gz; do
        zcat "$file" | awk '{gsub(/[\r\n]+$/, "", $0); print}' >> "$temp_file"
    done
    
    # Write the concatenated content to the .esd file
    mv "$temp_file" "$output_file_esd"
    
    echo "Concatenated chromosomes and saved to $output_file_esd"
    
    # rm the processed .gz files
    rm "$dir"/*.gz
    
    echo "Removed processed .gz files from $dir"
}

# loop through all .tar archives 
for tar_file in "$input_dir"/*.tar; do
    # create a temporary directory for this tar file
    tmp_dir=$(mktemp -d)
    
    # extract the .tar archive to the temporary directory
    tar -xvf "$tar_file" -C "$tmp_dir"
    
    # get the name of the directory created by the tar extraction
    extracted_dir=$(basename "$tar_file" .tar)
    extracted_path="$tmp_dir/$extracted_dir"
    
    # create a subdirectory for the processed files
    output_subdir="$output_dir/$extracted_dir"
    mkdir -p "$output_subdir"
    
    # check if there are any .gz files in the extracted directory
    gz_files=("$extracted_path"/*.gz)
    
    if [ -e "${gz_files[0]}" ]; then
        # process each .gz file in the extracted directory
        for gz_file in "${gz_files[@]}"; do
            process_gz_file "$gz_file" "$output_subdir"
        done
        
        # concatenate the processed files and add a header
       ##### REMOVING THIS PART FOR NOW concatenate_files "$output_subdir"
    else
        echo "No .gz files found in $tar_file"
    fi
    
    # clean up the temp directory for next archive
    rm -rf "$tmp_dir"
    
    echo "Processed all files in $tar_file"
done

echo "All tar archives processed."
