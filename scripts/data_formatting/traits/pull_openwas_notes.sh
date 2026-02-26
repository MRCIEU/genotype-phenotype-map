#!/bin/bash

source ../.env

# Output file
output_file=${home}/pre_steps/data/opengwas_notes.txt
study_list=${home}/pre_steps/data/opengwas_ukbstudies.txt
missing_ids=${home}/pre_steps/data/tmp_missing.txt

#echo '"id","trait","note"' > $output_file

# Loop through each ukb-* study
for study in `cat $missing_ids` ; do
  file=${RAWDATA_DIR}/opengwas/igd/${study}/${study}.json

  if [ -f $file ]; then
    # Extract first required key-values using jq
    line=$(jq -r '[.id, .trait, .note] | @csv' $file)
    echo $line >> $output_file
  else
    echo "File not found: $file"
    echo $study >> $missing_ids
  fi
done

echo "Done! Output saved to $output_file"
