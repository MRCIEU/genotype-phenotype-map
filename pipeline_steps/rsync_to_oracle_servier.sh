# We can't easily rsync due to too many subdirectories on the server.
# So, the idea is to flatten the directory heirarchy for both finemapped_studies.tsv column and rsync

#!/bin/bash
set -e

RSYNC_DIR=$DATA_DIR/rsync_to_server
mkdir -p $RSYNC_DIR

FINEMAPPED_FILES=$RSYNC_DIR/finemapped_original_files.txt
FLATTENED_FINEMAPPED_FILES=$RSYNC_DIR/finemapped_flattened_files.txt
FLATTENED_DIR=$RSYNC_DIR/finemapped_flattened/

cat /local-scratch/projects/genotype-phenotype-map/data/ld_blocks/EUR/*/*/finemapped_studies.tsv | awk '{print $3}' | grep -v file > $FINEMAPPED_FILES
sed "s|${DATA_DIR}study/||g; s|/finemapped/|_finemapped_|g" $FINEMAPPED_FILES > $FLATTENED_FINEMAPPED_FILES

echo "rsyncing finemapped directory to oracle server"
mkdir -p $FLATTENED_DIR
while IFS= read -r src <&3 && IFS= read -r dst <&4; do
    echo "rsyncing $src to $dst"
    rsync -a "$src" "opc@132.145.23.105:/oradiskvdb1/data/study/$dst"
    # if [ ! -e "$dst" ]; then
        # ln -s "$src" "$dst"
    # fi
done 3<$FINEMAPPED_FILES 4<$FLATTENED_FINEMAPPED_FILES

echo "preparing ld_blocks directory"
mkdir -p $RSYNC_DIR/ld_blocks
cd $DATA_DIR
find ld_blocks -type f -name "finemapped_studies.tsv" -exec cp --parents {} $RSYNC_DIR \;
find $RSYNC_DIR/ld_blocks -type f -name "finemapped_studies.tsv" -exec sed -i "s|$DATA_DIR|/oradiskvdb1/data/|g; s|/finemapped/|_finemapped_|g" {} \;

cd $RSYNC_DIR

echo "rsyncing ld_blocks to oracle server"
rsync -aR ld_blocks/ opc@132.145.23.105:/oradiskvdb1/data/