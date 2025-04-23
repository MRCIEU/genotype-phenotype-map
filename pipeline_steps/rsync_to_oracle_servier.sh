# We can't easily rsync due to too many subdirectories on the server.
# So, the idea is to flatten the directory heirarchy for both finemapped_studies.tsv column and rsync

#!/bin/bash
set -e

RSYNC_DIR=$DATA_DIR/rsync_to_server
mkdir -p $RSYNC_DIR

FINEMAPPED_FILES=$RSYNC_DIR/finemapped_original_files.txt
FLATTENED_FINEMAPPED_FILES=$RSYNC_DIR/finemapped_flattened_files.txt
FLATTENED_DIR=$RSYNC_DIR/study/
REMOTE_STUDY_CONTENTS=$RSYNC_DIR/remote_study_contents.txt

mkdir -p $FLATTENED_DIR

cat /local-scratch/projects/genotype-phenotype-map/data/ld_blocks/*/*/*/finemapped_studies.tsv | awk '{print $3}' | grep -v file | sort > $FINEMAPPED_FILES
echo "files to sync: $(wc -l $FINEMAPPED_FILES)"
sed "s|${DATA_DIR}study/|${FLATTENED_DIR}|g; s|/finemapped/|_finemapped_|g" $FINEMAPPED_FILES > $FLATTENED_FINEMAPPED_FILES

echo "saving remote directory listing"
ssh opc@132.145.23.105 "find /oradiskvdb1/data/study/ | sort" > $REMOTE_STUDY_CONTENTS

echo "removing already synced files from transfer list"
comm -23 $FLATTENED_FINEMAPPED_FILES $REMOTE_STUDY_CONTENTS > ${FLATTENED_FINEMAPPED_FILES}.tmp
mv ${FLATTENED_FINEMAPPED_FILES}.tmp $FLATTENED_FINEMAPPED_FILES

echo "preparing flattened finemapped directory"
mkdir -p $FLATTENED_DIR
while IFS= read -r src <&3 && IFS= read -r dst <&4; do
    if [ ! -e "$dst" ]; then
        ln -s "$src" "$dst"
    fi
done 3<$FINEMAPPED_FILES 4<$FLATTENED_FINEMAPPED_FILES

echo "preparing ld_blocks directory"
mkdir -p $RSYNC_DIR/ld_blocks
cd $DATA_DIR
find ld_blocks -type f -name "finemapped_studies.tsv" -exec cp --parents {} $RSYNC_DIR \;
find $RSYNC_DIR/ld_blocks -type f -name "finemapped_studies.tsv" -exec sed -i "s|$DATA_DIR|/oradiskvdb1/data/|g; s|/finemapped/|_finemapped_|g" {} \;

cd $RSYNC_DIR

echo "rsyncing ld_blocks to oracle server"
rsync -aR ld_blocks/ opc@132.145.23.105:/oradiskvdb1/data/

echo "rsyncing flattened finemapped directory to oracle server"
rsync -avRL study/ opc@132.145.23.105:/oradiskvdb1/data/

rm -rf $FLATTENED_DIR