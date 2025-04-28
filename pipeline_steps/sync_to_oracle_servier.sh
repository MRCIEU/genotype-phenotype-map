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

#STEP 1: Create a list of all the finemapped files to sync, squash the directory heirarchy

mkdir -p $FLATTENED_DIR

cat /local-scratch/projects/genotype-phenotype-map/data/ld_blocks/*/*/*/finemapped_studies.tsv | awk '{print $3}' | grep -v file | sort > $FINEMAPPED_FILES
echo "files to sync: $(wc -l $FINEMAPPED_FILES)"
sed "s|${DATA_DIR}study/|${FLATTENED_DIR}|g; s|/finemapped/|_finemapped_|g" $FINEMAPPED_FILES > $FLATTENED_FINEMAPPED_FILES

echo "saving remote directory listing"
ssh $ORACLE_SERVER "find /oradiskvdb1/data/study/ | sort" > $REMOTE_STUDY_CONTENTS

echo "removing already synced files from transfer list"
comm -23 $FLATTENED_FINEMAPPED_FILES $REMOTE_STUDY_CONTENTS > ${FLATTENED_FINEMAPPED_FILES}.tmp
mv ${FLATTENED_FINEMAPPED_FILES}.tmp $FLATTENED_FINEMAPPED_FILES

echo "preparing flattened finemapped directory"
mkdir -p $FLATTENED_DIR
while IFS= read -r src <&3 && IFS= read -r dst <&4; do
    ln -s "$src" "$dst"
done 3<$FINEMAPPED_FILES 4<$FLATTENED_FINEMAPPED_FILES

#STEP 2: Alter the ld block information accordingly and rsync

echo "preparing ld_blocks directory"
mkdir -p $RSYNC_DIR/ld_blocks
cd $DATA_DIR
find ld_blocks -type f -name "finemapped_studies.tsv" -exec cp --parents {} $RSYNC_DIR \;
find $RSYNC_DIR/ld_blocks -type f -name "finemapped_studies.tsv" -exec sed -i "s|$DATA_DIR|/oradiskvdb1/data/|g; s|/finemapped/|_finemapped_|g" {} \;

cd $RSYNC_DIR

echo "rsyncing ld_blocks to oracle server"
rsync -aR ld_blocks/ $ORACLE_SERVER:/oradiskvdb1/data/

#STEP 3: Rsync the flattened finemapped directory to the oracle server

echo "rsyncing flattened finemapped directory to oracle server"
rsync -avRL study/ $ORACLE_SERVER:/oradiskvdb1/data/

#STEP 4: Copy over the recently created db files to the oracle server.  Note: not automatically copying the gwas_upload.db file

echo "copying db files to oracle server"
rsync -av $RESULTS_DIR/latest/studies.db $ORACLE_SERVER:/oradiskvdb1/db/studies_new.db
rsync -av $RESULTS_DIR/latest/associations.db $ORACLE_SERVER:/oradiskvdb1/db/associations_new.db
rsync -av $RESULTS_DIR/latest/ld.db $ORACLE_SERVER:/oradiskvdb1/db/ld_new.db

#STEP 5: Clean up

rm -rf $FLATTENED_DIR

touch /tmp/sync_done