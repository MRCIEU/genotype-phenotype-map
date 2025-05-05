# We can't easily do a full rsync due to too many subdirectories on the server.
# So, the idea is to flatten the directory heirarchy for both finemapped_studies.tsv column and rsync

#!/bin/bash
set -e

RSYNC_DIR=$DATA_DIR/rsync_to_server
mkdir -p $RSYNC_DIR

REMOTE_SERVER_DATA_DIR=/oradiskvdb1/data
REMOTE_SERVER_STUDY_DIR=$REMOTE_SERVER_DATA_DIR/study
REMOTE_SERVER_DB_DIR=/oradiskvdb1/db
STUDY_DIR=${DATA_DIR}study

#STEP 1: Create a list of all the finemapped files to sync, squash the directory heirarchy

cd $STUDY_DIR

rsync -aRv --prune-empty-dirs --include='*/' --include='*finemapped/*' --exclude='*' . $ORACLE_SERVER:$REMOTE_SERVER_STUDY_DIR

#STEP 3: Alter the ld block information accordingly and rsync
echo "preparing ld_blocks directory"
cd $DATA_DIR/ld_blocks
mkdir -p $RSYNC_DIR/ld_blocks

cd $DATA_DIR
find ld_blocks -type f -name "finemapped_studies.tsv" -exec cp --parents {} $RSYNC_DIR \;
find $RSYNC_DIR/ld_blocks -type f -name "finemapped_studies.tsv" -exec sed -i "s|${STUDY_DIR}|$REMOTE_SERVER_STUDY_DIR|g;" {} \;

cd $RSYNC_DIR/ld_blocks

echo "rsyncing ld_blocks to oracle server"
rsync -aRv . $ORACLE_SERVER:$REMOTE_SERVER_DATA_DIR/ld_blocks_new/

#STEP 4: Copy over the recently created db files to the oracle server.  Note: not automatically copying the gwas_upload.db file

echo "copying db files to oracle server"
rsync -av $RESULTS_DIR/latest/studies.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/studies_new.db
rsync -av $RESULTS_DIR/latest/associations.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/associations_new.db
rsync -av $RESULTS_DIR/latest/ld.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/ld_new.db

touch /tmp/sync_done