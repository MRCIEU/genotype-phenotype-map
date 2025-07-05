#!/bin/bash
set -e

RSYNC_DIR=$DATA_DIR/rsync_to_server
mkdir -p $RSYNC_DIR

REMOTE_SERVER_DATA_DIR=/oradiskvdb1/data
REMOTE_SERVER_STUDY_DIR=$REMOTE_SERVER_DATA_DIR/study
REMOTE_SERVER_SVG_DIR=/oradiskvdb1/static/svgs
REMOTE_SERVER_DB_DIR=/oradiskvdb1/db
STUDY_DIR=${DATA_DIR}study

cd $STUDY_DIR

#Sync the finemapped GWAS results to the oracle server
#TODO: ucomment that back once we have deleted the old finemapped studies
# rsync -aRv --prune-empty-dirs --include='*/' --include='*finemapped/*' --exclude='*' . $ORACLE_SERVER:$REMOTE_SERVER_STUDY_DIR

#Sync the full svg files to the oracle server
# find . -type f -name "full.*" -exec cp --parents {} $RSYNC_DIR \;
# rsync -aRv --prune-empty-dirs --include='*/' --include='*svg/full*' --exclude='*' . $ORACLE_SERVER:$REMOTE_SERVER_SVG_DIR
find . -type f -regex '.*/[0-9]\+\.\(zip\|json\)$' -exec rsync -av --remove-source-files {} $ORACLE_SERVER:$REMOTE_SERVER_SVG_DIR/full \;

#Alter the ld block information accordingly and rsync
echo "preparing ld_blocks directory"
cd $DATA_DIR/ld_blocks
mkdir -p $RSYNC_DIR/ld_blocks

cd $DATA_DIR
find ld_blocks -type f -name "finemapped_studies.tsv" -exec cp --parents {} $RSYNC_DIR \;
find $RSYNC_DIR/ld_blocks -type f -name "finemapped_studies.tsv" -exec sed -i "s|${STUDY_DIR}|$REMOTE_SERVER_STUDY_DIR|g;" {} \;

cd $RSYNC_DIR/ld_blocks

echo "rsyncing ld_blocks to oracle server"
rsync -aRv . $ORACLE_SERVER:$REMOTE_SERVER_DATA_DIR/ld_blocks_new/

#Copy over the recently created db files to the oracle server.  Note: not automatically copying the gwas_upload.db file

echo "copying db files to oracle server"
rsync -av $RESULTS_DIR/latest/studies.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/studies_new.db
rsync -av $RESULTS_DIR/latest/associations.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/associations_new.db
rsync -av $RESULTS_DIR/latest/ld.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/ld_new.db

touch /tmp/sync_done