#!/bin/bash
set -e

SVG_DIR=$DATA_DIR/svgs
RSYNC_DIR=$DATA_DIR/rsync_to_server
mkdir -p $RSYNC_DIR

REMOTE_SERVER_DATA_DIR=/oradiskvdb1/data
REMOTE_SERVER_STUDY_DIR=$REMOTE_SERVER_DATA_DIR/study
REMOTE_SERVER_SVG_DIR=/oradiskvdb1/static/svgs_new
REMOTE_SERVER_DB_DIR=/oradiskvdb1/db
STUDY_DIR=${DATA_DIR}study

cd $STUDY_DIR

#Sync the finemapped GWAS results to the oracle server
# rsync -aRv --prune-empty-dirs --include='*/' --include='*finemapped/*_with_lbfs.tsv.gz' --exclude='*' . $ORACLE_SERVER:$REMOTE_SERVER_STUDY_DIR

#Sync the traits and groups svg files to the oracle server
rsync -av $SVG_DIR/traits/ $ORACLE_SERVER:$REMOTE_SERVER_SVG_DIR/traits
rsync -av $SVG_DIR/groups/ $ORACLE_SERVER:$REMOTE_SERVER_SVG_DIR/groups

#Alter the ld block information accordingly and rsync
echo "preparing ld_blocks directory"
mkdir -p $RSYNC_DIR/ld_blocks

cd $DATA_DIR
find ld_blocks -type f -name "finemapped_studies.tsv" -exec cp --parents {} $RSYNC_DIR \;
find ld_blocks -type f -name "coloc_pairwise_results.tsv.gz" -exec cp --parents {} $RSYNC_DIR \;
find $RSYNC_DIR/ld_blocks -type f -name "finemapped_studies.tsv" -exec sed -i "s|${STUDY_DIR}|$REMOTE_SERVER_STUDY_DIR|g;" {} \;
find $RSYNC_DIR/ld_blocks -type f -name "coloc_pairwise_results.tsv.gz" -exec sed -i "s|${STUDY_DIR}|$REMOTE_SERVER_STUDY_DIR|g;" {} \;

cd $RSYNC_DIR/ld_blocks
echo "rsyncing ld_blocks to oracle server"
rsync -aRv . $ORACLE_SERVER:$REMOTE_SERVER_DATA_DIR/ld_blocks_new/

#Copy over the recently created db files to the oracle server.  Note: not automatically copying the gwas_upload.db file
echo "copying db files to oracle server"
rsync -av $RESULTS_DIR/latest/studies.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/studies_new.db
rsync -av $RESULTS_DIR/latest/coloc_pairs.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/coloc_pairs_new.db
rsync -av $RESULTS_DIR/latest/associations.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/associations_new.db
rsync -av $RESULTS_DIR/latest/ld.db $ORACLE_SERVER:$REMOTE_SERVER_DB_DIR/ld_new.db
