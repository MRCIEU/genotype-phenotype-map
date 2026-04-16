#!/bin/bash
set -e

if [ -f .env ]; then
  export $(cat .env | xargs)
fi

# Defaults (env / .env can set these; CLI flags below override)
export PIPELINE_VERSION="${PIPELINE_VERSION:-current}"
export BLOCK_LIST="${BLOCK_LIST:-}"

# Parse CLI: --version / -v, --block-list / --block_list — removed from args passed to snakemake
SNAKEMAKE_ARGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --version|-v)
      if [ -z "${2:-}" ]; then
        echo "ERROR: $1 requires a value" >&2
        exit 1
      fi
      export PIPELINE_VERSION="$2"
      shift 2
      ;;
    --block-list|--block_list)
      if [ -z "${2:-}" ]; then
        echo "ERROR: $1 requires a value" >&2
        exit 1
      fi
      export BLOCK_LIST="$2"
      shift 2
      ;;
    *)
      SNAKEMAKE_ARGS+=("$1")
      shift
      ;;
  esac
done

snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake.log
mkdir -p $(dirname $snakemake_log)

export IMAGE=docker://mrcieu/genotype-phenotype-map:1.0.0
export APPTAINER_VARS="--nv -B /local-scratch -B /projects -B /local-scratch/tmp:/tmp -B /home/$(whoami) -B $(pwd):/home/pipeline --env PIPELINE_VERSION=$PIPELINE_VERSION --env BLOCK_LIST=$BLOCK_LIST --pwd /home/pipeline "

echo "Start time $(date)"
echo "PIPELINE_VERSION=$PIPELINE_VERSION BLOCK_LIST=${BLOCK_LIST:-<empty>}"
apptainer run $APPTAINER_VARS $IMAGE Rscript pipeline_steps/identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log

NUM_STUDIES=$(wc -l < ${DATA_DIR}/pipeline_metadata/studies_to_process.tsv)
if [[ $NUM_STUDIES == 0 ]]; then
  echo 'Nothing to process, exiting.'
  exit 0
elif [[ $NUM_STUDIES -gt 250000 ]]; then
  echo 'ERROR: too many studies to ingest at one time.  This will drastically slow down snakemake, consider splitting studies into smaller chunks'
  exit 0
fi

# Detect "clean" in any remaining snakemake arg (same idea as before when first arg was clean)
CLEAN_RUN=false
for a in "${SNAKEMAKE_ARGS[@]}"; do
  if [[ "$a" =~ clean ]]; then
    CLEAN_RUN=true
    break
  fi
done

if [[ "$CLEAN_RUN" == true ]]; then
  echo "Cleaning up leftover intermediate files from previous run"
  set +e
  rm $DATA_DIR/ld_blocks/*/*/*/*_complete
  set -e
  # Drop args that mention clean from snakemake invocation
  FILTERED=()
  for a in "${SNAKEMAKE_ARGS[@]}"; do
    if [[ ! "$a" =~ clean ]]; then
      FILTERED+=("$a")
    fi
  done
  SNAKEMAKE_ARGS=("${FILTERED[@]}")
fi

apptainer run $APPTAINER_VARS $IMAGE snakemake --profile ./ "${SNAKEMAKE_ARGS[@]}" &>> $snakemake_log

# apptainer run $APPTAINER_VARS $IMAGE Rscript pipeline_steps/post_pipeline_cleanup.R \
#   --current_results_dir $RESULTS_DIR/$PIPELINE_VERSION \
#   --pipeline_summary_file $RESULTS_DIR/$PIPELINE_VERSION/pipeline_summary.html &>> $snakemake_log

# rm $DATA_DIR/pipeline_metadata/studies_to_process.tsv
echo "End time $(date)"
