#rsync -avz $DATA_DIR/ld_blocks/ opc@141.147.82.95:/oradiskvdb1/data/ld_blocks/  --filter='+ */'  --filter='+ *finemapped_studies.tsv' --filter='- *'
#rsync -avz $DATA_DIR/study/ opc@141.147.82.95:/oradiskvdb1/data/study/ --filter='+ finemapped/*'

files_to_transfer=$(cat /local-scratch/projects/genotype-phenotype-map/data/ld_blocks/EUR/1/101384499-103762931/finemapped_studies.tsv | awk '{print $3}' | grep -v file)
for f in $files_to_transfer; do
  new_file=$(echo "$f" | sed 's#/local-scratch/projects/genotype-phenotype-map/#/oradiskvdb1/#');
  echo $new_file;
  new_dir=$(dirname $new_file)
  echo $new_dir
  ssh -i ~/.ssh/oracle_gpm_vm.key opc@141.147.82.95 "mkdir -p $new_dir"
  scp -i ~/.ssh/oracle_gpm_vm.key $f opc@141.147.82.95:$new_file;
done
