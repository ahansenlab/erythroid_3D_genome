#!/bin/bash

# test files 
basedir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/boundaries"
windows=("20kb" "50kb")

# databases 
databases=("/mnt/md0/varshini/Analysis/Motifs/H14CORE_meme_format.meme" "/mnt/md0/varshini/Analysis/Motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme")

for window in "${windows[@]}"; do

lost="lost_p12_boundary_2kb_${window}"
gained="gained_p13_boundary_2kb_${window}"
control="all_bounds_2kb_${window}"

for database in "${databases[@]}"; do

echo "$database"
echo "$control"
echo "$lost"

./prepare_ame_files.sh "$basedir/$control" "$basedir/$lost" "$database" "$basedir"
./prepare_ame_files.sh "$basedir/$control" "$basedir/$gained" "$database" "$basedir"


done

done
