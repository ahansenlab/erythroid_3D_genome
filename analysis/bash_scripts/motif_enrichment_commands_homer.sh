#!/bin/bash

# test files 
basedir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/boundaries"
windows=("20kb" "50kb")

# databases 
for window in "${windows[@]}"; do
lost="lost_p12_boundary_2kb_${window}"
gained="gained_p13_boundary_2kb_${window}"
control="all_bounds_2kb_${window}"


echo "$database"
echo "$control"
echo "$lost"

#./prepare_homer_files.sh "$basedir/$control" "$basedir/$lost" "$basedir"
#./prepare_homer_files.sh "$basedir/$control" "$basedir/$gained" "$basedir"

./prepare_homer_ame_files_summits.sh "$basedir/$control" "$basedir/$lost" "$basedir"
./prepare_homer_ame_files_summits.sh "$basedir/$control" "$basedir/$gained" "$basedir" 

done
