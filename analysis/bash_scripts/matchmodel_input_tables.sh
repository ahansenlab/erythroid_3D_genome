#!/bin/bash

ostr="$1"
match_dir="$2"
match_inputs=("${@:3}")

#track_input_files="/mnt/md0/varshini/Analysis/Blood/matchmakers/model_inputs/"
track_input_files="/mnt/md0/varshini/Analysis/Blood/complete_analyses/bigwigs/"
outdir="/mnt/md0/varshini/Analysis/Blood/matchmakers/model_input_tables"

mapfile -d $'\0' feature_files < <(find "$track_input_files" -maxdepth 1 -type f -name "*.bw" -print0)

for matchmaker_input in "${match_inputs[@]}"; do

  bn="${matchmaker_input%.*}"
  echo "${bn}_${ostr}.tab"
  computeMatrix reference-point -S "${feature_files[@]}" -R "${match_dir}/$matchmaker_input"\
               --outFileNameMatrix "${outdir}/${bn}_${ostr}.tab" -o "${outdir}/${bn}_${ostr}.gz" \
               --upstream 200_000 --downstream 200_000 --binSize 2000 -p 8 --averageTypeBins "mean" \
               --verbose --referencePoint center

done
