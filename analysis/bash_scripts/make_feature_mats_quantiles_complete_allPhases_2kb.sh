#!/bin/bash

dist="${1:-100kb}"
feature_file_dir="${2:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/model_inputs/}" # read in all files "
outdir="${3:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/metaplots/}"

phases=(1 2 3)
for i in "${phases[@]}"; do
mapfile -d $'\0' feature_files < <(
    find "$feature_file_dir" -maxdepth 1 -type f -name "P${i}*.bw" -print0
)

track_file_dir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/match_quantiles/quantiles_${dist}_clus/"

mapfile -d $'\0' track_files < <(find "$track_file_dir" -type f -name "*.bed" -print0)

echo "$feature_file_dir"
echo "${feature_files[@]}"
echo "$track_file_dir"

computeMatrix reference-point -S "${feature_files[@]}" -R "${track_files[@]}" -o "${outdir}quant_${dist}_v2_all_2kb_P${i}.gz" --upstream 200_000 --downstream 200_000 --binSize 2_000 -p 8 --averageTypeBins "mean" --verbose 
plotHeatmap -m "${outdir}quant_${dist}_v2_all_2kb_P${i}.gz" -o "${outdir}quant_${dist}_all_2kb_P${i}.svg" --colorMap Blues
done

