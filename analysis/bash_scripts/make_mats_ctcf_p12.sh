#!/bin/bash

dist="${1:-100kb}"
feature_file_dir="${2:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/model_inputs/}" # read in all files "
outdir="${3:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/metaplots/}"

# do just CTCF since its depleted
phases=(1 2 3)
for i in "${phases[@]}"; do

mapfile -d $'\0' feature_files < <(
    find "$feature_file_dir" -maxdepth 1 -type f -name "P${i}_CTCF.bw" -print0
)
track_file_dir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/match_quantiles/quantiles_${dist}_clus/"

mapfile -d $'\0' track_files < <(find "$track_file_dir" -type f -name "*.bed" -print0)

computeMatrix reference-point -S "${feature_files[@]}" -R "${track_files[@]}" -o "${outdir}quant_${dist}_v2_all_2kb_CTCF${i}.gz" --upstream 200_000 --downstream 200_000 --binSize 2000 -p 8 --averageTypeBins "mean" --verbose 
plotHeatmap -m "${outdir}quant_${dist}_v2_all_2kb_CTCF${i}.gz" -o "${outdir}quant_${dist}_all_2kb_CTCF${i}.svg" --perGroup --colorMap Blues


new_file="/mnt/md1/varshini/Blood/chip_remapping_all/phase_defs_1e-2/P${i}_CTCF.bed"
track_files+=("$new_file")
computeMatrix reference-point -S "${feature_files[@]}" -R "${track_files[@]}" -o "${outdir}quant_${dist}_v2_all_2kb_CTCF${i}_scale.gz" --upstream 200_000 --downstream 200_000 --binSize 2000 -p 8 --averageTypeBins "mean" --verbose 
plotHeatmap -m "${outdir}quant_${dist}_v2_all_2kb_CTCF${i}_scale.gz" -o "${outdir}quant_${dist}_all_2kb_CTCF${i}_scale.svg" --perGroup --colorMap Blues

done
