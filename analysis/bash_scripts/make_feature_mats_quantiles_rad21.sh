
#!/bin/bash

dist="${1:-100kb}"
feature_file_dir="${2:-/mnt/md1/varshini/Blood/georgiades_chip/}" # read in all files "
outdir="${3:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/metaplots/}"

track_file_dir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/match_quantiles/quantiles_${dist}_clus/"

phases=("1" "2" "3")
reps=("rep1" "rep2" "rep3")
for i in "${phases[@]}"; do 

for rep in "${reps[@]}"; do 
mapfile -d $'\0' feature_files < <(find "$feature_file_dir" -type f -name "*Don*0${i}*${rep}*.bw" -print0 -maxdepth 1)
mapfile -d $'\0' track_files < <(find "$track_file_dir" -type f -name "*clus.bed" -print0)
#track_files+=("/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs/erythroid_tracks/ATAC_not_CTCF.bed")

echo "$feature_file_dir"
echo "${feature_files[@]}"
echo "${track_files[@]}"

computeMatrix reference-point -S "${feature_files[@]}" -R "${track_files[@]}" -o "${outdir}quant_${dist}_rad21_Don${i}_reps_mean_${rep}.gz" --upstream 200000 --downstream 200000 --binSize 2000 -p 8 --averageTypeBins "mean" --verbose
plotHeatmap -m "${outdir}quant_${dist}_rad21_Don${i}_reps_mean_${rep}.gz" -o "${outdir}quant_${dist}_rad21_Don${i}_reps_mean_${rep}.svg" --colorMap "Blues"

#computeMatrix reference-point -S "${feature_files[@]}" -R "${track_files[@]}" -o "${outdir}quant_${dist}_rad21_small_allReps.gz" --upstream 20000 --downstream 20000 --binSize 200 --averageTypeBins "sum" --verbose -p 8
#plotHeatmap -m "${outdir}quant_${dist}_rad21_small_allReps.gz" -o "${outdir}quant_${dist}_rad21_small_allReps.svg" --colorMap "Blues"

done
done
