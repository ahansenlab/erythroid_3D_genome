#!/bin/bash

tag="$1"
loops="${2:-merged_calls_125kb_q0.01.bedpe}"
loop_dir="${3:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/loops_anchors}"
track_dir="${4:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs/erythroid_tracks}"
#outdir="${5:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/loops_anchors}"

# ATAC, Tf-sensitive, erythroid-sensitive, both
files_to_subset=("032524_combined_DAP_JAN_028_ATAC.bed" "acrs_unique.bed" "erythroid_associated_ACRs.bed" "erythroid_tf_acrs_unique.bed")

for file in "${files_to_subset[@]}"; do
   # make either and both pairtobeds
   bn="${file%.*}"
   pairToBed -a "${loop_dir}/${loops}" -b "${track_dir}/${file}" | awk '{OFS=FS="\t"} {print $1, $2, $3, $4, $5, $6}' | sort -k1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq > "${loop_dir}/${bn}_${tag}.bedpe"
   pairToBed -a "${loop_dir}/${loops}" -b "${track_dir}/${file}" -type both | awk '{OFS=FS="\t"} {print $1, $2, $3, $4, $5, $6}' | sort -k1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq > "${loop_dir}/${bn}_${tag}_both.bedpe"

done
