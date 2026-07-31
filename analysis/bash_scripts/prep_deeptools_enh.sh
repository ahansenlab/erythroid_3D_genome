#!/bin/bash 

peaks_dir="$1"
peaks_fname="$2"
bigwig_fname="$3"
tag="$4"

basedir="/mnt/md1/varshini/Blood/chip_remapping_all/peaks"
bigwig_dir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/bigwigs"
output_dir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/matchmakers/genomic_features"

input_peaks="${basedir}/${peaks_dir}/${peaks_fname}"
input_bigwig="${bigwig_dir}/${bigwig_fname}"

mkdir -p "$output_dir"
./run_deeptools_enh.sh "$input_peaks" "$input_bigwig" "$tag" "$output_dir"

