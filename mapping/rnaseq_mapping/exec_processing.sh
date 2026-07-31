#!/bin/bash
samples="$1"

basedir="/mnt/md1/varshini/Blood/rnaseq_v2" # assume fq is in fq_dir/sample_name/sample_name_1/2.fq.gz 
fq_dir="${basedir}/01.RawData"

align_dir="${basedir}/aligned"
counts_dir="feature_counts"

if [ ! -d "$align_dir" ]; then
mkdir "$align_dir"
fi

if [ ! -d "$counts_dir" ]; then 
mkdir "$counts_dir"
fi

mapfile -t fq_names < <(cut -f1 "$samples")
mapfile -t sample_names < <(cut -f2 "$samples")

# 2. Loop through the arrays simultaneously using the index
#for i in "${!fq_names[@]}"; do
#    
#    # Access the corresponding values from both arrays using the index
#    fq_name="${fq_names[$i]}"
#    sample_name="${sample_names[$i]}"
#    
#    echo "  og name: $fq_name"
#    echo "  saving to: $sample_name"

#   # alignment command 
#   ./align_sort.sh -f "$fq_dir" -n "$fq_name" -s "$sample_name" -o "$align_dir" 
#    
#done

./featurecount.sh "$counts_dir" "v2"
