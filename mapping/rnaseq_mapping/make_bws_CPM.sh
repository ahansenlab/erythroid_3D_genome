#!/bin/bash

mapfile -d $'\0' dedup < <(find aligned/ -name *sorted.rmdup.bam -print0)

for bamfile in "${dedup[@]}"; do
sample_name=$(echo "$bamfile" | cut -d'.' -f1)
TMPDIR="/mnt/md1/varshini/Blood" bamCoverage -b "$bamfile" -o "${sample_name}_rmdup_CPM.bw" --normalizeUsing "CPM" 
echo "$sample_name done"

done

