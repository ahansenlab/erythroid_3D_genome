#!/bin/bash 
outdir="$1"
tag="$2"

annotation="/mnt/md0/varshini/Gene_Annotations/gencode.v45.annotation.gtf"

mapfile -d $'\0' dedup < <(find aligned/ -name *sorted.rmdup.bam -print0)
mapfile -d $'\0' sorted < <(find aligned/ -name *sorted.bam -print0)

featureCounts -T 8 -p -s 2 -t exon -a "$annotation" -o "${outdir}/rmdup_counts_${tag}.txt" "${dedup[@]}"
featureCounts -T 8 -p -s 2 -t exon -a "$annotation" -o "${outdir}/counts_${tag}.txt" "${sorted[@]}"
