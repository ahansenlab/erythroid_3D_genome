#!/bin/bash
treat="$1"
thresh="$2"
genome="/mnt/md0/DataRepository/chromsizes/hg38/hg38.sorted.chrom.sizes"

# call peaks on each condition
bedtools bamtobed -bedpe -i bam/$treat.rmdup.sorted.bam > processed/$treat.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' processed/$treat.bed > processed/$treat.clean.bed
cut -f 1,2,6 processed/$treat.clean.bed | sort -k1,1 -k2,2n -k3,3n > processed/$treat.fragments.bed
bedtools genomecov -bg -i processed/$treat.fragments.bed -g "$genome" > processed/$treat.fragments.bedgraph

bash /home/varshini/anaconda3/envs/seacr/bin/SEACR_1.3.sh "processed/${treat}.fragments.bedgraph" "$thresh" norm stringent "processed/${treat}"
