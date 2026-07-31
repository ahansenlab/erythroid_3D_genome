#!/bin/bash

peaks_file="$1"
bw_file="$2"
ostr="$3"
outdir="$4"

bn="${peaks_file%.*}"
bnn=$(basename "$bn")

sort_options=("mean" "sum")
for sort_option in "${sort_options[@]}"; do
computeMatrix scale-regions -S "$bw_file" -R "$peaks_file"\
               --sortRegions descend --sortUsing "$sort_option" \
               --outFileSortedRegions "${outdir}/${bnn}_sortRegions_${sort_option}_${ostr}.bed"\
               --outFileNameMatrix "${outdir}/${bnn}_${sort_option}_${ostr}.tab" -o "${outdir}/${bnn}_${sort_option}_${ostr}.gz" \
               --binSize 10 -p 8 --averageTypeBins "mean" \
               --verbose

zcat "${outdir}/${bnn}_${sort_option}_${ostr}.gz" \
| awk -v sort="${sort_option}" '
BEGIN{OFS="\t"}
!/^@/ {
    n=0; s=0
    for(i=7;i<=NF;i++){
        if($i != "nan"){
            s += $i
            n++
        }
    }

    score = (sort=="mean" ? (n ? s/n : "NA") : s)

    print $1, $2, $3, score
}' > "${outdir}/${bnn}_${sort_option}_${ostr}_score.bed"
done




