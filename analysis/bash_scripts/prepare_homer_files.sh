#!/bin/bash
control="$1"
test="$2"
outdir="${3:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/boundaries/motif_enrichment}"

# create the fasta files 
fasta_filename="/mnt/coldstorage/Varshmallow/hg38_analysisSet/hg38.analysisSet.fa"
atac_filename="/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs/erythroid_tracks/032524_combined_DAP_JAN_028_ATAC.bed"

bedtools getfasta -fi "$fasta_filename" -bed "${control}.bed" > "${control}.fa"
bedtools getfasta -fi "$fasta_filename" -bed "${test}.bed" > "${test}.fa" 

# this will report multiple atac peaks that overlap a single boundary
bedtools intersect -a "$atac_filename" -b "${control}.bed" -wa > "${control}_atac.bed"
bedtools intersect -a "$atac_filename" -b "${test}.bed" -wa > "${test}_atac.bed"
 
bedtools getfasta -fi "$fasta_filename" -bed "${control}_atac.bed" > "${control}_atac.fa"
bedtools getfasta -fi "$fasta_filename" -bed "${test}_atac.bed" > "${test}_atac.fa"

btest=$(basename "$test")
bcontrol=$(basename "$control")

./run_homer.sh "${control}_atac.fa" "${test}_atac.fa" "${outdir}/motif_output_anchors_atac_${btest}_over_${bcontrol}/homer"
#./run_ame_merged.sh "${control}.fa" "${test}.fa" "${outdir}/motif_output_anchors_${btest}_over_${bcontrol}/homer" "$database"

