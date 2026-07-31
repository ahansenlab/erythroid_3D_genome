#!/bin/bash
control="$1"
test="$2"
outdir="${3:-/mnt/md0/varshini/Analysis/Blood/complete_analyses/boundaries/motif_enrichment}"

# create the fasta files 
fasta_filename="/mnt/coldstorage/Varshmallow/hg38_analysisSet/hg38.analysisSet.fa"
atac_filename="/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs/erythroid_tracks/032524_combined_DAP_JAN_028_ATAC.bed"

bedtools intersect -a "$atac_filename" -b "${control}.bed" -wa > "${control}_atac.bed"
bedtools intersect -a "$atac_filename" -b "${test}.bed" -wa > "${test}_atac.bed"

homer_dir="/mnt/md0/varshini/Analysis/Motifs/homer_summits/"
mkdir -p "$outdir"

btest=$(basename "$test")
bcontrol=$(basename "$control")

cd "/mnt/md0/varshini/Analysis/Motifs/homer/bin/"
pwd
findMotifsGenome.pl "${test}_atac.bed"  hg38 "${outdir}/motif_output_anchors_atac_${btest}_over_${bcontrol}/homer_summits" -bg "${control}_atac.bed" -size 100 -p 12

#./run_homer.sh "${control}_atac.fa" "${test}_atac.fa" "${outdir}/motif_output_anchors_atac_${btest}_over_${bcontrol}/homer"
#./run_ame_merged.sh "${control}.fa" "${test}.fa" "${outdir}/motif_output_anchors_${btest}_over_${bcontrol}/homer" "$database"


