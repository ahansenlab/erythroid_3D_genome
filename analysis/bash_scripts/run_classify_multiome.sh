basedir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs/erythroid_tracks"
trackdir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs"
script_dir="/mnt/md0/varshini/Analysis/Blood/scripts"

files=("032524_combined_DAP_JAN_028_ATAC.bed" "tf_acrs_sorted.bed" "erythroid_associated_ACRs.bed" "erythroid_tf_acrs_all.bed")
files=("erythroid_tf_acrs_all.bed" "tf_acrs_sorted.bed")

cd "$script_dir"
for file in "${files[@]}"; do
   ./classify_multiome_acrs_uniqueE.sh -b "$basedir/$file" -t "$trackdir" -e merged_enhancers_500bp.bed -c merged_CTCF_P123.bed -p /mnt/md0/varshini/Gene_Annotations/hg38.refGene.tss_strand.bed -o "${basedir}/classifications_uniqueE" 
done 

for file in "${files[@]}"; do
   ./classify_multiome_acrs.sh -b "$basedir/$file" -t "$trackdir" -e merged_enhancers_500bp.bed -c merged_CTCF_P123.bed -p /mnt/md0/varshini/Gene_Annotations/hg38.refGene.tss_strand.bed -o "${basedir}/classifications" 
done 

