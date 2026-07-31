trackdir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs"
script_dir="/mnt/md0/varshini/Analysis/Blood/scripts"

cd "$script_dir"

./classify_tracks.sh -t "$trackdir" -e merged_enhancers_500bp.bed -c merged_CTCF_P123.bed -p /mnt/md0/varshini/Gene_Annotations/hg38.refGene.tss_strand.bed -o "$trackdir/track_classifications"
