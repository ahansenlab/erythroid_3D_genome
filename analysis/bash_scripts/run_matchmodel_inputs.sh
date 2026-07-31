#!/bin/bash 

script_dir="/mnt/md0/varshini/Analysis/Blood/scripts"

match_track_dir="/mnt/md0/varshini/Analysis/Blood/complete_analyses/track_defs/erythroid_tracks/classifications"
track_names="erythroid_tf_acrs_all_enhancernotpromoternotctcf.bed"

cd "$script_dir"

./matchmodel_input_tables.sh "complete_v2" "$match_track_dir" "$track_names"


