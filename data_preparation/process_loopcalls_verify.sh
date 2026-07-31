#!/bin/bash

qval=0.01
anch=2500
anch_name=$(( anch * 2 / 1000 ))

# filter by FDR and format as bedpe 
./bedpe_qval.sh /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_5kb_st0.88.tsv "$qval"
./bedpe_qval.sh /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_2kb_st0.88.tsv "$qval"
./bedpe_qval.sh /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_1kb_st0.92.tsv "$qval"

# change binsizes to 5kb
./format_bedpe.sh /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_2kb_st0.88_q${qval}.bedpe "$anch" > /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_2kb_st0.88_q${qval}_${anch_name}kb.bedpe
./format_bedpe.sh /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_1kb_st0.92_q${qval}.bedpe "$anch"  > /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify/merged_calls_1kb_st0.92_q${qval}_${anch_name}kb.bedpe 

# get rid of duplicates
./exec_p2p_v2.sh merged_calls_1kb_st0.92_q${qval}_${anch_name}kb.bedpe merged_calls_2kb_st0.88_q${qval}_${anch_name}kb.bedpe merged_calls_5kb_st0.88_q${qval}.bedpe /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify 

# concatenate 
cd /mnt/md0/varshini/Analysis/Blood/Loop_Calls/mustache_calls/verify 
cat merged_calls_1kb_st0.92_q${qval}_${anch_name}kb.bedpe merged_calls_2kb_st0.88_q${qval}_${anch_name}kb_clipped.bedpe merged_calls_5kb_st0.88_q${qval}_clipped.bedpe > merged_calls_125kb_q0.01_verify.bedpe
