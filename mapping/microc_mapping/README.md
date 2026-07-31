# Alignment docs 

- The processing steps are a bit involved since there were two sequencing runs, and only some samples were included on both runs. 
- These steps briefly describe how I processed the data and point to the scripts used for each.

- Note: This snakefile is mostly for documentation purposes. It includes some cleanup functions
of the .pairs files that could be avoided with some different flags (originally I was including the sam string which makes the pairs files really large). 

If you would like to re-process this data, I would recommend running the pipeline with snakefile_clean on both sequencing runs, then following the sampling and merging steps as described below. Then you wouldn't have to drop those extra columns in a convoluted way and will save several days of processing time. Or even better, we have since uploaded a more user-friendly snakemake pipeline for Micro-C and RCMC processing (https://github.com/ahansenlab/MicroC_RCMC_analysis/tree/main).

# ORIGINAL RUN:
/mnt/md0/varshini/Analysis/genomics_general/

1. ./generate_microc_samples.sh FASTQ_PATH OUTPUT_PATH
2. snakefile: snakefile_dropsam 

 - dropsam scripts: dropsam_2.sh and dropsam_4.sh

4. sample_sankaran.sh for P1, Rep2

# RESEQ RUN:
/mnt/md0/varshini/Analysis/reseq

1. ./generate_microc_samples.sh FASTQ_PATH OUTPUT_PATH
2. snakefile: snakefile

# COMBINING SEQRUNS
3. use pairtools merge and dedup on Phase1_Rep1, Phase2_Rep2, Phase3_Rep2, Phase1_Rep1 (merge_sankaran.sh)
4. phase23_merge.sh 
5. phase1_merge.sh to merge the sampled Phase1_Rep2 with the merged Phase1_Rep1


