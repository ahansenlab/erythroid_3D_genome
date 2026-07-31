#!/bin/bash
### ENHANCERS ###

# Phase 1
awk '{OFS=FS="\t"} {print $1, $2, $3}' "peaks/q0.01/P1_H3K27ac_q0.01_peaks.narrowPeak" > "peaks/q0.01/P1_H3K27ac_q0.01_peaks.bed"
awk '{OFS=FS="\t"} {print $1, $2, $3}' "peaks/q0.01/P1_H3K4me1_q0.01_peaks.broadPeak" > "peaks/q0.01/P1_H3K4me1_q0.01_peaks.bed"

bedtools intersect -a "peaks/q0.01/P1_H3K4me1_q0.01_peaks.bed" -b "peaks/q0.01/P1_H3K27ac_q0.01_peaks.bed" > phase_defs/phase1_enhancers.bed

# Phase 2
bedtools intersect -a "peaks/stringent/H3K4me1_consensus_broad.bed" -b "peaks/stringent/H3K27ac_consensus_narrow.bed" > phase_defs/phase2_enhancers.bed

# Phase 3
awk '{OFS=FS="\t"} {print $1, $2, $3}' "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_H3K27ac_Rep1.stringent.bed" > "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_H3K27ac_Rep1.clean.bed"
awk '{OFS=FS="\t"} {print $1, $2, $3}' "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_H3K27ac_Rep2.stringent.bed" > "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_H3K27ac_Rep2.clean.bed"

bedtools intersect -a /mnt/md1/varshini/Blood/p3_cutrun/processed/P3_H3K27ac_Rep1.clean.bed -b /mnt/md1/varshini/Blood/p3_cutrun/processed/P3_H3K27ac_Rep2.clean.bed > phase_defs/phase3_k27ac.bed

cd phase_defs
cat phase1_enhancers.bed phase2_enhancers.bed phase3_k27ac.bed | bedtools sort | bedtools merge -d 500 > merged_enhancers_500bp.bed


#### CTCF ####
# at some point, I need to make the version that keeps the peak stats so that I can filter on that 

# Phase 1
# since there's only one rep for phase 1, use a more stringent fdr-cutoff
awk '{OFS=FS="\t"} {print $1, $2, $3}' "peaks/stringent/A5_CTCF_narrow_peaks.narrowPeak" > "phase_defs/P1_CTCF.bed"

# Phase 2
# lenient-2 had a lenient p-value and relies on the consensus to filter out FPs
cp "peaks/lenient2/CTCF_consensus_narrow.bed" "phase_defs/P2_CTCF.bed"

# Phase 3
awk '{OFS=FS="\t"} {print $1, $2, $3}' "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_CTCF_Rep1.stringent.bed" > "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_CTCF_Rep1.clean.bed"
awk '{OFS=FS="\t"} {print $1, $2, $3}' "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_CTCF_Rep2.stringent.bed" > "/mnt/md1/varshini/Blood/p3_cutrun/processed/P3_CTCF_Rep2.clean.bed"

bedtools intersect -a /mnt/md1/varshini/Blood/p3_cutrun/processed/P3_CTCF_Rep1.clean.bed -b /mnt/md1/varshini/Blood/p3_cutrun/processed/P3_CTCF_Rep2.clean.bed > phase_defs/P3_CTCF.bed

# after running CTCF_processing.sh
cd phase_defs
cat P1_CTCF.bed P2_CTCF.bed P3_CTCF.bed | bedtools sort | bedtools merge > merged_CTCF_P123_noMotifs.bed
cd CTCF
cat P1_CTCF.UNIQUE.peaks.bed P2_CTCF.UNIQUE.peaks.bed P3_CTCF.UNIQUE.peaks.bed | bedtools sort | bedtools merge > merged_CTCF_P123.bed
cat P1_CTCF.UNIQUE.peaks.bed P2_CTCF.UNIQUE.peaks.bed P3_CTCF.UNIQUE.peaks.bed | bedtools sort | bedtools merge -s -c 4,5,6 -o distinct,max,distinct > merged_CTCF_P123_strand.bed
cd /mnt/md1/varshini/Blood/chip_remapping_all
