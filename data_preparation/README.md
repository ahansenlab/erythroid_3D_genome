## Data Preparation

This folder contains scripts related to preparation of the main datasets used in this manuscript.

### related to loop preparation. 
1. Call loops with Mustache and q-value filter (process_loopcalls_verify.sh and bedpe_qval.sh)
2. Deduplicate loops with priority (p2p_notboth.sh)
3. Classify loops (classify_blood_loops_final.R) which runs loop_anchor_classification.R


### related to gene expression datasets. 
1. get_tpms_blood.R approximates TPMs from counts given in Ludwig and Lareau et al., Cell Rep (2019). 
2. perturb_multiome_bulkExp_markerGenes.R gets marker genes and differentially expressed genes for populations grouped according to the P1-P3 from this study, for the non-targeting and control cells from the single-cell experiments in Martin-Rufino and Caulier et al., Science (2025).
