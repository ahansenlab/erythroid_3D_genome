# erythroid_3D_genome
Code for performing analyses and simulations for the publication "Cohesin loading at regulatory elements shapes 3D genome folding during erythropoiesis" (BioRxiv, Aug 2026)

## This repo contains the following:
1. mapping/ contains code related to processing sequencing data (Micro-C, RNA-Seq, public datasets that were re-analyzed for this study).
2. data_preparation/ contains code related to making the erythroid epigenomic annotations and 1D tracks used in this study.
3. simulations/ contains polymer simulation code to model the effect of preferential cohesin loading and chromatin compaction, and some processing code to make maps and loop strength quantifications from simulated polymers.
4. figure_generation/ contains Jupyter notebooks split by main figure and associated supplementary figures. Each notebook yields at least one of the figures used in the manuscript exactly, barring aesthetic modifications.

Please see the indicated subfolders for further documentation of each of these specific portions. Each folder also contains the necessary environment(s); note that environments cannot be shared between each of these. For example, the mapping environment is incompatible with the simulation code. 
