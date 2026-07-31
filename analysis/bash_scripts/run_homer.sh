#!/bin/bash

ctrl="$1"
test="$2"
outdir="$3"

homer_dir="/mnt/md0/varshini/Analysis/Motifs/homer/"

mkdir -p "$outdir"

/mnt/md0/varshini/Analysis/Motifs/homer/bin/findMotifs.pl "$test" fasta "$outdir" -fastaBg "$ctrl"

