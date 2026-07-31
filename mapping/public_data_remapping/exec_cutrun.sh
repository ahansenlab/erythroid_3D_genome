#!/bin/bash

# done already: prefetch/dump 

# read sample info with sample name and condition in other columns 
while IFS=$'\t' read -r srr sample_name condition; do
    echo "$srr" "$sample_name"
    #./cutrun_processing.sh "$srr" "$sample_name"
    ./cutrun_seact.sh "$sample_name" 0.01
done < srr_accessions.tsv


