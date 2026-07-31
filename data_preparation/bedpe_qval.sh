#!/bin/bash
name="$1"
qval="$2"

tail -n +2 "$name" | awk -v q="$qval" '{OFS="\t"} ($7 < q) {print $1, $2, $3, $4, $5, $6, $7}' > "${name%.*}_q${qval}.bedpe" 
