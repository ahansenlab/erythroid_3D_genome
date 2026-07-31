#!/bin/bash
input_file="$1"
output_file="$2"

zcat "$input_file" | awk -v last="$(zgrep '^#' "$input_file" | tail -n 1)" 'BEGIN {OFS=FS} {if ($0 == last) {NF-=4} print}' | pbgzip -c > "$output_file"
