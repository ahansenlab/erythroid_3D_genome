#!/bin/bash

base_path=$1
output_file=$2
ext=${3:-'.fastq.gz'}
# Write header to the output file
echo -e "sample_id\tsample\trep\tlane\tfastq_r1\tfastq_r2" > temp_samples.txt

# Read all files matching the pattern
for file in "$base_path"/*${ext}; do
    # Extract the filename without path
    filename=$(basename "$file")

    echo $filename
    sample_2=$(echo "$filename" | awk -F'_' '{sub(/_S[0-9]+.*/, ""); print substr($0, index($0,$4))}')
    echo "s2:" $sample_2

    # Extract relevant parts
    read lane <<< $(echo "$filename" | awk -F'_S[0-9]+' '{split($2, arr, "_"); lane=arr[2]; print lane}')
    #read rep <<< $(echo "$filename" | awk -F'_S[0-9]+' '{split($1, arr, "_"); print arr[length(arr)]}')
    IFS='_' read proj day rep <<< "$sample_2"
    echo "fields:" $lane, $proj, $day, $rep

    # Determine R1 and R2 filenames
    if [[ "$filename" == *_R1_* ]]; then
        fastq_r1="${base_path}/${filename}"
        fastq_r2="${base_path}/${filename/_R1_/_R2_}"
        if [[ ! -f "$fastq_r2" ]]; then
            fastq_r2=""
        fi
    elif [[ "$filename" == *_R2_* ]]; then
        fastq_r2="${base_path}/${filename}"
        fastq_r1="${base_path}/${filename/_R2_/_R1_}"
        if [[ ! -f "$fastq_r1" ]]; then
            fastq_r1=""
        fi
    fi

    # Write to output
    if [[ "$sample_2" != "Undetermined"* ]]; then
        echo -e "${sample_2}_${lane}\t${proj}_${day}\t$rep\t$lane\t$fastq_r1\t$fastq_r2" | uniq >> temp_samples.txt
    fi
done
uniq temp_samples.txt > "$output_file"
rm temp_samples.txt
echo "Metadata extraction complete. Output written to $output_file."
