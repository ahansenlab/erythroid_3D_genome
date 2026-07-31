#!/bin/bash 
ctrl="$1"
test="$2"
outdir="$3"
database="${4:-/mnt/md0/varshini/Analysis/Motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme}"

databasename=$(basename "$database")

mkdir -p "${outdir}/${databasename}"
echo "outputting to ${outdir}/${databasename}"

if [[ -d $database ]]; then
    files=( $(find "$database" -maxdepth 1 -type f) )
    ame --oc "${outdir}/${databasename}" --control "$ctrl" "$test" "${files[@]}"
elif [[ -f $database ]]; then
    ame --oc "${outdir}/${databasename}" --control "$ctrl" "$test" "$database"
else
    echo "$database is not valid"
    exit 1
fi



