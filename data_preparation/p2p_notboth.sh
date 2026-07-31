
#!/bin/bash

f1="$1"
f2="$2"
f3="$3"
odir="$4"
tag="$5"


# generate basenames for all 3 files 
bn1=$(basename "$f1" ".bedpe")
bn2=$(basename "$f2" ".bedpe")
bn3=$(basename "$f3" ".bedpe")

# compare f2 against f1 
pairToPair -a "${odir}/$f2" -b "${odir}/$f1" -type "notboth" -slop 2500 > "f2_not_f1.bedpe"
pairToPair -a "${odir}/$f3" -b "${odir}/$f1" -type "notboth" -slop 2500 > "f3_not_f1.bedpe"
pairToPair -a "${odir}/$f3" -b "${odir}/$f2" -type "notboth" -slop 2500 > "f3_not_f2.bedpe"
pairToPair -a "${odir}/f3_not_f2.bedpe" -b "${odir}/f3_not_f1.bedpe" -type "both" > "f3_only.bedpe"


cat "$f1" "f2_not_f1.bedpe" "f3_only.bedpe" > merged_loops.bedpe
