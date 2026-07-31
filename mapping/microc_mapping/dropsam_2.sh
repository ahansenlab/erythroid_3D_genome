#!/bin/bash

input_file="$1"
output_file="${input_file%.gz}.dropped.gz"

zcat "$input_file" | awk '
    BEGIN { in_comments = 1; }
    /^#/ {
        if (in_comments) {
            print;
        } else {
            split($0, a, FS, seps);
            for (i = 1; i <= NF - 4; i++) {
                printf "%s%s", a[i], seps[i];
            }
            print "";
        }
    }
    !/^#/ {
        in_comments = 0;
        split($0, a, FS, seps);
        for (i = 1; i <= NF - 4; i++) {
            printf "%s%s", a[i], seps[i];
        }
        print "";
    }
' | gzip > "$output_file"
