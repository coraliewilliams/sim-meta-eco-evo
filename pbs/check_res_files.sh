#!/bin/bash

# Set the directory path
dir="/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results/raw"

# Set the range of numbers to check
start=100000
end=160000

# Iterate over the sequence and check for missing files
for ((i=start; i<=end; i++)); do
    filename="res_${i}.RDATA"
    filepath="${dir}/${filename}"
    
    if [ ! -f "$filepath" ]; then
        echo "Missing file: $filename"
    fi
done
