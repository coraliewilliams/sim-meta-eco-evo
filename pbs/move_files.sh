#!/bin/bash

# Define the source directory
source_dir="/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results/raw"

# Define the destination directories
dest_dir1="/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results/set1"
dest_dir2="/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results/set2"
dest_dir3="/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results/set3"

# Move files to set1 (res_1.RDATA to res_100000.RDATA)
for i in $(seq 1 100000); do
  mv "$source_dir/res_${i}.RDATA" "$dest_dir1/"
done

# Move files to set2 (res_100001.RDATA to res_200000.RDATA)
for i in $(seq 100001 200000); do
  mv "$source_dir/res_${i}.RDATA" "$dest_dir2/"
done

# Move files to set3 (res_200001.RDATA to res_360000.RDATA)
for i in $(seq 200001 360000); do
  mv "$source_dir/res_${i}.RDATA" "$dest_dir3/"
done

echo "Files have been moved successfully!"
