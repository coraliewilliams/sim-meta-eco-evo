#!/bin/bash

echo "Enter missed job IDs (one per line). Type ';' on a new line to finish:"
missed_jobs=()
while IFS= read -r line; do
  [[ "$line" == ";" ]] && break
  [[ -n "$line" ]] && missed_jobs+=("$line")
done

read -p "Enter full path to the PBS script: " pbs_script
[[ ! -f "$pbs_script" ]] && echo "PBS script not found." && exit 1

for job_id in "${missed_jobs[@]}"; do
  echo "Submitting job for index $job_id"
  qsub -v PBS_ARRAY_INDEX=$job_id "$pbs_script"
done
