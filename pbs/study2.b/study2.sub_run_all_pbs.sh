#!/usr/bin/env bash
# submit_all_pbs.sh
# Generate and qsub 36 PBS array jobs to cover 1–360,000 in chunks of 10,000

# (1) adjust this to your actual working directory:
WORKDIR="/srv/scratch/z5394590/phylo_meta_sandwich/main/study2.sub"

# (2) create a place to store all the little PBS scripts
PBS_DIR="${WORKDIR}/pbs_scripts"
mkdir -p "${PBS_DIR}"

# (3) for each of the 36 blocks:
for (( chunk=15; chunk<16; chunk++ )); do
  start=$(( chunk*10000 + 1 ))
  end=$(( start + 10000 - 1 ))
  jobname="sim_meta_${start}_${end}"
  script="${PBS_DIR}/${jobname}.pbs"

  cat > "${script}" <<EOF
#!/bin/bash
#PBS -N ${jobname}
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:08:00
#PBS -J ${start}-${end}

cd ${WORKDIR}

module purge
module add r/4.3.1

Rscript sim_meta_study2.sub.R
EOF

  # submit it
  qsub "${script}"
done
