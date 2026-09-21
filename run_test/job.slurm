#!/bin/bash
#SBATCH -A stf218
#SBATCH -J alloy_mc_smoke
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:10:00

#SBATCH -o alloy_mc_smoke-%j.out
#SBATCH -e alloy_mc_smoke-%j.err

# --- environment -----------------------------------------------------------
module load PrgEnv-cray cce cray-mpich

set -euo pipefail
cd /lustre/orion/stf218/scratch/junqi/alloy-thermo-skill/run_test

# Number of MPI ranks == number of parallel-tempering temperature replicas.
# The engine reads composition.input and coupling.input from this directory.
echo "Starting alloy-mc: 8 replicas on 1 node(s)"
srun -n 8 /lustre/orion/stf218/scratch/junqi/alloy-thermo-skill/engine/alloy_mc control.input 12345
echo "alloy-mc finished with exit code $?"
