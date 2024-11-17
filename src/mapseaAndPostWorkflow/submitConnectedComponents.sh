#!/bin/bash

#SBATCH --job-name=connectedcomponents
#SBATCH --requeue
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=skm6640@psu.edu

# #SBATCH --time=1-00:00:00
# #SBATCH --account=kdm16_sc
# #SBATCH --partition=sla-prio

echo "Making graph for hsa"

parallel -v -j $SLURM_NTASKS 'python3 makeConnectedComponents.py output/hsa{}/hsa{}.egs > output/hsa{}/hsa{}.graph' ::: {1..22} X Y
