#!/bin/bash

#SBATCH --job-name=rmredundantG4s
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=skm6640@psu.edu

##SBATCH --time=1-00:00:00
##SBATCH --account=kdm16_sc
##SBATCH --partition=burst
##SBATCH --qos=burst4x

# Use the pairs folder data, with paralle invoking for each hsa
# $1 = hsa #
# $2 = species 1 name
# $3 = species 2 name
# $4 = species 1 chr no.
# $5 = species 2 chr no.

# this is a placeholder directory
# change to the directory (as needed) where the data is located
dataDir=../../mapping

mkdir -p ${dataDir}/datasets/hsa$1

echo "Removing redundant G4s for hsa$1 -> $2:chr$4 and $3:chr$5"

merge_factor=3

python3 ../../mapsea/src/refiner.py \
        -d ${dataDir}/datasets/hsa$1/chr$4_$2_vs_chr$5_$3.mappedG4s.dat \
        -f $merge_factor \
        -c $SLURM_CPUS_PER_TASK > ${dataDir}/datasets/hsa$1/chr$4_$2_vs_chr$5_$3.mappedG4s.rmredundant.df
