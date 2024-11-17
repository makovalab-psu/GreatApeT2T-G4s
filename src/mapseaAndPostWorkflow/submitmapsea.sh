#!/bin/bash

#SBATCH --job-name=mapG4s
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=skm6640@psu.edu

#SBATCH --time=1-00:00:00
#SBATCH --account=kdm16_sc
#SBATCH --partition=burst
#SBATCH --qos=burst4x

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

echo "Mapping for hsa$1 -> $2:chr$4 and $3:chr$5"

python3 ../../mapsea/src/mapsea.py \
        -m ${dataDir}/lastzPairwiseAlignments/mafs/hsa$1/chr$4_$2_vs_chr$5_$3.maf \
        -b ../../datasets/pG4s \
        -o ${dataDir}/datasets/hsa$1/chr$4_$2_vs_chr$5_$3.mappedG4s.dat \
        -t tmp/hsa$1 \
        -r 1.0 \
        -d ${dataDir}/speciesDict.json