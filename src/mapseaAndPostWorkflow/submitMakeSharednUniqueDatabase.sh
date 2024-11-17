#!/bin/bash

#SBATCH --job-name=includeUnique
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=skm6640@psu.edu

#SBATCH --time=1-00:00:00
#SBATCH --account=kdm16_sc
#SBATCH --partition=burst
#SBATCH --qos=burst4x

python3 makeSharednUniqueDatabase.py
