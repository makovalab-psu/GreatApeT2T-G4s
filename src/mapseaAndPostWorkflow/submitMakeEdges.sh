#!/bin/bash

#SBATCH --job-name=makeedges
#SBATCH --requeue
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=skm6640@psu.edu

##SBATCH --time=1-00:00:00
##SBATCH --account=kdm16_sc
##SBATCH --partition=burst
##SBATCH --qos=burst4x

# $1 = hsa #

mkdir -p output/hsa$1

echo "Making edges for hsa$1"
echo "Identifying aligned unique nodes for hsa$1"

if [ -f output/hsa$1/hsa$1.egs ]; then
	rm -f output/hsa$1/hsa$1.egs
fi

if [ -f output/hsa$1/hsa$1_alignedUnique.nds ]; then
        rm -f output/hsa$1/hsa$1_alignedUnique.nds
fi

parallel -v -j $SLURM_NTASKS 'python3 makeEdges.py {} output/hsa'"$1"'/hsa'"$1"'_alignedUnique.nds >> output/hsa'"$1"'/hsa'"$1"'.egs' ::: output/hsa$1/*rmredundant.df

mv output/hsa$1/hsa$1_alignedUnique.nds output/hsa$1/hsa$1_alignedUnique.nds.bk
cat output/hsa$1/hsa$1_alignedUnique.nds.bk | sort | uniq > output/hsa$1/hsa$1_alignedUnique.nds
rm -f output/hsa$1/hsa$1_alignedUnique.nds.bk
