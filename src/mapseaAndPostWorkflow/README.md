##executing the mapsea pipeline
#repeat the steps from hsa 1..22 to Y
hsa=1; cat lastzPairwiseAlignments/pairs/pairs_hsa${hsa}.dat | parallel --colsep "\t" sbatch submitmapsea.sh ${hsa} {1} {2} {3} {4}

##for removing the redundant G4s (flank can be adjusted inside the job submit file)
#repeat the steps from hsa 1..22 to Y
hsa=1; cat lastzPairwiseAlignments/pairs/pairs_hsa${hsa}.dat | parallel --colsep "\t" sbatch submitrefiner.sh ${hsa} {1} {2} {3} {4}


##for making edges
#repeat the steps from hsa 1..22 and X,Y
hsa=1; sbatch submitMakeEdges.sh ${hsa}

##make connected components
sbatch submitConnectedComponents.sh

##determine unique G4s and make the final dataframe readable output
sbatch submitMakeSharednUniqueDatabase.sh
