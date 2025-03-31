This directory contains all the Jupyter notebooks required to analyze data and generate the plots for this paper. To ensure proper execution, run the notebooks in the sequence listed below, as they may have dependencies on prior notebooks or generated datasets. The order is as follows:

1. `upsetPlots_Phylogeny_CorrDivergenceTimes.ipynb`: Generates upset plots, phylogenetic analyses, and correlation with divergence times.
2. `G4sOverlappingNewT2Tregions.ipynb`: Analyzes pG4s overlapping newly resolved T2T (telomere-to-telomere) regions.
3. These notebooks require running the Python scripts (`../src/functional*.py`) beforehand to generate the necessary datasets:  
   - `functionalEnrichmentnMethylHuman.ipynb`: Analyzes functional enrichment and methylation data for humans.  
   - `functionalEnrichmentnMethylNHPs.ipynb`: Analyzes functional enrichment and methylation data for non-human great apes.
4. `alignedvsNonAlignedpG4s.ipynb`: Locates aligned and non-aligned species-specific pG4s across datasets.
5. `density.ipynb`: Calculates the G4 density and perfoms regression on variables-G4 density, gene density and GC content of chromosomes.
