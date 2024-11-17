This directory contains the database of pG4s generated using the mapsea tool followed by the graph-based pipeline, as described in the paper. 

**Files in This Directory**
- `allhsaG.graph.df.gz`: Contains information about pG4s that are shared or unique across the whole genome for all great apes.
- `uniquehsaN.egs.gz`: Includes markers of G4s that are unique to any species. Format: `species|chr|start|strand|length|score`.
- `alignedUniquehsaG.egs.gz`: Includes markers of G4s unique to any species that are also present in alignments. Format: `species|chr|start|strand|length|score`.

**Subdirectories (HSA 1–22, X, Y)**
Each subdirectory corresponds to a specific human homolog chromosome group (HSA 1–22, X, Y) and contains the following files:
- `allhsaN.graph.df.gz`: Contains chromosome-specific information about pG4s that are shared or unique.
- `uniquehsaN.egs.gz`: Includes markers of G4s unique to any species within the specific HSA. Format: `species|chr|start|strand|length|score`.