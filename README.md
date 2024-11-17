# Evolutionary Dynamics of G-Quadruplexes in Human and Other Great Ape Telomere-to-Telomere Genomes  
#### This repository contains all the code and data generated for the paper.

Saswat K. Mohanty, Francesca Chiaromonte, Kateryna D. Makova*
 
Department of Biology, Penn State University, University Park, PA 16802  
Department of Statistics, Penn State University, University Park, PA 16802  
*Correspondence to Kateryna D. Makova ([kdm16@psu.edu](mailto:kdm16@psu.edu))  
 
## Directory Structure

This repository contains the following directories:

- **`src`**: Includes all Python and Bash scripts necessary for analyzing data, generating datasets, or creating plots as required for this manuscript.

- **`jupyterNotebooks`**: Contains Python-based Jupyter notebooks for analyzing data and generating plots or datasets relevant to the manuscript.

- **`datasets`**: Hosts all datasets used or generated during the analyses described in this manuscript.

- **`plots`**:  Stores all plots generated as part of the manuscript’s analyses.

**Some key dependencies required for the analyses in this manuscript are:**

- *`mapsea`*: A general-purpose repository containing scripts for mapping short elements stored in `.BED` files (e.g., G-quadruplex intervals) to `.MAF` (alignment) files.

- *`g4Discovery`*: A general-purpose repository with scripts for predicting and annotating G-quadruplexes (G4s) in genome sequences, combining `pqsfinder` and `G4Hunter`.

## **Predicted G4s with Evolutionary Information**

For downloading pG4s across the six great ape species—*Homo sapiens*, *Pan paniscus*, *Pan troglodytes*, *Gorilla gorilla*, *Pongo abelii*, and *Pongo pygmaeus*—navigate to the `datasets/pG4s/` directory. Each species subdirectory contains a file named: **`wholeGenome.pG4s.withSharingInfo.bed`**. This file uses a custom `.BED` format that includes additional evolutionary and sequence information. The columns in this file are as follows:

```markdown
- CHR: chromosome
- START: pG4 start
- END: pG4 end
- ID: #unique 10-digit ID assigned to this pG4, other species sharing this pG4 will have same ID
- PQSSCORE: pqsfinder score
- STRAND: strand in which pG4 is present
- G4HUNTERSCORE: G4Hunter score
- SHARED: shared pG4(S), aligned (ASS) or unaligned species-specific (USS) pG4
- WITH: with which species it is shared 
		{B: Bonobo, C: Chimpanzee, H: Human, G: Gorilla, 
		Bo: Bornean Orangutan, So: Sumatran Orangutan}
- DUPLICATED: If it is shared with the same chromosome in the same species (>1 indicates duplicated)
- SEQ: pG4 sequence in positive-strand
- MOTIF: standard (S) or bulged (B) pG4 motif
```

**Using the file as a standard BED file**: To convert the file into a standard 6-column `.BED` format (removing evolutionary and sequence details), use the following command:

```bash
cut -f 1-6 wholeGenome.pG4s.withSharingInfo.bed > wholeGenome.pG4s.bed
```

## Citation
If you use any of this data or tool in your research, please cite the following paper:

> Mohanty, S. K., Chiaromonte, F., & Makova, K. (2024). "[Evolutionary Dynamics of G-Quadruplexes in Human and Other Great Ape Telomere-to-Telomere Genomes](https://www.biorxiv.org/content/10.1101/2024.11.05.621973v1). *bioRxiv*, 2024-11. `doi: 10.1101/2024.11.05.621973`
