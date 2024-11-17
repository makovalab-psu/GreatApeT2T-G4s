### Files name with: wholeGenome.pG4s.withSharingInfo.bed 

These files are per species BED files, containing the information of sharing, as well as their location and sequence information. The columns in the `wholeGenome.pG4s.withSharingInfo.bed` are as follows in the sequential manner:

- CHR: chromosome
- START: pG4 start
- END: pG4 end
- ID: #unique 10-digit random ID assigned to this pG4, other species sharing this pG4 will have same ID
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

**NOTE**: *To extract as a bed file, please use the first 6 columns, since the last 6 columns do not follow the specified BED file guidelines, they are used for providing non-canonical information about pG4 sharing. To extract the first 6 cols use:*
`cut -f 1-6 wholeGenome.pG4s.withSharingInfo.bed`

### Files name with: chrG.pqsfinder.filtered.bed

These files are needed for analysis. The columns in the `chrG.pqsfinder.filtered.bed` are as follows in the sequential manner:

- CHR: chromosome
- START: pG4 start
- END: pG4 end
- PQSSCORE: pqsfinder score
- LENGTH: pG4 length
- STRAND: strand in whch pG4 is present
- G4HUNTERSCORE: G4Hunter score

#### InflectionTest folder in Homo_sapiens

This folder contains, BED file of pG4s generated using a pqsfinder score threshold of 30 and G4Hunter score thresgold of 0, to look at inflection points. 
