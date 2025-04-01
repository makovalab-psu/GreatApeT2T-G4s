## To find multiple annotation intersections
bedtools multiinter -header -names promoter utr5 cds intron utr3 enhancer npcgene or cpgi repeats ngnr -i <(gzcat promoter_regions.bed.gz) <(gzcat utr5_regions.bed.gz) <(gzcat cds_regions.bed.gz) <(gzcat intron_regions.bed.gz) <(gzcat utr3_regions.bed.gz) <(gzcat enhancer_regions.bed.gz) <(gzcat npcgene_regions.bed.gz) <(gzcat ../../dataForAnalysis/geneAnnotations/Homo_sapiens/originsofReplication/core_stochastic_origins_chm13v2.0_hglft.ucsc.bed.gz | bedtools sort -i - -faidx ../../dataForAnalysis/genomeLength/Homo_sapiens/chr.genome) <(gzcat ../../dataForAnalysis/geneAnnotations/Homo_sapiens/cpgislands/chrG_cpgislands.bed.gz) <(gzcat repeats_regions.bed.gz) <(gzcat ngnr_regions.bed.gz) -g ../../dataForAnalysis/genomeLength/Homo_sapiens/chr.genome -empty > humanFunctionalAnnotationMultiInterval.bed
gzip humanFunctionalAnnotationMultiInterval.bed

## To get the number of bases for each type of intersection
awk '{len = $3 - $2; bases[$5] += len} END {for (type in bases) print bases[type]"\t"type}' <(gzcat humanFunctionalAnnotationMultiInterval.bed.gz) | sort -nr > typeSorted.dat

## To get the top 20 hits
less typeSorted.dat | grep "," | head -n 20
