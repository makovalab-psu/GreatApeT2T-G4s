# DO NOT RUN THIS SCRIPT #
#!/bin/bash
exit 1

# This is an example script of how the G4 Enrichment vs GC content was generated in Humans for 5Mb window size
set -uex

sci=Homo_sapiens
spe=human
genome=3117275501

## window in kb (only integer)
winSize=5000

gzip -d ../genomes/Homo_sapiens/T2T.chm13v2.0.chr.fasta.gz
genomefa=../genomes/Homo_sapiens/T2T.chm13v2.0.chr.fasta

## to get windows
bedtools makewindows -g ../genomeLength/${sci}/chr.genome -w ${winSize}000 > ${sci}/genomeWindows${winSize}kb.${sci}.bed

## to intersect g4s and get coverage
bedtools coverage -a ${sci}/genomeWindows${winSize}kb.${sci}.bed -b ../../pG4s/${sci}/chrG.pqsfinder.filtered.bed -F 1.0 > ${sci}/all_g4s_in_Windows.${sci}.bed

## to get gc_content
for type in all; do
    bedtools nuc -fi ${genomefa} -bed ${sci}/${type}_g4s_in_Windows.${sci}.bed > ${sci}/${type}_g4s_in_Windows.withGC.${sci}.bed; 
done

## to get enrichment vs gc
types=(all)
total=$(wc -l < ../../pG4s/${sci}/chrG.pqsfinder.filtered.bed)
length=${#types[@]}
 
for (( i=0; i<length; i++ )); do
    tail -n +2 "${sci}/${types[i]}_g4s_in_Windows.withGC.${sci}.bed" | cut -f 4,6,9 | \
    awk -v genome=$genome -v totalg4s="${total[i]}" '{ enrich = ($1 / totalg4s) / ($2 / genome); print enrich "\t" $3 }' > "${sci}/${types[i]}_g4sEnrichmentvsGCcontent.${winSize}kb_Wind.${sci}.dat"
done

## remove intermediate files
rm ${sci}/genomeWindows${winSize}kb.${sci}.bed ${sci}/all_g4s_in_Windows.${sci}.bed ${sci}/all_g4s_in_Windows.withGC.${sci}.bed

gzip ../genomes/Homo_sapiens/T2T.chm13v2.0.chr.fasta