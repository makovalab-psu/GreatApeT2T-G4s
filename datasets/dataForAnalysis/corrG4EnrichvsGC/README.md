set -uex

sci=Homo_sapiens
spe=human
genome=3117275501
# window in kb (only integer)
winSize=5000

genomefa=../../../../mapG4ToAlignments/pqsfinder/fasta_chr/T2T.chm13v2.0.chr.fasta
# genomefa=../../../genomes/${sci}/chrG_${spe}.processed.fa

# to get windows
bedtools makewindows -g ../../../gene_annotations/${sci}/chr.genome -w ${winSize}000 > genomeWindows${winSize}kb.${sci}.bed

# to intersect g4s and get coverage
bedtools coverage -a genomeWindows${winSize}kb.${sci}.bed -b ../../../g4set_fork/${sci}/chrG.pqsfinder.filtered.bed -F 1.0 > all_g4s_in_Windows.${sci}.bed

# to get gc_content
for type in all; do
    bedtools nuc -fi ${genomefa} -bed ${type}_g4s_in_Windows.${sci}.bed > ${type}_g4s_in_Windows.withGC.${sci}.bed; 
done

# to get enrichment vs gc
types=(all)
total=$(wc -l < ../../../g4set_fork/${sci}/chrG.pqsfinder.filtered.bed)

length=${#types[@]}

for (( i=0; i<length; i++ )); do
    tail -n +2 "${types[i]}_g4s_in_Windows.withGC.${sci}.bed" | cut -f 4,6,9 | \
    awk -v genome=$genome -v totalg4s="${total[i]}" '{ enrich = ($1 / totalg4s) / ($2 / genome); print enrich "\t" $3 }' > "${types[i]}_g4sEnrichmentvsGCcontent.${winSize}kb_Wind.${sci}.dat"
done

for (( i=0; i<length; i++ )); do
   tail -n +2 "${types[i]}_g4s_in_Windows.withGC.${sci}.bed" | cut -f 4,6,9 | \
   awk -v genome=3117275501 -v totalg4s="${total[i]}" '{ enrich = $1 / $2; print enrich "\t" $3 }' > "${types[i]}_g4sDensityvsGCcontent.${winSize}kb_Wind.${sci}.dat"
done