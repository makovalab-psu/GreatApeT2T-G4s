# DO NOT RUN THIS SCRIPT #
#!/bin/bash
exit 1

## This is an example script to generate G4s in the newly resolved regions in T2T genome
## The below steps are valid for human only
## This file `g4s_inNewRegionsT2T.bed` contains the g4s in the newly resolved regions in human T2T genome
## 1. To get the G4s in the new regions
for i in {1..22} X Y; do cat g4s_inNewRegionsT2T.bed | grep "^chr$i[^0-9]" | wc -l | awk -v chr=$i '{print "chr"chr "\t" $1}' >> newG4sinnewRegions.dat; done;

## 2. to get the new regions
## The file `chm13v2-unique_to_hg38.bed` is downloaded form human-T2T github repository
for i in {1..23} X Y; do cat chm13v2-unique_to_hg38.bed | grep "^chr$i[^0-9]" | awk -v chr=$i '{sum += $3-$2; } END {print "chr"chr "\t" sum;}' >> newRegions.dat; done;

### 3*. To intersect the new G4s in human to repeats
bedtools intersect -a g4s_inNewRegionsT2T.bed -b <(gzcat ../../repeatAnnotations/Homo_sapiens/genomic.repeats.bed.gz) -v > g4s_inNewregions_NOTinRepeats.bed
bedtools intersect -a ../../../pG4s/Homo_sapiens/chrG.pqsfinder.filtered.bed -b g4s_inNewRegionsT2T.bed -v > g4s_inOldRegionsT2T.bed
bedtools intersect -a g4s_inOldRegionsT2T.bed -b <(gzcat ../../repeatAnnotations/Homo_sapiens/genomic.repeats.bed.gz) -v > g4s_inOldregions_NOTinRepeats.bed

## The below steps are valid for NHPs (Non-human great apes) only
## This file in NHPs however, `g4s_inNewRegionsT2T.bed` contains the newly resolved regions in T2T genome
## 1. For getting the desired file type columns, with only the chr information
cat g4s_inNewRegionsT2T.bed | awk '{split ($1, a, "_"); print a[1] "\t" $2 "\t" $3}' | bedtools merge -i - > g4s_inNewRegionsT2T.modified.bed

for i in {1..23} X Y; do 
    cat ../../../pG4s/Pan_paniscus/chrG.pqsfinder.filtered.bed | \
    grep "^chr$i[^0-9]" > chr${i}.pqsfinder.filtered.bed;
done

## 2. For intersecting it to g4 annotations
for i in {1..23} X Y; do 
    bedtools intersect -a chr${i}.pqsfinder.filtered.bed -b g4s_inNewRegionsT2T.modified.bed -wa -wb > chr${i}.g4sinNewRegions.bed; 
done

## 3. To get the total G4s
for i in {1..23} X Y; do 
    cat chr${i}.g4sinNewRegions.bed | \
    wc -l | \
    awk -v chr=$i '{print "chr"chr "\t" $1}' >> newG4sinnewRegions.dat; 
done

## 4. To get the sum of total newly resolved regions
for i in {1..23} X Y; do 
    cat g4s_inNewRegionsT2T.modified.bed | \
    grep "^chr$i[^0-9]" | \
    awk -v chr=$i '{sum += $3-$2; } END {print "chr"chr "\t" sum;}' >> newRegions.dat; 
done

## 5. Remove the chromosome bed files
rm -f chr*.pqsfinder.filtered.bed chr*.g4sinNewRegions.bed

## 6. See the total number of new G4s
cut -f 2 newG4sinnewRegions.dat | awk '{sum+=$1} END {print sum}'
