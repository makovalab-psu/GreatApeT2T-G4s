# DO NOT RUN THIS SCRIPT #
#!/bin/bash
exit 1

# This is an example script of how original gene annotation gff files were cnverted to the required format
## Example: Pongo_pygmaeus (Bornean orangutan)

## 1. To get the only single transcript per protein coding gene

python3 ../../../src/extractProteinCodingGenes.py GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff

## 2. Converting the accession IDs to chromosome names

awk 'BEGIN {FS=OFS="\t"} FNR==NR {if (NR>1) map[$2] = "chr" $1; next} /^#/ {print; next} {if ($1 in map) $1=map[$1]; print}' chr2acc GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.protcoding_singletrans.gff | \
grep "^chr" > genomic.protcoding_singletrans.chr.nocomments.gff && \
awk 'BEGIN {FS=OFS="\t"} FNR==NR {if (NR>1) map[$2] = "chr" $1; next} /^#/ {print; next} {if ($1 in map) $1=map[$1]; print}' chr2acc GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff | \
grep "^chr" > genomic.chr.nocomments.gff

## 3. To get non-protein coding genes

cat genomic.chr.nocomments.gff | \
grep 'gene_biotype=' | \
grep -v 'gene_biotype=protein_coding' | \
grep -v 'pseudogene' > genomic.chr.nocomments.nonprotcoding.gff

## 4. Remove the raw files, and zip the files for optimal storage

rm GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.protcoding_singletrans.gff GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff && \
gzip *.gff



# This is an exmaple of how enhancer annotations for NHPs (e.g. Gorilla) were generated using chain files 
# derived from T2TCHM13v2.0 to the respective NHPs T2T genome

## 1. To generate liftover
## the file `enhancer_regions.liftNames.human.bed` contains the enhancer regions in human genome with chromsome names like, chm13#1#chr1
liftOver enhancer_regions.liftNames.human.bed chainFiles/primates16.20231205_wfmash-v0.12.5.chain/chm13#1.aln.chain enhancer_regions.liftNames.NHPs.bed unMapped -multiple

## 2. To extract for species, e.g. Gorilla_gorilla (gorilla)
cat enhancer_regions.liftNames.NHPs.bed | \
grep "^mGorGor" > enhancer_regions.liftNames.gorilla.bed

## 3. To change # to _
sed -i -e 's/#/_/g' enhancer_regions.liftNames.gorilla.bed

## 4. Only get the primary chromosomes, as mentioned in the file `mGorGor1_haplo_pri`
grep -Ff mGorGor1_haplo_pri enhancer_regions.liftNames.gorilla.bed > enhancer_regions.liftNames.gorilla.pri.bed

## 5. Get the final file as required by removing chainfile names
awk 'BEGIN{FS=OFS="\t"} {split($1, arr, "_"); $1 = arr[3]; $5="enhancer"; print}' enhancer_regions.liftNames.gorilla.pri.bed | \
bedtools sort -i - > enhancer_regions.bed

## 6. Remove the intermediate files and zip the final file
rm enhancer_regions.liftNames.gorilla.bed enhancer_regions.liftNames.gorilla.pri.bed && \
gzip enhancer_regions.bed



# To generate cpg island annotations for the respective non-human great ape species
# The script can be obtained by `wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/cpg_lh`

## 1. This command was invoked with different chromsomes seprataed into fasta files `chr${chr}_${speciesName}.fa.gz`
for chr in $(seq 1 ${maxAUTOSOMES}) X Y; do
        ./cpg_lh <(zcat ../genomes/Gorilla_gorilla/GenomeSplit/chr${chr}_${speciesName}.fa.gz) | awk '{$2 = $2 - 1; width = $3 - $2; \
        printf("%s\t%d\t%s\tCpGI\t%s\t.\n", $1, $2, $3, width);}' \
        | sort -k1,1 -k2,2n > chr${chr}_cpgislands.org.bed;

        awk '{split($1, a, "."); $1=a[2]; print}' OFS='\t' chr${chr}_cpgislands.org.bed > chr${chr}_cpgislands.bed;
        rm -f chr${chr}_cpgislands.org.bed;
done;

## 2. Concatenate all chromsome files into a single file
cat chr*_cpgislands.bed >> chrG_cpgislands.bed

## 3. Remove the intermediate files and zip the final file
rm -f chr*_cpgislands.bed && \
gzip chrG_cpgislands.bed