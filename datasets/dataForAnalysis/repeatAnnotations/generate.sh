# DO NOT RUN THIS SCRIPT #
#!/bin/bash
exit 1

# This is an example of how repeat annotations were converted to desired format
## 1. For human get the desired columns
less -S chm13v2.0.rmsk.bed | \
cut -f 1,7,8,4,5,6 | \
awk '{print $1,$5,$6,$2,$3,$4}' OFS='\t' > genomic.repeats.bed

## For non-huamn great apes
## 1. Remove the hsa identifier
awk 'BEGIN{FS=OFS="\t"} {split($1, arr, "_"); $1 = arr[1]; print}' mPanPan1.pri.cur.20231122.combo3.bed > genomic.repeats.bed

## 2. Remove the raw files, and zip the final file
rm chm13v2.0.rmsk.bed mPanPan1.pri.cur.20231122.combo3.bed && \
gzip genomic.repeats.bed