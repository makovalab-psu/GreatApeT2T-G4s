# May 3 2024, creating subset annotation file
#
# Script will read through the annotation file and collect 'gene', 'mRNA', 'CDS',
# and 'exon' entries for one (the longest) transcript 
# 
# @karolpal-jr

from BCBio import GFF
from BCBio.GFF import GFFExaminer
examiner = GFFExaminer()
import copy 
import sys


def main():

    # you should be able to find the input file by googling the 'GCF' identifier 
    # and download it from the NCBI website  
    annotation_file =  sys.argv[1]
    out_file = annotation_file.replace(".gff",".protcoding_singletrans.gff")
    handle = open(annotation_file)

    records = GFF.parse(handle)
    recs = []
    for rec in records:
        
        rec_out = copy.copy(rec)
        rec_out.features = []
        
        for feature in rec.features:

            if feature.type != "gene":
                continue

            if not "gene_biotype" in feature.qualifiers:
                continue

            if not "protein_coding" in feature.qualifiers["gene_biotype"]:
                continue

            feature_out = copy.copy(feature)
            feature_out.sub_features = []
            has_mrna = False
            subf_out = None
            
            max_len = 0
            for subf in feature.sub_features:

                if subf.type == "mRNA":
                    
                    length = 0
                    this_out = copy.copy(subf)

                    for subsubf in subf.sub_features:
                        if subsubf.type != "CDS":
                            continue
                        length += len(subsubf)

                    if length > max_len:
                        max_len = length
                        subf_out = copy.copy(this_out)
                        has_mrna = True

            if (has_mrna ):
                feature_out.sub_features.append(subf_out)
                rec_out.features.append(feature_out)
        recs.append(rec_out)    
    with open(out_file, "w") as out_handle:
        GFF.write(recs, out_handle)
    handle.close()

    
if __name__ == "__main__":
    main()
