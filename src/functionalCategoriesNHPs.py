import gff2bed
import gzip
import numpy as np
import pandas as pd
from pybedtools import BedTool
import warnings

from pandas.errors import SettingWithCopyWarning
warnings.simplefilter("ignore", category=SettingWithCopyWarning)

''' Functions defined '''
def generate_filtered_species(df, presList, species):
    df_grouped = df.groupby("ID")
    indexInclude = []
    for _, group in df_grouped:
        species_present = list(np.sort(group["SPECIES"].unique()))
        if species_present == presList:
            indexInclude.extend((group[group["SPECIES"] == species]).index)
    return df.loc[indexInclude].set_index("ID")

def mutx_to_bed(dataframe):
    bed = pd.DataFrame()
    dataframe = dataframe.reset_index()
    bed["CHROM"] = "chr" + dataframe["CHR"]#[f"chr{str(hsa)}"]*len(dataframe)
    bed["START"] = dataframe["START"].astype(int)
    bed["END"] = dataframe["START"].astype(int) + dataframe["LENGTH"].astype(int)
    bed["SCORE"] = dataframe["SCORE"].astype(float)
    bed["LENGTH"] = dataframe["LENGTH"].astype(int)
    bed["STRAND"] = dataframe["STRAND"]
    return bed

''' Non-human great apes'''
speciesSciName = ["Pan_paniscus","Pan_troglodytes","Gorilla_gorilla","Pongo_abelii","Pongo_pygmaeus"]
speciesName = ["bonobo","chimp","gorilla","sorang","borang"]

''' Directories and files '''
geneAnnot = "../datasets/dataForAnalysis/geneAnnotations"
repeatAnnot = "../datasets/dataForAnalysis/repeatAnnotations"
faidxLoc = "../datasets/dataForAnalysis/genomeLength"
outDir= "../datasets/functionalOutputs"

''' Annotations and pG4s overlap'''
hsa = "G"
df_ann_mutx = pd.read_csv(f"../datasets/pipelineOutputs/allhsa{hsa}.graph.df.gz", header=0, sep="\t", low_memory=False, compression='gzip') #read the file

for nos, name in enumerate(speciesSciName):

    print(f"Starting analysis for {speciesName[nos]} ...")
    whole_gff = gzip.open(f"{geneAnnot}/{name}/genomic.chr.nocomments.gff.gz",'rt') # read the gff file
    whole_gff_gene = pd.DataFrame(gff2bed.convert(gff2bed.parse(whole_gff, type="gene"))) # convert the gff to bed format using the gff2bed module
    whole_gff_gene = BedTool.from_dataframe(whole_gff_gene) #convert to a pybedtool object

    protCode_singleTrans_gff_file = gzip.open(f"{geneAnnot}/{name}/genomic.chr.nocomments.protcoding_singletrans.gff.gz",'rt') # read the gff file
    protCode_singleTrans_gff = pd.DataFrame(gff2bed.convert(gff2bed.parse(protCode_singleTrans_gff_file, type=None))) # convert the gff to bed format using the gff2bed module
    protCode_singleTrans_gff = BedTool.from_dataframe(protCode_singleTrans_gff) #convert to a pybedtool object
    protCode_singleTrans_gff_df = protCode_singleTrans_gff.to_dataframe()

    print("Curating different annotation files ...")

    # to get the gene gff annotations from gene data
    gene_gff = BedTool.from_dataframe(protCode_singleTrans_gff_df[protCode_singleTrans_gff_df["name"].str.startswith("gene-")])
    BedTool.sort(gene_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/gene_regions.bed.gz")

    # to get the non-protein coding gene gff annotations from gene data
    npcgene_gff = gzip.open(f"{geneAnnot}/{name}/genomic.chr.nocomments.nonprotcoding.gff.gz",'rt')
    npcgene_gff = pd.DataFrame(gff2bed.convert(gff2bed.parse(npcgene_gff)))
    npcgene_gff = BedTool.from_dataframe(npcgene_gff)
    BedTool.sort(npcgene_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/npcgene_regions.bed.gz")

    # to get the rna gff annotations from gene data
    rna_gff = BedTool.from_dataframe(protCode_singleTrans_gff_df[protCode_singleTrans_gff_df["name"].str.startswith("rna-")])
    BedTool.sort(rna_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/rna_regions.bed.gz")

    # to get the cds gff annotations from gene data
    cds_gff = BedTool.from_dataframe(protCode_singleTrans_gff_df[protCode_singleTrans_gff_df["name"].str.startswith("cds-")])
    BedTool.sort(cds_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/cds_regions.bed.gz")

    # to get the exon gff annotations from gene data
    exon_gff = BedTool.from_dataframe(protCode_singleTrans_gff_df[protCode_singleTrans_gff_df["name"].str.startswith("exon-")])
    BedTool.sort(exon_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/exon_regions.bed.gz")

    # to get the intron gff annotations from gene data
    for_intron = protCode_singleTrans_gff_df[(protCode_singleTrans_gff_df["name"].str.startswith("rna-")) | (protCode_singleTrans_gff_df["name"].str.startswith("exon-"))] # for all kinds of RNA
    intron_list = []
    for idx, row in for_intron.iterrows():
        nameRow = row["name"].split("-")
        nameElem = nameRow[0]
        if nameElem == "rna":
            namernaId = "-".join(nameRow[1:])
            rnaStrand = row["strand"]
            exon_count = 0
        elif nameElem == "exon" and "namernaId" in locals():
            if "-".join(nameRow[1:-1]) == namernaId and rnaStrand == row["strand"]:
                exon_count += 1
                if rnaStrand == "+":
                    if exon_count == 1:
                        start = row["end"]
                    elif exon_count > 1:
                        end = row["start"]
                        intron_list.append([row["chrom"], start, end, f"intron-{namernaId}-{exon_count-1}", 0, rnaStrand])
                        start = row["end"]
                elif rnaStrand == "-":
                    if exon_count == 1:
                        end = row["start"]
                    elif exon_count > 1:
                        start = row["end"]
                        intron_list.append([row["chrom"], start, end, f"intron-{namernaId}-{exon_count-1}", 0, rnaStrand])
                        end = row["start"]
    intron_gff = pd.DataFrame(intron_list, columns=["chrom","start","end","name","score","strand"])
    intron_gff = BedTool.from_dataframe(intron_gff)
    BedTool.sort(intron_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/intron_regions.bed.gz")

    # to get the UTR gff annotations from gene data
    for_utr = protCode_singleTrans_gff_df[(protCode_singleTrans_gff_df["name"].str.startswith("exon-")) | (protCode_singleTrans_gff_df["name"].str.startswith("cds-"))].reset_index(drop=True) # for all kinds of RNA
    utr5_list = []
    utr3_list = []
    for idx, row in for_utr.iterrows():
        exon_start = 0
        exon_end = 0
        cds_start = 0
        cds_end = 0
        nameRow = row["name"].split("-")
        if row["name"].startswith("exon-"):
            if nameRow[-1] == "1":
                current_exon_id = "-".join(nameRow[1:-1])
                first_exon_loc = idx
                cdsFound = 0
        elif row["name"].startswith("cds-"):
            nameRow = row["name"].split("-")
            if cdsFound == 0:
                current_cds_id = "-".join(nameRow[1:])
                last_exon_loc = idx - 1
                first_cds_loc = idx
            cdsFound = 1
            if len(for_utr) == idx + 1 or for_utr.loc[idx+1]["name"].startswith("exon-"):
                last_cds_loc = idx
                if row["strand"] == "+":
                    utr5_start = for_utr.loc[first_exon_loc]["start"]
                    utr5_end = for_utr.loc[first_cds_loc]["start"]
                    utr3_start = for_utr.loc[last_cds_loc]["end"]
                    utr3_end = for_utr.loc[last_exon_loc]["end"]
                    if utr5_start != utr5_end:
                        utr5_list.append([row["chrom"], utr5_start, utr5_end, f"utr-5-{current_exon_id}", 0, row["strand"]])
                    if utr3_start != utr3_end:
                        utr3_list.append([row["chrom"], utr3_start, utr3_end, f"utr-3-{current_exon_id}", 0, row["strand"]])
                elif row["strand"] == "-":
                    utr5_start = for_utr.loc[first_cds_loc]["end"]
                    utr5_end = for_utr.loc[first_exon_loc]["end"]
                    utr3_start = for_utr.loc[last_exon_loc]["start"]
                    utr3_end = for_utr.loc[last_cds_loc]["start"]
                    if utr5_start != utr5_end:
                        utr5_list.append([row["chrom"], utr5_start, utr5_end, f"utr-5-{current_exon_id}", 0, row["strand"]])
                    if utr3_start != utr3_end:
                        utr3_list.append([row["chrom"], utr3_start, utr3_end, f"utr-3-{current_exon_id}", 0, row["strand"]])
    utr5_gff = BedTool.from_dataframe(pd.DataFrame(utr5_list, columns=["chrom","start","end","name","score","strand"]))
    utr5_gff = utr5_gff.subtract(intron_gff,s=True) #remove the intronic regions
    BedTool.sort(utr5_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/utr5_regions.bed.gz")
    utr3_gff = BedTool.from_dataframe(pd.DataFrame(utr3_list, columns=["chrom","start","end","name","score","strand"]))
    utr3_gff = utr3_gff.subtract(intron_gff,s=True) #remove the intronic regions
    BedTool.sort(utr3_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/utr3_regions.bed.gz")

    # to get the promoter gff annotations from gene data, as 1kb upstream of the gene
    protCode_singleTrans_gff_df = protCode_singleTrans_gff.to_dataframe()
    for_promoter = protCode_singleTrans_gff_df[protCode_singleTrans_gff_df["name"].str.contains("gene-")] #contains both genes and pseudogenes
    #For genes on the positive strand
    promoter_start = for_promoter['start'] - 1000
    promoter_end = for_promoter['start']
    # For genes on the negative strand
    promoter_start[for_promoter['strand'] == '-'] = for_promoter['end']  
    promoter_end[for_promoter['strand'] == '-'] = for_promoter['end'] + 1000
    promoter_gff = pd.DataFrame({
        'chrom': for_promoter['chrom'],
        'start': promoter_start,
        'end': promoter_end,
        'name': for_promoter["name"].str.replace("gene-","promoter-"),
        'score': for_promoter['score'],
        'strand': for_promoter['strand']
    })
    promoter_protCode_singleTrans_gff = BedTool.from_dataframe(promoter_gff)
    promoter_gff = promoter_protCode_singleTrans_gff.subtract(gene_gff, s=True) #remove the gene regions
    BedTool.sort(promoter_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/promoter_regions.bed.gz")

    # # to get the enhancer gff annotations from gene data
    enhancer_gff = BedTool(f"{outDir}/{name}/enhancer_regions.bed.gz")

    # to get the repeats annotations from gene data, as 1kb upstream
    repeats = BedTool(f"{repeatAnnot}/{name}/genomic.repeats.bed.gz")
    repeats = repeats.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    repeats["score"] = "repeats"
    repeats["name"] = "repeats-" + repeats["name"]
    repeats_gff = BedTool.sort(BedTool.from_dataframe(repeats)) #named gff to match other annotation variables
    BedTool.sort(repeats_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/repeats_regions.bed.gz")

    #get cpg annotations
    cpgi_gff = pd.read_csv(f"{geneAnnot}/{name}/cpgIslands/chrG_cpgislands.bed.gz",header=None,sep="\t", compression='gzip')
    cpgi_gff[3] = "cpgi-." # to account for generalisation
    cpgi_gff = BedTool.from_dataframe(cpgi_gff)
     
    # to get the regions which are non-(genes, enhancers, promoters, repeats, non-protein coding genes, ors, cpgis)
    gr_gff = BedTool.sort(BedTool.cat(gene_gff, promoter_gff, enhancer_gff, repeats_gff, npcgene_gff, cpgi_gff), faidx=f"{faidxLoc}/{name}/chr.genome")
    ngnr_gff = gr_gff.complement(g=f"{faidxLoc}/{name}/chr.genome", L=True).to_dataframe()
    ngnr_gff["name"] = ngnr_gff.apply(lambda row: f"id-{row.name + 1}", axis=1)
    ngnr_gff["score"] = "ngnr"
    ngnr_gff["strand"] = "."
    ngnr_gff = BedTool.sort(BedTool.from_dataframe(ngnr_gff))
    BedTool.sort(ngnr_gff,faidx=f"{faidxLoc}/{name}/chr.genome").saveas(f"{outDir}/{name}/ngnr_regions.bed.gz")

    print("Generating shared and unique G4 datasets ...")

    speciesofInterest = nos + 2
    speciesSpecificSpec = [speciesofInterest]
    speciesSpecificBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, speciesSpecificSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
    speciesSpecificBed.saveas(f"{outDir}/{name}/speciesSpecificG4s{speciesName[nos]}.bed.gz")

    hominidSpec = [1,2,3,4,5,6]
    hominidBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, hominidSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
    hominidBed.saveas(f"{outDir}/{name}/HominidG4s{speciesName[nos]}.bed.gz")

    if speciesofInterest == 2 or speciesofInterest == 3: # for bonobo and chimp
        paniniSpec = [2,3]
        paniniBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, paniniSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
        paniniBed.saveas(f"{outDir}/{name}/PaniniG4s{speciesName[nos]}.bed.gz")

        homininiSpec = [1,2,3]
        homininiBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, homininiSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
        homininiBed.saveas(f"{outDir}/{name}/HomininiG4s{speciesName[nos]}.bed.gz")
    
        homininaeSpec = [1,2,3,4]
        homininaeBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, homininaeSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
        homininaeBed.saveas(f"{outDir}/{name}/HomininaeG4s{speciesName[nos]}.bed.gz")
        g4types = ["speciesSpecific", "hominid", "homininae", "hominini", "panini"]

    elif speciesofInterest == 4: #for gorilla
        homininaeSpec = [1,2,3,4]
        homininaeBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, homininaeSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
        homininaeBed.saveas(f"{outDir}/{name}/HomininaeG4s{speciesName[nos]}.bed.gz")
        g4types = ["speciesSpecific", "hominid", "homininae"]

    elif speciesofInterest == 5 or speciesofInterest == 6: # for orangutans
        ponginiSpec = [5,6]
        ponginiBed = BedTool.from_dataframe(mutx_to_bed(generate_filtered_species(df_ann_mutx, ponginiSpec, speciesofInterest)).sort_values(by=["CHROM","START","END"]).drop_duplicates().reset_index(drop=True))
        ponginiBed.saveas(f"{outDir}/{name}/PonginiG4s{speciesName[nos]}.bed.gz")
        g4types = ["speciesSpecific", "hominid", "pongini"]

    print("Intersecting shared and unique G4s with functional annotations ...")

    elements = ["gene", "rna", "exon", "intron", "utr5", "utr3", "cds", "npcgene", "promoter", "enhancer", "cpgi", "repeats", "ngnr"]

    for g4type in g4types:
        for element in elements:

            if element == "gene" or element == "npcgene" or element == "rna" or element == "cds" or element == "exon" or element == "intron" or element.startswith("utr"): #elements corresponding to a gene
                intersect = BedTool.to_dataframe(eval(f"{g4type}Bed").intersect(eval(f"{element}_gff"), f=1.0, s=True, wb=True)) #intersect the bed file with the gff file for same strand
                intersect_opp = BedTool.to_dataframe(eval(f"{g4type}Bed").intersect(eval(f"{element}_gff"), f=1.0, S=True, wb=True)) #intersect the bed file with the gff file for opposite strand
            else:
                intersect = BedTool.to_dataframe(eval(f"{g4type}Bed").intersect(eval(f"{element}_gff"), f=1.0, wb=True)) #intersect the bed file with the gff file 

            #self-defined types
            if element != None and element != "intron" and element != "promoter" and element != "utr5" and element != "utr3" and element != "enhancer" and element != "repeats" and element != "ngnr" and element != "npcgene":
                intersect = intersect[intersect["blockCount"].str.startswith(element.lower())]
                if element == "gene" or element == "npcgene" or element == "rna" or element == "cds" or element == "exon" or element == "intron" or element.startswith("utr"):
                    intersect_opp = intersect_opp[intersect_opp["blockCount"].str.startswith(element.lower())]

            try:
                new_rows = intersect["blockCount"].str.split("-", n = 1, expand = True)
                new_rows.columns = ["type","typename"]

                intersect = pd.concat([intersect, new_rows], axis=1)
                intersect.drop(columns=['blockCount'], inplace=True)
                intersect.drop_duplicates(inplace=True)

                if element == "gene" or element == "npcgene" or element == "rna" or element == "cds" or element == "exon" or element == "intron" or element.startswith("utr"):
                    new_rows_opp = intersect_opp["blockCount"].str.split("-", n = 1, expand = True)
                    new_rows_opp.columns = ["type","typename"]

                    intersect_opp = pd.concat([intersect_opp, new_rows_opp], axis=1)
                    intersect_opp.drop(columns=['blockCount'], inplace=True)
                    intersect_opp.drop_duplicates(inplace=True)
            except:
                pass

            intersect.drop_duplicates(inplace=True)
            intersect.columns = ["chrom","start","end","score","length","strand","typeChrom","typeStart","typeEnd","type","typeStrand","typeTag","typeName"]

            if element == "gene" or element == "npcgene" or element == "rna" or element == "cds" or element == "exon" or element == "intron" or element.startswith("utr"):
                intersect_opp.drop_duplicates(inplace=True)
                intersect_opp.columns = ["chrom","start","end","score","length","strand","typeChrom","typeStart","typeEnd","type","typeStrand","typeTag","typeName"]

            # self-defined types
            if element == "utr5" or element == "utr3" or element == "intron":
                intersect["type"] = intersect["typeTag"]
                intersect_opp["type"] = intersect_opp["typeTag"]

            intersect.reset_index(drop=True)
            if element == "gene" or element == "npcgene" or element == "rna" or element == "cds" or element == "exon" or element == "intron" or element.startswith("utr"):
                intersect_opp.reset_index(drop=True)

            intersect.to_csv(f"{outDir}/{name}/allhsa{hsa}.intersected.betn.{speciesName[nos]}_{element if element != None else 'all'}.{g4type}G4s.bed.gz", sep="\t", index=False, header=None, compression='gzip') #save the intersected file
            if element == "gene" or element == "npcgene" or element == "rna" or element == "cds" or element == "exon" or element == "intron" or element.startswith("utr"):
                intersect_opp.to_csv(f"{outDir}/{name}/allhsa{hsa}.intersected.betn.{speciesName[nos]}_{element if element != None else 'all'}.opp.{g4type}G4s.bed.gz", sep="\t", index=False, header=None, compression='gzip') #save the intersected file


    length = []
    for element in elements:
        length_gff = BedTool.merge(BedTool.sort(eval(f"{element}_gff"))).to_dataframe()
        length_gff["length"] = length_gff["end"] - length_gff["start"]
        length.append(sum(length_gff['length']))
    print(f"\n{elements}:\n{length}\n")