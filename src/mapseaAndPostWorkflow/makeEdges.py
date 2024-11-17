import pandas as pd
import sys

''' 
This script takes a dataframe with the output of MAP-SEA: STEP2 
and returns the edges of the graph, as well as the G4s which are
aligned, but not shared within the alignment.
'''

inputfile = sys.argv[1]
outputfile = sys.argv[2]

# read the dataframe given by the output of mapsea in the format .df
df = pd.read_csv(inputfile, header=0, sep="\t")
df = df[df["TYPE"] == "fullMAF"] # filter the dataframe to only fullMAF entries

# to account the real start prior to merging
df["REAL_START"] = df["START"] + df["INDEX"] 
# give each g4 a unique label
df['LABEL'] = df['SPECIES'].astype(str) + '|' + df['CHR'].astype(str) + '|' + df['REAL_START'].astype(int).astype(str) + '|' + df['STRAND'] + '|' + df['QUERY_LENGTH'].astype(int).astype(str) + '|' + df['SCORE'].astype(float).astype(str)

# get each pair of G4s that are aligned/unique by grouping based on IDs and G4s
grouped = df.groupby(["ID","G4"])

labels = []
aligned_labels = []
for idx, group in grouped:
    if len(group) > 1: # i.e. 2 for pairwise alignments
        sublabels = [i for i in group["LABEL"]]
        labels.append(sublabels) #append the two g4s to a list
    elif len(group) == 1: # for G4s that are unique but aligned
        aligned_labels.append(group.index[0])

#to account for pairwise shared G4s
labels = pd.DataFrame(labels).sort_values(by=[0])
labels.to_csv(sys.stdout, header=False, index=False, sep="\t")

#to account for unique G4s which are aligned
aligned_labels = df['LABEL'].loc[aligned_labels] 
aligned_labels.to_csv(outputfile, header=False, index=False, mode="a")