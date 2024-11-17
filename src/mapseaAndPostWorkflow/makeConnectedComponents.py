import networkx as nx
import pandas as pd
import sys

''' 
This script takes a file with two columns representing 
the edges of a graph and returns the connected components 
of the graph.
'''

edgefile = sys.argv[1]

edges = pd.read_csv(edgefile, header=None, sep="\t")

G = nx.Graph()
for _, row in edges.iterrows():
    G.add_edge(row[0],row[1])

# Finding connected components in the graph
connected_components = list(nx.connected_components(G)) 

graph = pd.DataFrame(connected_components) 
graph.to_csv(sys.stdout, header=False, index=False, sep="\t", na_rep="NA")