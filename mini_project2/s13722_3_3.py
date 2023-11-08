"""
Author : Ayesha Sanahari
date : 01/Dec/2020
Find the degree of DREB1A protein using networkx package
Input: string_interactions.tsv file containing information about protein interactions
Output: degree of DREB1A protein
"""

import networkx as nx;

# Define an empty dictionary to store interactions
dict={}
# define the key of dictionary
num=1
#open/read file and store in a variable
with open('string_interactions.tsv', 'r') as file:
    for line in file:
        if '\n' != line and 'node1' not in line:
            interaction = line.strip().split('\t')
            dict[num]= interaction[0:2]
            num +=1
print(dict)

G = nx.Graph()
for key, values in dict.items():
    G.add_edge(values[0], values[1])

print("number_of_edges: " ,G.number_of_edges())
print("degree of the rice DREB1A protein: ", G.degree['ERF24'])









