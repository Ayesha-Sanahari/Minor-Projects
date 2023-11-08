"""
Author : Ayesha Sanahari
date : 15/Jan/2020
Predict the majority voting scores for all the unknown protein members
for the stress tolerance biological process of the ATDREB2A network

Input:
•	text file containing a list of known Arabidopsis thaliana proteins
annotated to the particular function – “AT_stress_proteins.txt”
•	unknown protein members for the stress tolerance biological process of the
ATDREB2A network (“string_interactions_1.tsv” file in tabular format)

Output: list of unknown proteins with the predicted majority voting score
"""
import networkx as nx;

# Read the file containing the list of known stress proteins
stress_proteins=[]
with open('AT_stress_proteins.txt', 'r') as file:
    for line in file:
        if '\n' != line :
            stress_proteins_all = line.strip().split('\t')
            # Extract the relevant information & make a list called stress_proteins
            stress_proteins += stress_proteins_all[1:2]

print("stress_proteins:", stress_proteins)


# Read the tsv file which contains the protein-protein interactions of unknown proteins
tsv_data=[]
with open('string_interactions_1.tsv', 'r') as file:
    for line in file:
        if '\n' != line and 'node1' not in line:
            tsv_data_all = line.strip().split('\t')
            # Extract the relavant information & save it into a list called tsv_data
            tsv_data += tsv_data_all[0:2]

print("tsv_data:", tsv_data)

# Define an empty dictionary to store interactions
tsv_dict={}
# define the key of dictionary
num=1
#open/read file and store in a variable
with open('string_interactions_1.tsv', 'r') as file:
    for line in file:
        if '\n' != line and 'node1' not in line:
            interaction = line.strip().split('\t')
            tsv_dict[num]= interaction[0:2]
            num +=1

print(tsv_dict)

# Convert the names of proteins into  UPPERCASE
Capilalize_tsv_dict = {key: [x.upper() for x in tsv_dict[key]] for key in tsv_dict}
print("Capilalize_tsv_dict:", str(Capilalize_tsv_dict))


# Create an interaction graph of the unknown proteins from the tsv_dict
G = nx.Graph()
for key, values in Capilalize_tsv_dict.items():
    G.add_edge(values[0], values[1])

# Convert the names of proteins into  UPPERCASE
unique_tsv_data = set(x.upper() for x in tsv_data)
unique_stress_proteins = set(x.upper() for x in stress_proteins)

print("unique_stress_proteins:", unique_stress_proteins)
print("unique_tsv_data:",unique_tsv_data)

print("length unique_tsv_data:",  len(unique_tsv_data))
print("length tsv_data:",  len(tsv_data))
print("length unique_stress_proteins:",  len(unique_stress_proteins))
print("length stress_proteins:",  len(stress_proteins))

# Seperate the proteins in the tsv_data whose functions aren’t already know (which are not belong to stress_proteins ) as unknown_proteins
unknown_proteins = unique_tsv_data - unique_stress_proteins
print("\n Difference, unknown_proteins:", unknown_proteins)


# For each of these unknown proteins, determine how many of their neighbors are known proteins
unknown_dict ={}
score =0
for protein in unknown_proteins:
    for neighbour in G.neighbors(protein):
        # Assign them a score based on the number of known neighbours. ( Predict majority voting scores)
        if neighbour in unique_stress_proteins:
            score += 1
            unknown_dict[protein] = score

print("unknown_dict")
print(unknown_dict)

# Sort the unknown proteins in descending order based on the scores, with proteins with high scores at the top.
from collections import OrderedDict

import operator
desending = sorted(unknown_dict.items(),key=operator.itemgetter(1), reverse=True)


print()
print("desending:", desending)
print("len desending:", len(desending))


# Write the ordered list to an output file.
ordered_list = str(desending)
with open('Ordered_list.txt', 'w') as output:
        output = output.write(ordered_list)

print("Descending ordered list:", ordered_list)

print("\n No. of unknown proteins in the network for stress tolerance: " ,len(unknown_proteins))
print("degree of the ATDREB2A protein: ", G.degree['DREB2A'])





