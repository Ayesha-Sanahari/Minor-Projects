#!/usr/bin/python3

"""
Author: Ayesha Sanahari
Date: 22/Jan/2021
Retrieve the GenBank records for the given accession Numbers
"""
from Bio import Entrez

seq = ["AAK43967.1", "AED90870.1", "NP_567720.1", "AAK59861.1"]
Entrez.email = "A.N.Other@example.com"
for sequences in seq:
    handle = Entrez.efetch(db="Protein", id=sequences, rettype="fasta", retmode="text")
    file = str(sequences) + ".fasta"
    with open(file, 'w') as newfile:
        newfile.writelines(handle)
