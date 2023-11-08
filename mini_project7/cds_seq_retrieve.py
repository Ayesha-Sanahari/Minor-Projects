#!/usr/bin/python3

"""
Author: Ayesha Sanahari
Date: 22/Jan/2021
Retrieve the GenBank record for an accession number (with version)
of a protein coding DNA sequence or a reverse-transcribed mRNA complement
given by the user and save the sequence in FASTA format (cds_seq.fasta).

Input: Accession number (with version)
Output: cds_seq.fasta file
"""
from Bio import Entrez

#getting Accesion number
AcNum = input("enter the accession number with version")
#retrive fasta file
Entrez.email = "Name@example.org"
handle = Entrez.efetch(db="Nucleotide", id=AcNum, rettype="fasta", retmode="text")

#write sequence into new fasta file
with open("cds_seq.fasta",'w') as newfile :
		newfile.writelines(handle)

