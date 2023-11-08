#!/usr/bin/python3
"""
Author: Ayesha Sanahari
Date: 22/Jan/2021
Translate the given sequence in the FASTA file and save the translated amino acid
sequence in another FASTA file.
The FASTA header must contain the added word “translated” at the end.

Input: “mRNA_seq.fasta” file.
Output: aa_seq..fasta file
"""

from Bio import SeqIO
for record in SeqIO.parse("mRNA_seq.fasta", "fasta"):
	seq=record.seq
	# To perform biopython translate function perfectly, added Ns until the final sequence is a multiple of 3
	seq = seq + (len(seq) % 3 - 1) * 'N'
	aa=seq.translate()
	title=str(record.description)
	# remove the transcribed word at the end of header
	header=str(title[:-12])+"_translated"

with open("aa_seq.fasta", 'w') as newfile:
	newfile.write(">"+str(header)+"\n"+str(aa))




