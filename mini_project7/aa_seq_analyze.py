#!/usr/bin/python3
"""
Author: Ayesha Sanahari
Date: 22/Jan/2021
Analyze the aa_seq.fasta file and calculate the length, molecular weight, alanine percentage,
and glycine percentage of the sequence. Save the calculated parameters in a new text file called “aa_stats.txt”.

Input: aa_seq.fasta file
Output: “aa_stats.txt” file
"""

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

for seq_record in SeqIO.parse("aa_seq.fasta", "fasta"):
	myseq=str(seq_record.seq)
	
	#remove unwanted letters for make protein unambigous else it gives a error
	protein=myseq.replace('*','')
	protein=protein.replace('X','')
new = ProteinAnalysis(str(protein))
header=seq_record.id
weight = new.molecular_weight()
newfile=open("aa_stats.txt" , 'w') 
newfile.write(str(header)+"\n"+"length of amino acid :"+str(len(seq_record))+"\n"
		+"molecular weight of amino acid :"+str(weight)+"\n"+"alanine percentage :"+str(new.get_amino_acids_percent()['A'])+"\n"
		+"glycine percentage :"+str(new.get_amino_acids_percent()['G']))
newfile.close()

