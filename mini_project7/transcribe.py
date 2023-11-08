#!/usr/bin/python3
"""
Author: Ayesha Sanahari
Date: 22/Jan/2021
Transcribe the given sequences in the FASTA file and save the
transcribed mRNA sequence in another FASTA file.
The FASTA header must contain the added word “transcribed” at the end.

Input: “cds_seq.fasta” file.
Output: mRNA_seq.fasta file
"""
from Bio import SeqIO

for record in SeqIO.parse( "cds_seq.fasta" , 'fasta'):
	seq=record.seq
	seq=seq.transcribe()
	title= record.description
	header=str(title)+"_transcribed"

# write transcribed sequence into a new fasta file
newfile = open( "mRNA_seq.fasta" , 'w')
newfile.write(">"+str(header)+"\n"+str(seq))

newfile.close()


