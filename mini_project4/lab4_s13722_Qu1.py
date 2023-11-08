"""
Author : Ayesha Sanahari
date : 15/Jan/2021
Use built-in or third-party Python modules and packages (Biopython and re module)
to solve biological questions

Input: ATdreb2a.fasta file containing DREB2A gene sequence
Output: blast result as “dreb2a_blast.xml”
        hits that are below the threshold E value 0.05
        number of blast hits with ABRE element
"""
from Bio import SeqIO
import re

# Load the downloaded FASTA file as a sequence record
for seq_record in SeqIO.parse("ATdreb2a.fasta", "fasta"):
    print("sequence ID :", seq_record.id)
    print("Description :", seq_record.description)
    print(repr(seq_record.seq))
    print("Seq length: ", len(seq_record))

from Bio.Blast import NCBIWWW

record = SeqIO.read("ATdreb2a.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

print('aye')

# Write the BLAST result into an xml file
with open("dreb2a_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

result_handle.close()
result_handle = open("dreb2a_blast.xml")

from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)
#blast_records = NCBIXML.parse(result_handle)

#Define a counter for count the number of blast hits with ABRE element
hit_count =0
E_VALUE_THRESH = 0.05

# check for each blast hit sequence for the presence of the ABRE element
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:

            print("\n ****Alignment****")
            print("blast hit title:", alignment.title)
            print("alignment length:", alignment.length)
            print("E value:", hsp.expect)
            print("score:", hsp.score)
            print("hit/subject sequence:", hsp.sbjct)
            print("hit sequence length:", len(hsp.sbjct))

            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")

            # Finding the ABRE elements present in the sequence
            # print the detected sequence fragment (e.g., ACGTGC or ACGTTC) along with the sequence location of each finding.
            # (Please note that the ABRE element can be found in multiple locations in the same sequence) So used re.finditer
            pattern = re.compile("ACGT[GT]C")
            matches = re.finditer(pattern, hsp.sbjct)

            ay = "ABRE elements : Not present"
            for item in matches:
                hit_count += 1
                ay= "ABRE elements : " , item.group() , item.span()

            print(ay)
print("\n number of blast hits with ABRE element: ", hit_count)


