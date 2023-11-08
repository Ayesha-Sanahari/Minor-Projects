"""
Author : Ayesha Sanahari
date : 27/Nov/2020
Translate mRNA sequence into an Amino acid sequence
Input: transcribed mRNA sequence fasta file and codon_table
Output: Amino acid sequence
"""
#Define an empty dictionary to store codons & amino acid letters
Codon_map = {}
# Open/read codon_table and store it in a variable
with open('codon_table.txt', 'r') as Codon_table:
    for line in Codon_table:

        # Remove header and empty lines in Codon table
        if '#' not in line and line != '\n':
            (codon, amino_acid, Letter, FullName) = line.strip().split('\t')
            # make a dictionary from codons as key and amino acid letters as values
            Codon_map[codon] = Letter
print(Codon_map)
protein = ''
seq = ''
Header = ''

# open/read mRNA and store in a variable
with open('OSDREB1A_mRNA.fasta', 'r') as sequence:
    for line in sequence:
        # Remove empty lines
        if line != '\n':
            line = line.strip()
            # remove header
            if '>' not in line:
                seq = seq + line
        # get header to write in protein file
        if '>' in line:
            Header = Header + line
            Header = Header.replace('mRNA transcribed','translated amino acid sequence')

print(seq)
print("Length of transcribed mRNA seq: ", len(seq),"bp")

# Define the position of the sequence
position = 0
# iterate through the sequence
while position  <= len(seq):
    #divide the sequence into codons by three bases assuming 1st reading frame
    codon1 = seq[position:position+3]
    #print(codon1)

    # go through the dictionary and find relevant amino acid letters
    if codon1 in Codon_map.keys():
        protein += Codon_map[codon1]

    position = position + 3

print(protein)
print("Length of amino acid seq: ",len(protein),'aa')

# Find the position of aa relevant to start codon and end codon
AA_No_BeforeMeth= protein.find('M')
AA_No_Before_O= protein.find('O')
print("No. of amino acids before Methionine: ", AA_No_BeforeMeth)
print("No. of amino acids before the first stop codon: ", AA_No_Before_O)

# Assuming translation is starting with Methionine, get the translated protein seq. and length
protein1 = protein[AA_No_BeforeMeth:AA_No_Before_O]
print(protein1)
print("Assuming translation is starting with Methionine",'\n' ,"Length of translated amino acid sequence is: ", len(protein1), "aa")


# Assuming translation is starting with first letter, get the translated protein seq. and length
print(Header)
print(protein[:AA_No_Before_O])
print("Assuming translation is starting with first letter",'\n' ,"Length of translated amino acid sequence is: ", len(protein[:AA_No_Before_O]), "aa")

# write the translated protein into a fasta file
with open("Translated_mRNA.fasta", 'w') as output:
    output = output.write(Header + '\n' + protein[:AA_No_Before_O])

