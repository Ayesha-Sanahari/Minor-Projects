"""
Author : Ayesha Sanahari
date : 16/Nov/2020
Calculating the length of a given DNA sequence
Input: DNA sequence in FASTA format
Output: sequence length
"""
sequence =''
# open/read sequence & store it in a variable
sequenceFile = open('sequence.fasta', 'r')
for line in sequenceFile:
    #To remove blank lines
    if line != '\n':
        line = line.strip()
        # Remove the fasta header
        if('>') not in line:
            sequence = sequence + line
print(sequence)

# Define counter
baseCounter = 0

# For each letter/base in the sequence
for base in sequence:
    # Increase the counter by 1
    baseCounter = baseCounter + 1

# Return the counter/ length
print('The length of your sequence :', baseCounter)

sequenceFile.close()
