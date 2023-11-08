"""
Author : Ayesha Sanahari
date : 16/Nov/2020
Calculating the length of a given DNA sequence by len() function
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


# Return the length
length = len(sequence)
print('The length of your DNA sequence  : ' ,length)

#close the file
sequenceFile.close()
