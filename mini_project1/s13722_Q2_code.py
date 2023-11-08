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
    # Remove blank lines
    if line != '\n':
        line = line.strip()
        # Remove the fasta header
        if('>') not in line:
            sequence = sequence + line
print(sequence)

# Define counters for each type of base
A_Counter = 0
T_Counter = 0
G_Counter = 0
C_Counter = 0

# For each letter/base in the sequence
for base in sequence:
    # if A was found increase the no of A bases by 1
    if base == 'A':
        A_Counter = A_Counter + 1
    # else if T was found increase the no of T bases by 1
    elif base == 'T':
        T_Counter = T_Counter + 1
    # else if G was found increase the no of G bases by 1
    elif base == 'G':
        G_Counter = G_Counter + 1
    # else if C was found increase the no of C bases by 1
    elif base == 'C':
        C_Counter = C_Counter + 1
    # else mention another letter found except A, T, G, C
    else:
        print('There is a character not belongs to A,T,G,C')


# Return the count of each base
print('No. of A bases in your sequence :', A_Counter)
print('No. of T bases in your sequence :', T_Counter)
print('No. of G bases in your sequence :', G_Counter)
print('No. of C bases in your sequence :', C_Counter)

# Return the total no. of bases
total = A_Counter + T_Counter + G_Counter + C_Counter
print('Total No. of bases in your DNA sequence :', total)

sequenceFile.close()
