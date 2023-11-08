"""
Author : Ayesha Sanahari
date : 27/Nov/2020
Differentiate two FASTA records & replace 'T' with 'U' in the mRNA sequence
Input: A FASTA file with two fasta records
Output: FASTA file containing only transcribed mRNA sequence
"""

# Define an empty dictionary to store the fasta sequences
seq_dict = {}
# Define an empty string variable to store the sequence
sequence =''

# open/read fasta file (with 2 fasta records) & store it in a variable
with open('OSDREB1A.fasta','r') as sequenceFile:
    for line in sequenceFile:

        # Remove blank lines
        if line != '\n':
            line = line.strip()

            # seperate the headers with '>' as the keys in the dictionary
            if ('>') in line:
                header = line
                #print(header)
                seq_dict[header]=''

            else:
                sequence = sequence + line

                # select the protein sequence as the value
                if 'M' in line:
                    seq_dict[header] = sequence
                # select only the mRNA sequence as the value and replace 'T' to 'U' only in mRNA sequence
                else:
                    seq_dict[header] += line.strip().replace('T','U')

print(seq_dict)

# open a file to write the transcribed mRNA sequence
with open("OSDREB1A_mRNA.fasta",'w') as output:
    # select only the transcribed mRNA sequence
    for item in seq_dict:
         if 'mRNA' in item :
             # add transcribed to the end of the mRNA header
             header1= item +' transcribed'

    # write the transcribed mRNA sequence to a new fasta file
    output = output.write(header1 + '\n' +seq_dict[item])
    # print the content to be written in the file (header with transcribed mRNA sequence)
    print(header1 + '\n' +seq_dict[item])


