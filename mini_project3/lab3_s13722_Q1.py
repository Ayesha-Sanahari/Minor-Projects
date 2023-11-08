"""
Author : Ayesha Sanahari
date : 06/Dec/2020
Writing custom Python methods to analyze DNA sequences

Calculating the AT content of a given DNA or mRNA sequence
Input: DNA/mRNA sequences in FASTA format
Output: AT content of the given DNA or mRNA sequence
"""
class Sequence:

    def get_AT_content(dna):
        length = len(dna)
        A_count = dna.upper().count('A')
        T_count = dna.upper().count('T')
        AT_content = (A_count + T_count) / length
        return round(AT_content, 3)

    @staticmethod
    def splitingFastaFile(file_name):

        # Define an empty dictionary to store the fasta sequences
        seq_dict = {}

        # open/read fasta file (with many fasta records) & store it in a variable
        with open(file_name,'r') as sequenceFile:
            for line in sequenceFile:
                # Remove blank lines
                if line != '\n':
                    line = line.strip()

                    # seperate the headers with '>' as the keys in the dictionary
                    if ('>') in line:
                        header = line
                        #print(header)
                        seq_dict[header]=''
                    # if line not begin with '>', take the sequences as the values of dictionary
                    else:
                        seq_dict[header] += line.strip()
        return seq_dict

    def getType(sequence):
        if 'M' in sequence:
            return "Amino acid"
        elif 'U' in sequence:
            return "mRNA"
        else:
            return "DNA"


dic = Sequence.splitingFastaFile("OSDREB_sequences.fasta")
print(dic)

for key,value in dic.items():
    print(key)
    print('Type : ',Sequence.getType(value))

    if 'M' not in value:
        print('AT_content : ', Sequence.get_AT_content(value),'\n')
    else:
        print('AT_content : Not Found \n')











