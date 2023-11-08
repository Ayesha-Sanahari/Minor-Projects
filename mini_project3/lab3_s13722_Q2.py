"""
Author : Ayesha Sanahari
date : 20/Dec/2020
Familiarizing with Python Object-oriented programming (OOP) techniques.
Writing custom Python methods to analyze DNA sequences

Eg: Calculating the AT content of a given DNA or mRNA sequence
     Input: DNA/mRNA sequences in FASTA format
     Output: AT content of the given DNA or mRNA sequence
"""
class Sequence:
    sequence_count = 0

    def __init__(self, sequence, Gene_name, Gene_ID, Species_name, Subsp_name):
        self.Gene_ID = Gene_ID
        self.Gene_name = Gene_name
        self.Seq_Type = Sequence.get_Seq_Type(sequence)
        self.Seq_length = len(sequence)
        self.Species_name = Species_name
        self.Subsp_name = Subsp_name
        Sequence.sequence_count += 1

    @staticmethod
    def fasta_Split(file_name):
        """
        A static method to split multiple FASTA sequences in a single text file
        and return a dictionary containing the Gene name
        (first item in the hyphen-separated list) as the key and
        a list containing hyphen-separated fields in the header plus the sequence as the value.

        Input: Argument - file_name
        Output: dictionary containing keys and values
        """

        seq_dict = {}
        with open(file_name, 'r') as sequenceFile:

            for line in sequenceFile:
                # Remove blank lines
                if line != '\n':
                    line = line.strip()

                    # seperate the headers with '>' as the keys in the dictionary
                    if ('>' in line) :
                        #the fasta header includes the ">" sign too, so to remove it
                        key_list = line.strip('>').split('-')
                        key = key_list[0]
                        seq_dict[key] = []
                        sequence = ''

                    else:
                        sequence += line
                        seq_dict[key] = [sequence]

                        for i in key_list:
                            # seq_dict[key].insert(0,i)
                            seq_dict[key].append(i)
        return (seq_dict)

    def get_Seq_Type(sequence):
        """
        A method to check the sequence type to
        distinguish between all 3 sequence types: DNA, mRNA, amino acid.

        Input: Argument - sequence
        Output: sequence type (DNA or mRNA or amino acid)
        """

        amino_acid = ['M','N','I','R','K','Q','E','S','P','L','O','H','T','V','W','D','Y']

        if 'U' in sequence:
            return 'mRNA'
        else:
            for letter in sequence:
                if letter in amino_acid:
                    return 'Amino acid'
                else:
                    return 'DNA'


    def get_Character_Count(sequence):
        """
        this method should return a dictionary of character counts with
        each character as the key and count as the value.
        A character can be a nucleobase or an amino acid.

        Input: Argument - sequence
        Output: dictionary of character counts
        """
        character_dic ={}               # initialize a dictionary

        for ch in sequence:             #read through each character in sequence
            if ch in character_dic:
                character_dic[ch] += 1   #if it's been seen before, increment counter
            else:
                character_dic[ch] = 1    #otherwise, insert it into the dictionary

        return character_dic


class DNAseq(Sequence):
    def __init__(self, sequence, Gene_name, Gene_ID, Species_name, Subsp_name):
        super().__init__(sequence, Gene_name, Gene_ID, Species_name, Subsp_name)
        self.AT_content = self.get_ATcontent(sequence)
        self.Transcribed_sequence = self.transcribe_Sequence(sequence)


    def get_ATcontent(self, sequence):
        length = len(sequence)
        A_count = sequence.upper().count('A')
        T_count = sequence.upper().count('T')
        AT_content = (A_count + T_count) / length
        return round(AT_content, 3)

    def transcribe_Sequence(self, sequence):
    # transcribe the given DNA sequence into its mRNA sequence and
    # store it in the Transcribed sequence instance variable.
        Transcribed_sequence = sequence.replace('T','U')
        return Transcribed_sequence


class MRNAseq(Sequence):
    Amino_acid_codons = ''

    def __init__(self, sequence, Gene_name, Gene_ID, Species_name, Subsp_name ):
        super().__init__(sequence, Gene_name, Gene_ID, Species_name, Subsp_name)
        self.AT_content = self.get_ATcontent(sequence)
        self.Translated_sequence = self.translate_Sequence(sequence)

    def get_ATcontent(self,sequence):
        length = len(sequence)
        A_count = sequence.upper().count('A')
        T_count = sequence.upper().count('U')
        AT_content = (A_count + T_count) / length
        return round(AT_content, 3)

    @ classmethod
    def upload_Codons(cls, text_file):
        """
        This class method should create and return a dictionary to store codon-amino acid pairs from a text file
        Input: Argument- text_file
        Output: dictionary which stores codon-amino acid pairs
        """
        cls.Amino_acid_codons = text_file

        #Define an empty dictionary to store codons & amino acid letters
        Codon_map = {}
        # Open/read codon_table and store it in a variable
        with open(text_file, 'r') as Codon_table:
            for line in Codon_table:
                # Remove header and empty lines in Codon table
                if '#' not in line and line != '\n':
                    (codon, amino_acid, Letter, FullName) = line.strip().split('\t')
                    # make a dictionary from codons as key and amino acid letters as values
                    Codon_map[codon] = Letter
        return (Codon_map)

    def translate_Sequence(self, sequence):
        """
        This method is translating a given mRNA sequence into its amino acid sequence.
        (Consider only its first reading frame)
        Input: argument- mRNA sequence
        Output: translated Amino acid sequence
        """
        protein = ''
        #Get codon map from the defined method upload_Codons
        Codon_map = MRNAseq.upload_Codons('codon_table.txt')

        # Define the position of the sequence
        position = 0
        # iterate through the sequence
        while position <= len(sequence)-2:
            # divide the sequence into codons by three bases assuming 1st reading frame
            codon1 = sequence[position:position + 3]

            # go through the dictionary and find relevant amino acid letters
            if codon1 in Codon_map.keys():
                protein += Codon_map[codon1]

            position = position + 3

        # Find the position of aa relevant to start codon and end codon
        AA_No_BeforeMeth = protein.find('M')   #No. of amino acids before Methionine
        AA_No_Before_O = protein.find('O')     #No. of amino acids before the first stop codon

        # Assuming translation is starting with Methionine, get the translated protein seq. and length
        protein1 = protein[AA_No_BeforeMeth:AA_No_Before_O]
        return (protein1)


class Proteinseq(Sequence):

    def __init__(self, sequence, Gene_name, Gene_ID, Species_name, Subsp_name, Uniprot_ID, Reviewed_status):
        super().__init__(sequence, Gene_name, Gene_ID, Species_name, Subsp_name)

        self.Uniprot_ID = Uniprot_ID
        self.Reviewed_status = Reviewed_status
        self.Hydrophobicity = self.get_Hydrophobicity(sequence)

    def get_Hydrophobicity(self, sequence):
        """
        this method should return the percentage of the total hydrophobic amino acid residues
        (A, I, L, M, F, W, Y, V) in the sequence & update the Hydrophobicity object variable

        Input: argument- protein sequence
        Output: percentage of the total hydrophobic AA residues
        """
        count = 0
        hydrophobic_AA = ['A','I','L','M','F','W','Y','V']

        for letter in sequence:
            if letter in hydrophobic_AA:
                count += 1
                hydrophobicity = (count / len(sequence))*100
            else:
                hydrophobicity = 0

        return round(hydrophobicity,2)





