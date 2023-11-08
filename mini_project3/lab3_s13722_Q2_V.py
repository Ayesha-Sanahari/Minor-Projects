from lab3_s13722_Q2 import *

seq_set = Sequence.fasta_Split("OSDREB_sequences.fasta")
print(seq_set)

DREB1A_P = seq_set['DREB1A P']
DREB1B_P = seq_set['DREB1B P']
DREB2A_P = seq_set['DREB2A P']
DREB2B_P = seq_set['DREB2B P']
DREB1A_CDS = seq_set['DREB1A CDS']
DREB1B_CDS = seq_set['DREB1B CDS']
DREB2A_CDS = seq_set['DREB2A CDS']
DREB2B_CDS = seq_set['DREB2B CDS']

#Creating objects manually from the seq_set dictionary
obj_DREB1A_P =Proteinseq( DREB1A_P[0],DREB1A_P[1],DREB1A_P[2],DREB1A_P[3],DREB1A_P[4],DREB1A_P[5],DREB1A_P[6])
obj_DREB1B_P =Proteinseq( DREB1B_P[0],DREB1B_P[1],DREB1B_P[2],DREB1B_P[3],DREB1B_P[4],DREB1B_P[5],DREB1B_P[6])
obj_DREB2A_P =Proteinseq( DREB2A_P[0],DREB2A_P[1],DREB2A_P[2],DREB2A_P[3],DREB2A_P[4],DREB2A_P[5],DREB2A_P[6])
obj_DREB2B_P =Proteinseq( DREB2B_P[0],DREB2B_P[1],DREB2B_P[2],DREB2B_P[3],DREB2B_P[4],DREB2B_P[5],DREB2B_P[6])

obj_DREB1A_CDS = DNAseq (DREB1A_CDS[0],DREB1A_CDS[1],DREB1A_CDS[2],DREB1A_CDS[3],DREB1A_CDS[4])
obj_DREB2B_CDS = DNAseq (DREB2B_CDS[0],DREB2B_CDS[1],DREB2B_CDS[2],DREB2B_CDS[3],DREB2B_CDS[4])



print('i.')
print('Gene ID :', obj_DREB1A_CDS.Gene_ID)
print('sequence length :', obj_DREB1A_CDS.Seq_length)
print('sequence type :', obj_DREB1A_CDS.Seq_Type)
print('AT content :',obj_DREB1A_CDS.AT_content ,'\n' )

print('ii.')
mRNA_sequence = obj_DREB2B_CDS.Transcribed_sequence
obj_DREB2B_mRNA = MRNAseq( mRNA_sequence, DREB2B_CDS[1],DREB2B_CDS[2],DREB2B_CDS[3],DREB2B_CDS[4])

print('sequence length :', obj_DREB2B_mRNA.Seq_length)
print('sequence type :', obj_DREB2B_mRNA.Seq_Type)
print('AT content :',obj_DREB2B_mRNA.AT_content)
print('Transcribed_mRNA sequence of OSDREB2B :',obj_DREB2B_CDS.Transcribed_sequence, '\n')

print('iii.')
Protein = obj_DREB2B_mRNA.Translated_sequence
print('Translated_sequence of OSDREB2B mRNA :',Protein)
print('length :', len(Protein), 'aa', '\n')

print('iv.')
print('Uniprot_ID :', obj_DREB2A_P.Uniprot_ID)
print('Reviewed status :', obj_DREB2A_P.Reviewed_status)
print('Type :',obj_DREB2A_P.Seq_Type)
print('amino acid composition :', Sequence.get_Character_Count(DREB2A_P[0]))
print('Hydrophobicity ;', obj_DREB2A_P.Hydrophobicity , '\n')

print('v.')
print('number of sequences created :', Sequence.sequence_count)
