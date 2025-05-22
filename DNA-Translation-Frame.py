#Modify the DNA translation program to translate in each forward frame (1,2,3)

#Modify the DNA translation program to translate in each reverse (complement)
#translation frame too.

#Modify the translation program to handle 'N' symbols in the third position of a codon
#If all four codons represented correspond to the same amino-acid, then output that amino-acid.
#Otherwise, output 'X'.

import sys
table = sys.argv[1]
anthrax = sys.argv[2]

table = open('standard.code')
data = {}
for l in table:
    sl = l.split()
    key = sl[0]
    value = sl[2]
    data[key] = value    
table.close()

b1 = data['Base1']
b2 = data['Base2']
b3 = data['Base3']
aa = data['AAs']
st = data['Starts']

codons = {}
init = {}
n = len(aa)
for i in range(n):
    codon = b1[i] + b2[i] + b3[i]
    codons[codon] = aa[i]
    init[codon] = (st[i] == 'M')

def translate(seq, index):
    aaseq = []                                                                #empty list to store translated amino acids
    for i in range(index, len(seq), 3):                                       #iterate from the start of the index every 3
        if i + 2 < len(seq):                                                  #check if theres 3 bases for a full codon
            codon = seq[i:i+3]                                                #get codon
        if codon[2] == 'N':                                                   #checks if third nucleotide is 'N'
            possible_codons = [codon[0] + codon[1] + base for base in 'ACTG'] #if the third base is 'N', creates list of all possible codons by combining first 2 bases of codon with 'ACTG'
            amino_acids = set(codons.get(c, 'X') for c in possible_codons)    #hold amino acids corresponding to possible codons, if codon isn't found, will be 'X'
            if len(amino_acids) == 1:                                         #check if 1 unique amino acid in set
                aaseq.append(list(amino_acids)[0])                            #if it is, add it to aaseq
            else:                                                             #if there are many or no amino acids
                aaseq.append('X')                                             #assign 'X' to unrecognized codons
        else:                                                                 #if codon doesn't contain 'N'
            aa = codons.get(codon, 'X')                                       #looks up amino acid in dict
            aaseq.append(aa)                                                  #add amino acid to aaseq
    return ''.join(aaseq)

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def reverse_comp(seq):
    return ''.join(complement.get(base) for base in reversed(seq))

anthrax = open('anthrax.txt')
seq = ''.join(anthrax.read().split())
anthrax.close()

# Translate in all 3 forward frames
for frame in range(3):
    translated_seq = translate(seq, frame)
    print("Frame:", frame + 1, "Translation in Forward Frame:", translated_seq)

# Translate in all 3 reverse complement frames
rev_comp_seq = reverse_comp(seq)
for frame in range(3):
    translated_seq = translate(rev_comp_seq, frame)
    print("Frame:", frame + 1, "Translation in Reverse Complement Frame:", translated_seq)

start_codon = seq[:3] # iterate over every three items find start codon
if init.get(start_codon): # use built in function get to get the start codon
    print("True. The initial codon of the anthrax SASP gene is a valid translation start site.")
else:
    print("False. The initial codon of the anthrax SASP gene is NOT a valid translation start site")
