'''
How to use the BioSTAR library
'''

from BioStar.nucleic_acids import DNA, RNA, OpenReadFrame
from BioStar.protein import Protein
from BioStar.tools.compare_sequence import CompareNucleotideSequence

'''
Easily transform DNA to RNA and to Protein
DNA, Protein and RNA already have some built-in physical-chemical properties.

DNA, RNA and Protein constructors share one optional argument "sequence"
They must be passed as the object is instantiated (recommended approach)

'''
dna_seq = "ATGCACGAGGGGGAGGCTTACTTATGGACCATTTAAAGGGATCTTTAG"
dna = DNA(sequence=dna_seq)
# print(dna.sequence)
rna = dna.to_rna()
# print(rna.sequence)
protein = dna.to_protein()
# print(protein.sequence)
# print(protein.get_aromacity())

'''
Get an complete mapping of all possible reading frames
for any given dna sequence.
The positive (5' -> 3') frames are marked as "+" "++" and "+++"
The negative (3' -> 5') frames are marked as "-" "--" and "---"

Additionaly you can set up a threshold to consider only sequences of a certain peptide count.
'''
orf_map = dna.get_orf_map()
# print(orf_map)

orf_map_with_threshold = dna.get_orf_map(length_threshold=25)
# print(orf_map_with_threshold)

'''
Get a mutation map when comparing between two dna sequences too.
Additionaly you can use the property "DNA().sequence"
if you don't want to pass the raw string.
Although both approaches will work,
I'd recommend using object properties for clarity.

The show_mutations_only parameter is by default 

'''
dna_seq2 = "ATGCACGAGGGGGAGGCTTACTTATGGACCATTTAGAGGGATCTTTAG"
dna2 = DNA(dna_seq2)

compare_sequence = CompareNucleotideSequence(dna.sequence, dna2.sequence)
print(compare_sequence.compare(show_only_mutations=True))


