How to use the BioSTAR library

from BioStar.nucleic_acids import DNA, RNA, OpenReadFrame
from BioStar.protein import Protein
from BioStar.tools.compare_sequence import CompareNucleotideSequence

Easily transform DNA to RNA and to Protein
DNA, Protein and RNA already have some built-in physical-chemical properties.

DNA, RNA and Protein constructors share one optional argument "sequence"
They must be passed as the object is instantiated (recommended approach)

There is a how_to_use.py on the directory which can be helpful, too!

There will be more changes added latter.
- Molecular mass and pH calculation
- Optimize sequence for expression in model organisms
- Identify possible signal peptide sequences
