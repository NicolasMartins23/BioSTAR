BioSTAR (Bioinformatics Software for Targeted Analysis and Research) is my passion project.
It kickstarted my journey to learn how to program and transition from a Biochemistry Major to working with software.
The library is inteded to be as simple as possible to use and also complete complex tasks effectively.

I'm setting up an API for usage on this library to host it on my python website, and I decided to share and update the code on github too.

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
