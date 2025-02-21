The web version of this project can be accessed here: https://softwaremartins.pythonanywhere.com/biostar

1 - PROJECT NAME
BioSTAR (Bioinformatics System for Targeted Analysis and Research)

2 - PROJECT DESCRIPTION
It kickstarted my journey to learn how to program and transition from a Biochemistry Major to working with software.
The library is inteded to be as simple as possible to use and also complete complex tasks effectively.

3 - MAIN USAGE
The 3 main classes on this library as DNA, RNA and Protein

All of those classes will have a property named sequence, which stores the genetic/peptide sequence as string. I also added a how_to_use.txt to make things more practical

DNA and RNA both inherit from the NucleicAcid parent class
get_peptide_sequence()
    Returns a string which represents the peptide sequence of a given NucleicAcid when translated intto a protein.
    If you wish to return a protein object instead, use the method to_protein()

to_protein()
    Returns a protein object of a given NucleicAcid when translated into a protein.
    If you wish to return a protein object instead, use the method to_protein()

When creating an instance of DNA("SEQUENCE") or RNA("SEQUENCE") you must pass the sequence when instantiating the new variable as a DNA object.
The recommended approach is to create a string variable first then pass it as an property, although it is only a matter of preference.

This is a summary of all functionality
DNA: peptide_sequence(), to_protein(), rna_sequence(), to_rna(), at_skew(), at_content(), gc_content(), gc_skew(), template_strand(), get_orf_map()
RNA: peptide_sequence(), to_protein(), dna_sequence(), to_dna()
Protein: aromacity(), charge_at_pH(), composition_ratio(), extinction_coefficient(), hydrophobic_index(), isoelectric_point(), molecular_weight(), secondary_structure_propensity()

There will be more changes added latter.
- Protein: pI and charge at pH, identify possible signal peptide sequences
- Nucleic Acids: Optimize sequence for expression in model organisms
