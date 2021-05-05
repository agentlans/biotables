# Biotables

Scripts for generating fixed lookup tables for bioinformatics and molecular biology

Requires

- Biopython
- Biostrings, tidyr, readr packages for R

Run `Make.sh` to generate the following CSV files.
The lines in each file have the formats listed below.

- `Complement.csv` - nucleotide base, its complement
- `GeneticCode/Translation.csv` - codon, amino acid
- `GeneticCode/ReverseTranslation.csv` - amino acid, codon
- `Matrices/*.csv` - amino acid 1, amino acid 2, score

Note that the "reverse translation" is one way and the converse may not be true.
For example, L can only be translated from a YTN codon ({C,T}T{A,T,G,C}) but not
all YTN codons translate to L, such as TTC which translates to F.

## Author, License

Copyright 2021 Alan Tseng

GNU General Public License v3

The generated CSV files are public domain licensed under Creative Commons Zero v1.0 Universal

