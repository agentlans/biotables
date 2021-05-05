from Bio.Seq import Seq

dna = Seq("ATGCYRWSKMDVHBN")
dna_comp = dna.complement()

comple_file = open('Complement.csv', 'w')
for base, comple in zip(dna, dna_comp):
  comple_file.write(base + ',' + comple + '\n')

comple_file.close()

