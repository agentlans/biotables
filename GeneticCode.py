from Bio.Seq import translate

amb_nt = dict() # Maps ambiguous symbol to nucleotides
nt_amb = dict() # Maps nucleotides to ambiguous symbol
amb_nt_file = open('AmbiguousNucleotides.txt', 'r')
for line in amb_nt_file:
  (symb, nts) = line.strip().split('\t')
  nt_list = frozenset(nts.split(','))
  # Update the dictionaries
  nt_amb[nt_list] = symb
  amb_nt[symb] = nt_list

amb_nt_file.close()

def cross(lst_of_lst):
  "Returns list of combinations of characters from lists."
  if len(lst_of_lst) == 1:
    return lst_of_lst
  elif len(lst_of_lst) == 2:
    return [a + b for a in lst_of_lst[0] for b in lst_of_lst[1]]
  else:
    return cross([lst_of_lst[0], cross(lst_of_lst[1:])])

# cross([['a','b'], ['2'], ['c','d']])

def summarize_nucleotides(nt_list):
  "Summarizes a list of nucleotides into one letter."
  return nt_amb[frozenset(nt_list)]

def expand_nucleotides(amb_str):
  "Expands a letter to a list of nucleotides."
  return cross([amb_nt[nt] for nt in amb_str])

# expand_nucleotides('ATN')

def summarize_aminoacid(aa_list):
  "Summarizes a list of amino acids into one letter."
  if len(aa_list) == 1:
    return aa_list[0]
  elif frozenset(aa_list) == frozenset(['D','N']):
    return 'B'
  elif frozenset(aa_list) == frozenset(['E','Q']):
    return 'Z'
  else:
    return 'X'

def letters_at(seqs, pos):
  "Returns the letters at given position in sequences."
  return [x[pos] for x in seqs]

def summarize_nt_seqs(seqs):
  "Summarizes a list of nucleotide sequences to one string."
  if seqs == []:
    return ""
  l = len(seqs[0])
  return ''.join([summarize_nucleotides(letters_at(seqs,i)) for i in range(l)])

#summarize_nt_seqs(reverse_trans['K'])

# Generate a list of possible codons
# Possible bases
#bases = ['A', 'T', 'G', 'C']
bases = list(amb_nt.keys())
codons = [b1 + b2 + b3 for b1 in bases for b2 in bases for b3 in bases]

# Translation
trans_file = open("GeneticCode/Translation.csv", "w")
for codon in codons:
  aa = translate(codon)
  if aa != 'X':
    trans_file.write(codon + ',' + aa + '\n')

trans_file.close()

# "Reverse translation"
reverse_trans = dict()
for codon in codons:
  aa = translate(codon)
  if aa in reverse_trans:
    reverse_trans[aa].append(codon)
  else:
    reverse_trans[aa] = [codon]

def is_nonambiguous(nts):
  for base in nts:
    if base not in ['A','T','G','C']:
      return False
  return True

def filter_nonambiguous(seqs):
  return [x for x in seqs if is_nonambiguous(x)]

def reverse_translate_ambiguous(aa1, aa2):
  return reverse_trans[aa1] + reverse_trans[aa2]

reverse_trans['Z'] = reverse_translate_ambiguous('E','Q')
reverse_trans['J'] = reverse_translate_ambiguous('L','I')
reverse_trans['B'] = reverse_translate_ambiguous('D','N')

rev_trans_file = open("GeneticCode/ReverseTranslation.csv", "w")
for aa, cd in reverse_trans.items():
  if aa == 'X':
    rev_trans_file.write('X,NNN\n')
  else:
    rev_trans_file.write(aa + ',' + summarize_nt_seqs(filter_nonambiguous(cd)) + '\n')

rev_trans_file.close()

