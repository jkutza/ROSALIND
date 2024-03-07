file_path = "rosalind_prot.txt"
with open(file_path, "r") as bioinf7:
  bioinf7_contents = bioinf7.read()

# Split original sequence into chunks
seq_length = len(bioinf7_contents)

split_sequence = []

for i in range(0, seq_length, 3): # Third variable specifies spacing between iterations
  split_sequence += [bioinf7_contents[i:i+3]]

# Translate to amino acid sequence by iterating over protein_sequence
protein_sequence = []

for codon in split_sequence:
  # U _ _
  if codon == "UUU" or codon == "UUC":
    protein_sequence += "F"
  if codon == "UUA" or codon == "UUG":
    protein_sequence += "L"
  if codon == "UCU" or codon == "UCC" or codon == "UCA" or codon == "UCG":
    protein_sequence += "S"
  if codon == "UAU" or codon == "UAC":
    protein_sequence += "Y"
  # if codon == "UAA" or codon == "UAG" or codon == "UGA":
    # protein_sequence += "stop"
  if codon == "UGU" or codon == "UGC":
    protein_sequence += "C"
  if codon == "UGG":
    protein_sequence += "W"

  # C _ _
  if codon == "CUU" or codon == "CUC" or codon == "CUA" or codon == "CUG":
    protein_sequence += "L"
  if codon == "CCU" or codon == "CCC" or codon == "CCA" or codon == "CCG":
    protein_sequence += "P"
  if codon == "CAU" or codon == "CAC":
    protein_sequence += "H"
  if codon == "CAA" or codon == "CAG":
    protein_sequence += "Q"
  if codon == "CGU" or codon == "CGC" or codon == "CGA" or codon == "CGG":
    protein_sequence += "R"

  # A _ _
  if codon == "AUU" or codon == "AUC" or codon == "AUA":
    protein_sequence += "I"
  if codon == "AUG":
    protein_sequence += "M"
  if codon == "ACU" or codon == "ACC" or codon == "ACA" or codon == "ACG":
    protein_sequence += "T"
  if codon == "AAU" or codon == "AAC":
    protein_sequence += "N"
  if codon == "AAA" or codon == "AAG":
    protein_sequence += "K"
  if codon == "AGU" or codon == "AGC":
    protein_sequence += "S"
  if codon == "AGA" or codon == "AGG":
    protein_sequence += "R"

  # G _ _
  if codon == "GUU" or codon == "GUC" or codon =="GUA" or codon =="GUG":
    protein_sequence += "V"
  if codon == "GCU" or codon == "GCC" or codon =="GCA" or codon == "GCG":
    protein_sequence += "A"
  if codon == "GAU" or codon == "GAC":
    protein_sequence += "D"
  if codon == "GAA" or codon == "GAG":
    protein_sequence += "E"
  if codon == "GGU" or codon == "GGC" or codon == "GGA" or codon == "GGG":
    protein_sequence += "G"

# Polish and print final protein string
final_string = ''.join(protein_sequence)
print(final_string)
