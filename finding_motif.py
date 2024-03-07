file_path = "rosalind_subs.txt"
with open(file_path, "r") as subs:
  subs_contents = subs.read()

# Split by newline characters
subs_list = subs_contents.split("\n")

s = subs_list[0] # Full sequence
t = subs_list[1] # Motif
motif_indices = []

# Iterate over the full sequence, at each index seeing if the same length of the motif at that point matches the motif
for nt in range(len(s)):   # nt = numerical index, s[nt] = nucleotide base
  if s[nt:nt+len(t)] == t:
    motif_indices.append(nt+1)

# Format into final output with spaces instead of a list with commas. Since values are integers, .join() doesn't work
b = [str(motif_indices) for motif_indices in motif_indices]
print(" ".join(b))
