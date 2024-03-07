file_path = "rosalind_gc.txt"
with open(file_path, "r") as bioinf5:
  bioinf5_contents = bioinf5.read()

############# Getting names and Sequences Into Their Own Lists #################

# Empty lists for names and sequences
all_names = []
all_sequences = []

new_entries = bioinf5_contents.splitlines()

temp_seq = []

# iterate over new list and add to either all_names or all_sequences
for i in range(len(new_entries)-1):

  # If the first character is ">", just add item to all_names
  if new_entries[i][0] == ">":
    all_names += [new_entries[i]] # Use brackets to add a string to a list, instead of breaking it up into letters

  # If the first character is NOT ">"
  if new_entries[i][0] != ">":

    # If the previous first character in the next item down is also not ">"
    if new_entries[i+1][0] != ">":
      temp_seq += [new_entries[i]]

    # If the first character in the next item down IS ">"
    if new_entries[i+1][0] == ">":
      temp_seq += [new_entries[i]]
      joined_string = ''.join(temp_seq)
      all_sequences += [joined_string]
      temp_seq = []

# Add the last line to temp_seq, then add that to all_sequences
temp_seq += [new_entries[len(new_entries)-1]]
last_joined_string = ''.join(temp_seq)
all_sequences += [last_joined_string]

print("All Names: " + str(all_names))
print("All Sequences: " + str(all_sequences))
print("------------------------------")

#################### Calculating GC Content ####################################
# First getting GC percentages into their own list
gc_percentages = []

# iterate over all sequences
for seq in range(len(all_sequences)):
  counter = 0

  # iterate over each nt per sequence
  for nt in all_sequences[seq]:
    if nt == "C" or nt == "G":
      counter += 1

  gc_percent = (counter / len(all_sequences[seq])) * 100
  gc_percentages += [gc_percent]

print("All GC Percentages: " + str(gc_percentages))

# Find the max and index it
max_gc_percent = max(gc_percentages)
index_of_max = gc_percentages.index(max_gc_percent)

print("Max GC Percentage is: " + str(max_gc_percent))
print("The Index of the Max is: " + str(index_of_max))
print("------------------------------")

# Final output
print(all_names[index_of_max][1:len(new_entries[0])])
rounded_max = "{:.6f}".format(max_gc_percent)         # Rounds to 6 decimal places
print(rounded_max)
