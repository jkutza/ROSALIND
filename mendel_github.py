# k = number of homozygous dominant members in population
# m = number of heterozygous members in population
# n = number of homozygoues recessive members in population
k = 17
m = 19
n = 23

total = k + m + n

# Probability of mm
mm = (m/total) * ((m-1)/(total-1)) * 0.25 # Tt x Tt = 25% chance of tt offspring

# Probability of mn
mn = (m/total) * (n/(total-1)) * 0.50 # Tt x tt = 50% chance of tt offspring

# Probability of nm
nm = (n/total) * (m/(total-1)) * 0.50 # tt x Tt = 50% chance of tt offspring

# Probability of nn
nn = (n/total) * ((n-1)/(total-1)) # tt x tt = 100% change of tt offspring

# Sum of probabilities
sum = mm + mn + nm + nn

# Final output for probability of an offspring having a dominant allele
final_prob = "{:.5f}".format(1 - sum)

print(final_prob)