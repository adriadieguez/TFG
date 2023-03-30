# Import modules
import sys
from collections import Counter
from Bio import SeqIO
import itertools

# Read FASTA file
path = sys.argv[1]
sequences = [i for i in SeqIO.parse(path, 'fasta')]

# Read sequences
n_sequences = len(sequences) # Number of sequences
len_sequence = len(sequences[0].seq) # Length of sequences
v = []
for i in range (len_sequence):
    s = ''
    for k in range (n_sequences):
        s += sequences[k].seq[i]
    v.append(s)

def normalize(x, len_sequence):
    return round(x/len_sequence, 4)

# Compute all possible alignments
alphabet = [0, 1, 2, 3]
combinations = itertools.product(alphabet, repeat=n_sequences)
sorted_combinations = sorted(combinations)
for i in range(len(sorted_combinations)):
    sorted_combinations[i] = ''.join([str(s) for s in sorted_combinations[i]])


counter = Counter(v) # Count ocurrences
for i in range(len(sorted_combinations)):
    if sorted_combinations[i] not in counter:
        counter[sorted_combinations[i]] = 0
freqs = dict(map(lambda x: (x[0], normalize(x[1],len_sequence)),
            counter.items())) # From ocurrences to frequences

# Display sorted results
sorted_freqs = dict(sorted(freqs.items(), key=lambda x:x[1], reverse=True))
[print(key,':',value) for key, value in sorted_freqs.items()]
