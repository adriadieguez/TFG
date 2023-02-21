# Import modules
import sys
from collections import Counter
from Bio import SeqIO

# Read FASTA file
path = sys.argv[1]
sequences = [i for i in SeqIO.parse(path, 'fasta')]

# Read sequences
n_sequences = len(sequences) # Number of sequences
len_sequence = len(sequences[0].seq) # Length of sequences
v = []
for i in range (len_sequence):
    str = ''
    for k in range (n_sequences):
        str += sequences[k].seq[i]
    v.append(str)

def normalize(x, len_sequence):
    return round(x/len_sequence, 4)

counter = Counter(v) # Count ocurrences
freqs = dict(map(lambda x: (x[0], normalize(x[1],len_sequence)),
            counter.items())) # From ocurrences to frequences

# Display sorted results
sorted_freqs = dict(sorted(freqs.items(), key=lambda x:x[1], reverse=True))
[print(key,':',value) for key, value in sorted_freqs.items()]
