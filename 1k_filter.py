from Bio import SeqIO

import sys
inFile = sys.argv[1]
outFile = sys.argv[2]

long_sequences = []
with open (inFile, 'r') as f:

	for record in SeqIO.parse(f, "fasta"):
	    if len(record.seq) >= 1000 :
	        long_sequences.append(record)

SeqIO.write(long_sequences, outFile, "fasta")

