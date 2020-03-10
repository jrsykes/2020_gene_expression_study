
from Bio.SeqIO.FastaIO import SimpleFastaParser

count = 0
total_len = 0
with open("test.fastq") as in_handle:
    for title, seq in SimpleFastaParser(in_handle):
        count += 1
        total_len += len(seq)

print (count)