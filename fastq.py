#!/usr/bin/env python

import sys
import string
from itertools import islice


def read_fastq(fn):
    with open(fn) as infile:
        while True:
            lines_gen = islice(infile, 4)
            lines=list(lines_gen)
            if len(lines) == 4:
                yield lines
            else:
                return



a_seq=list()
a_annot=list()
for line4 in read_fastq(sys.argv[1]):
    annot=line4[0].strip()
    seq=line4[1].strip()
    a_seq.append(seq.translate(string.maketrans('acgt','ACGT')))
    a_annot.append(annot)

# kmers

'''
def __init__(self, ...):
    self.trans = string.maketrans('TAGCtagc', 'ATCGATCG')

def reverse_complement(seq):
    return seq[::-1].translate(self.trans)

k=30
kmers=dict()
s=0
for seq in a_seq:
    s+=len(seq)
    for i in range(len(seq)-k):
       kmer=seq[i:i+k]
       if kmer in kmers:
           kmers[kmer]+=1
       else:
           kmers[kmer]=1

print 'total kmers:', len(kmers)
print 'total nucl:', s
'''
k=30
for annot, seq in zip(a_annot, a_seq):
    kmers=set()
    for i in range(len(seq)-k):
       kmer=seq[i:i+k]
       kmers.add(kmer)
    print annot
    print len(kmers)
    for kmer in kmers:
       print kmer
