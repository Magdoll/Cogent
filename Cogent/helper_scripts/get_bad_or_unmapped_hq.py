import os, sys
from Bio import SeqIO

input  = sys.argv[1] # ex: hq.fasta
prefix = sys.argv[2] # ex: hq.no5merge

group_file = prefix + '.collapsed.group.txt'
out_file   = prefix + '.badmapped.fasta'

if not os.path.exists(group_file):
    print >> sys.stderr, "Collapsed group file {0} does not exist! Abort!".format(group_file)

good = set()
for line in open(group_file):
    a,b=line.strip().split()
    for x in b.split(','): good.add(x)

fa_or_fq = 'fasta' if input.upper().split('.')[-1] in ('FA', 'FASTA') else 'fastq'
f = open(out_file, 'w')
print >> sys.stderr, "Reading transcript file {0}....".format(input)
for r in SeqIO.parse(open(input),fa_or_fq):
    if r.id not in good: 
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
f.close()

print >> sys.stderr, "Unmapped or badly mapped transcripts written to: {0}".format(out_file)
