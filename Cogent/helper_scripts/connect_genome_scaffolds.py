# IPython log file

import os, sys
from collections import defaultdict
import GFF


def process_contig(dirname, f):
    d = defaultdict(lambda: []) # path_number --> (start, end, scaffold)
    reader = GFF.gmapGFFReader(os.path.join(dirname,'aloha2.fa.cuttlefish.gff'))
    for r in reader:
        if r.strand == '+': s,e = r.seq_exons[0].start, r.seq_exons[-1].end
        else: s,e = r.seq_exons[-1].start, r.seq_exons[0].end
        d[r.seqid].append((s, e, r.chr, r.start, r.end))

    if len(d) == 0: return
    
    for path_i, x in d.iteritems():
        x.sort(key=lambda x: x[0])
        xx = [_chr for s,e,_chr,_chr_s,_chr_e in x]
        f.write("{0}\t{1}\t{2}\n".format(dirname, path_i, ",".join(xx)))
        for s,e,_chr,_chr_s,_chr_e in x:
            f.write("#{0}:{1}-{2}\t{3}:{4}-{5}\n".format(path_i,s,e,_chr,_chr_s,_chr_e))

