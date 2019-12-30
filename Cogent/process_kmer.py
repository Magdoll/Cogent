import sys
from collections import namedtuple

import networkx as nx


MashDist = namedtuple('MashDist', 'id1 id2 pval sim')

def mash_distance_reader(dist_filename):
    """
    Output from mash dist should have format:
    id1   id2    p-value    ???     <shared>/<total>
    """
    for line in open(dist_filename):
        raw = line.strip().split()
        a, b = list(map(int, raw[4].split('/')))
        if raw[0] != raw[1]:
            yield MashDist(id1=raw[0], id2=raw[1], pval=float(raw[2]), sim=a*1./b)


def make_weighted_graph_from_mash_dist(nodelist, dist_filename, threshold):
    G = nx.Graph()

    for r in mash_distance_reader(dist_filename):
        if r.sim >= threshold:
            if r.id1 not in nodelist or r.id2 not in nodelist:
                print("id1 {0} or id2 {1} not in nodelist. Ignore.".format(r.id1, r.id2), file=sys.stderr)
                continue
            G.add_edge(nodelist[r.id1], nodelist[r.id2], weight=r.sim)

    return G


## OBSOLETE now that I am using mash dist
# def process_kmer(kmer_file):
#     """
#     format:
#     @<seqid>
#     <# of kmers to follow>
#     ...
#     ...
#     each line is a kmer
#     """
#     kmer_index_max = 0
#     kmer_d = {}  # kmer --> kmer index
#     seq_d  = {}  # seq --> set of kmer index
#     f = open(kmer_file)
#     while True:
#         line = f.readline().strip()
#         if len(line) == 0: break
#         assert line.startswith('@')
#         seqid = line[1:].split()[0]
#         seq_d[seqid] = set()
#         num_kmers = int(f.readline().strip())
#         for i in xrange(num_kmers):
#             kmer = f.readline().strip()
#             if kmer not in kmer_d:
#                 kmer_d[kmer] = kmer_index_max
#                 kmer_index_max += 1
#             seq_d[seqid].add(kmer_d[kmer])
#     return seq_d

## OBSOLETE now that I am using mash dist
# def make_weighted_graph_from_kmer(seq_d, threshold):
#     G = nx.Graph()
#
#     keys = seq_d.keys()
#     n = len(keys)
#     for i in xrange(n-1):
#         x1 = seq_d[keys[i]]
#         for j in xrange(i+1, n):
#             assert i!=j and keys[i]!=keys[j]
#             x2 = seq_d[keys[j]]
#             sim = len(x1.intersection(x2))*1. / min(len(x1), len(x2))
#             if sim >= threshold:
#                 G.add_edge(i, j, weight=sim)
#     return keys, G







