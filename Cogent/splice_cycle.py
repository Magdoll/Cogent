import pdb
import logging
from collections import defaultdict
from Bio import SeqIO
from Cogent.process_path import stitch_string_from_path

log = logging.getLogger('Cogent.splice_cycle')

def max_common_sequence_length(seq, indices, cur_size):
    """
    Return the length of the longest commmon sequence for all subseq starting at
    index indices[0], indices[1], ....etc
    """
    seqlen = len(seq)
    j = cur_size+1
    while j < seqlen-indices[-1]:
        firstnt = seq[indices[0]+j]
        if not all(seq[i+j]==firstnt for i in indices[1:]):
            return j
        j += 1
    return j

def precycle_kmer_adjustment(kmer_size):
    """
    Even with detect_and_replace_cycle (which only looks at repeat kmers within a seq),
    there can be cycles in the graph. So instead, try to detect k-mer re-usage and do simple
    tests using networkx.simple_cycles to see if we can eliminate them as much as possible.
    """
    seqdict = SeqIO.to_dict(SeqIO.parse(open('in.trimmed.fa'), 'fasta'))
    kmer_usage = defaultdict(lambda: defaultdict(lambda: [])) # kmer -> (seqid, i)
    for r in SeqIO.parse(open('in.trimmed.fa'),'fasta'):
        for i in xrange(len(r.seq)-kmer_size):
            kmer_usage[r.seq.tostring()[i:i+kmer_size]][r.id].append(i)

    max_kmer_needed = []
    for kmer,v in kmer_usage.iteritems():
        for seqid, indices in v.iteritems():
            if len(indices) > 1: # kmer appeared more than once in the same sequence
                max_kmer_needed.append(max_common_sequence_length(seqdict[seqid], indices, kmer_size))
    if len(max_kmer_needed) == 0:
        log.info("K-mer in-seq cycle detection: None found at k={0}".format(kmer_size))
        return kmer_size
    else:
        min_adjusted_k = min(max_kmer_needed)
        max_adjusted_k = max(max_kmer_needed)
        mean_adjusted_k = sum(max_kmer_needed)*1./len(max_kmer_needed)
        log.info("K-mer in-seq cycle detection: {0} in-seq cycle found with max common \
        sequence analysis showing min: {1}, max: {2}, mean: {3}".format(\
            len(max_kmer_needed), min_adjusted_k, max_adjusted_k, mean_adjusted_k))

        return min(100, max_adjusted_k)


def detect_and_replace_cycle(G, path_d, weight_d, mermap, max_node_index, kmer_size, debug=False):
    """
    For each path in path_d, if the same node occurs multiple times,
    replace the whole occurrence of <n> -> ... -> <n> with a new node

    UGLY but temporary fix :(

    FOr each cycle replacement:
    (a) add prev -> newnode -> succ in G
    (b) reduce one weight in the path and if falls to zero, remove nodes
    (c) update path
    """
    for seqid, path in path_d.iteritems():
        has_cycle = True
        while has_cycle:
            has_cycle = False
            for n in path:
                if path.count(n) > 1:
                    # add in sanity check
                    oldseq = stitch_string_from_path(path, mermap)
                    has_cycle = True
                    if debug:
                        pdb.set_trace()
                    # find first and last occurrence of n
                    # NOTE: need to handle special case when i=0 and j=last
                    i = path.index(n)
                    j = len(path) - path[::-1].index(n)  # is actually one beyond the last occurrence
                    newnode = max_node_index + 1
                    max_node_index += 1
                    newmer = mermap[path[i]]
                    for p in path[i+1:j]: newmer += mermap[p][kmer_size-1:]
                    mermap[newnode] = newmer
                    # update G before update path
                    if i > 0: # if i==0, then nothing to connect with prev
                        G.add_edge(path[i-1], newnode, weight=weight_d[seqid])
                    if j < len(path): # if j=last, then nothing to connect with after
                        G.add_edge(newnode, path[j], weight=weight_d[seqid])
                    for k in xrange(max(0,i-1), min(j, len(path)-1)):
                        s, t = path[k], path[k+1]
                        G[s][t]['weight'] -= weight_d[seqid]
                    # now we can update the path
                    path_d[seqid] = path[:i] + [newnode] + path[j:]
                    path = path_d[seqid]


                    assert stitch_string_from_path(path, mermap) == oldseq
                    break

    # clean up zero-weight edges and zero degree nodes
    # TODO here
    for s,t,d in G.edges(data=True):
        if d['weight']==0:
            G.remove_edge(s, t)
    bad = []
    for n in G.nodes_iter():
        if G.degree(n) == 0: bad.append(n)
    G.remove_nodes_from(bad)

