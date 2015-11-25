import sys
import numpy as np
import matplotlib.pyplot as plt
import math
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
        a, b = map(int, raw[4].split('/'))
        if raw[0] != raw[1]:
            yield MashDist(id1=raw[0], id2=raw[1], pval=float(raw[2]), sim=a*1./b)


def make_weighted_graph_from_mash_dist(nodelist, dist_filename, threshold):
    G = nx.Graph()

    for r in mash_distance_reader(dist_filename):
        if r.sim >= threshold:
            if r.id1 not in nodelist or r.id2 not in nodelist:
                print >> sys.stderr, "id1 {0} or id2 {1} not in nodelist. Ignore.".format(r.id1, r.id2)
                continue
            G.add_edge(nodelist[r.id1], nodelist[r.id2], weight=r.sim)

    return G


## OBSOLETE if using mash dist!!!
def process_kmer(kmer_file):
    """
    format:
    @<seqid>
    <# of kmers to follow>
    ...
    ...
    each line is a kmer
    """
    kmer_index_max = 0
    kmer_d = {}  # kmer --> kmer index
    seq_d  = {}  # seq --> set of kmer index
    f = open(kmer_file)
    while True:
        line = f.readline().strip()
        if len(line) == 0: break
        assert line.startswith('@')
        seqid = line[1:].split()[0]
        seq_d[seqid] = set()
        num_kmers = int(f.readline().strip())
        for i in xrange(num_kmers):
            kmer = f.readline().strip()
            if kmer not in kmer_d:
                kmer_d[kmer] = kmer_index_max
                kmer_index_max += 1
            seq_d[seqid].add(kmer_d[kmer])
    return seq_d

def make_weighted_graph_from_kmer(seq_d, threshold):
    G = nx.Graph()

    keys = seq_d.keys()
    n = len(keys)
    for i in xrange(n-1):
        x1 = seq_d[keys[i]]
        for j in xrange(i+1, n):
            assert i!=j and keys[i]!=keys[j]
            x2 = seq_d[keys[j]]
            sim = len(x1.intersection(x2))*1. / min(len(x1), len(x2))
            if sim >= threshold:
                G.add_edge(i, j, weight=sim)
    return keys, G


def plot_proportion_related_vs_unrelated(seq_d, output_filename):
    f = open(output_filename, 'w')
    f.write("id1\tid2\ttag\tsim\n")

    n = len(seq_d)
    keys = seq_d.keys()
    for i in xrange(n-1):
        x1 = seq_d[keys[i]]
        for j in xrange(i+1, n):
            x2 = seq_d[keys[j]]
            sim = len(x1.intersection(x2))*1. / min(len(x1), len(x2))
            if keys[i].split('.')[1] == keys[j].split('.')[1]: # ex: PB.1.3 vs PB.1.4
                tag = 'related'
            else:
                tag = 'unrelated'
            f.write("{0}\t{1}\t{2}\t{3}\n".format(keys[i], keys[j], tag, sim))
    f.close()


def draw_reduced_graph(G):
    plt.figure(figsize=(14, 14))
    plt.clf()
    pos = nx.spring_layout(G)

    # generate N points along the same line


def draw_assignment(G, pos, labels, labels2, answer_d, output_prefix):
    plt.figure(figsize=(14, 14))
    plt.clf()
    groups = np.unique(labels2)
    m = len(groups)
    for i in xrange(m):
        nx.draw_networkx_nodes(G, pos, nodelist=list(labels[labels2==groups[i]]), \
                node_color=plt.cm.jet(i*1./m), node_size=300, alpha=.5)

    ans_labels = {}
    for k,v in answer_d.iteritems():
        for x in v: ans_labels[x] = k
    nx.draw_networkx_labels(G, pos, ans_labels, font_size=12)

    plt.axis('off')
    plt.savefig(output_prefix +'.png')
    plt.show()

    # show only edges that are > 0.9
    pp = filter(lambda (a,b,d): d['weight']>=.9, G.edges(data=True))
    nx.draw_networkx_edges(G, pos, edgelist=pp, width=1, alpha=.5)
    plt.savefig(output_prefix+'.with_edges.png')
    plt.show()

def assign_pos_by_answer(G, pos, answer_d):
    m = math.sqrt(len(answer_d)) + 1
    keys = answer_d.keys()
    for i in xrange(len(keys)):
        _row, _col = i/m, i%m
        x_min = _row * 1. / m
        x_max = x_min + 1./m
        y_min = _col * 1. / m
        y_max = y_min + 1./m
        pts = generate_points_in_limit(len(answer_d[keys[i]]), x_min, x_max, y_min, y_max)
        for j,v in enumerate(answer_d[keys[i]]): pos[v] = pts[j,:]

def generate_points_in_limit(n, x_min, x_max, y_min, y_max):
    p = np.random.random((n*2, 2))
    func = np.bitwise_and
    good = p[func(func(x_min <= p[:,0], p[:,0] < x_max), func(y_min <= p[:,1], p[:,1] < y_max))]
    if len(good) >= n: return good
    while True:
        p = np.random.random((n, 2))
        tmp = p[func(func(x_min <= p[:,0], p[:,0] < x_max), func(y_min <= p[:,1], p[:,1] < y_max))]
        good = np.vstack([good, tmp])
        if len(good) >= n: break
    return good[:n, :]