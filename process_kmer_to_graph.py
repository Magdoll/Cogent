# IPython log file
import os, sys, time
import itertools
import networkx as nx
import numpy  as np
import matplotlib.pyplot as plt
import process_kmer as sp
from skimage.future import graph
from cPickle import *
from collections import defaultdict
from Bio import SeqIO
from collections import namedtuple

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

def run_ncut(G, labels2_map, ncut_map, nodelist):
    start_t = time.time()
    labels = np.array([i for i in G.nodes_iter()])
    labels2 = graph.cut_normalized(labels, G, thresh=0.2)
    for i1, i2 in itertools.izip(labels, labels2):
        labels2_map[i2].append(nodelist[i1])
        ncut_map[i1] = i2

    print >> sys.stderr, "(subgraph) has {0} nodes. ncut down to {1} partitions in {2} sec.".format(\
        G.number_of_nodes(), len(set(labels2)), time.time()-start_t)

def main(fastq_filename, dist_filename, output_prefix):

    # nodelist: dict of seqid --> index
    seqdict = dict((r.id.split()[0], r) for r in SeqIO.parse(open(fastq_filename),'fastq'))
    nodelist = seqdict.keys()
    nodelist = dict((x,i) for i,x in enumerate(nodelist))

    G = sp.make_weighted_graph_from_mash_dist(nodelist, dist_filename, threshold=0.05)
    for i in G.nodes_iter():
        G.node[i]['labels'] = [i]

    # now we convert nodelist back to index --> seqid
    nodelist = dict((i,x) for (x,i) in nodelist.iteritems())

    ncut_map = {} # label1/node id --> ncut label
    labels2_map = defaultdict(lambda: []) # ncut label --> list of seqids in that cut
    for g in nx.connected_component_subgraphs(G):
        run_ncut(g, labels2_map, ncut_map, nodelist)

    with open(output_prefix + '.log', 'w') as f:
        for k,v in labels2_map.iteritems():
            print k,v
            f.write("{0}\t{1}\t{2}\n".format(k, len(v), v))


# ----- version where there is no pbid (no genome answer)
    for n in G:
        G.node[n]['label'] = str(ncut_map[n])  # label is the assignment of ncut
        G.node[n]['gene'] =  str(nodelist[n])   # since we don't know the gene, just let it be seqid

# ----- version where there is pbid
#    nodelist = np.array(nodelist)
#    xx = defaultdict(lambda: [])
#    for i,x in enumerate(nodelist):
#        pbid = x.split('.')[1]
#        xx[pbid].append(i)
#        if i in G:
#            G.node[i]['gene'] = str(pbid) # for graphml drawing
#            G.node[i]['label'] = str(ncut_map[i])
#
    nx.write_graphml(G, output_prefix + '.graphml')

#    labels = np.array([i for i in G.nodes_iter()])
#    labels2 = np.array([ncut_map[i] for i in G.nodes_iter()])
#    pos = nx.random_layout(G)
#    sp.assign_pos_by_answer(G, pos, xx)
#    sp.draw_assignment(G, pos, labels, labels2, xx, output_prefix)

    # make the directory and the subdirs
    os.makedirs(output_prefix)
    for ncut_label, members in labels2_map.iteritems():
        d2 = os.path.join(output_prefix, str(ncut_label))
        os.makedirs(d2)
        with open(os.path.join(d2, 'in.fa'), 'w') as f:
            for seqid in members:
                f.write(">{0}\n{1}\n".format(seqid, seqdict[seqid].seq))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename")
    parser.add_argument("dist_filename")
    parser.add_argument("output_prefix")

    args = parser.parse_args()

    main(args.fastq_filename, args.dist_filename, args.output_prefix)