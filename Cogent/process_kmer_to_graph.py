import os, sys, time, itertools
from collections import defaultdict, namedtuple
import networkx as nx
import numpy  as np
from Bio import SeqIO
from skimage.future import graph
from settings import process_kmer as pk
from settings import draw_kmer_graphs as drawk


MashDist = namedtuple('MashDist', 'id1 id2 pval sim')


def run_ncut(G, labels2_map, ncut_map, nodelist):
    start_t = time.time()
    labels = np.array([i for i in G.nodes_iter()])
    labels2 = graph.cut_normalized(labels, G, thresh=0.2)
    for i1, i2 in itertools.izip(labels, labels2):
        labels2_map[i2].append(nodelist[i1])
        ncut_map[i1] = i2

    print >> sys.stderr, "(subgraph) has {0} nodes. ncut down to {1} partitions in {2} sec.".format(\
        G.number_of_nodes(), len(set(labels2)), time.time()-start_t)


def write_output_dirs(labels2_map, seqdict, weightdict, output_prefix):
    """
    For each partition, create <output_prefix>/<partition>/in.fa and in.weights
    """
    # make the directory and the subdirs
    os.makedirs(output_prefix)
    for ncut_label, members in labels2_map.iteritems():
        d2 = os.path.join(output_prefix, str(ncut_label))
        os.makedirs(d2)
        with open(os.path.join(d2, 'in.fa'), 'w') as f:
            for seqid in members:
                f.write(">{0}\n{1}\n".format(seqid, seqdict[seqid].seq))
        with open(os.path.join(d2, 'in.weights'), 'w') as f:
            for seqid in members:
                f.write("{0}\t{1}\n".format(seqid, weightdict[seqid]))


def family_finding(dist_filename, seqdict, output_prefix, has_pbid=False):
    '''
    Make a weighted (undirected) graph where each node is a sequence, each edge is k-mer similarity
    Then do normalized cut to find the family partitions

    For each partition, make <output_prefix>/<partition_number>/in.fa

    If the IDs are in PB id format (has genome answer), like PB.1.3
    then write that out as the "gene" (ground truth) label in graphml output for visualization
    '''

    # nodelist: dict of seqid --> index
    nodelist = seqdict.keys()
    nodelist = dict((x,i) for i,x in enumerate(nodelist))  # seqid --> index

    G = pk.make_weighted_graph_from_mash_dist(nodelist, dist_filename, threshold=0.05)
    for n in G:
        G.node[n]['labels'] = [n]

    # now we convert nodelist back to index --> seqid
    nodelist = dict((i,x) for (x,i) in nodelist.iteritems())

    ncut_map = {} # label1/node id --> ncut label
    labels2_map = defaultdict(lambda: []) # ncut label --> list of seqids in that cut
    for g in nx.connected_component_subgraphs(G):
        run_ncut(g, labels2_map, ncut_map, nodelist)

    seqid_unassigned = set(seqdict.keys())
    with open(output_prefix + '.partition.txt', 'w') as f:
        f.write("Partition\tSize\tMembers\n")
        for k,v in labels2_map.iteritems():
            print >> sys.stderr, k,v
            f.write("{0}\t{1}\t{2}\n".format(k, len(v), ",".join(v)))
            for seqid in v:
                seqid_unassigned.remove(seqid)
        f.write("#unassigned:{0}\n".format(",".join(seqid_unassigned)))

    if not has_pbid:
        # ----- version where there is no pbid (no genome answer)
        for n in G:
            G.node[n]['label'] = str(ncut_map[n])  # label is the assignment of ncut
            G.node[n]['gene'] =  str(nodelist[n])   # since we don't know the gene, just let it be seqid
            G.node[n]['labels'] = str(G.node[n]['labels'])  # make it string so can write to graphml
    else:
        # ----- version where there is pbid
        nodelist = np.array(nodelist)
        gene_answer = defaultdict(lambda: [])  # gene index --> list of nodes in that gene (based on truth)
        for i,seqid in enumerate(nodelist):
            pbid = seqid.split('.')[1]
            gene_answer[pbid].append(i)
            if i in G:
                G.node[i]['gene'] = str(pbid) # this is the ground truth
                G.node[i]['label'] = str(ncut_map[i]) # this is the ncut label
        labels = np.array([i for i in G.nodes_iter()])
        labels2 = np.array([ncut_map[i] for i in G.nodes_iter()])
        pos = nx.random_layout(G)
        drawk.assign_pos_by_answer(G, pos, gene_answer)
        drawk.draw_assignment(G, pos, labels, labels2, gene_answer, output_prefix)

    nx.write_graphml(G, output_prefix + '.graphml')

    return labels2_map



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename")
    parser.add_argument("dist_filename")
    parser.add_argument("output_prefix")
    parser.add_argument("-c", "--count_filename", help="Count filename (if not given, assume all weight is 1)")

    args = parser.parse_args()

    if not os.path.exists(args.fastq_filename):
        print >> sys.stderr, "Fastq filename {0} does not exist. Abort.".format(args.fastq_filename)
        sys.exit(-1)

    if not os.path.exists(args.dist_filename):
        print >> sys.stderr, "Mash dist filename {0} does not exist. Abort.".format(args.dist_filename)
        sys.exit(-1)

    if os.path.exists(args.output_prefix):
        print >> sys.stderr, "Output directory {0} already exists. Not OK! Abort!".format(args.output_prefix)
        sys.exit(-1)


    seqdict = dict((r.id.split()[0], r) for r in SeqIO.parse(open(args.fastq_filename),'fastq'))

    if args.count_filename is not None:
        weightdict = {}
        for line in open(args.count_filename):
            seqid, weight = line.strip().split('\t')
            weightdict[seqid] = int(weight)
    else:
        weightdict = dict((seqid, 1) for seqid in seqdict)

    labels2_map = family_finding(args.dist_filename, seqdict, args.output_prefix, has_pbid=False)
    write_output_dirs(labels2_map, seqdict, weightdict, args.output_prefix)
