__author__ = 'etseng@pacb.com'

import os, sys, time
import logging
from Bio import SeqIO
import splice_graph as sp
import splice_cycle
import GFF
import networkx as nx
from collections import defaultdict
sys.setrecursionlimit(999999)

def trim_ends(seq):
    for i in xrange(len(seq)):
        if str.isupper(seq[i]): break
    for j in xrange(len(seq)-1,-1,-1):
        if str.isupper(seq[j]): break
    return seq[i:j+1]


def split_files(input_filename='in.fa', split_size=20):
    i = 0
    count = 0
    d = "split_0"
    if not os.path.exists(d):
        os.makedirs(d)
    f = open(os.path.join(d, 'in.fa'), 'w')
    split_dirs = [d]
    os.system("cp in.weights {0}/in.weights".format(d))

    for r in SeqIO.parse(open(input_filename), 'fasta'):
        if count >= split_size:
            f.close()
            count = 0
            i += 1
            d = "split_" + str(i)
            if not os.path.exists(d):
                os.makedirs(d)
            split_dirs.append(d)
            os.system("cp in.weights {0}/in.weights".format(d))
            f = open(os.path.join(d, 'in.fa'), 'w')
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        count += 1
    f.close()
    return split_dirs

def run_Cogent_on_split_files(split_dirs):
    time1 = time.time()
    olddir = os.getcwd()
    for d in split_dirs:
        os.chdir(d)
        run_Cogent_on_input()
        os.chdir(olddir)

    if os.path.exists('combined'):
        os.system("rm -rf combined")
    os.makedirs('combined')
    # now combine all the aloha2 results and run LP again
    f = open('combined/aloha.fa', 'w')
    i = 0
    for d in split_dirs:
        for r in SeqIO.parse(open(os.path.join(d, 'aloha2.fa')), 'fasta'):
            f.write(">path{0}\n{1}\n".format(i, r.seq))
            i += 1
    f.close()

    f = open('in.trimmed.fa', 'w')
    for r in SeqIO.parse(open('in.fa'),'fasta'):
        f.write(">{0}\n{1}\n".format(r.id, trim_ends(r.seq.tostring())))
    f.close()

    os.chdir('combined')
    os.system("ln -s ../in.weights in.weights")
    os.system("ln -s ../in.trimmed.fa in.trimmed.fa")
    sp.run_gmap()
    post_gmap_processing()
    os.chdir('../')

    # now the result we want is in combined/aloha2.fa, do postprocessing on it with the full in.fa

    os.system("ln -f -s combined/aloha2.fa aloha2.fa")
    sp.run_gmap(dbname='aloha2', infile='in.trimmed.fa')
    #post_gmap_processing()

    time4 = time.time()
    logger.info("[RUNTIME] Total time in run_Cogent: {0}".format(time4-time1))


def run_Cogent_on_input():
    time1 = time.time()
    # first trim in.fa away all lower case
    f = open('in.trimmed.fa', 'w')
    for r in SeqIO.parse(open('in.fa'),'fasta'):
        f.write(">{0}\n{1}\n".format(r.id, trim_ends(r.seq.tostring())))
    f.close()

    seqweights = {}
    # read in the weights for each sequence
    with open('in.weights') as f:
        for line in f:
            seqid, weight = line.strip().split('\t')
            seqweights[seqid] = int(weight)

    seqdict = SeqIO.to_dict(SeqIO.parse(open('in.trimmed.fa'),'fasta'))
    # setting up the DiGraph
    G = nx.DiGraph()
    node_d = {None: -1}
    path_d = {}
    reader = SeqIO.parse(open('in.trimmed.fa'),'fasta')
    for r in reader: sp.initfunc(G, node_d, path_d, r.seq.tostring(), r.id, seqweights[r.id])
    del node_d[None]
    mermap = dict((v,k) for k,v in node_d.iteritems())

    # resolve all homopolymers
    homo_nodes = filter(lambda n: G.has_edge(n, n), G.nodes_iter())
    for n in homo_nodes:
        sp.untangle_homopolymer_helper(G, path_d, mermap, n)

    splice_cycle.detect_and_replace_cycle(G, path_d, seqweights, mermap, max(G.nodes()), sp.KMER_SIZE)

    visited = {}
    sp.reachability(G, mermap, visited, path_d)

    # cycle detection and abort if detected
    for k,v in path_d.iteritems():
        for x in v:
            if v.count(x) > 1:
                logger.info("CYCLE detected! Abort!")
                os.system("touch CYCLE_DETECTED")
                sys.exit(-1)
    #iter = nx.simple_cycles(G)
    #for it in iter:
    #    print >> sys.stderr, "CYCLE detected! Abort!"
    #    os.system("touch CYCLE_DETECTED")
    #    sys.exit(-1)

    nx.write_graphml(G, 'in.0.graphml')

    logger.info("Initial Graph Size: {0} nodes, {1} edges".format(G.number_of_nodes(), G.number_of_edges()))

    ## sanity check: confirm that all sequences can be reconstructed via the collapsed graph
    ## also check that all nodes are visited
    #for n in G.nodes_iter(): assert n in visited
    #for k,v in path_d.iteritems():
    #    s = sp.stitch_string_from_path(v, mermap)
    #    s2 = seqdict[k].seq.tostring().upper()
    #    assert s.find(s2) >= 0

    while True:
        cur_num_nodes = G.number_of_nodes()
        sp.find_source_bubbles(G, path_d, mermap)
        sp.reachability(G, mermap, {}, path_d)
        sp.find_bubbles(G, path_d, mermap)
        sp.reachability(G, mermap, {}, path_d)
        sp.contract_sinks(G, path_d, mermap)
        sp.find_dangling_sinks(G, path_d, mermap)
        sp.reachability(G, mermap, {}, path_d)
        if G.number_of_nodes() == cur_num_nodes: break

    nx.write_graphml(G, 'in.1.graphml')

    logger.info("Post-Reduction Graph Size: {0} nodes, {1} edges".format(G.number_of_nodes(), G.number_of_edges()))

    time2 = time.time()

    keys = path_d.keys()
    keys.sort()
    good_for, paths = sp.find_minimal_path_needed_to_explain_pathd(G, path_d, keys)
    sp.solve_with_lp_and_reduce(good_for, paths, mermap)

    time3 = time.time()

    sp.run_gmap()
    post_gmap_processing()

    time4 = time.time()

    logger.info("[RUNTIME] for graph construction and reduction: {0}".format(time2-time1))
    logger.info("[RUNTIME] for path finding and LP solving: {0}".format(time3-time2))
    logger.info("[RUNTIME] for GMAP and post-processing: {0}".format(time4-time3))
    logger.info("[RUNTIME] Total time in run_Cogent: {0}".format(time4-time1))

def post_gmap_processing():
    good_for = defaultdict(lambda: [])
    reader = GFF.gmapGFFReader('in.trimmed.fa.gff')
    for r in reader:
        if r.coverage >= 98.: good_for[r.seqid].append(int(r.chr[4:]))

    N = max(max(v) for v in good_for.itervalues())+1
    prob = sp.make_into_lp_problem(good_for.items(), N)
    prob.solve()
    touse = []
    for v in prob.variables():
        logger.debug("{0} = {1}".format(v.name, v.varValue))
        if v.varValue == 1: touse.append(int(v.name))

    with open('aloha2.fa', 'w') as f:
        for r in SeqIO.parse(open('aloha.fa'),'fasta'):
            if int(r.id[4:]) in touse:
                f.write(">{0}\n{1}\n".format(r.id, r.seq))


def run_gmap_for_final_GFFs():
    cmd = '~/bin/gmapl -D ~/share/gmap_db_new/ -d cuttlefish -f gff3_gene -n 0 aloha2.fa > aloha2.fa.cuttlefish.gff'
    os.system(cmd)
    cmd = '~/bin/gmapl -D ~/share/gmap_db_new/ -d cuttlefish -f gff3_gene -n 0 in.trimmed.fa > in.trimmed.fa.cuttlefish.gff'
    os.system(cmd)
    cmd = "~/bin/gmap_build -D . -d aloha2 aloha2.fa"
    os.system(cmd)
    cmd = '~/bin/gmap -D . -d aloha2 -f gff3_gene -n 0 in.trimmed.fa > in.trimmed.fa.aloha2.gff'
    os.system(cmd)

    os.system("gff3_to_collapsed.py aloha2.fa.cuttlefish.gff")
    os.system("gff3_to_collapsed.py in.trimmed.fa.cuttlefish.gff")
    os.system("gff3_to_collapsed.py in.trimmed.fa.aloha2.gff")


def main():
    assert os.path.exists('in.fa')

    if not os.path.exists('in.weights'):
        # make weight file
        with open('in.weights', 'w') as f:
            for r in SeqIO.parse(open('in.fa'), 'fasta'):
                # ex: i0HQ|c96271/f2p7/1352
                f.write("{0}\t{1}\n".format(r.id, int(r.id.split('/')[1].split('p')[0][1:])))

    num_size = int(os.popen("grep -c \">\" in.fa").read().strip())

    if num_size <= 20:
        run_Cogent_on_input()
    else:
        dirs = split_files()
        run_Cogent_on_split_files(dirs)


if __name__ == "__main__":
    logger = logging.getLogger('Cogent')
    logger.setLevel(logging.INFO)

    # create a file handler
    handler = logging.FileHandler(os.path.join(sys.argv[1], 'hello.log'))
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    os.chdir(sys.argv[1])
    main()
    os.system("touch COGENT.DONE")
    run_gmap_for_final_GFFs()