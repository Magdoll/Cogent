__author__ = 'etseng@pacb.com'

import os
import sys
import time
import logging
from Bio import SeqIO
import networkx as nx

from Cogent.__init__ import get_version
from Cogent import settings as cc_settings
from Cogent import sanity_checks, splice_cycle
from Cogent import splice_graph as sp
from Cogent.Utils import trim_ends, run_external_call, run_gmap, cleanup_gmap, post_gmap_processing, run_gmap_for_final_GFFs
from Cogent.process_path import solve_with_lp_and_reduce, find_minimal_path_needed_to_explain_pathd
from Cogent.sanity_checks import sanity_check_gmapl_exists, sanity_check_path_all_valid


class CycleDetectedException(Exception):
    pass

sys.setrecursionlimit(999999)


def split_files(input_filename='in.fa', split_size=20):
    """
    Split input files into split_0/in.fa, split_1/in.fa...
    Return the list of split directories
    """
    i = 0
    count = 0
    d = "split_0"
    if not os.path.exists(d):
        os.makedirs(d)
    f = open(os.path.join(d, 'in.fa'), 'w')
    split_dirs = [d]
    run_external_call("cp in.weights {0}/in.weights".format(d))

    for r in SeqIO.parse(open(input_filename), 'fasta'):
        if count >= split_size:
            f.close()
            count = 0
            i += 1
            d = "split_" + str(i)
            if not os.path.exists(d):
                os.makedirs(d)
            split_dirs.append(d)
            run_external_call("cp in.weights {0}/in.weights".format(d))
            f = open(os.path.join(d, 'in.fa'), 'w')
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        count += 1
    f.close()
    return split_dirs


def run_Cogent_on_split_files(split_dirs):
    """
    1. run Cogent individually on each split directory
    2. combine all cogent2.fa from split directories, pretend they are the "INPUT", run Cogent on it

    """
    time1 = time.time()
    olddir = os.getcwd()
    for d in split_dirs:
        os.chdir(d)
        if os.path.exists('cogent2.fa'):
            print >> sys.stderr, "skipping {0} because done already".format(d)
            os.chdir(olddir)
            continue
        run_Cogent_on_input()
        # clean up cogent in the split dir
        if os.path.exists('cogent') and os.path.isdir('cogent'):
            cleanup_gmap('cogent')
        if os.path.exists('cogent2') and os.path.isdir('cogent2'):
            cleanup_gmap('cogent2')
        os.chdir(olddir)

    if os.path.exists('combined'):
        run_external_call("rm -rf combined")
    os.makedirs('combined')
    # now combine all the cogent2 results and pretend they are the "INPUT"
    f = open('combined/in.fa', 'w')
    f2 = open('combined/in.weights', 'w')
    i = 0
    for d in split_dirs:
        for r in SeqIO.parse(open(os.path.join(d, 'cogent2.fa')), 'fasta'):
            f.write(">fake_input_path{0}\n{1}\n".format(i, r.seq))
            f2.write("fake_input_path{0}\t1\n".format(i))
            i += 1
    f.close()
    f2.close()

    os.chdir('combined')
    if i > cc_settings.MAX_POST_SPLIT_IN_SIZE:
        dirs = split_files(input_filename='in.fa', split_size=cc_settings.MAX_POST_SPLIT_IN_SIZE)
        run_Cogent_on_split_files(dirs)
    run_Cogent_on_input()
    os.chdir('../')

    # now take the output from combined and run LP against it,
    # using the real input this time

    with open('in.trimmed.fa', 'w') as f:
        for r in SeqIO.parse(open('in.fa'), 'fasta'):
            f.write(">{0}\n{1}\n".format(r.id, trim_ends(str(r.seq))))

    if os.path.exists('post_combined'):
        run_external_call("rm -rf post_combined")
    os.makedirs('post_combined')
    os.chdir('post_combined')
    run_external_call("ln -s ../combined/cogent2.fa cogent.fa")
    run_external_call("ln -s ../in.weights in.weights")
    run_external_call("ln -s ../in.trimmed.fa in.trimmed.fa")
    run_gmap()
    post_gmap_processing(seqrecs=[r for r in SeqIO.parse(open('in.trimmed.fa'), 'fasta')])
    os.chdir('../')

    # now the result we want is in combined/cogent2.fa, do postprocessing on it with the full in.fa

    run_external_call("ln -f -s post_combined/cogent2.fa cogent2.fa")
    run_gmap(dbname='cogent2', infile='in.trimmed.fa')
    #post_gmap_processing()

    time4 = time.time()
    log.info("[RUNTIME] Total time in run_Cogent: {0}".format(time4-time1))


def run_Cogent_on_input():
    """
    The main reconstruction function.

    Homopolymers and repeated nodes in path must be resolved first.
    (however, it's possible the graph contains cycles not manifested in path,
     this is a bug that will result in failure to *explain* the sequences later,
     right now I catch the bug by using the sequence pth itself but this should be fixed eventually)

    Graph reduction is iteratively done until cannot be further reduced

    Two points of failure:
    (1) graph is not reduced to small enough, too many paths, mem explosion
        cur soln: fall back to using own paths
    (2) cycle in graph
        cur soln: fall back to using own paths (still wrong)
    """
    time1 = time.time()
    # first trim in.fa away all lower case
    f = open('in.trimmed.fa', 'w')
    for r in SeqIO.parse(open('in.fa'),'fasta'):
        f.write(">{0}\n{1}\n".format(r.id, trim_ends(str(r.seq))))
    f.close()

    seqweights = {}
    # read in the weights for each sequence
    with open('in.weights') as f:
        for line in f:
            seqid, weight = line.strip().split('\t')
            seqweights[seqid] = int(weight)

    adjusted_kmer = splice_cycle.precycle_kmer_adjustment(cc_settings.KMER_SIZE)
    if adjusted_kmer != cc_settings.KMER_SIZE:
        log.info("Adjusting k-mer size to: {0}".format(adjusted_kmer))
        cc_settings.KMER_SIZE = adjusted_kmer

    # setting up the DiGraph
    G = nx.DiGraph()
    node_d = {None: -1}  # this is just used to initialize the graph, delete it later
    path_d = {}
    reader = SeqIO.parse(open('in.trimmed.fa'),'fasta')
    seqrecs = []
    for r in reader:
        sp.add_seq_to_graph(G, node_d, path_d, str(r.seq), r.id, seqweights[r.id])
        seqrecs.append(r)
    del node_d[None]
    mermap = dict((v,k) for k,v in node_d.iteritems())



    # resolve all homopolymers
    homo_nodes = filter(lambda n: G.has_edge(n, n), G.nodes_iter())
    for n in homo_nodes:
        sp.untangle_homopolymer_helper(G, path_d, mermap, seqweights, n)

    splice_cycle.detect_and_replace_cycle(G, path_d, seqweights, mermap, max(G.nodes()), cc_settings.KMER_SIZE)

    visited = {}
    sp.reachability(G, mermap, visited, path_d)

    # cycle detection and abort if detected
    # (this should not happen with splice_cycle.detect_and_replace_cycle run)
    for k,v in path_d.iteritems():
        for x in v:
            if v.count(x) > 1:
                log.info("CYCLE detected through path analysis! Raise CycleDetectedException!")
                os.system("touch CYCLE_DETECTED")
                raise CycleDetectedException

    if cc_settings.NX_CYCLE_DETECTION:
        log.info("Doing nx.cycle_detection....")
        iter = nx.simple_cycles(G)
        for _it in iter:
            print >> sys.stderr, "CYCLE detected through simple_cycles! Raise CycleDetectedException!"
            os.system("touch CYCLE_DETECTED")
            raise CycleDetectedException

    nx.write_graphml(G, 'in.0.graphml')

    log.info("Initial Graph Size: {0} nodes, {1} edges".format(G.number_of_nodes(), G.number_of_edges()))

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
        #assert sanity_check_path_all_valid(path_d, G)
        if G.number_of_nodes() == cur_num_nodes:
            break

    nx.write_graphml(G, 'in.1.graphml')

    log.info("Post-Reduction Graph Size: {0} nodes, {1} edges".format(G.number_of_nodes(), G.number_of_edges()))

    time2 = time.time()

    keys = path_d.keys()
    keys.sort()
    good_for, paths = find_minimal_path_needed_to_explain_pathd(G, path_d, keys)
    solve_with_lp_and_reduce(good_for, paths, mermap)

    time3 = time.time()

    run_gmap()
    post_gmap_processing(seqrecs=seqrecs)

    time4 = time.time()

    log.info("[RUNTIME] for graph construction and reduction: {0}".format(time2-time1))
    log.info("[RUNTIME] for path finding and LP solving: {0}".format(time3-time2))
    log.info("[RUNTIME] for GMAP and post-processing: {0}".format(time4-time3))
    log.info("[RUNTIME] Total time in run_Cogent: {0}".format(time4-time1))



def main():
    assert os.path.exists('in.fa')
    assert os.path.exists('in.weights')

    sanity_checks.sanity_check_fasta('in.fa')

    num_size = int(os.popen("grep -c \">\" in.fa").read().strip())

    if num_size <= cc_settings.MAX_SPLIT_IN_SIZE:
        run_Cogent_on_input()
    else:
        dirs = split_files(input_filename='in.fa', split_size=cc_settings.MAX_SPLIT_IN_SIZE)
        run_Cogent_on_split_files(dirs)

    # clean up GMAP db files
    if os.path.exists('cogent') and os.path.isdir('cogent'):
        cleanup_gmap('cogent')
    if os.path.exists('cogent2') and os.path.isdir('cogent2'):
        cleanup_gmap('cogent2')

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("dirname")
    parser.add_argument("-e", "--expected_error_rate", type=int, default=1, help="Expected error rate (default: 1%)")
    parser.add_argument("--nx_cycle_detection", default=False, action="store_true", help="Cycle detection using networkx (default: off), will increase run-time. Recommend for debugging failed cases only.")
    parser.add_argument("-k", "--kmer_size", type=int, default=30, help="kmer size (default: 30)")
    parser.add_argument("-D", "--gmap_db_path", help="GMAP database location (optional)")
    parser.add_argument("-d", "--gmap_species", help="GMAP species name (optional)")
    parser.add_argument("--small_genome", action="store_true", default=False, help="Genome size is smaller than 3GB (use gmap instead of gmapl)")
    parser.add_argument('--version', action='version', version='%(prog)s ' + str(get_version()))
    parser.add_argument("--debug", action="store_true", default=False)

    args = parser.parse_args()

    cc_settings.KMER_SIZE = args.kmer_size
    assert sp.cc_settings.KMER_SIZE == args.kmer_size

    cc_settings.EXPECTED_ERR_RATE = args.expected_error_rate
    cc_settings.NX_CYCLE_DETECTION = args.nx_cycle_detection



    if not args.small_genome:
        sanity_check_gmapl_exists()

    log = logging.getLogger('Cogent')
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    # clean out the log and DONE file if there was a previous run
    if os.path.exists("COGENT.DONE"):
        os.remove("COGENT.DONE")
    if os.path.exists("hello.log"):
        os.remove('hello.log')


    # create a file handler
    handler = logging.FileHandler(os.path.join(args.dirname, 'hello.log'))
    if args.debug:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    log.addHandler(handler)

    log.info("Setting k-mer size to: {0}".format(args.kmer_size))
    log.info("Setting expected error rate to: {0}%".format(args.expected_error_rate))

    os.chdir(args.dirname)
    while cc_settings.KMER_SIZE <= 200:
        try:
            main()
            break
        except CycleDetectedException:
            cc_settings.KMER_SIZE += 20
            log.info("Trying to raise k={0}. Rerun attempt.".format(cc_settings.KMER_SIZE))

    os.system("touch COGENT.DONE")

    if args.gmap_db_path is not None and args.gmap_species is not None:
        run_gmap_for_final_GFFs(small_genome=args.small_genome, gmap_db_path=args.gmap_db_path, species_db=args.gmap_species)

