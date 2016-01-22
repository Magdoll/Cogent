from Cogent import splice_align
from Cogent.process_path import stitch_string_from_path

__author__ = 'etseng@pacb.com'

import sys
from Bio import SeqIO

def sanity_check_fasta(fasta_filename):
    """
    Check that all sequences are A/T/C/G
    """
    passed = True
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        seq_passed = True
        seq = str(r.seq).upper()
        for s in seq:
            if s not in ('A','T','C','G'):
                print r.id, s
                seq_passed = False
                break
        if not seq_passed:
            print >> sys.stderr, "sequence {0} contains non A/T/C/G characters! Not OK!".format(r.id)
            passed = False

    if not passed:
        print >> sys.stderr, "Please fix the offending sequences first. Abort."
        sys.exit(-1)


def sanity_check_is_chain(G, chain):
    """
    chain definition:
     u -> n1 -> n2 -> ... -> v

    u has only one outgoing (but can have many incoming)
    v has only one incoming (but can have many outgoing)
    n1, n2, ... have exactly one in coming and one outgoing
    """
    assert G.out_degree(chain[0]) == 1
    for i,n1 in enumerate(chain[1:-1]):
        assert G.in_degree(n1) == 1 and G.out_degree(n1) == 1 and G.predecessors(n1)[0] == chain[i] \
                and G.successors(n1)[0] == chain[i+2]


def sanity_check_reconstruction(path_d, mermap, seqdict):
    for seqid, path in path_d.iteritems():
        seq = stitch_string_from_path(path, mermap)
        orig = seqdict[seqid].seq.tostring().upper()

        flag, cigar = splice_align.validate_reconstructed_seq(seq, orig)
        if not flag:
            print "not validated:", seqid, cigar
            raw_input()


def sanity_check_is_subset(path_d, mermap, seqdict):
    for seqid, path in path_d.iteritems():
        seq = stitch_string_from_path(path, mermap)
        orig = seqdict[seqid].seq.tostring().upper()
        if seq!=orig:
            print "not valid subset for:", seqid
            return False
    return True


def sanity_check_path_all_valid(path_d, G):
    """
    Make sure all the paths are still a valid path in G
    """
    check_passed = True
    for seqid, path in path_d.iteritems():
        for i in xrange(len(path)-1):
            if not G.has_edge(path[i], path[i+1]):
                print "missing path {0}->{1}".format(path[i], path[i+1])
                check_passed = False
    return check_passed