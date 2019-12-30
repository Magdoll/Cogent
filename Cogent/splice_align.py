__author__ = 'etseng@pacb.com'

import os, sys, subprocess
from Bio import SeqIO
from .settings import EXPECTED_ERR_RATE
import parasail

def iter_cigar_string(cigar_string):
    num = cigar_string[0]
    for s in cigar_string[1:]:
        if str.isalpha(s):
            yield int(num), s
            num = ''
        else:
            num += s

def node_is_similar(seq1, seq2):
    l1 = len(seq1)
    l2 = len(seq2)
    if l1 == 0 or l2 == 0: return False
    if l1 <= 2 and l2 <= 2: return True
    if l1 < l2:
        l1, l2 = l2, l1
        seq1, seq2 = seq2, seq1
    # always make seq1 the longer one
    o1 = parasail.sg_qx_trace(seq1, seq2, 3, 1, parasail.matrix_create("ACGT", 2, -5))
    # require the the whole (shorter) seq2 must be aligned
    # and set min score to approx 90% accuracy

    if EXPECTED_ERR_RATE == 0:
        return o1.score > l2*2*1.0
    elif EXPECTED_ERR_RATE < 2:
        return o1.score > l1*2*0.8
    else:
        raise Exception("Expected error rate not implemented for {0}% and above".format(EXPECTED_ERR_RATE))
    return res is not None

def node_is_skipping(seq1, seq2, kmer_size):
    """
    Let seq1 always be the longer one. If seq2 is the exon skipped one, then

    seq_next is the next node for both seq1 and seq2

    we already know that seq1 and seq2 share a common prefix and suffix

    seq1 = [common prefix w/ seq2] + <extra exon here> + [common prefix w/ seq_next]
    seq2 = [common prefix w/ seq1] + [common prefix w/ seq_next]

    seq1 and seq2 is given as trimmed away of the common prefix & suffix
    so if it's already skipped, then we expect one of them to be of length 0 or close
    in practice we allow the seq that is supposed to be 0-length to have a little insertions (5%)
    """
    if not node_is_similar(seq1[:kmer_size-1], seq2[:kmer_size-1]) or not node_is_similar(seq1[-kmer_size+1:], seq2[-kmer_size+1:]):
        return False, False

    seq1 = seq1[kmer_size-1:-(kmer_size-1)]
    seq2 = seq2[kmer_size-1:-(kmer_size-1)]
    l1 = len(seq1)
    l2 = len(seq2)
    if l1 < l2: return 'SEQ2', l1 <= max(1, .05 * l2)
    else: return 'SEQ1', l2 <= max(1, .05 * l1)


def get_consensus_through_voting(seq1, weight1, seq2, weight2):
    """
    Use majority vote. This is a temporary replacement for get_consensus_through_pbdagcon()
    since I don't want to deal with dependency right now.
    """
    if weight1 > weight2: return seq1
    else: return seq2


# NOT SUPPORTED unless ToFU is present -_0
def get_consensus_through_pbdagcon(seq1, weight1, seq2, weight2):
    fname = os.tempnam("/tmp") + '.fasta'
    with open(fname, 'w') as f:
        for i in range(min(10, weight1)): # max cutoff at 10
            f.write(">test1_{0}\n{1}\n".format(i, seq1))
        for i in range(min(10, weight2)): # max cutoff at 10
            f.write(">test2_{0}\n{1}\n".format(i, seq2))

    cmd = "ice_pbdagcon.py --maxScore 10 {0} {0}.g_con g_con".format(fname)
    if subprocess.check_call(cmd, shell=True) == -1:
        raise Exception("Trouble running cmd: ").with_traceback(cmd)

    gcon_filename = fname + '.g_con.fasta'
    gref_filename = fname + '.g_con_ref.fasta'
    if os.path.exists(gcon_filename) and os.stat(gcon_filename).st_size > 0:
        return str(SeqIO.parse(open(gcon_filename), 'fasta').next().seq)
    elif os.path.exists(gref_filename):
        return str(SeqIO.parse(open(gref_filename), 'fasta').next().seq)
    else:
        raise Exception("Could not get results from running cmd:").with_traceback(cmd)


def validate_reconstructed_seq(seq, orig):
    """
    seq --- the sequence that is reconstructed
    orig --- the original sequence

    because the reconstructed seq can be longer, we don't care about deletions
      (deletions w.r.t could just be exon skipping or minor base errors)
    we only care that there is NOT a lot of insertions (which would indicate error in my bubble solution)
    """
    o1 = parasail.sg_qx_trace(seq, orig, 3, 1, parasail.matrix_create("ACGT", 2, -5))
    if o1.score < l2*2*.90: return False, o1.cigar.decode
    for num, type in iter_cigar_string(o1.cigar.decode):
        if type == 'I' and num > 5:
            return False, o1.cigar.decode
    return True, o1.cigar.decode

