import os, sys
import cupcake.io.GFF as GFF
from Bio import SeqIO
from collections import defaultdict

def read_cogent2_aligned_to_genome_gff(filename):
    """
    Read cogent2 mapped to a genome.

    Return: dict of {cogent path} --> list of gmapRecord; set of mapped genome contigs

    NOTE: (gmap was run with -n 0 so if multiple must be chimeric)
    """
    d = defaultdict(lambda: [])
    contigs_seen = set()

    if not os.path.exists(filename):
        return {}, set()

    try:
        for r in GFF.gmapGFFReader(filename):
            d[r.seqid].append(r)
            contigs_seen.add(r.chr)
    except IndexError:
        pass
    return dict(d), contigs_seen

def is_true_gmap_chimeric(records):
    """
    Given a list of gmap records of a single input, if it truly chimeric if:
    (a) at least two of the records overlap (same chromosome, overlapping location) by 100 bp
    OR
    (b) at least two of the records are on different chromsomes

    In other words, if all the records are on the same chromosome and do not overlap, they are
    NOT truly chimeric.
    """
    by_chr = defaultdict(lambda: [])
    for r in records: by_chr[r.chr].append(r)
    if len(by_chr) > 1: return True # on multiple chromosomes, really chimeric
    else: # all on same chromosome
        flag = False
        records.sort(key=lambda r: r.start)
        for i in xrange(len(records)-1):
            if records[i].end - records[i+1].start >= 100: # overlap by more than 100 bp, true chimeric
                flag = True
                break
        return flag

def calculate_cov_acc(d):
    """
    Given dict of {cogent path} --> list of gmapRecord
    (gmap was run with -n 0 so if multiple must be chimeric)

    If a Cogent contig was mapped chimerically (even if it's on the same contig)
    possible with "switchbacks", then it's considered bad and we want to report it.
    But we also want to rule out any bad mapping. So, check that for any GMAP multimapping,
    that at least two of the mapped loci overlap. Otherwise it is considered "OK".


    """
    worst_cov, worst_acc, has_chimeric = 100, 100, False

    for v in d.itervalues():
        if len(v) > 1 and is_true_gmap_chimeric(v): # is truly chimeric
            has_chimeric = True
            cov = sum(x.coverage for x in v)
            acc = sum(x.identity*x.coverage for x in v)*1./cov
        else: # not chimeric
            cov = v[0].coverage
            acc = v[0].identity
        if cov < worst_cov:
            worst_cov, worst_acc = cov, acc
    return worst_cov, worst_acc, has_chimeric

def tally_for_a_Cogent_dir(dirname, f1, f2, genome1, genome2):
    """
    1. read input mapped to cogent2 (in.trimmed.fa.cogent2.gff)
    2. read cogent2 mapped to genome1
    3. read cogent2 mapped to genome2 (if genome2 does not exist, just repeat genome1)
    """
    if not os.path.exists(os.path.join(dirname, 'COGENT.DONE')):
        return
    seq_info = defaultdict(lambda: [])
    contigs_seen = set()
    # input mapped to Cogent contigs
    filename = os.path.join(dirname, 'in.trimmed.fa.cogent2.gff')
    reader = GFF.gmapGFFReader(filename)
    for r in reader:
        seq_info[r.seqid].append(r)
        contigs_seen.add(r.chr)
    # sanity check that all sequences in in.fa are mapped to cogent2.fa
    for r in SeqIO.parse(open(os.path.join(dirname, 'in.fa')), 'fasta'):
        assert r.id in seq_info

    d_genome1, contig_genome1 = read_cogent2_aligned_to_genome_gff(os.path.join(dirname,'cogent2.fa.'+genome1+'.gff'))
    d_genome2, contig_genome2 = read_cogent2_aligned_to_genome_gff(os.path.join(dirname,'cogent2.fa.'+genome2+'.gff'))

    # write:
    # dirname, # of input, # of cogent contig, # of pacbio_contig, total pacbio cov, pacbio iden
    f1.write("{0}\t{1}\t{2}\t".format(dirname, len(seq_info), len(contigs_seen)))
    cov1, acc1, has_chimeric1 = calculate_cov_acc(d_genome1)
    f1.write("{0}\t{1:.2f}\t{2:.2f}\t{3}\t{4}\t".format(len(contig_genome1), cov1, acc1, has_chimeric1, ",".join(contig_genome1)))
    # (for genome2), # of contig, total worst cov, iden, is_chimeric, comma-separated list of contigs
    cov2, acc2, has_chimeric2 = calculate_cov_acc(d_genome2)
    f1.write("{0}\t{1:.2f}\t{2:.2f}\t{3}\t{4}\n".format(len(contig_genome2), cov2, acc2, has_chimeric2, ",".join(contig_genome2)))


    in_aligned_to_genome1 = os.path.join(dirname, 'in.trimmed.fa.'+genome1+'.gff')
    if os.path.exists(in_aligned_to_genome1):
        d3, junk = read_cogent2_aligned_to_genome_gff(in_aligned_to_genome1)
    else:
        d3 = {}

    for seqid, v in seq_info.iteritems():
        contigs = [x.chr for x in v]
        acc = sum(x.identity*x.coverage for x in v)/sum(x.coverage for x in v)
        f2.write("{0}\t{1}\t{2}\t{3}\t".format(seqid, dirname, ",".join(contigs), acc))

        if not seqid in d3:
            f2.write("NA\t0\tNA\tNA\n")
        else:
            scaffolds = [x.chr for x in d3[seqid]]
            cov = sum(x.coverage for x in d3[seqid])
            acc = sum(x.identity*x.coverage for x in d3[seqid])/cov
            f2.write("{0}\t{1}\t{2}\t{3}\n".format(",".join(scaffolds), len(scaffolds), cov, acc))
