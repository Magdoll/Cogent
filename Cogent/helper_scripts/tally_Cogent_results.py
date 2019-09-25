import os, re, sys
from cupcake.io import GFF
from Bio import SeqIO
from collections import defaultdict
from bx.intervals.cluster import ClusterTree
from Cogent import BioReaders

"""
# BLASTN 2.6.0+
# Query: Domino_testis_i0_HQ_sampled07afe|c338324/f2p10/426
# Database: /pbi/dept/secondary/siv/gconcepcion/db/ncbi/nt
# Fields: query acc.ver, subject acc.ver, subject title, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 1 hits found
Domino_testis_i0_HQ_sampled07afe|c338324/f2p10/426      NR_146154.1     Homo sapiens RNA, 28S ribosomal (LOC109910382), ribosomal RNA   99.750  400     1       0   1
        400     2557    2956    0.0     734
"""

rex = re.compile('#\sFields:[\S\s]+')
def read_blastn(filename, qlen_dict):
    """
    Read a BLASTn output file and return just "subject title"
    """
    f = open(filename) 
    best_of = defaultdict(lambda: (100, 'NA'))
    f.readline()
    assert f.readline().strip().split()[1]=='Query:'
    assert f.readline().strip().split()[1]=='Database:'
    line = f.readline().strip()
    if line.find('0 hits found') > 0: return best_of
    fields = line.strip().split(',')
    try:
        i = fields.index (' subject title')
        j = fields.index (' evalue')
        k = fields.index(' q. start')
        l = fields.index(' q. end')
    except ValueError:
        print >> sys.stderr, "Unable to find fields 'evalue' and 'subject title' in {0}. Abort!".format (filename)
    f.readline()
    for line in f:
        if line.startswith('#'): continue
        raw = line.strip().split('\t')
        seqid, e, name, qstart, qend = raw[0], float(raw[j]), raw[i], int(raw[k]), int(raw[l])
        if qend-qstart >= .8*qlen_dict[seqid]:  #BLASTn result have to cover 80% of the sequence length
            if e < best_of[seqid]: best_of[seqid] = (e, name)
    return best_of


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


def read_cogent2_aligned_to_genome_sam(input, filename):
    """
    Read cogent2 mapped to a genome.

    Return: dict of {cogent path} --> list of SAM Record; set of mapped genome contigs

    NOTE: (minimap2 was run with --secondary=no so if multiple must be chimeric)
    """
    d = defaultdict(lambda: [])
    contigs_seen = set()

    if not os.path.exists(filename):
        return {}, set()

    try:
        for r in BioReaders.GMAPSAMReader(filename, True, query_len_dict=dict((r.id, len(r.seq)) for r in SeqIO.parse(open(input),'fasta'))):
            if r.sID == '*': continue # unmapped
            d[r.qID].append(r)
            contigs_seen.add(r.sID)
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

def is_true_minimap2_chimeric(records):
    """
    Given a list of minimap2 records of a single input, if it truly chimeric if:
    (a) at least two of the records overlap (same chromosome, overlapping location) by 100 bp
    OR
    (b) at least two of the records are on different chromsomes

    In other words, if all the records are on the same chromosome and do not overlap, they are
    NOT truly chimeric.
    """
    by_chr = defaultdict(lambda: [])
    for r in records: by_chr[r.sID].append(r)
    if len(by_chr) > 1: return True # on multiple chromosomes, really chimeric
    else: # all on same chromosome
        flag = False
        records.sort(key=lambda r: r.sStart)
        for i in xrange(len(records)-1):
            if records[i].sEnd - records[i+1].sStart >= 100: # overlap by more than 100 bp, true chimeric
                flag = True
                break
        return flag

def calculate_cov_acc(d):
    """
    Given dict of {cogent path} --> list of minimap2 record
    (minimap2 was run with --secondary=no so if multiple must be chimeric)

    If a Cogent contig was mapped chimerically (even if it's on the same contig)
    possible with "switchbacks", then it's considered bad and we want to report it.
    But we also want to rule out any bad mapping. So, check that for any multimapping,
    that at least two of the mapped loci overlap. Otherwise it is considered "OK".
    """
    worst_cov, worst_acc, has_chimeric = 100, 100, False
    if len(d) == 0: return 0, 0, False

    for v in d.itervalues():
        c = ClusterTree(0,0)
        for x in v:
            qlen = x.qLen
            c.insert(x.qStart, x.qEnd, -1)
        cov = sum(_e-_s for _s,_e,_junk in c.getregions())*100./qlen
        acc = sum(x.identity*x.qCoverage for x in v)*1./sum(x.qCoverage for x in v)
        if len(v) > 1 and is_true_minimap2_chimeric(v): # is truly chimeric
            has_chimeric = True
        if cov < worst_cov:
            worst_cov, worst_acc = cov, acc
    return worst_cov, worst_acc, has_chimeric

def tally_for_a_Cogent_dir(dirname, writer1, writer2, genome1, genome2=None, blastn_filename=None):
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
    filename = os.path.join(dirname, 'in.trimmed.fa.cogent2.sam')
    reader = BioReaders.GMAPSAMReader(filename, True, \
                                      query_len_dict=dict((r.id, len(r.seq)) for r in SeqIO.parse(open(os.path.join(dirname, 'in.trimmed.fa')), 'fasta')))
    for r in reader:
        seq_info[r.qID].append(r)
        contigs_seen.add(r.sID)
    # sanity check that all sequences in in.fa are mapped to cogent2.fa
    for r in SeqIO.parse(open(os.path.join(dirname, 'in.fa')), 'fasta'):
        assert r.id in seq_info

    d_genome1, contig_genome1 = read_cogent2_aligned_to_genome_sam(os.path.join(dirname, 'cogent2.fa'), os.path.join(dirname,'cogent2.fa.'+genome1+'.sam'))
    if genome2 is not None:
        d_genome2, contig_genome2 = read_cogent2_aligned_to_genome_sam(os.path.join(dirname, 'cogent2.fa'), os.path.join(dirname,'cogent2.fa.'+genome2+'.sam'))

    if blastn_filename is not None:
        qlen_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(os.path.join(dirname, 'in.trimmed.fa')),'fasta'))
        best_of = read_blastn(os.path.join(dirname, blastn_filename), qlen_dict)

    # write:
    # dirname, # of input, # of cogent contig, # of pacbio_contig, total pacbio cov, pacbio iden
    cov1, acc1, has_chimeric1 = calculate_cov_acc(d_genome1)
    rec1 = {'gene_family': dirname,
            'input_size': len(seq_info),
            'num_Cogent_contigs': len(contigs_seen),
            'num_genome_contig': len(contig_genome1),
            'genome_cov': cov1,
            'genome_acc': acc1,
            'genome_chimeric': has_chimeric1,
            'genome_contigs': ",".join(contig_genome1)}


    # (for genome2), # of contig, total worst cov, iden, is_chimeric, comma-separated list of contigs
    if genome2 is not None:
        cov2, acc2, has_chimeric2 = calculate_cov_acc(d_genome2)
        rec1['num_genome2_contig'] = len(contig_genome2)
        rec1['genome2_cov'] = cov2
        rec1['genome2_acc'] = acc2
        rec1['genome2_chimeric'] = has_chimeric2
        rec1['genome2_contigs'] = ",".join(contig_genome2)
    # (for blastn, optional) best name with best e-value
    if blastn_filename is not None:
        if len(best_of) == 0:
            rec1['num_blastn'] = 0
            rec1['blastn_best'] = 'NA'
        else:
            stuff = best_of.values() # list of (e-value, name)
            stuff.sort()
            rec1['num_blastn'] = sum(_n!='NA' for _e,_n in best_of.values())
            rec1['blastn_best'] = stuff[0][1]
    writer1.writerow(rec1)

    in_aligned_to_genome1 = os.path.join(dirname, 'in.trimmed.fa.'+genome1+'.sam')
    if os.path.exists(in_aligned_to_genome1):
        d3, junk = read_cogent2_aligned_to_genome_sam(os.path.join(dirname, 'in.trimmed.fa'), in_aligned_to_genome1)
    else:
        d3 = {}

    for seqid, v in seq_info.iteritems():
        contigs = [x.sID for x in v]
        acc = sum(x.identity*x.qCoverage for x in v)/sum(x.qCoverage for x in v)

        rec2 = {'seqid': seqid,
                'gene_family': dirname,
                'Cogent_contig': ",".join(contigs),
                'Cogent_contig_acc': acc}

        if not seqid in d3:
            rec2['scaffold'] = 'NA'
            rec2['num_scaffold'] = 0
            rec2['scaffold_coverage'] = 'NA'
            rec2['scaffold_acc'] = 'NA'
            if blastn_filename is not None:
                rec2['blastn_best'] = 'NA'
        else:
            scaffolds = [x.sID for x in d3[seqid]]
            # calculate cov and acc
            c = ClusterTree(0,0)
            for x in d3[seqid]:
                qlen = x.qLen
                c.insert(x.qStart, x.qEnd, -1)
            cov = sum(_e-_s for _s,_e,_junk in c.getregions())*100./qlen
            acc = sum(x.identity*x.qCoverage for x in d3[seqid])*1./sum(x.qCoverage for x in d3[seqid])
            rec2['scaffold'] = ",".join(scaffolds)
            rec2['num_scaffold'] = len(scaffolds)
            rec2['scaffold_coverage'] = cov
            rec2['scaffold_acc'] = acc
            if blastn_filename is not None:
                rec2['blastn_best'] = best_of[seqid][1]
        writer2.writerow(rec2)

