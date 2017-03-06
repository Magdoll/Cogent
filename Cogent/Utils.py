import os, subprocess, logging
from collections import defaultdict
from Bio import SeqIO

from Cogent import GFF
from Cogent.process_path import make_into_lp_problem

log = logging.getLogger('Cogent.Utils')

def trim_ends(seq, dun_trim_if_all_lower=True):
    """
    Remove ends that are lower case
    However if whole sequence is lower case???? return everything =___=
    """
    for i in xrange(len(seq)):
        if str.isupper(seq[i]): break
    for j in xrange(len(seq)-1,-1,-1):
        if str.isupper(seq[j]): break

    if i>=j+1 and dun_trim_if_all_lower:
        return seq.upper()
    return seq[i:j+1]

def run_external_call(cmd):
    if subprocess.check_call(cmd, shell=True) != 0:
        raise Exception, "Failed to run: {0}".format(cmd)


def run_gmap(dbname='cogent', infile='in.trimmed.fa'):
    run_external_call("gmap_build -D . -d {0} {0}.fa".format(dbname))
    run_external_call("gmap -D . -d {0} -n 100 -f gff3_gene -t 1 {1} > {1}.gff".format(dbname, infile))

def cleanup_gmap(dbname='cogent'):
    run_external_call("rm -rf " + dbname)

def post_gmap_processing(db_name='cogent', gff_filename='in.trimmed.fa.gff', output_prefix='cogent2', seqrecs=[]):
    good_for = defaultdict(lambda: [])
    reader = GFF.gmapGFFReader(gff_filename)
    for r in reader:
        assert r.chr.startswith('path')  # chr should be path0, path1, etc
        if r.coverage >= 98.: good_for[r.seqid].append(int(r.chr[4:]))

    touse = []
    if len(good_for) == 0:
        log.warning("[BUG] good_for in post_gmap_processing is empty. Probably from cycles. CHECK!")
    else:
        N = max(max(v) for v in good_for.itervalues())+1
        prob = make_into_lp_problem(good_for.items(), N)
        prob.solve()
        for v in prob.variables():
            log.debug("{0} = {1}".format(v.name, v.varValue))
            if v.varValue == 1: touse.append(int(v.name))


    with open(output_prefix + '.fa', 'w') as f:
        for r in SeqIO.parse(open(db_name + '.fa'),'fasta'):
            if int(r.id[4:]) in touse:
                f.write(">{0}\n{1}\n".format(r.id, r.seq))
        # if there are some sequences that didn't map (possibly from cycles)
        # then just use THEMSELVES
        fake_path_i = max(touse)+1 if len(touse) >= 1 else 0
        for r in seqrecs:
            if r.id not in good_for:
                log.warning("[BUG] {0} is not fully mapped to cogent in GMAP. \
                Likely cycle issues. Use itself in output.".format(r.id))
                f.write(">path{0}\n{1}\n".format(fake_path_i, r.seq))
                fake_path_i += 1

def run_gmap_for_final_GFFs(small_genome=False, gmap_db_path='~/share/gmap_db_new/', species_db='cuttlefish',\
                            input='in.trimmed.fa', output='cogent2'):

    if small_genome:
        prog = 'gmap'
    else:
        prog = 'gmapl'
    cmd = prog + " -D {p} -d {sp} -f gff3_gene -n 0 {o}.fa > {o}.fa.{sp}.gff".format(\
        p=gmap_db_path, sp=species_db, o=output)
    run_external_call(cmd)
    cmd = prog + " -D {p} -d {sp} -f gff3_gene -n 0 {i} > {i}.{sp}.gff".format(\
        p=gmap_db_path, sp=species_db, i=input)
    run_external_call(cmd)
    cmd = "gmap_build -D . -d {o} {o}.fa".format(o=output)
    run_external_call(cmd)
    cmd = "gmap -D . -d {o} -f gff3_gene -n 0 {i} > {i}.{o}.gff".format(\
        o=output, i=input)
    run_external_call(cmd)

    run_external_call("gff3_to_collapsed.py {o}.fa.{sp}.gff".format(o=output, sp=species_db))
    run_external_call("gff3_to_collapsed.py {i}.{sp}.gff".format(i=input, sp=species_db))
    run_external_call("gff3_to_collapsed.py {i}.{o}.gff".format(i=input, o=output))