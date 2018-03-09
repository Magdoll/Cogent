import os, subprocess, logging
from collections import defaultdict
from Bio import SeqIO
from Cogent import BioReaders
#from Cogent import MinimapIO
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


def run_minimap2(ref='cogent.fa', infile='in.trimmed.fa', format='PAF'):
    """
    Map input to cogent contig using minimap2. Replacement of run_gmap().
    NOTE: output will be in PAF format!
    """
    if format=='PAF':
        outfile = infile + '.paf'
        run_external_call("minimap2 -x splice -t 1 {d} {i} > {o}".format(d=ref, i=infile, o=outfile))
    elif format=='SAM':
        outfile = infile + '.sam'
        run_external_call("minimap2 -ax splice -t 1 {d} {i} > {o}".format(d=ref, i=infile, o=outfile))
    else:
        raise Exception, "Unrecognized minimap2 output format: {0}. Abort!".format(format)
    return outfile

def post_minimap2_processing(ref='cogent.fa', sam='in.trimmed.fa.sam', output_prefix='cogent2', seqrecs=[]):
    good_for = defaultdict(lambda: [])
    reader = BioReaders.GMAPSAMReader(sam, True, query_len_dict=dict(((r.id, len(r.seq)) for r in seqrecs)))
    for r in reader:
        assert r.sID.startswith('path')  # chr should be path0, path1, etc
        assert 0 < r.qCoverage <= 1
        assert 0 < r.identity <= 1
        if r.qCoverage >= 0.98 and r.identity >= 0.98: good_for[r.qID].append(int(r.sID[4:]))

    touse = []
    if len(good_for) == 0:
        log.warning("[BUG] good_for in post_minimap2_processing is empty. Probably from cycles. CHECK!")
    else:
        N = max(max(v) for v in good_for.itervalues())+1
        try:
            prob = make_into_lp_problem(good_for.items(), N, add_noise=False)
            prob.solve()
        except:
            prob = make_into_lp_problem(good_for.items(), N, add_noise=True)
            prob.solve()
        for v in prob.variables():
            log.debug("{0} = {1}".format(v.name, v.varValue))
            if v.varValue == 1: touse.append(int(v.name))


    with open(output_prefix + '.fa', 'w') as f:
        for r in SeqIO.parse(open(ref),'fasta'):
            if int(r.id[4:]) in touse:
                f.write(">{0}\n{1}\n".format(r.id, r.seq))
        # if there are some sequences that didn't map (possibly from cycles)
        # then just use THEMSELVES
        fake_path_i = max(touse)+1 if len(touse) >= 1 else 0
        for r in seqrecs:
            if r.id not in good_for:
                log.warning("[BUG] {0} is not fully mapped to cogent in minimap2. \
                Likely cycle issues. Use itself in output.".format(r.id))
                f.write(">path{0}\n{1}\n".format(fake_path_i, r.seq))
                fake_path_i += 1

def run_minimap2_for_final_SAM(input='in.trimmed.fa', output='cogent2.fa', ref='genome.fasta', species_name='speciesX'):
    # map in.trimmed.fa to genome
    run_external_call("minimap2 -t 1 -ax splice -uf --secondary=no {r} {i} > {i}.{sp}.sam".format(r=ref, i=input, sp=species_name))
    # map cogent2.fa to genome
    run_external_call("minimap2 -t 1 -ax splice -uf --secondary=no {r} {o} > {o}.{sp}.sam".format(r=ref, o=output, sp=species_name))
    # map in.trimmed.fa to cogent2
    run_external_call("minimap2 -t 1 -ax splice -uf --secondary=no {o} {i} > {i}.cogent2.sam".format(o=output, i=input))

#
# def run_gmap(dbname='cogent', infile='in.trimmed.fa'):
#     run_external_call("gmap_build -D . -d {0} {0}.fa".format(dbname))
#     run_external_call("gmap -D . -d {0} -n 100 -f gff3_gene -t 1 {1} > {1}.gff".format(dbname, infile))
#
# def cleanup_gmap(dbname='cogent'):
#     run_external_call("rm -rf " + dbname)
#
# def post_gmap_processing(db_name='cogent', gff_filename='in.trimmed.fa.gff', output_prefix='cogent2', seqrecs=[]):
#     good_for = defaultdict(lambda: [])
#     reader = GFF.gmapGFFReader(gff_filename)
#     for r in reader:
#         assert r.chr.startswith('path')  # chr should be path0, path1, etc
#         if r.coverage >= 98.: good_for[r.seqid].append(int(r.chr[4:]))
#
#     touse = []
#     if len(good_for) == 0:
#         log.warning("[BUG] good_for in post_gmap_processing is empty. Probably from cycles. CHECK!")
#     else:
#         N = max(max(v) for v in good_for.itervalues())+1
#         try:
#             prob = make_into_lp_problem(good_for.items(), N, add_noise=False)
#             prob.solve()
#         except:
#             prob = make_into_lp_problem(good_for.items(), N, add_noise=True)
#             prob.solve()
#         for v in prob.variables():
#             log.debug("{0} = {1}".format(v.name, v.varValue))
#             if v.varValue == 1: touse.append(int(v.name))
#
#
#     with open(output_prefix + '.fa', 'w') as f:
#         for r in SeqIO.parse(open(db_name + '.fa'),'fasta'):
#             if int(r.id[4:]) in touse:
#                 f.write(">{0}\n{1}\n".format(r.id, r.seq))
#         # if there are some sequences that didn't map (possibly from cycles)
#         # then just use THEMSELVES
#         fake_path_i = max(touse)+1 if len(touse) >= 1 else 0
#         for r in seqrecs:
#             if r.id not in good_for:
#                 log.warning("[BUG] {0} is not fully mapped to cogent in GMAP. \
#                 Likely cycle issues. Use itself in output.".format(r.id))
#                 f.write(">path{0}\n{1}\n".format(fake_path_i, r.seq))
#                 fake_path_i += 1
#
# def run_gmap_for_final_GFFs(small_genome=False, gmap_db_path='~/share/gmap_db_new/', species_db='cuttlefish',\
#                             input='in.trimmed.fa', output='cogent2'):
#
#     if small_genome:
#         prog = 'gmap'
#     else:
#         prog = 'gmapl'
#     cmd = prog + " -D {p} -d {sp} -f gff3_gene --cross-species --max-intronlength-ends 200000 -n 0 {o}.fa > {o}.fa.{sp}.gff".format(\
#         p=gmap_db_path, sp=species_db, o=output)
#     run_external_call(cmd)
#     cmd = prog + " -D {p} -d {sp} -f gff3_gene --cross-species --max-intronlength-ends 200000 -n 0 {i} > {i}.{sp}.gff".format(\
#         p=gmap_db_path, sp=species_db, i=input)
#     run_external_call(cmd)
#     cmd = "gmap_build -D . -d {o} {o}.fa".format(o=output)
#     run_external_call(cmd)
#     cmd = "gmap -D . -d {o} -f gff3_gene -n 0 {i} > {i}.{o}.gff".format(\
#         o=output, i=input)
#     run_external_call(cmd)
#
#     run_external_call("gff3_to_collapsed.py {o}.fa.{sp}.gff".format(o=output, sp=species_db))
#     run_external_call("gff3_to_collapsed.py {i}.{sp}.gff".format(i=input, sp=species_db))
#     run_external_call("gff3_to_collapsed.py {i}.{o}.gff".format(i=input, o=output))