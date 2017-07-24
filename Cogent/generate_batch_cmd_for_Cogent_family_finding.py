__author__ = 'etseng@pacb.com'


import os, sys
from csv import DictReader
from Bio import SeqIO

"""
USAGE:

generate_batch_cmd_for_Cogent_family_finding.py [csv] [dir] [cmd]

where [csv] is the preCluster_out.cluster_info.csv and
      [dir] is the preCluster_out/ directory with all the bins.

[csv] should contain <cluster> and <size>.

Generates a series of commands that can be run either locally or split for qsub.
"""

def generate_batch_cmds(csv_filename, dirname, cmd_filename, output_dir, cpus):
    cmd_f = open(cmd_filename, 'w')
    for r in DictReader(open(csv_filename), delimiter=','):
        cid = r['cluster']
        d2 = os.path.join(dirname, cid)
        if not os.path.exists(d2):
            print >> sys.stderr, "Directory {0} does not exist! Abort!".format(d2)
            sys.exit(-1)

        cmd_f.write("cd {0}\n".format(os.path.abspath(d2)))

        cmd_f.write("run_mash.py -k 30 --cpus={0} {1}/isoseq_flnc.fasta\n".format(cpus, os.path.abspath(d2)))

        cmd_f.write("process_kmer_to_graph.py " + \
                    "{0}/isoseq_flnc.fasta {0}/isoseq_flnc.fasta.s1000k30.dist {1} {2}\n".format(\
                        os.path.abspath(d2), output_dir, cid))


    cmd_f.close()




if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Generate batch commands for running Cogent family finding for each preCluster output bin")
    parser.add_argument("precluster_csv", help="Cluster CSV file (ex: preCluster.cluster_info.csv)")
    parser.add_argument("precluster_dir", help="preCluster out directory (ex: preCluster_out/)")
    parser.add_argument("output_dir", help="Output directory (of family finding)")
    parser.add_argument("--cpus", default=20, type=int, help="Number of CPUs (default: 20)")
    parser.add_argument("--cmd_filename", default='cmds', help="Output command filename (default: cmds)")
    args = parser.parse_args()

    generate_batch_cmds(args.precluster_csv, os.path.abspath(args.precluster_dir),
                        args.cmd_filename,
                        os.path.abspath(args.output_dir),
                        args.cpus)
