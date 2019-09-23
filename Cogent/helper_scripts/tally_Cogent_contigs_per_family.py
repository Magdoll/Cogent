__author__ = 'etseng@pacb.com'

import os, sys, glob

import cupcake.io.GFF as GFF
import tally_Cogent_results as sp


def main(cogent_dir, genome1, genome2, output_prefix, blastn_filename=None):

    dirs = filter(lambda x: os.path.isdir(x) and os.path.exists(os.path.join(x,'COGENT.DONE')), glob.glob(cogent_dir+'/*'))

    f1 = open(output_prefix+'.family_summary.txt', 'w')
    f2 = open(output_prefix+'.seq_summary.txt', 'w')
    f3 = open(output_prefix+'.err.log', 'w')


    f1.write("gene_family\tinput_size\tnum_Cogent_contigs\t")
    f1.write("num_genome1_contig\tgenome1_cov\tgenome1_acc\tgenome1_chimeric\tgenome1_contigs\t")
    f1.write("num_genome2_contig\tgenome2_cov\tgenome2_acc\tgenome2_chimeric\tgenome2_contigs")
    if blastn_filename is not None: f1.write("\tnum_blastn\tblastn_best\n")
    else: f1.write("\n")


    f2.write("seqid\tgene_family\tCogent_contig\tCogent_contig_acc\tscaffold\tnum_scaffold\tscaffold_coverage\tscaffold_acc")
    if blastn_filename is not None: f2.write("\blastn_best\n")
    else: f2.write("\n")
    for d in dirs:
        try:
            sp.tally_for_a_Cogent_dir(d, f1, f2, genome1, genome2, blastn_filename)
            print d
        except:
            f3.write("error with {0}. ignore\n".format(d))

    f1.close()
    f2.close()
    f3.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("cogent_dir")
    parser.add_argument("genome1")
    parser.add_argument("genome2")
    parser.add_argument("output_prefix")
    parser.add_argument("--blastn_filename", default=None, help="BLASTN output filename in each directory. If provided, will be parsed in evaluation summary file.")

    args = parser.parse_args()
    main(args.cogent_dir, args.genome1, args.genome2, args.output_prefix, args.blastn_filename)
