__author__ = 'etseng@pacb.com'

import os, sys, glob
from csv import DictWriter

import cupcake.io.GFF as GFF
import tally_Cogent_results as sp



def main(cogent_dir, genome1, output_prefix, genome2=None, blastn_filename=None):

    dirs = [x for x in glob.glob(cogent_dir+'/*') if os.path.isdir(x) and os.path.exists(os.path.join(x,'COGENT.DONE'))]

    FIELDS_fam = ['gene_family', 'input_size', 'num_Cogent_contigs',
              'num_genome_contig', 'genome_cov', 'genome_acc',
              'genome_chimeric', 'genome_contigs']

    if genome2 is not None:
        FIELDS_fam += ['num_genome2_contig', 'genome2_cov', 'genome2_acc',
                'genome2_chimeric', 'genome2_contigs']

    if blastn_filename is not None:
        FIELDS_fam += ['num_blastn', 'blastn_best']

    f1 = open(output_prefix+'.family_summary.txt', 'w')
    f2 = open(output_prefix+'.seq_summary.txt', 'w')
    f3 = open(output_prefix+'.err.log', 'w')

    writer1 = DictWriter(f1, FIELDS_fam, delimiter='\t')
    writer1.writeheader()

    FIELDS_seq = ['seqid', 'gene_family', 'Cogent_contig', 'Cogent_contig_acc', 'scaffold',
                  'num_scaffold', 'scaffold_coverage', 'scaffold_acc']
    if blastn_filename: FIELDS_seq += ['blastn_best']

    writer2 = DictWriter(f2, FIELDS_seq, delimiter='\t')
    writer2.writeheader()
    for d in dirs:
        #try:
        sp.tally_for_a_Cogent_dir(d, writer1, writer2, genome1, genome2, blastn_filename)
        print(d)
        #except:
        #    f3.write("error with {0}. ignore\n".format(d))

    f1.close()
    f2.close()
    f3.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("cogent_dir", help="Cogent directory name")
    parser.add_argument("genome", help="Genome name")
    parser.add_argument("output_prefix")
    parser.add_argument("--genome2", help="Optional second genome name")
    parser.add_argument("--blastn_filename", default=None, help="(optional) BLASTN output filename in each directory. If provided, will be parsed in evaluation summary file.")

    args = parser.parse_args()
    main(args.cogent_dir, args.genome, args.output_prefix, args.genome2, args.blastn_filename)
