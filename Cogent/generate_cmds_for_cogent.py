__author__ = 'etseng@pacb.com'

import os, sys

def main(args):
    for d in os.listdir(args.dirname):
        d2 = os.path.join(args.dirname, d)
        if not os.path.isdir(d2): continue
        fa = os.path.join(d2, 'in.fa')
        we = os.path.join(d2, 'in.weights')
        if not os.path.exists(fa):
            print >> sys.stderr, "{0} is missing! Please fix this first.".format(fa)
            sys.exit(-1)
        if not os.path.exists(we):
            print >> sys.stderr, "{0} is missing! Please fix this first.".format(we)
            sys.exit(-1)

        cmd = "reconstruct_contig.py " + d2
        if args.gmap_db_path is not None and args.gmap_species is not None:
            cmd += " -D {0} -d {1}".format(\
            os.path.join(args.dirname, d), args.gmap_db_path, args.gmap_species)
        if args.small_genome:
            cmd += " --small_genome"
        print cmd


if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("dirname")
    parser.add_argument("-D", "--gmap_db_path", help="GMAP database location, optional (ex: ~/share/gmap_db_new)")
    parser.add_argument("-d", "--gmap_species", help="GMAP species name, optinal (ex: hg19)")
    parser.add_argument("--small_genome", action="store_true", default=False, help="Genome size is smaller than 3GB (use gmap instead of gmapl)")

    args = parser.parse_args()

    main(args)