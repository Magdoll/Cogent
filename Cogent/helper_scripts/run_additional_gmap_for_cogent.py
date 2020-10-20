
# ------------------------------------------------- 
# Run additional GMAP jobs after Cogent is done
# NOTE: must first confirm that all Cogent jobs are done! check with other scripts
# -------------------------------------------------

import os, sys, glob

GMAP_CMD_GFF = "-n 0 -t 1 -f gff3_gene --max-intronlength-ends 200000 --max-intronlength-middle 200000 --cross-species "
GMAP_CMD_SAM = "-n 0 -t 1 -f samse --max-intronlength-ends 200000 --max-intronlength-middle 200000 --cross-species "


def main(args):
    if args.small_genome: prog = 'gmap'
    else: prog = 'gmapl'
    files = filter(lambda x: os.path.isdir(x), glob.glob(args.cogent_dir+'/*'))

    for dd in files:
        if args.force_rerun or not os.path.exists(os.path.join(dd, "in.trimmed.fa.{0}.gff".format(args.gmap_species))):
            if args.sam:
                print(prog + " {cmd} -D {db} -d {name} {dd}/in.trimmed.fa > {dd}/in.trimmed.fa.{name}.sam".format(\
                    cmd=GMAP_CMD_SAM, db=args.gmap_db_path, dd=dd, name=args.gmap_species))
                print(prog + " {cmd} -D {db} -d {name} {dd}/cogent2.fa > {dd}/cogent2.fa.{name}.sam".format(\
                    cmd=GMAP_CMD_SAM, db=args.gmap_db_path, dd=dd, name=args.gmap_species))
            print(prog + " {cmd} -D {db} -d {name} {dd}/in.trimmed.fa > {dd}/in.trimmed.fa.{name}.gff; gff3_to_collapsed.py {dd}/in.trimmed.fa.{name}.gff".format(\
                cmd=GMAP_CMD_GFF, db=args.gmap_db_path, dd=dd, name=args.gmap_species))
            print(prog + " {cmd} -D {db} -d {name} {dd}/cogent2.fa > {dd}/cogent2.fa.{name}.gff; gff3_to_collapsed.py {dd}/cogent2.fa.{name}.gff".format( \
                cmd=GMAP_CMD_GFF, db=args.gmap_db_path, dd=dd, name=args.gmap_species))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Run additional GMAP alignment for each Cogent directory.")
    parser.add_argument("cogent_dir")
    parser.add_argument("-D", "--gmap_db_path", help="GMAP database location")
    parser.add_argument("-d", "--gmap_species", help="GMAP species name")
    parser.add_argument("--small_genome", action="store_true", default=False, help="Genome size is smaller than 3GB (use gmap instead of gmapl)")
    parser.add_argument("--force_rerun", action="store_true", default=False, help="Re-run GMAP even if output file already exists. (default: off)")
    parser.add_argument("--sam", action="store_true", default=False, help="Also output SAM (default: off)")

    args = parser.parse_args()
    main(args)