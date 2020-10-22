#!/usr/bin/env python3
__author__ = 'etseng@pacb.com'
import os, sys, subprocess
from multiprocessing import Pool, TimeoutError
from Bio import SeqIO
from Cogent.__init__ import get_version
from Cogent.sanity_checks import sanity_check_mash_exists

"""
Wrapper for running mash (https://github.com/marbl/mash)

1. Split the input files into sizable chunks
2. Make sketches for each chunk
3. Do pairwise dist for each chunk where chunk_i < chunk_j
4. Combine the results
"""


def split_input(fasta_filename, chunk_size):
    i = 0
    split_i = 0
    f = open(fasta_filename + '.' + str(split_i), 'w')
    files = [f.name]
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        print(i)
        i += 1
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        if i >= chunk_size:
            i = 0
            split_i += 1
            f.close()
            f = open(fasta_filename + '.' + str(split_i), 'w')
            files.append(f.name)
    f.close()
    return [x for x in files if os.stat(x).st_size > 0]

def run_sketch(fasta_filename, kmer_size, sketch_size):
    cmd = "mash sketch -i -k {k} -s {s} -o {i}.s{s}k{k} {i}".format(i=fasta_filename, k=kmer_size, s=sketch_size)
    if subprocess.check_call(cmd, shell=True)!=0:
        print("ERROR running:", cmd, file=sys.stderr)
        sys.exit(-1)
    return "{i}.s{s}k{k}.msh".format(i=fasta_filename, k=kmer_size, s=sketch_size)


def run_dist(sketch1, sketch2, min_dist):
    cmd = "mash dist -d {0} {1} {2} > {1}.{2}.dist".format(min_dist, sketch1, sketch2)
    if subprocess.check_call(cmd, shell=True)!=0:
        print("ERROR running:", cmd, file=sys.stderr)
        sys.exit(-1)
    print("{0}.{1}.dist completed".format(sketch1, sketch2), file=sys.stderr)
    return "{0}.{1}.dist".format(sketch1, sketch2)

def main(fasta_filename, kmer_size, sketch_size, min_dist, chunk_size, cpus):

    olddir = os.getcwd()
    dirname = os.path.dirname(os.path.abspath(fasta_filename))
    fasta_filename = os.path.basename(fasta_filename)

    os.chdir(dirname)

    sanity_check_mash_exists()

    o = "{i}.s{s}k{k}.dist".format(i=fasta_filename, k=kmer_size, s=sketch_size)
    if os.path.exists(o):
        print("Expected output {0} already exists. Just use it.".format(o), file=sys.stderr)
        os.chdir(olddir)
        return os.path.join(dirname, o)

    inputs = split_input(fasta_filename, chunk_size)
    N = len(inputs)

    pool = Pool(processes=cpus)

    # get the sketch files
    results = [ pool.apply_async(run_sketch, args=(inputs[i], kmer_size, sketch_size,)) for i in range(N) ]
    files = [ r.get() for r in results ]
    pool.close()
    pool.join()
    # get the dist files
    pool = Pool(processes=cpus)
    results = []
    for i in range(N):
        for j in range(i, N):
            results.append(pool.apply_async(run_dist, args=(files[i], files[j], min_dist,)))
            #run_dist(files[i], files[j], min_dist)
    pool.close()
    pool.join()
    outputs = [ r.get() for r in results ]


    # combine the dist files
    with open(o, 'w') as f:
        for file in outputs:
            for line in open(file):
                f.write(line)

    # clean up the split files
    for file in inputs: os.remove(file)
    for file in files: os.remove(file)
    for file in outputs: os.remove(file)
    print("Output written to:", os.path.join(dirname, o), file=sys.stderr)

    os.chdir(olddir)
    return os.path.join(dirname, o)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fasta_filename")
    parser.add_argument("-k", "--kmer_size", default=30, help="K-mer size (default: 30)")
    parser.add_argument("-s", "--sketch_size", default=1000, help="Sketch size (default: 1000)")
    parser.add_argument("-d", "--min_dist", default=0.95, help="Minimum distance (default: 0.95)")
    parser.add_argument("--chunk_size", default=1000, type=int, help="Chunk size")
    parser.add_argument("--cpus", default=1, type=int, help="# of CPUs (default: 1)")
    parser.add_argument('--version', action='version', version='%(prog)s ' + str(get_version()))

    args = parser.parse_args()
    main(args.fasta_filename, args.kmer_size, args.sketch_size, args.min_dist, args.chunk_size, args.cpus)
