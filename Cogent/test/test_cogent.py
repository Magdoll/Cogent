__author__ = 'etseng@pacb.com'

import os, sys, logging, subprocess, shutil
import unittest
import tempfile
from Cogent.run_mash import main as run_mash_main
from Cogent.process_kmer_to_graph import family_finding, write_output_dirs
from Cogent.reconstruct_contig import run_Cogent_on_input
import Cogent.test.test_data as data
from Bio import SeqIO


class TestCogent(unittest.TestCase):
    def test_cogent(self):
        d = '/home/UNIXHOME/etseng/TEST/aloha/test'
        #d = tempfile.mkdtemp()
        fname = os.path.join(d ,'human_test.fa')

        with open(fname, 'w') as f:
            f.write(data.input + '\n')

        seqdict = dict((r.id.split()[0], r) for r in SeqIO.parse(open(fname),'fasta'))
        weightdict = dict((seqid, 1) for seqid in seqdict)
        
        dist_filename = run_mash_main(fname, kmer_size=30, sketch_size=1000, min_dist=0.95, chunk_size=100, cpus=1)

        output_prefix = os.path.join(d, 'human_test.k30')

        labels2_map = family_finding(\
            dist_filename, seqdict, output_prefix, \
            has_pbid=False, weight_threshold=0.05, \
            ncut_threshold=0.2)

        # check that the ncut result is consistent
        for k,v in labels2_map.iteritems():
            assert set(v) == set(data.labels2_map[k])
        output_dirs = write_output_dirs(labels2_map, seqdict, weightdict, output_prefix)


        for o in output_dirs:
            os.chdir(o)
            i = int(os.path.basename(o))
            subprocess.check_call("reconstruct_contig.py .", shell=True)
            assert os.path.exists('cogent2.fa')
            result = SeqIO.to_dict(SeqIO.parse(open('cogent2.fa'), 'fasta'))
            for _path, _seq in data.cogent2[i].iteritems():
                if str(_seq) != str(result[_path].seq):
                    raise Exception, "TEST FAILED! Results don't agree in cogent2.fa for {0}".format(o)
            os.chdir(d)

        # DONE remove test directory
        shutil.rmtree(d)

if __name__ == '__main__':
    unittest.main()





