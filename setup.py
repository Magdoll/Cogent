#!/usr/bin/env python

from setuptools import setup, find_packages

version = '1.9'

setup(name='Cogent',
      version=version,
      description='COding GENome reconstruction Tool using transcript sequences',
      author='Elizabeth Tseng',
      author_email='etseng@pacb.com',
#      install_requires=[
#          'numpy',
#          'networkx==1.10',
#          'scikit-image>=0.11.3',
#          'pulp',
#          'biopython',
#          'bx-python'
#      ],
      packages=['Cogent', 'Cogent.test'],
      package_dir = {'Cogent':'Cogent'},
      zip_safe=False,
      test_suite = "Cogent.test.test_cogent",
      scripts=['Cogent/run_mash.py',
               'Cogent/process_kmer_to_graph.py',
               'Cogent/gff3_to_collapsed.py',
               'Cogent/reconstruct_contig.py',
               'Cogent/generate_cmds_for_cogent.py',],
     )
