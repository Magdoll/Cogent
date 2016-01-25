#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='Cogent',
      version='1.0',
      description='COding GENome reconstruction Tool using transcript sequences',
      author='Elizabeth Tseng',
      author_email='etseng@pacb.com',
      install_requires=[
          'numpy',
          'networkx==1.10',
          'scikit-image>=0.11.3',
          'pulp',
          'biopython',
          'bx-python'
      ],
      packages=['Cogent'],
      package_dir = {'Cogent':'Cogent'},
      zip_safe=False,
      scripts=['Cogent/process_kmer_to_graph.py',
               'Cogent/gff3_to_collapsed.py',
               'Cogent/reconstruct_contig.py'],
     )
