#!/usr/bin/env python
import os, sys
from Cogent import GFF

input = sys.argv[1]
output = input[:input.rfind('.')] + '.collapsed.gff'

f = open(output, 'w')
reader = GFF.gmapGFFReader(input)
for r in reader: GFF.write_collapseGFF_format(f, r)
f.close()
