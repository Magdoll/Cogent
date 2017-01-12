# Cogent: COding GENome reconstruction Tool

Current version: 1.6

Cogent is a tool for reconstructing the coding genome using high-quality full-length transcriptome sequences. It is designed to be used on [Iso-Seq data](https://github.com/PacificBiosciences/cDNA_primer/wiki) and in cases where there is no reference genome or the ref genome is highly incomplete. 

See a [recent presentation](https://www.dropbox.com/s/mn6hwhguh0pqceu/20160106_Cogent_developers_conference_slides_Cuttlefish.pdf?dl=0) on Cogent being applied to the Cuttlefish Iso-Seq data. 

A [preliminary draft paper](https://www.dropbox.com/s/gmndqsihsv7i4gt/20151023_Cogent_RECOMB2016_draft_LizTseng_v4.pdf?dl=0) is available. Note that there have been minor implementation changes since the draft was written, however the concept and overall methodology remains.

Please see [wiki](https://github.com/Magdoll/Cogent/wiki) for details on usage.


## version updates
2017.01.11  updated to version 1.6. Fixed the changes in test.

2016.11.22  updated to version 1.5. Added auto k-mer size increment (up to k=200) in cycle detection.

2016.10.20  updated to version 1.4. Added features that detects cycles and tries to update k-mers to larger sizes to accommodate for that. Also found a bug in bubble detection that caused errors. Both are fixed.

2016.03.30  updated to version 1.3. Important bug in `splice_graph.contract_sinks` fixed where it previously was accidentally contracting sinks with predecessors that had multiple outgoing edges, which causes incorrect reconstructions in cases where there are a lot of isoforms with alternative 3' ends.

2016.03.21  updated to version 1.2. Fixed more minor bugs related to edge cases esp caused by non-obvious graph cycles.

2016.03.02  updated to version 1.1. Found several cases where border cases caused program crash. Fixed.

