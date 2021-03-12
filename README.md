<a href="https://programs.pacificbiosciences.com/l/1652/2021-03-12/42p2d7">![](https://github.com/Magdoll/images_public/blob/master/banners/test.png?raw=true)</a>

# Cogent: COding GENome reconstruction Tool

Current version: 8.0.0

Cogent is a tool for reconstructing the coding genome using high-quality full-length transcriptome sequences. It is designed to be used on [Iso-Seq data](https://github.com/PacificBiosciences/cDNA_primer/wiki) and in cases where there is no reference genome or the ref genome is highly incomplete. 

See a [recent presentation](https://www.dropbox.com/s/mn6hwhguh0pqceu/20160106_Cogent_developers_conference_slides_Cuttlefish.pdf?dl=0) on Cogent being applied to the Cuttlefish Iso-Seq data. 

[Cogent preliminary draft paper (updated 2016Dec version)](https://www.dropbox.com/s/kz0gi7qg0w82k9a/20161026_Cogent_manuscript_forGitHub.pdf?dl=0), [Supplementary](https://www.dropbox.com/s/37412o8glvnfhf9/20161026_Cogent_ManuscriptPlusSupplement_forGitHub.pdf?dl=0)

Please see [wiki](https://github.com/Magdoll/Cogent/wiki) for details on usage.


## version updates

2020.10.21 updated to version 8.0.0. Added `--randomly_resolve_ambiguous_dangling` (`-R`) option to `reconstruct_contigs.py`

2020.09.23 updated to version 7.0.0. networkx now on v2.5.

2020.08.28 updated to version 6.1.0. fixed `os.chdir` bug in `reconstruct_contigs.py`

2019.12.30 updated to version 6.0.0. Python 3.7 version! Switched off SSW to using parasail instead.

2019.09.25 updated to version 4.0.0. Final 4.0.x is for Py2. Added tallying scripts to helper scripts.

<details>
   <summary>Click to see older version logs</summary>
   
    2018.10.30 updated to version 3.9. added `--dun_trim_sequence` option to reconstruct and also fixed GFF.py
    
    2018.10.15 updated to version 3.8. fix `sam_to_gff3.py` indentation
    
    2018.10.12 updated to version 3.7. reconstruct now will accept kmer size `-k` greater than 200.
    
    2018.05.21 updated to version 3.5. Fixed bug for reachability.
    
    2018.05.14 updated to version 3.4. Fixed bug for not adding weights in `find_bubbles()`.
    
    2018.04.16 updated to version 3.2. Fixed edge case in untangling homopolymer.
    
    2018.03.09 updated to version 3.1. Replaced all GMAP with minimap2!! Rest of the algorithm remains same, but selection of final cogent2 (from cogent.fa) now can be different. Increased cogent -> cogent2 selection stringency to >= 98% cov AND >= 98% identity.
    
    2017.10.23 updated to version 2.1. Fixed test. 
    
    2017.07.24 updated to version 2.0. MAJOR CHANGE in adding preclustering as an option to speed up family finding.
    
    2017.06.21 updated to version 1.9. Fixed LP solver bug from multiple optimal solutions. Added `--output_prefix` to append to cogent2 IDs. 
    
    2017.03.06  updated to version 1.7. Automatically removed GMAP DBs to reduce space usage. Recursive handling of large inputs (combine/post-combine).
    
    2017.01.11  updated to version 1.6. Fixed the changes in test.
    
    2016.11.22  updated to version 1.5. Added auto k-mer size increment (up to k=200) in cycle detection.
    
    2016.10.20  updated to version 1.4. Added features that detects cycles and tries to update k-mers to larger sizes to accommodate for that. Also found a bug in bubble detection that caused errors. Both are fixed.
    
    2016.03.30  updated to version 1.3. Important bug in `splice_graph.contract_sinks` fixed where it previously was accidentally contracting sinks with predecessors that had multiple outgoing edges, which causes incorrect reconstructions in cases where there are a lot of isoforms with alternative 3' ends.
    
    2016.03.21  updated to version 1.2. Fixed more minor bugs related to edge cases esp caused by non-obvious graph cycles.
    
    2016.03.02  updated to version 1.1. Found several cases where border cases caused program crash. Fixed.

</details>
