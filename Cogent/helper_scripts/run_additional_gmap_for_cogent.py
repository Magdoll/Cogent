
# ------------------------------------------------- 
# Run additional GMAP jobs after Cogent is done
# NOTE: must first confirm that all Cogent jobs are done! check with other scripts
# -------------------------------------------------

import os, sys, glob

cogent_dir = "oneal_k30"
gmap_name = "blueberry_jersey_Falcon"
gmap_cmd = "~/bin/gmap -D /home/UNIXHOME/etseng/share/gmap_db_new/ -n 0 -t 6 -f gff3_gene --max-intronlength-ends 200000 --max-intronlength-middle 200000 --cross-species "
force_rerun = False

files = filter(lambda x: os.path.isdir(x), glob.glob(cogent_dir+'/*'))
for dd in files:
    if force_rerun or not os.path.exists(os.path.join(dd, "in.trimmed.fa.{0}.gff".format(gmap_name))):
        print "cd {0}; {1} -d {2} in.trimmed.fa > in.trimmed.fa.{2}.gff; gff3_to_collapsed.py in.trimmed.fa.{2}.gff; {1} -d {2} cogent2.fa > cogent2.fa.{2}.gff; gff3_to_collapsed.py cogent2.fa.{2}.gff; cd ../../".format(dd, gmap_cmd, gmap_name)
