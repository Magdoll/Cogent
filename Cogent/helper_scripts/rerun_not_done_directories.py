
# ==========================================================
# Rerun any directories that didn't have COGENT.DONE flag
# ==========================================================

import os, sys, glob

cogent_dir = sys.argv[1]
additional_cmd = sys.argv[2:] # ex: -D /home/UNIXHOME/etseng/share/gmap_db_new/ -d Calypte_falcon --small_genome

dirs = filter(lambda x: os.path.isdir(x) and not os.path.exists(os.path.join(x,'COGENT.DONE')), glob.glob(cogent_dir+'/*'))

with open('cmd_rerun', 'w') as f:
    for d in dirs: f.write("reconstruct_contig.py {0} {1}\n".format(d, " ".join(additional_cmd)))
