
import os, sys
import run_Cogent as x
os.chdir(sys.argv[1])
x.run_gmap_for_final_GFFs()
os.chdir('../')
