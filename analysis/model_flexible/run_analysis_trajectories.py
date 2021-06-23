import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 20
top_dir =  sys.argv[1] 
analys_dir = '../results/results_flexible/'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')
print('out', out_dirs)

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence

XLs_cutoffs = {'DSSO':30.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          analysis_dir = analys_dir,
                          nproc=nproc,
                          nskip=100)

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs)
AT.set_analyze_Binding_restraint()
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
#AT.set_analyze_Distance_restraint()
AT.set_analyze_Occams_restraint()

# Read stat files
AT.read_stat_files()
#AT.get_Psi_stats()
AT.write_models_info()
AT.hdbscan_clustering(['EV_sum', 'XLs_noA3G', 'XLs_A3G', 'Struct_Equiv_sum', 'BR_sum'],
                       min_cluster_size=10000,
                       min_samples=2,
                       skip=2)
AT.summarize_XLs_info()

exit()


