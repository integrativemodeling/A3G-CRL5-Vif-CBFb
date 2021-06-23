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


nproc = 1
top_dir =  sys.argv[1] 
cluster = sys.argv[2]
analys_dir = '../results/results_flexible/'
print(top_dir)

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')
#out_dirs = out_dirs[0:6]
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
                          nproc=nproc)

# Create dir
gsms_A_dir = f'{analys_dir}GSMs_cl{cluster}/sample_A'
gsms_B_dir = f'{analys_dir}GSMs_cl{cluster}/sample_B'

AT.create_gsms_dir(gsms_A_dir)
AT.create_gsms_dir(gsms_B_dir)


HA = AT.get_models_to_extract(f'{analys_dir}/selected_models_A_cluster{cluster}_detailed_random.csv')
HB = AT.get_models_to_extract(f'{analys_dir}/selected_models_B_cluster{cluster}_detailed_random.csv')

rmf_file_out_A = f'A_models_clust{cluster}.rmf3'
rmf_file_out_B = f'B_models_clust{cluster}.rmf3'

# Arguments for do_extract_models_single_rmf:
# HA :: Dataframe object from AT.get_models_to_extract()
# file_out :: The output rmf file
AT.do_extract_models_single_rmf(HA, 
                            rmf_file_out_A, # RMF file outputted for Sample A
                            top_dir,        # Top directory containing the PMI output folders
                            analys_dir,     # Analysis directory to write RMF and score files
                            scores_prefix = f'A_models_clust{cluster}')  # Prefix for the scores file


AT.do_extract_models_single_rmf(HB, rmf_file_out_B, top_dir, analys_dir, scores_prefix = f'B_models_clust{cluster}')

exit()


