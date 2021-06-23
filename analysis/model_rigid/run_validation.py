from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import random
import numpy as np
import pandas as pd
import sys
import os
import multiprocessing as mp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from validation import *

analysis_dir = sys.argv[1]
cluster = sys.argv[2]

clustering_dir = os.path.join(
    analysis_dir,'clustering_cl'+str(cluster))
scores_sample_A = os.path.join(
    analysis_dir,'selected_models_A_cluster'+str(cluster)+'_detailed.csv')
scores_sample_B = os.path.join(
    analysis_dir,'selected_models_B_cluster'+str(cluster)+'_detailed.csv')

print('clustering dir: ', clustering_dir)

#analysis_dir = '/home/ignacia/Research/pE-MAP/mod_RNAP_XLs_prob_04/analys_test/'
#clustering_dir = '/home/ignacia/Research/pE-MAP/mod_RNAP_XLs_prob_04/analys_test/clustering_cl3/'
#scores_sample_A = analysis_dir + 'selected_models_A_cluster3_detailed.csv'
#scores_sample_B = analysis_dir + 'selected_models_B_cluster3_detailed.csv'

XLs_cutoffs = {'DSSO':30.0}

V = ValidationModels(analysis_dir,
                     clustering_dir,
                     scores_sample_A,
                     scores_sample_B,
                     XLs_cutoffs)


V.get_excluded_volume_satisfaction()
