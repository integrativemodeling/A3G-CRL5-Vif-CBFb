#!/usr/bin/env python

'''
Different distance metrics to compare
structures in a ensemble
'''

import os
import random
import math
import itertools
#import numpy as np

import sys

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from compute_distance_metrics import *

sel = int(sys.argv[1])
print(sel)

alignment = [('Vif',6,175),
             ('EloB',1,105),
             ('EloC',7,112),
             ('CBFB',1,156),
             ('A3G',6,194),
             ('A3G',200,380),
             ('CUL5',11,302),
             ('CUL5',308,382),
             ('CUL5',405,515),
             ('CUL5',521,568),
             ('CUL5',574,687),
             ('CUL5',695,780),
             ('Rbx2',27,113)]



alignment_sel = [alignment[sel]]
print('Running (alignment): ', alignment_sel)

D = get_distance_metrics('../../results/model_flexible/clustering_flexible/',
                         0,
                         rmf_A='../../results/model_flexible/A_models_clust-1.rmf3',
                         rmf_B='../../results/model_flexible/B_models_clust-1.rmf3',
                         align_to = alignment_sel,
                         out_dir_name = 'RMSDs_5K_%s_%s'%(alignment_sel[0][0],alignment_sel[0][1]),
                         number_of_models = 5000)

D.compute_RMSD_all_versus_centroid()
del(D)

# nohup sh -c 'for i in {0..12}; do python get_RMSD_all.py $i > rmsd_$i.log; done > log' & 



