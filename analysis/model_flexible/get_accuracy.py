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
import multiprocessing as mp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from accuracy import *

selection_dictionary={"CBFB":[(1,156,"CBFB")],
                       "Vif":[(6,175,"Vif")],
                       "EloB":[(1,105,"EloB")],
                       "EloC":[(17,112,"EloC")],
                       "all_sel":[(1,156,"CBFB"),(6,175,"Vif"),(1,98,"EloB"),(17,112,"EloC")]}

nproc = 20
refrmf = '../../modeling/run_1/ini_all.rmf3'
clustering_dir = sys.argv[1]

rmf_A = '../../results/model_flexible/A_models_clust-1.rmf3' 
rmf_B = '../../results/model_flexible/B_models_clust-1.rmf3'
file_frames_A = '../../results/model_flexible/clustering_flexible/cluster.0.sample_A.txt'
file_frames_B = '../../results/model_flexible/clustering_flexible/cluster.0.sample_B.txt'

frames_A = []
frames_B = []
for line in open(file_frames_A,'r'):
    vals = line.split()
    frames_A.append(int(vals[0]))
for line in open(file_frames_B,'r'):
    vals = line.split()
    frames_B.append(int(vals[0])-len(frames_A))

print(len(frames_A))
print(len(frames_B))

frames_all = []

#frames=[0]*len(rmfs)
model=IMP.Model()
pr=IMP.pmi.analysis.Precision(model,
                            resolution=1,
                            selection_dictionary=selection_dictionary)
pr.set_precision_style('pairwise_rmsd')
pr.set_reference_structure(refrmf,0)
for fr in frames_A:
    pr.add_structure(rmf_name = rmf_A,rmf_frame_index = fr, structure_set_name='set0')
    frames_all.append(['A', fr])
for fr in frames_B:
    print(fr)
    try:
        pr.add_structure(rmf_name = rmf_B,rmf_frame_index = fr, structure_set_name='set0')
        frames_all.append(['B', fr])
    except:
        print('Missing frame', fr)
pr.set_reference_structure(refrmf,0)
vals = pr.get_rmsd_wrt_reference_structure_with_alignment('set0',('all_sel'))
print(vals['all_sel'],np.mean(vals['all_sel']['all_distances']) )
frames_all = np.array(frames_all)

all_vals = pd.DataFrame(columns = ['Sample', 'frame', 'RMSD'])
all_vals['Sample'] = frames_all[:,0]
all_vals['frame'] = frames_all[:,1]
all_vals['RMSD'] = np.array(vals['all_sel']['all_distances'])

print(all_vals.head())

print(all_vals[all_vals.RMSD == all_vals.RMSD.min()])

all_vals.to_csv('./../results/model_flexible/accuracy_cl-1_VCBC.dat',index=False)

exit()


