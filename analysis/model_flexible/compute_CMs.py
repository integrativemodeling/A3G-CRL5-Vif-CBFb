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
import itertools
import scipy
import scipy.spatial
import scipy.spatial.distance
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
#from accuracy import *

def get_resi_dict():
    resi_dict = {}
    model = IMP.Model()
    hier = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf_A)[0]
    prot_dictionary = get_particles_at_lowest_resolution(hier)

    for prot in prot_dictionary.keys():
        resi = []
        for p in prot_dictionary[prot]:
            pp = IMP.atom.Selection(p).get_selected_particles()[0]
            r = IMP.atom.Residue(pp)
            if 'bead' not in r.get_name():
                resi.append([r.get_index(),str(r.get_residue_type())])
            else:
                rn = r.get_name().split('-')
                rr = np.arange(int(rn[0]),int(rn[1].split('_')[0])+1,1)
                for ri in rr:
                    resi.append([ri,'BEA'])

        resi_dict[prot] = resi
    return resi_dict

def get_close_contacts(threshold = 0.2):
    resi_dict = get_resi_dict()

    cont_dict = {}
    for name1, name2 in itertools.combinations(resi_dict,2):
        contacts = []
        # Check if proteins have any contact
        if name1+'-'+name2 in cm_all.keys():
            p1 = name1
            p2 = name2
            mat = cm_all[name1+'-'+name2]
        elif name2+'-'+name1 in cm_all.keys():
            p1 = name2
            p2 = name1
            mat = cm_all[name2+'-'+name1]
        else:
            raise TypeError("No contact matrix for "+ name1, name2)

        if np.sum(mat) > 0:
            loc = np.where(mat> threshold)
            frq = mat[loc]
            
            k = 0
            for i in np.transpose(loc):
                contacts.append(resi_dict[p1][i[1]]+resi_dict[p2][i[0]]+[frq[k]])
                k += 1
            # Sort array
            contacts = np.array(contacts)
            sort_contacts = contacts[contacts[:,4].astype(float).argsort()]
            np.savetxt(os.path.join(
                out_dir,'contacts_%s_%s.dat'%(p1,p2)),sort_contacts[::-1],fmt='%s')
            cont_dict[p1+'-'+p2] = contacts

def get_coords_array(particles_list):
    '''
    Get all beads coordinates and radii
    '''
    coords = []
    radii = []
    for p in particles_list:
        residue_indexes = IMP.pmi.tools.get_residue_indexes(p)
        if len(residue_indexes) !=0 :
            for res in range(min(residue_indexes), max(residue_indexes) + 1):
                d = IMP.core.XYZR(p)
                crd = np.array([d.get_x(), d.get_y(), d.get_z()])
                coords.append(crd)
                radii.append(d.get_radius())
                
    return np.array(coords), np.array(radii)

def get_particles_at_lowest_resolution(hier):
    '''
        Read rmf3 file and return only coordinates of beads at
        lowest resolution
    '''
    
    particles_dict = {}
    for mol in IMP.atom.get_by_type(hier,
                                    IMP.atom.MOLECULE_TYPE):
         
        copy = IMP.atom.Copy(mol).get_copy_index()
        if (len(selection) == 0) or (mol.get_name() in selection):
            sel = IMP.atom.Selection(mol,resolution=1)
            particles_dict[mol.get_name()+'.'+str(copy)] = sel.get_selected_particles()
            

    return particles_dict

def get_contactmap_pair( particles_1, particles_2):
    '''
        Given two proteins, computes the contact map
    '''
    
    coords1, radii1 = get_coords_array(particles_1)
    coords2, radii2  = get_coords_array(particles_2)
    distances = scipy.spatial.distance.cdist(coords1, coords2)
    distances = (distances - radii2).T - radii1
    #contact_map = np.where((distances>0) & (distances <= self.cutoff), 1.0, 0)
    contact_map = np.where((distances <= cutoff), 1.0, 0)
    return contact_map

def update_CMs(rmf_file, frames):
    model=IMP.Model()
    for fr in frames:
        hier = IMP.pmi.analysis.get_hiers_from_rmf(model,fr,rmf_file)[0]
        prot_dictionary = get_particles_at_lowest_resolution(hier)
        for name1, name2 in itertools.combinations_with_replacement(prot_dictionary.keys(),2):

            if name1 == name2:
                particles_1 = prot_dictionary[name1]
                particles_2 = prot_dictionary[name2]
                cm = get_contactmap_pair(particles_1,particles_2)
                if str(name1+'-'+name2) not in cm_all.keys():
                    cm_all[name1+'-'+name2] = cm
                else:
                    cm_all[name1+'-'+name2] = cm_all[name1+'-'+name2] + cm
            else:
                particles_1 = prot_dictionary[name1]
                particles_2 = prot_dictionary[name2]
                cm = get_contactmap_pair(particles_1, particles_2)
                if str(name1+'-'+name2) not in cm_all.keys():
                    cm_all[name1+'-'+name2] = cm
                else:
                    cm_all[name1+'-'+name2] = cm_all[name1+'-'+name2] + cm
    del model

######################################################
# MAIN
######################################################

rmf_A = '../../results/model_flexible/A_models_clust-1.rmf3'
rmf_B = '../../results/model_flexible/B_models_clust-1.rmf3'
file_frames_A = '../../results/model_flexible/clustering_flexible/cluster.0.sample_A.txt'
file_frames_B = './../results/model_flexible/clustering_flexible/cluster.0.sample_B.txt'

selection = ['A3G','Vif']
cutoff = 14.0
nproc = 8
nproc2 = int(nproc/2)

frames_A = []
frames_B = []
for line in open(file_frames_A,'r'):
    vals = line.split()
    frames_A.append(int(vals[0]))
for line in open(file_frames_B,'r'):
    vals = line.split()
    frames_B.append(int(vals[0])-len(frames_A))

print('Frames A', len(frames_A))
print('Frames B', len(frames_B))


manager = mp.Manager()
cm_all = manager.dict()


# Define an output queue
output = mp.Queue()

# Divide frames into groups
ND_A = int(np.ceil(len(frames_A)/float(nproc2)))
ND_B = int(np.ceil(len(frames_B)/float(nproc2)))
frames_A_dict = {}
frames_B_dict = {}

print(ND_A, ND_B)

for k in range(nproc2-1):
    frames_A_dict[k] = frames_A[(k*ND_A):(k*ND_A+ND_A)]
frames_A_dict[nproc2-1] = frames_A[((nproc2-1)*ND_A):(len(frames_A))]

for k in range(nproc2-1):
    frames_B_dict[k] = frames_B[(k*ND_B):(k*ND_B+ND_B)]
frames_B_dict[nproc2-1] = frames_B[((nproc2-1)*ND_B):(len(frames_B))]

# Setup a list of processes that we want to run
processes = [mp.Process(target=update_CMs,
                        args=(rmf_A, frames_A_dict[x])) for x in range(nproc2)] + \
             [mp.Process(target=update_CMs,
                        args=(rmf_B, frames_B_dict[x])) for x in range(nproc2)]



# Run processes
for p in processes:
    p.start()

# Exit the completed processes
for p in processes:
    p.join()

# Normalize matrices
for key in cm_all.keys():
    cm_all[key] = cm_all[key]/float(len(frames_A)+len(frames_B))
            
# Save matrices
out_dir = '../../results/model_flexible/clustering_flexible/CMs/'
if os.path.exists(out_dir):
    print('Overwriting '+ out_dir)
else:
    os.makedirs(out_dir)

for key in cm_all.keys():
    np.savetxt(os.path.join(
        out_dir,'ContMap_%s.dat'%(key)),np.array(cm_all[key]),fmt='%s')
print(cm_all)

get_close_contacts()
