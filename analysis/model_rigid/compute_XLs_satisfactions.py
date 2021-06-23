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

def read_XLs_database(file):
    D = pd.read_csv(file)
    return D

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
        elif name2+'-'+name1 in self.cm_all.keys():
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
    print(coords1, coords2)
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
        print(fr, prot_dictionary)
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


def update_XLs_distances(rmf_file, frames, XLs_database, dist_all, XLs_not=[], dist_not=None):
    model=IMP.Model()
    for fr in frames:
        coor_A = []
        coor_B = []
        try:
            hier = IMP.pmi.analysis.get_hiers_from_rmf(model,fr,rmf_file)[0]
        except:
            continue
        for i, row in XLs_database.iterrows():
            s1 = IMP.atom.Selection(hier, molecule=row['Protein1'],residue_index=row['AbsPos1']).get_selected_particles()
            s2 = IMP.atom.Selection(hier, molecule=row['Protein2'],residue_index=row['AbsPos2']).get_selected_particles()
            coor_A.append(IMP.core.XYZ(s1[0]).get_coordinates())
            coor_B.append(IMP.core.XYZ(s2[0]).get_coordinates())
        coor_A = np.array(coor_A)
        coor_B = np.array(coor_B)
        dist_all[fr] = np.sqrt(np.sum((coor_A-coor_B)**2, axis=1))

        if len(XLs_not)>0:
            ncoor_A = []
            ncoor_B = []
            for i, row in XLs_not.iterrows():
                s1 = IMP.atom.Selection(hier, molecule=row['Protein1'],residue_index=row['AbsPos1']).get_selected_particles()
                s2 = IMP.atom.Selection(hier, molecule=row['Protein2'],residue_index=row['AbsPos2']).get_selected_particles()
                ncoor_A.append(IMP.core.XYZ(s1[0]).get_coordinates())
                ncoor_B.append(IMP.core.XYZ(s2[0]).get_coordinates())
            ncoor_A = np.array(ncoor_A)
            ncoor_B = np.array(ncoor_B)
            dist_not[fr] = np.sqrt(np.sum((ncoor_A-ncoor_B)**2, axis=1))
            
    del model
    
######################################################
# MAIN
######################################################

rmf_A = 'analys/A_models_clust-1.rmf3'
rmf_B = 'analys/B_models_clust-1.rmf3'
file_frames_A = 'analys/clustering_cl-1/cluster.0.sample_A.txt'
file_frames_B = 'analys/clustering_cl-1/cluster.0.sample_B.txt'

nproc = 16
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

#frames_A = random.sample(frames_A, 500)
#frames_B = random.sample(frames_B, 500)

manager = mp.Manager()
dist_all_A = manager.dict()
dist_all_B = manager.dict()
#dist_not_A = manager.dict()
#dist_not_B = manager.dict()

XLs_sampled = read_XLs_database('Interlinks_A3G_Vif_CRL5_unique_modeling.csv')

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

#######################
# Used for modeling
#######################
out = open('XLs_satisfaction.dat','w')


# Setup a list of processes that we want to run
processes = [mp.Process(target=update_XLs_distances,
                        args=(rmf_A, frames_A_dict[x], XLs_sampled, dist_all_A)) for x in range(nproc2)] + \
             [mp.Process(target=update_XLs_distances,
                        args=(rmf_B, frames_B_dict[x], XLs_sampled, dist_all_B)) for x in range(nproc2)]

# Run processes
for p in processes:
    p.start()

# Exit the completed processes
for p in processes:
    p.join()

#update_XLs_distances(rmf_A, frames_A, XLs_sampled, dist_all_A)
#update_XLs_distances(rmf_B, frames_B, XLs_sampled, dist_all_B)

# Compute satisfaction, used for modeling
DA = [v for k,v in dist_all_A.items()]
DA = np.array(DA)
DB = [v for k,v in dist_all_B.items()]
DB = np.array(DB)
DD = np.vstack((DA, DB))
print('all distances', np.shape(DD))
Dmin = np.min(DD,axis=0)
exit()
#print(np.shape(Dmin), len(Dmin), np.shape(DD[0]))
satif = np.sum([1 for d in Dmin if d<35.])/np.shape(DD[0])

print(np.sum([1 for d in Dmin if d<35.])/np.shape(DD[0]))

# Compute satisfaction, not used for modeling
#nDA = [v for k,v in dist_not_A.items()]
#nDA = np.array(nDA)
#nDB = [v for k,v in dist_not_B.items()]
#nDB = np.array(nDB)
#nDD = np.vstack((nDA, nDB))
#nDmin = np.min(nDD,axis=0)
#print('all distances', np.shape(nDD))
#satif_not = np.sum([1 for d in nDmin if d<35.])/np.shape(nDD[0])

out.write(f'XLs satisfaction, used for modeling: {satif[0]}\n')
#out.write(f'XLs satisfaction, not used for modeling: {satif_not[0]}\n')
out.close()

exit()



