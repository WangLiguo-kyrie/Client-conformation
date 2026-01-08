#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os

    
def contactsMatrix_within_cutoff(group_a, group_b,box_dim, radius=10):
    '''
    BB distance cutoff 1nm instead of 0.6nm distance for whole residue 
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions,box=box_dim)
    # determine which distances <= radius
    contacts_matrix = contacts.contact_matrix(dist, radius)*1
    return contacts_matrix     

#trajectory pbc and protein extraction
trajectory_dir='../'

u = mda.Universe('production-protein.tpr',trajectory_dir+'production-pbc-center.xtc')
# virtual stie in Martini couldn't be recognized as part of protein chain,so must use BB to avoid the fragments of virtual site bead. 
N49 = u.select_atoms('name BB and index 13356-13429')
print(len(N49.atoms))
print(N49.atoms)
Condensate_BB=u.select_atoms('name BB and index 0-13355')
print(len(Condensate_BB.atoms))
fragments = Condensate_BB.atoms.fragments
print(len(fragments))

n_res_client=36
n_res_condensate=163   #FUSLCD 163res
copy=len(fragments)

#discard initial 1us 
start = 2000
end = -1
#every 1ns
step = 2

inter_matrix_collect=np.zeros((n_res_client, n_res_condensate))
timeseries = []
#loop over each fragment
for i, frag in enumerate(fragments): 
    
    print('***** Condensate scaffold Monomer {} *****'.format(i))      
    monomer_A=frag.select_atoms('name BB')
    print(monomer_A.atoms)
    inter_matrix_monomer=np.zeros((n_res_client, n_res_condensate))   
        
    #directly construct the contact matrix between monomers, and then add all matrix together
    for index, ts in enumerate(u.trajectory[start:end:step]):    
        inter_matrix=contactsMatrix_within_cutoff(N49, monomer_A,ts.dimensions)
        inter_matrix_monomer=inter_matrix_monomer+inter_matrix
        timeseries.append(ts.time)
        print(' Time {}ps'.format(ts.time))
            
    inter_matrix_collect=inter_matrix_collect+inter_matrix_monomer
frames=len(timeseries)  
print(frames)
contacts_2D=inter_matrix_collect*copy/frames  
contacts_strength=contacts_2D.sum()   
print( contacts_strength)
print(contacts_2D.shape)
np.savetxt('2D_interchain_contact.xvg',contacts_2D,delimiter='	')

    


