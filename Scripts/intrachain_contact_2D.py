#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os
    
def contactsMatrix_within_cutoff(group_a, group_b, radius=10):
    '''
    BB distance cutoff 1nm instead of 0.6nm distance for whole residue 
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions)
    # determine which distances <= radius
    contacts_matrix = contacts.contact_matrix(dist, radius)*1
    return contacts_matrix    

#trajectory pbc and protein extraction
trajectory_dir='../'

u = mda.Universe(trajectory_dir+'production-pbc-center.gro',trajectory_dir+'production-pbc-center.xtc')
# virtual stie in Martini couldn't be recognized as part of protein chain,so must use BB to avoid the fragments of virtual site bead. 
N49 = u.select_atoms('name BB and index 13356-13429')
print(len(N49.atoms))


n_res=36

#discard initial 1us 
start = 2000
end = -1
#every 1ns
step = 2

    
intra_matrix_collect=np.zeros((n_res, n_res))
timeseries = []
BB=N49.select_atoms('name BB')
print(BB.atoms)    

for index, ts in enumerate(u.trajectory[start:end:step]):    
    intra_matrix=contactsMatrix_within_cutoff(BB, BB)
    intra_matrix_collect=intra_matrix_collect+intra_matrix
    timeseries.append(ts.time)
    print(' Time {}ps'.format(ts.time))

frames=len(timeseries)  
print(frames)
contacts_2D=intra_matrix_collect/frames  
contacts_strength=contacts_2D.sum()   
print( contacts_strength)
print(contacts_2D.shape)
np.savetxt('2D_intrachain_contact.xvg',contacts_2D,delimiter='	')
  

    


