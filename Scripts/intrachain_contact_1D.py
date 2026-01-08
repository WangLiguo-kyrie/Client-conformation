#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os
    
def contacts_within_cutoff(group_a, group_b, radius=10):
    '''
    BB distance cutoff 1nm instead of 0.6nm distance for whole residue 
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions)
    # determine which distances <= radius
    n_contacts = contacts.contact_matrix(dist, radius).sum()
    return n_contacts    

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

     
results_contact=[]
for k in range(1,n_res+1):
    res_A=k	
    # neighboring i-3,i-2,i-1,i, i+1,i+2,i+3 residues
    neighbor=[x for x in range(res_A-3,res_A+4)]
    neighbor_filter=[x for x in neighbor if (37> x >0)]
        
    resA=N49.select_atoms('resid {} and name BB'.format(res_A))
    neighbor_list=''
    for j in neighbor_filter:
        neighbor_list=neighbor_list+' '+str(j)
    print(neighbor_list)
    resB=N49.select_atoms('not resid {} and name BB '.format(neighbor_list))
    print(resA.atoms)
    print(resB.atoms)
    timeseries = []

    for index, ts in enumerate(u.trajectory[start:end:step]):    
        ca=contacts_within_cutoff(resA, resB)
        timeseries.append([ts.time, ca])
        print('Res {}, Time {}ps, Contacts {}'.format(res_A,ts.time,ca)) 
    timeseries=np.array(timeseries)

    ca_contact=timeseries[:,1].sum()  
    frames=timeseries.shape[0] 
    print(frames)
    contacts_strength=ca_contact/frames
    results_contact.append([res_A,contacts_strength])    
results_contact=np.array(results_contact)    
np.savetxt('intrachain_contact_1D.xvg',results_contact,delimiter='	')
  

    


