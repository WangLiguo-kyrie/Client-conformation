#! /usr/bin/python

#calculate the client COG Z-position evolution
import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis import *
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

trajectory_dir='./'
u = mda.Universe(trajectory_dir+'production-pbc-center.gro',trajectory_dir+'production-pbc-center.xtc')
# virtual stie in Martini couldn't be recognized as part of protein chain,so must use BB to avoid the fragments of virtual site bead. 
client = u.select_atoms('name BB and index 13356-13429')
print(len(client.atoms))
condensate=u.select_atoms('name BB and index 0-13355')

#discard initial 1us 
start = 2000
end = -1
#every 100ns
step = 200

dist_result=[]
for index, ts in enumerate(u.trajectory[start:end:step]):  
    #use the whole COG to anchor the position, consistent to the setting in MDP condensate (both folded domains and IDR)
        
    dist_z = client.center_of_geometry()[2]-condensate.center_of_geometry()[2]
    print(dist_z)
    dist_result.append([ts.time,dist_z])
dist_result=np.array(dist_result)
print(dist_result.shape)
np.savetxt('client_distr.xvg',dist_result,delimiter='	')


    
    



