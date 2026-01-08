#! /usr/bin/python
# -*- coding: utf-8 -*-
# calculate gyrate distribution in different trajectory interval to analyze convergence and conformation ensemble
import sys
import os
import numpy as np
import argparse
import pandas as pd
import seaborn as sns
#sns.set(color_codes=True)
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from statistics import mean
from scipy.stats import pearsonr,spearmanr

fig, ax = plt.subplots(figsize=(8,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
tick_prop = font_manager.FontProperties(fname=font_path, size=15)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)  
cm = plt.cm.get_cmap('Accent')
data = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/monomer/Martini3-IDP2_Bonded-BBP5-Dih7-noscfix-SCdih15/N49/intrachain_contact/intrachain_contact_1D.xvg')
contact = data[:,1]
contact_stremgth=contact.sum()
print(contact_stremgth)
resid = data[:,0]      
data2 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/FUSLCD/Client/N49/interchain_contact/interchain_contact_1D.xvg')
contact2 = data2[:,1]
contact2_stremgth=contact2.sum()
print(contact2_stremgth)
resid2 = data2[:,0]  
data3 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/WGR5/Client/N49/interchain_contact/interchain_contact_1D.xvg')
contact3 = data3[:,1]
contact3_stremgth=contact3.sum()
print(contact3_stremgth)
resid3 = data3[:,0]  
data4 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/mfp3S/Client/N49/interchain_contact/interchain_contact_1D.xvg')
contact4 = data4[:,1]
contact4_stremgth=contact4.sum()
print(contact4_stremgth)
resid4 = data4[:,0]    
    
corr_pearson2, _ = pearsonr(contact2,contact3)
corr_pearson3, _ = pearsonr(contact4,contact3)
corr_pearson4, _ = pearsonr(contact2,contact4)
print(corr_pearson2)
print(corr_pearson3)    
print(corr_pearson4)    
#ax.plot(resid,contact/contact_stremgth, color=cm.colors[0],label='Solution')  
ax.plot(resid2,contact2/contact2_stremgth, color=cm.colors[1],label='FUSLCD Condensate',linewidth=2)  
ax.plot(resid3,contact3/contact3_stremgth, color=cm.colors[2],label='WGR5 Condensate',linewidth=2) 
ax.plot(resid4,contact4/contact4_stremgth, color=cm.colors[4],label='mfp-3S Condensate',linewidth=2) 
#ax.text( 0.5,0.052,'FUSLCD vs WGR5: {:.2f}'.format(corr_pearson2),horizontalalignment='left',fontproperties=legend_prop,color=cm.colors[1])
#ax.text( 0.5,0.055,'mfp-3S vs WGR5: {:.2f}'.format(corr_pearson3),horizontalalignment='left',fontproperties=legend_prop, color=cm.colors[2])
#ax.text( 0.5,0.058,'FUSLCD vs mfp-3S: {:.2f}'.format(corr_pearson4),horizontalalignment='left',fontproperties=legend_prop, color=cm.colors[4])
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('Normalized Contact Number',fontproperties=font_prop)
ax.set_xlabel('N49 Residue',fontproperties=font_prop)
plt.ylim(0.015,0.04)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
plt.savefig('Normalized_interchain_contact_1D.png',dpi=600,bbox_inches='tight')
plt.show()  


fig, ax = plt.subplots(figsize=(8,6))
ax.plot(resid2,contact2, color=cm.colors[1],label='FUSLCD Condensate')  
ax.plot(resid3,contact3, color=cm.colors[2],label='WGR5 Condensate')  
ax.plot(resid4,contact4, color=cm.colors[4],label='mfp-3S Condensate') 
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('Contact Number',fontproperties=font_prop)
ax.set_xlabel('N49 Residue',fontproperties=font_prop)
#plt.xlim(0.2,0.45)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
plt.savefig('interchain_contact_1D.png',dpi=600,bbox_inches='tight')
plt.show() 

