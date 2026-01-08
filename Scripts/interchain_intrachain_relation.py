#! /usr/bin/python
# -*- coding: utf-8 -*-
# analyze interchain contact VS intrachain contact change
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
from scipy.stats import linregress

fig, ax = plt.subplots(figsize=(6,5))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
tick_prop = font_manager.FontProperties(fname=font_path, size=15)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)    
cm = plt.cm.get_cmap('Set2')

data_inter2 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/FUSLCD/Client/N49/interchain_contact/interchain_contact_1D.xvg')
contact_inter2 = data_inter2[:,1]
contact_inter2_stremgth=contact_inter2.sum()
print(contact_inter2_stremgth)
data_inter3 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/WGR5/Client/N49/interchain_contact/interchain_contact_1D.xvg')
contact_inter3 = data_inter3[:,1]
contact_inter3_stremgth=contact_inter3.sum()
print(contact_inter3_stremgth)
data_inter4 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/mfp3S/Client/N49/interchain_contact/interchain_contact_1D.xvg')
contact_inter4 = data_inter4[:,1]
contact_inter4_stremgth=contact_inter4.sum()
print(contact_inter4_stremgth)

data_intra = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/monomer/Martini3-IDP2_Bonded-BBP5-Dih7-noscfix-SCdih15/N49/intrachain_contact/intrachain_contact_1D.xvg')
contact_intra = data_intra[:,1]
contact_intra_stremgth=contact_intra.sum()
print(contact_intra_stremgth)    
data_intra2 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/FUSLCD/Client/N49/intrachain_contact/intrachain_contact_1D.xvg')
contact_intra2 = data_intra2[:,1]
contact_intra2_stremgth=contact_intra2.sum()
print(contact_intra2_stremgth)
data_intra3 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/WGR5/Client/N49/intrachain_contact/intrachain_contact_1D.xvg')
contact_intra3 = data_intra3[:,1]
contact_intra3_stremgth=contact_intra3.sum()
print(contact_intra3_stremgth)
data_intra4 = np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project3_Client/Condensate/mfp3S/Client/N49/intrachain_contact/intrachain_contact_1D.xvg')
contact_intra4 = data_intra4[:,1]
contact_intra4_stremgth=contact_intra4.sum()
print(contact_intra4_stremgth)

corr_pearson2, _ = pearsonr(contact_intra2/contact_intra2_stremgth-contact_intra/contact_intra_stremgth,contact_inter2/contact_inter2_stremgth)
print(corr_pearson2)
corr_pearson3, _ = pearsonr(contact_intra3/contact_intra3_stremgth-contact_intra/contact_intra_stremgth,contact_inter3/contact_inter3_stremgth)
print(corr_pearson3)
corr_pearson4, _ = pearsonr(contact_intra4/contact_intra4_stremgth-contact_intra/contact_intra_stremgth,contact_inter4/contact_inter4_stremgth)
print(corr_pearson4)
    
#ax.plot(contact_inter2/contact_inter2_stremgth, contact_intra2/contact_intra2_stremgth-contact_intra/contact_intra_stremgth, 'o', color=cm.colors[0])    

slope, intercept, r_value, p_value, std_err = linregress(contact_inter2/contact_inter2_stremgth, contact_intra2/contact_intra2_stremgth-contact_intra/contact_intra_stremgth)
# Create the fitted line
print(p_value)
#fitted_line = slope *(contact_inter2/contact_inter2_stremgth) + intercept
#plt.plot(contact_inter2/contact_inter2_stremgth, fitted_line, color=cm.colors[0])

df = pd.DataFrame({'x_column': contact_inter2/contact_inter2_stremgth, 'y_column': contact_intra2/contact_intra2_stremgth-contact_intra/contact_intra_stremgth})
sns.regplot(data=df, x="x_column", y="y_column",marker='o', color=cm.colors[0])

ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel(r'$\Delta$ Intrachain contact',fontproperties=font_prop)
ax.set_xlabel('Residue interchain contact',fontproperties=font_prop)
ax.set_title('N49-FUSLCD',fontproperties=font_prop)
ax.text(.05, .12, r'$corr$={:.2f}'.format(r_value),horizontalalignment='left',fontproperties=font_prop,transform=ax.transAxes)
ax.text(.05, .05, r'$p$={:.2g}'.format(p_value),horizontalalignment='left',fontproperties=font_prop,transform=ax.transAxes)
#plt.xlim(0.2,0.45)
plt.xticks([0.015,0.020,0.025,0.030,0.035,0.040],fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
plt.savefig('contact_coupling_1D.png',dpi=600,bbox_inches='tight')
plt.show()  

