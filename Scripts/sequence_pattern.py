#! /usr/bin/python
# used for aromatic pattern parameter calculation in '2020-science-Valence and patterning of aromatic residues determine the phase behavior of prion-like domains'
import sys
import os
import numpy as np
import pandas as pd

sequence_fasta='GYGYDLGYNAPWPYNNGYYGYNGYNGYHGRYGWNKGWNNGPWGGY'
stickers=['F','W','Y'] #sticke definition, 'K/R' could also be included
#transform sticke to 1, others to 0
def seq_encode(fasta):
    code=[]
    for letter in fasta:
        if letter in stickers:
            code.append(1)
        else:code.append(0)
    return code

sequence_code=seq_encode(sequence_fasta)
print(sequence_code)

#total sequence aromatic asymmetry
def asymmetry(binary_code):
    fraction=sum(binary_code)/len(binary_code)
    asymm=(2*fraction-1)**2
    return asymm

global_asymmetry=asymmetry(sequence_code)
print(global_asymmetry)
    
#‘average’ deviation from the total aromatic asymmetry in the sequence interatively
blob_size=5

deviation=0
for blob in range(len(sequence_code)-blob_size+1):
    blob_code=sequence_code[blob:blob+blob_size]
    local_asymmetry=asymmetry(blob_code)
    deviation+=(local_asymmetry-global_asymmetry)**2
deviation_ave=deviation/(len(sequence_code)-blob_size+1)
print(deviation_ave)
