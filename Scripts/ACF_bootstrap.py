#! /usr/bin/python
# -*- coding: utf-8 -*-
# calculate the correlation time uncertainty by bootstrap
import numpy as np
import argparse
import os
from statistics import mean, stdev

def simple_bootstrap(data):
    """
    Simple bootstrap with replacement.
    - Randomly sample n indices (with replacement) from 0..n-1.
    - Sort indices to preserve time order.
    - Return the corresponding rows.
    """
    n = len(data)
    idx = np.random.choice(n, size=n, replace=True)
    idx.sort()
    return data[idx]
    
def extract_sum_values(filename):
    """
    extract correlation time form gmx analyze -ac output
    """
    sum_values = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Look for lines that start with "COR:"
            if line.startswith("COR:"):
                tokens = line.split()
                # Ensure we have at least 6 tokens (COR: + 5 numbers)
                if len(tokens) >= 6:
                    try:
                        # Check if the first token after "COR:" is a float
                        float(tokens[1])
                        # If successful, this is a data line; the 4th number is at index 4
                        sum_val = float(tokens[4])
                        sum_values.append(sum_val)
                    except ValueError:
                        # The line contains words, skip it
                        continue
    return sum_values

data=np.loadtxt('end2end_dist.xvg')
#100 bootstrap sampling
for i in range(100):
    resample=simple_bootstrap(data)
    np.savetxt('bootstrap{}.xvg'.format(i),resample,delimiter='       ')
    os.system('gmx analyze -f bootstrap{}.xvg  -ac auto_bootstrp{}.xvg -fitfn exp >> autocorr_bootstrp.txt'.format(i,i)) 
    os.system('rm bootstrap{}.xvg'.format(i))
    os.system('rm auto_bootstrp{}.xvg'.format(i))
corr_time = extract_sum_values('autocorr_bootstrp.txt')
print(corr_time)
print(mean(corr_time))
print(stdev(corr_time))


