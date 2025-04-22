#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:33:24 2025

@author: Rita O
"""
import os
import sys
import matplotlib.pyplot as plt

plt.close('all');
os.system('clear')
os.system('cls')

import sys

import pandas as pd
import re
import matplotlib.pyplot as plt
import distinctipy
from distinctipy import colorblind
from distinctipy import examples
import numpy as np

plt.close('all')

## Create diffusion vecs - use matlab instead
bvecs_file      = '/home/localadmin/Documents/Rita/Data/common/bvecs_multishell_20240124.txt'
bvecs_file      = '/home/localadmin/Documents/Rita/Data/common/bvecs_multishell_20250411.txt'
bvals           = [1000, 2000, 3500, 5000, 7000]
nshells         = len(bvals)
ndir_per_shell  = [12, 24, 45, 64, 64]

color_list = distinctipy.get_colors(len(bvals),pastel_factor=0.5)


## Load bvecs
with open(bvecs_file, 'r') as file:
    data = file.read()
pattern = r"Vector\[\d+\] = \( ([^,]+), ([^,]+), ([^\)]+) \)"
matches = re.findall(pattern, data)
columns = ['X', 'Y', 'Z']
bvecs = pd.DataFrame(matches, columns=columns, dtype=float).to_numpy().T

expanded_bvals = [bvals[i] for i in range(len(bvals)) for _ in range(ndir_per_shell[i])]
expanded_bvals = np.multiply(expanded_bvals,1/max(bvals))
bvecs_plot = np.multiply(bvecs,np.sqrt(expanded_bvals))
 
## Plot
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111, projection='3d')
for b in range(len(bvals)):
    for v in range(ndir_per_shell[b]):
        
        if b==0:
            start_v = v
        else:
            start_v = v + sum(ndir_per_shell[0:b])
        
        ax.scatter(bvecs_plot[0, start_v], bvecs_plot[1, start_v], bvecs_plot[2,start_v], color=color_list[b], marker='*')
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])
ax.set_box_aspect([1, 1, 1])  


## Transform into bruker format
bruker_file = bvecs_file.replace('.txt','_bruker.txt')

f = open(bruker_file, "w")
f.write("[shells=" + str(nshells) + "]\n\n")
f.close()


fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111, projection='3d')
for b in range(len(bvals)):
    f = open(bruker_file, "a")
    f.write("[bvalue = " + str(bvals[ b])+ ".000000]\n")
    f.write("[directions=" + str(ndir_per_shell[b])+ "]\n")
    f.write("CoordinateSystem = xyz\n")
    f.write("Normalisation = unity\n")
    for v in range(ndir_per_shell[b]):
        
        if b==0:
            start_v = v
        else:
            start_v = v + sum(ndir_per_shell[0:b])
            
        pattern = r"Vector\[" + str(start_v) + r"\] = \(([^,]+), ([^,]+), ([^\)]+)\)"
        match = re.search(pattern, data)

        # If a match is found, write it to the file
        if match:
            x, y, z = match.groups()
            f.write(f"Vector[{v}] = ( {x}, {y}, {z} )\n")
            ax.scatter(float(x), float(y), float(z), color=color_list[b], marker='*')

f.close()
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])
ax.set_box_aspect([1, 1, 1])  


## Randomize order and save in bruker format
bvals_random        = [1000, 7000, 2000, 5000, 3500]
bvals_order         = (0, 4, 1, 3, 2)

bruker_file_random = bruker_file.replace('.txt','_random.txt')

with open(bruker_file, 'r') as file:
    data = file.read()
pattern = r"\[bvalue = [^\]]+\](?:\n(?!\[bvalue = ).*)*"

matches = re.findall(pattern, data, re.MULTILINE)
reordered_matches = [matches[i] for i in bvals_order if i < len(matches)]

f = open(bruker_file_random, "w")
f.write("[shells=" + str(nshells) + "]\n\n")
f.close()

with open(bruker_file_random, 'a') as file:
    for match in reordered_matches:
        file.write(match.strip() + "\n")  
