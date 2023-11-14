#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  21 17:53:25 2022

@author: amansharma
REMEMBER TO HAVE YeaZ ACTIVATED IN CONDA ENV FOR USING read_roi
All relevant scripts in Documents/Data/test/Aman_scripts
SEGMENTED CELLS IN ROI ZIP MUST 
"""


#%% libraries
from tifffile import imwrite
from tifffile import imread
from PIL import Image
from all_funcs import crop_using_roi_tuple
from all_funcs import extract_roi_coordinates

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from os import walk
from os import listdir
from os.path import isfile, join



#from read_roi import read_roi_file
#from read_roi import read_roi_zip

import trap_ex_func as trapf
import cell_size_func as cellf

import shutil
#%%

exp_folder = "/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/Best_examples/1/fluor";
path = "/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/Best_examples/1/";
seg = path+"segm/pics/";
bg = imread("/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/Best_examples/bg.tif");
bg = np.array(bg);
bg_avg = np.average(bg,axis=0);
print(np.shape(bg_avg));
org_v = np.zeros((18,1));
exp_folderl = listdir(exp_folder);
exp_folderl.sort();
vlr=[[11.64,50.35],[14.73,56.87],[18.66,68.84],[29.55,80.31],[42.9,79.9],[16.2,63.9],[28.56,72.74],[40.71,82.42],[12.24,66.8],[15.58,66.07],[26.73,73.89],[0,67.18],[21.38,88.27],[29.58,91.86],[35.45,97.18],[39.1,98.39],[0,77.09],[14.35,96.34]];
vl = [vlr[i][0]+vlr[i][1] for i in range(18)];
max_prj =[];

st=0;

seg_img =[];
seg_img_pths = listdir(seg);
seg_img_pths.sort();
for imgg in seg_img_pths:
	if(imgg!='.DS_Store'): seg_img.append(imread(join(seg,imgg)));


for fov_flds in (exp_folderl):
	if(fov_flds!=".DS_Store" and not(isfile(join(exp_folder,fov_flds)))):
		img_pth = join(exp_folder,fov_flds);#FOV folders
		i = 0;
		ls = listdir(img_pth);
		ls.sort();
		for imgs in ls:
			if(imgs!='.DS_Store'):
				
				img = imread(join(img_pth,imgs));
				if(st==0): max_prj.append(img);
				else: max_prj[i]+=img;

				img = np.array(img) - bg_avg[2:,:];
				img[img<0] =0;
				si = seg_img[i]
				si[si==2] =0;
				si[si==3] =0;
				nz = np.nonzero(si);
				si[nz]=1; 
				img = np.multiply(img,si);
				org_v[i] += np.sum(img);
				


				i+=1;
		st+=1;

    	
for i in range(len(max_prj)):
	imwrite(path+"/max_prj/"+str(i+1)+".tif",max_prj[i] );

frac = [(org_v[i]/10**7)/(vl[i]) for i in range(18)];
frac = [frac[i]/frac[0] for i in range(18)];

fig,ax = plt.subplots();

ax.plot(np.array(range(len(org_v)))*15,frac,label='volume fraction fold change');
ax.scatter(np.array(range(len(org_v)))*15,frac,color='red');
ax.set_ylabel('AU');
ax.set_xlabel('Time(mins)');
ax.legend();
plt.savefig(path+"/frac.svg",format='svg');

fig.clf();
plt.close();

org_v = [org_v[i]/(org_v[0]*10**7) for i in range(18)];

fig,ax = plt.subplots();

ax.plot(np.array(range(len(org_v)))*15,org_v,label='organelle fold change');
ax.scatter(np.array(range(len(org_v)))*15,org_v
	,color='red');
ax.set_ylabel('AU');
ax.set_xlabel('Time(mins)');
ax.legend();
plt.savefig(path+"/orgn.svg",format='svg');
