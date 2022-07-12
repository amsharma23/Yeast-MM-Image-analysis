#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:10:08 2022

@author: amansharma
REMEMBER TO HAVE YeaZ ACTIVATED IN CONDA ENV FOR USING read_roi
All relevant scripts in Documents/Data/test/Aman_scripts
SEGMENTED CELLS IN ROI ZIP MUST 
"""


#%% libraries
from tifffile import imwrite
from all_funcs import crop_using_roi_tuple
from all_funcs import extract_roi_coordinates
from nd2reader import ND2Reader
import os
import numpy as np
import matplotlib.pyplot as plt

from os import walk
from os import listdir
from os.path import isfile, join


from read_roi import read_roi_file
from read_roi import read_roi_zip

import trap_ex_func as trapf
import cell_size_func as cellf

import shutil


#%% Get fov.zip file names and roi coordinates

dirl = "/Users/amansharma/Documents/Data/test/";

zip_dir = "/Users/amansharma/Documents/Data/test/rois_zip/";
zip_files = [f for f in listdir(zip_dir) if isfile(join(zip_dir, f))];
seg_zip_files = [f for f in listdir(zip_dir) if (isfile(join(zip_dir, f))  and ("seg" in f))];


zip_files = zip_files[1:];
zfs = zip_files.sort();


files =[];
for s in zip_files:
    files.append(s.replace('.zip', '')); 
    

if(not os.path.exists(dirl+"plots/")):
    os.mkdir(dirl+"plots");


#%%
roi_d = read_roi_zip("/Users/amansharma/Documents/Data/test/FOVs-20220531_ylb128_gal2p/FOV2_zip/FOV2.zip");

roi_traps_arr = trapf.gettraparr(roi_d); #gives array of the name of segmented cells that are traps

roi_traps={};

for el in range(len(roi_traps_arr)):
    for el1 in (roi_traps_arr[el].keys()):
       
        roi_traps[el1] = roi_d[el1];

roi_seg_in_mul_f = trapf.getmtimeframes(roi_d);

roi_seg_cells = cellf.getycells(roi_seg_in_mul_f,roi_traps);


cells = cellf.getuniquecells(roi_seg_cells);

cells = cellf.addarea(cells); #add an value of cell area to each cell - key is 'area'

#keys of cells are the name for a unique cell; it has array with the segmentation properties stored for each time frame
#so cells[key]['time_f'] will give you number of time frames the cell is in; cells[key][int] will give segmentation element; cells[key][int]['area'] will give the calculated area for that time frame


cellf.plot_AvT(cells,dirl); #plots area vs time for cells with more than 4 time point measures
cellf.savekym(cells,dirl); #save kymograph for cells 

    

    