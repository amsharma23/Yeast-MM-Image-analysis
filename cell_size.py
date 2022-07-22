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
from tifffile import imread
from PIL import Image
from all_funcs import crop_using_roi_tuple
from all_funcs import extract_roi_coordinates
from nd2reader import ND2Reader
import os
import pandas as pd
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
#%%
def get_seg_images(path):#outputs array of tif files read as numpy arrays

    seg_images={};

    fovs_folder = [join(path,f) for f in listdir(path) if not(isfile(join(path, f)))];

    
    for folders in fovs_folder:
        if(("seg" in folders) and ("im" in folders)):
            for imgs in listdir(folders):
                    if (isfile(join(folders,imgs)) and (".tif" in imgs)):
                        im_ar = imread(join(folders,imgs));
                        seg_images[imgs] = im_ar;
                        
                    
    return(seg_images);


#%% Get fov.zip file names and roi coordinates 
#Give directory of where all the zip files are - zip file of segmentation from YeastMate Fiji Plugin

exp_folder = "/Users/amansharma/Documents/Data/test/temp_fovs";
trajcs_a ={};
ct = 1;
for fov_flds in listdir(exp_folder):
    if(fov_flds!=".DS_Store" and not(".tif" in fov_flds)):
        pth = exp_folder+"/"+fov_flds;#FOV folders
        
        seg_imgs = get_seg_images(pth); #store all segmented images       
        
        roi_d = read_roi_zip(pth+"/"+fov_flds+"_zip/"+fov_flds+".zip"); #get the segmentation data

        trapf.remtraps(roi_d,seg_imgs,pth); #makes and saves cropped and trap removed tif images in a Crops/ folder in the path given

        trajcs = cellf.el_vol(pth+"/"); #the function will plot the estimated vol vs time as fitted by prolate ellipsoids, and spits back the trajectories
        
        trajcs_a[ct] = [trajcs[key] for key in trajcs.keys()];
        ct+=1;
        print(ct);

cellf.plt_trj(trajcs_a,exp_folder); #plots the volume trajectories of every trap in a FOV
print("Plotted trajectories");
cellf.plt_dist(trajcs_a,exp_folder); #plots a distribution for the volumes of the cells





    

    