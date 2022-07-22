#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 12:25:17 2022

@author: amansharma
EXTRACTS CROPPED FOVs TIME SERIES INTO SEPERATE FOLDERS
Make sure you have all_funcs added to path
"""
#%% libraries
from tifffile import imwrite
from all_funcs import crop_using_roi_tuple
from all_funcs import extract_roi_coordinates
from nd2reader import ND2Reader
from pims import ND2Reader_SDK
import os
import numpy as np
import matplotlib.pyplot as plt


from os import walk
from os import listdir
from os.path import isfile, join


from read_roi import read_roi_file
from read_roi import read_roi_zip

import shutil

#%%
path_dir = "/Volumes/GodardSSD/Microscope Images/Saransh"; #where the ND2 files are
output_path = "/Users/amansharma/Documents/Data/test/";#where you want the ouput



ND2_files = [f for f in listdir(path_dir) if (isfile(join(path_dir, f)) and (".nd2" in f))];

for file in ND2_files:
    if (file[0] == '.'): file = file[2:];
    
    fil = str.replace(file,'.nd2','');
    
    
    if (not os.path.exists(output_path+fil)):
        os.mkdir(output_path+fil);#make a folder where all files will be saved
        
    if (not os.path.exists(output_path+fil+"/FOVs/")):
        os.mkdir(output_path+fil+"/FOVs/"); #make a folder where all the FOV files will be saved
        

    with ND2Reader_SDK(path_dir+"/"+file) as frames:
        fovs = (frames.sizes);
        fvs = range(fovs['m']);
        
        for fov in fvs: #goes over all FoVs
            #print(frames.sizes);# --> xy: frame xy pixels to bundle; c: channels; t: time; z: z-slices; v: FOVs
            frames.iter_axes = 't';  # iterates through times
            frames.bundle_axes = 'yx';  # stitching the whole image in x and y dimension
            frames.default_coords['m'] = fov;  #picking the FOV
            frames.default_coords['z'] = 7; #outof 7 choose the middle z-slice
            #frames.default_coords['c'] = 0;    
            t=1;
     
            if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)):
                os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov));#makes folder with all the time frames for a particular fov
            if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_ts/"):
                os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_ts/");#makes folder with all the time frames for a particular fov                if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_seg_im/"):
                os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_seg_im/");#makes folder for the segmented images   
            if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_zip/"):
                os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_zip/");#makes folder for the segemented RoIs
                

            try:
                for frame in frames[:]: #iterates through time
                    imwrite(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_ts/"+str(t)+".tif",frame,photometric='minisblack');
                    t+=1;
            except KeyError:
                print("Will be missing last frame; due to ND2Reader error");