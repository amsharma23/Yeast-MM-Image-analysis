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
path_dir = "/Users/amansharma/Documents/Data/test/"; #where the ND2 files are
output_path = "/Users/amansharma/Documents/Data/test/";#where you want the ouput

ND2_files = [f for f in listdir(path_dir) if (isfile(join(path_dir, f)) and (".nd2" in f))];


for file in ND2_files:
    if (file[0] != '.'):
    
        fil = str.replace(file,'.nd2','');
    
    
        if (not os.path.exists(output_path+fil)):
            os.mkdir(output_path+fil);#make a folder where all files will be saved
        
        if (not os.path.exists(output_path+fil+"/FOVs/")):
            os.mkdir(output_path+fil+"/FOVs/"); #make a folder where all the FOV files will be saved
        

        with ND2Reader_SDK(path_dir+"/"+file) as frames:
            fovs = (frames.sizes);
            print(frames.sizes);
            print(type(frames));
            Fovs = range(fovs['m']);
# --> xy: frame xy pixels to bundle; c: channels; t: time; z: z-slices; v: FOVs            
            for fov in Fovs: #goes over all FoVs
                frames.default_coords['m'] = fov;  #picking the FOV
                if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)):
                     os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov));#makes folder with all the time frames and channels for a particular fov
                
                Cs = range(fovs['c']);
                for c in Cs: #goes over all channels
                    frames.default_coords['c'] = c;
                    if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/Chan"+str(c)):
                        os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/Chan"+str(c));#makes folder with all the time frames for a particular channel    

                    Ts =range(fovs['t']);
                    for t in Ts:
                        frames.default_coords['t'] = t;
                        if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/Chan"+str(c)+"/Time"+str(t)):
                            os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/Chan"+str(c)+"/Time"+str(t));#makes folder with all the time frames for a particular channel
                        
                        frames.iter_axes = ['z'];  # iterates through times
                        frames.bundle_axes = 'yx';  # stitching the whole image in x and y dimension

                    
            
    
                        z=0;
                        try:
                            for frame in frames[:]: #iterates through time
                                imwrite(output_path+fil+"/FOVs/FOV"+str(fov)+"/Chan"+str(c)+"/Time"+str(t)+"/Z_pos"+str(z)+".tif",frame,photometric='minisblack');
                                z=z+1;
        
                        except KeyError:
                            print("Will be missing last frame; due to ND2Reader error");