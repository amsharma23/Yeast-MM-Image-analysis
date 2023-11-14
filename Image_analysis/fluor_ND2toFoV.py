#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 12:25:17 2022

@author: amansharma
EXTRACTS FOVs TIME SERIES and ΤΙΜΕ SERIES FROM FLUROSENCE CHANNEL INTO SEPERATE FOLDERS
Make sure you have all_funcs added to path & activate YeaZ for read-roi
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
import shutil

from os import walk
from os import listdir
from os.path import isfile, join
import tqdm

#from read_roi import read_roi_file
#from read_roi import read_roi_zip

import shutil

#%%
path_dir = "/Volumes/Saransh_Toshiba_SSD/MM_aging"; #where the ND2 files are
output_path = "/Users/amansharma/Documents/Data/test";#where you want the ouput

#ND2_files = [f for f in listdir(path_dir) if (isfile(join(path_dir, f)) and (".nd2" in f))];

files = ["2023.02.05_ylb128_aging_glurich_swain_mito_delayedcapture.nd2"];

for file in files:
    if (file[0] != '.'):
    
        fil = str.replace(file,'.nd2','');
    
        if (not os.path.exists(output_path+"/"+fil)):
            os.mkdir(output_path+"/"+fil);#make a folder where all files will be saved
        
#        if (not os.path.exists(output_path+"/"+fil+"/FOVs_fl/")):
#            os.mkdir(output_path+"/"+fil+"/FOVs_fl/"); #make a folder where all the FOV files will be saved
        

        with ND2Reader_SDK(path_dir+"/"+file) as frames:
            fovs = (frames.sizes);
            fvs = range(fovs['m']);
        
            print((fvs));
            c = 1;
            for fov in [10]: #goes over all FoVs
                if(c==0):

                    frames.iter_axes = 't';  # iterates through times
                    frames.bundle_axes = 'yx';  # stitching the whole image in x and y dimension
                    frames.default_coords['m'] = fov;  #picking the FOV
                    frames.default_coords['z'] = 7; #outof 7 choose the middle z-slice
                    frames.default_coords['c'] = c; 
                    t=1;
         
                    if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)):
                        os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov));#makes folder with all the time frames for a particular fov
                
                    if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_ts_ph/"):
                        os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_ts_ph/");#makes folder with all the time frames for a particular fov                
                

                    if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_seg_im/"):
                        os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_seg_im/");#makes folder for the segmented images   
                    if not os.path.exists(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_zip/"):
                        os.mkdir(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_zip/");#makes folder for the segemented RoIs
                    

                    try:
                        for frame in frames[:]: #iterates through time
                            imwrite(output_path+fil+"/FOVs/FOV"+str(fov)+"/FOV"+str(fov)+"_ts_ph/"+str(t)+".tif",frame,photometric='minisblack');
                            t+=1;
                    except KeyError:
                        print("Will be missing last frame; due to ND2Reader error");
                
                elif(c==1):
                    frames.iter_axes = 't';  # iterates through times
                    frames.bundle_axes = 'yx';  # stitching the whole image in x and y dimension
                    frames.default_coords['m'] = fov;  #picking the FOV
                    frames.default_coords['c'] = c;
                    ts = [24*k for k in range(11)]; 
                    for t in ts:

                        if not os.path.exists(output_path+"/"+fil+"/FOV"+str(fov)+"_ts_fl/"):
                            os.mkdir(output_path+"/"+fil+"/FOV"+str(fov)+"_ts_fl/");#makes folder with all the time frames for a particular fov                

                        if not os.path.exists(output_path+"/"+fil+"/FOV"+str(fov)+"_ts_fl/"+str(t)+"/"):
                            os.mkdir(output_path+"/"+fil+"/FOV"+str(fov)+"_ts_fl/"+str(t)+"/");#makes folder with all the time frames for a particular fov                

                        try:
                            for i in range(13): #iterates through time
                                frames.default_coords['z'] = i;
                                frame = frames[t];
                                imwrite(output_path+"/"+fil+"/FOV"+str(fov)+"_ts_fl/"+str(t)+"/"+str(i)+".tif",frame,photometric='minisblack');
                                
                        except KeyError:
                            print("Will be missing last frame; due to ND2Reader error");

                #c=0;
                print("Fov:"+str(fov)+" of "+str(len(fvs)));