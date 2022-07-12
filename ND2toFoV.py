#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 12:25:17 2022

@author: amansharma
EXTRACTS CROPPED FOVs TIME SERIES INTO SEPERATE FOLDERS
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

import shutil

#%%

def extract_roi_coordinates(roizipfile):
    coordinates_list=[] #[(left,top,right,bottom)] row and column coordinates for the rectangular roi
    roi = read_roi_zip(roizipfile); #read the roi contaning zip file
    print(roi)
    a=list(roi.items()); #extract the items (individual rois) in the roi file
    l=len(a); #find the number of ROIs in the zip file
    for i in range(0,l):
        top=a[i][1]['top'];
        left=a[i][1]['left'];
        bottom=a[i][1]['top']+a[i][1]['height'];
        right=a[i][1]['left']+a[i][1]['width'];
        coordinates_list.append([(left,top,right,bottom)]);
        
    #print(np.shape(coordinates_list))
    return coordinates_list



#%%
path_dir = "/Volumes/GodardSSD/Microscope Images/Saransh"
output_path = "/Users/amansharma/Documents/Data/test/rois_im"
ND2_files = [f for f in listdir(path_dir) if isfile(join(path_dir, f))];
ND2_files = ND2_files[1:]; #first element is .DS_Store file
location = os.getcwd();

for file in ND2_files:
    
    #fil = path_dir+"/"+file;
    fil = str.replace(file,'.nd2','');
    
    if not os.path.exists(output_path+"/FOVs-"+fil):
        os.mkdir(output_path+"/FOVs-"+fil); #make a folder where all the FOV files will be saved
    
    with ND2Reader(path_dir+"/"+file) as frames:
        for fov in [2]: #goes over all FoVs
            #print(frames.sizes); --> xy: frame xy pixels to bundle; c: channels; t: time; z: z-slices; v: FOVs
            frames.iter_axes = 't';  # iterates through times
            frames.bundle_axes = 'yx';  # stitching the whole image in x and y dimension
            frames.default_coords['v'] = fov;  #picking the FOV
            frames.default_coords['z'] = 7; #outof 7 choose the middle z-slice
            #frames.default_coords['c'] = 0;    
            t =1;
            #os.chdir(path_dir+"/FOVs"+file);
            if not os.path.exists(output_path+"/FOVs-"+fil+"/FOV"+str(fov)):
                os.mkdir(output_path+"/FOVs-"+fil+"/FOV"+str(fov));#makes folder with all the time frames for a particular fov
            
            for frame in frames[:]:#iterates through time
            
                imwrite(output_path+"/FOVs-"+fil+"/FOV"+str(fov)+"/"+str(t)+".tif",frame,photometric='minisblack');
                t+=1;
                