#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:08:20 2022

@author: amansharma
"""
#%%
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



#%% Get fov.zip file names and roi coordinates
zip_dir = "/Users/amansharma/Documents/Data/test/rois_zip/";
zip_files = [f for f in listdir(zip_dir) if isfile(join(zip_dir, f))];

zip_files = zip_files[1:];
zfs = zip_files.sort();
files =[]
for s in zip_files:
    files.append(s.replace('.zip', '')); 


st = "/Users/amansharma/Documents/Data/test/rois_im";
shutil.rmtree(st);
os.mkdir(st);

fld = [st+'/'+s for s in files];
for s in fld: os.mkdir(s);



#%%

def crop_using_roi_tuple(coordinates_tuple,image_array):

    return image_array[coordinates_tuple[1]:coordinates_tuple[3], coordinates_tuple[0]:coordinates_tuple[2]]


def nd2_roi_extractor(file,output_path,roi_coordinates,fov,roi): 
    
    with ND2Reader(file) as frames:
        
        #print(frames.sizes); --> xy: frame xy pixels to bundle; c: channels; t: time; z: z-slices; v: FOVs
        frames.iter_axes = 't';  # iterates through times
        frames.bundle_axes = 'yx';  # stitching the whole image in x and y dimension
        frames.default_coords['v'] = fov-1;  #picking the FOV
        frames.default_coords['z'] = 7;
        frames.default_coords['c'] = 0;
        i=1;
        
        
        
        #cropped = crop_using_roi_tuple(roi_coordinates,frames[0]);
        
        os.mkdir(output_path+'/');
        
        for frame in frames[:]:
            cropped = crop_using_roi_tuple(roi_coordinates,frame);
            imwrite(output_path+'/'+str(i)+'.tif', cropped, photometric='minisblack');
            i+=1;
        
        
    frames.close()




def extract_roi_coordinates(roizipfile):
    coordinates_list=[] #[(left,top,right,bottom)] row and column coordinates for the rectangular roi
    roi = read_roi_zip(roizipfile); #read the roi contaning zip file
    a=list(roi.items()); #extract the items (individual rois) in the roi file
    l=len(a); #find the number of ROIs in the zip file
    for i in range(0,l):
        top=a[i][1]['top'];
        left=a[i][1]['left'];
        bottom=a[i][1]['top']+a[i][1]['height'];
        right=a[i][1]['left']+a[i][1]['width'];
        coordinates_list.append([(left,top,right,bottom)]);

    return coordinates_list




for fls in files:

    roi_tuples = extract_roi_coordinates(zip_dir+fls+'.zip');    
    fov = int(fls[3]);
    stf = fld[files.index(fls)];    
    for j in range(0,len(roi_tuples)):
        roi_coordinates = roi_tuples[j][0];
        st = stf+"/roi_"+str(j+1);
        nd2_roi_extractor('/Volumes/GodardSSD/Microscope Images/Saransh/20220531_ylb128_gal2p.nd2',st,roi_coordinates,fov,j);
            
        
        