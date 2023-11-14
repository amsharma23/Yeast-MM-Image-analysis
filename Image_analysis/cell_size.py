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

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from os import walk
from os import listdir
from os.path import isfile, join
from tqdm import tqdm


#from read_roi import read_roi_file
#from read_roi import read_roi_zip

import trap_ex_func as trapf
import cell_size_func as cellf

import shutil
#%%
def get_seg_images(path,traps):#outputs array of tif files read as numpy arrays

    seg_images={};
    if(not traps):
        fovs_folder = [join(path,f) for f in listdir(path) if not(isfile(join(path, f)))];

    
        for folders in fovs_folder:
            if(("seg" in folders) and ("im" in folders)):
                for imgs in listdir(folders):
                    if (isfile(join(folders,imgs)) and (".tif" in imgs)):
                        im_ar = imread(join(folders,imgs));
                        seg_images[imgs] = im_ar;
    else:
        #fovs_folder = [join(path,f) for f in listdir(path) if not(isfile(join(path, f)))];

    
        #for folders in fovs_folder:
        #    if("seg" in folders):
                #for trap_fld in listdir(folders):
        arr = listdir(path);
        #print(arr)
        arr = sorted(arr);      
        for img_name in arr:
            if (isfile(join(path,img_name)) and (".tif" in img_name)):
                
                try:
                    im_ar = imread(join(path,img_name));
                except:
                    print('');
                
                seg_images[img_name] = im_ar;



                    
    return(seg_images);

def get_bf_images(path):#outputs array of tif files read as numpy arrays

    bf_images={};
    arr = listdir(path);
    arr = sorted(arr);      
    for img_name in arr:
        if (isfile(join(path,img_name)) and (".tif" in img_name)):
            
            try:
                im_ar = imread(join(path,img_name));
            except:
                print('');
            
            bf_images[img_name] = im_ar;

    return(bf_images);


n = [2,3];

for i in (n):
    
    if(i==2):
        k = [0,1,2,3,4];

    elif(i==3):
        k = [0,1,2,3,4];

    for j in (k):

        masks_pth = "/Users/amansharma/Documents/Data/test/20220617_ylb128_glu_minimal/z_4/fov"+str(i)+"/set0/trap"+str(j)+"/masks/";
        bf_pth = "/Users/amansharma/Documents/Data/test/20220617_ylb128_glu_minimal/z_4/fov"+str(i)+"/set0/trap"+str(j)+"/BF/";
    
    
    
        pth = "/Users/amansharma/Documents/Data/test/20220617_ylb128_glu_minimal/z_4/fov"+str(i)+"/set0/trap"+str(j)+"/";
        
        traps = True;
        seg_imgs = get_seg_images(masks_pth,traps);
        bf_imgs = get_bf_images(bf_pth);
        
        time_scale = 1200; #1200sec for every time step
        trajcs_mth,trajcs_dau=cellf.ar_plot(seg_imgs,len(seg_imgs.keys()),time_scale,pth);
        
        cellf.make_kym(seg_imgs,bf_imgs,trajcs_mth,trajcs_dau,pth,time_scale);



            

