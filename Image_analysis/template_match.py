#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:54:49 2022

@author: amansharma
"""

from os import walk
from os import listdir
from os.path import isfile, join

import math
import tifffile as tf
import skimage.io as skio
import matplotlib.pyplot as plt 
import numpy as np
from PIL import Image 
from scipy import signal as sig 
import os
import cv2 

def getmindist(per_pts,pt):
    md = 0;
    for ppts in per_pts:
        dis = math.sqrt((ppts[0]-pt[0])**2 + (ppts[1]-pt[1])**2);
        if (dis<md or md==0): md = dis; 

    return(md);

def trap_loc(templ,img,ths):
      
    op = cv2.matchTemplate(img,templ,cv2.TM_CCORR_NORMED);
    locs = np.where(op >= ths)
    pts_p = [(0,0)];
    ct = 0;
    pos=[];
    imgg =0;
    for pts in zip(*locs[::-1]):
    
        min_d = getmindist(pts_p,pts);
        if(min_d>150):    
            imgg = cv2.rectangle(img,pts,(pts[0]+h,pts[1]+w),255,2)
            ct+=1;
            pts_p.append(pts);
    
    return(pts_p[1:],imgg,ct)        




trap = cv2.imread('/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/fluor_temp.tif',0);
#trap1 = trap.astype('uint8')
w, h = trap.shape[::-1];  
print(w,h);
fovs_loc = "/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/FOVs";
t_trps = 0;

for FOV in os.listdir(fovs_loc):
    if(FOV!='.DS_Store'):
        op_pth = os.path.join(fovs_loc,FOV);
        fl_pth = os.path.join(fovs_loc,FOV,FOV+"_ts_fl");
        ph_pth = os.path.join(fovs_loc,FOV,FOV+"_ts_ph");
        seg_pth = os.path.join(fovs_loc,FOV,FOV+"_seg_im");
        zip_pth = os.path.join(fovs_loc,FOV,FOV+"_zip");

        int_img = cv2.imread(fl_pth+'/1.tif',0);
        #int_img1 = int_img.astype('uint8');
        


        if not os.path.exists(op_pth+"/Crops/"):
            os.mkdir(op_pth+"/Crops/");
            
        th = 0.7;
        trapl, auto_ph, count = trap_loc(trap,int_img,th);
        if(count!=0):
            trps =1;
            tf.imwrite(op_pth+"/Crops/matched.tif",auto_ph);

            for pos in trapl:
                if not os.path.exists(op_pth+"/Crops/"+str(trps)):
                    os.mkdir(op_pth+"/Crops/"+str(trps));
                if not os.path.exists(seg_pth+"/"+str(trps)):
                    os.mkdir(seg_pth+"/"+str(trps));                    
                if not os.path.exists(zip_pth+"/"+str(trps)):
                    os.mkdir(zip_pth+"/"+str(trps));                    
                
                for ts_ph in os.listdir(ph_pth):
                    ph = tf.imread(os.path.join(ph_pth,ts_ph));
                    w_m,h_m = ph.shape[::-1];
                    if(pos[0]-20 <0): pl = 0;
                    else: pl = pos[0]-20;

                    if(pos[1]-20 <0): pb = 0;
                    else: pb = pos[1]-20;

                    pr = pos[0]+w+20;
                    pt = pos[1]+h+20;            
                    cr_im = ph[pb:pt,pl:pr];

                    tf.imwrite(op_pth+"/Crops/"+str(trps)+"/"+ts_ph,np.array(cr_im));

                trps+=1;
                t_trps+=1;


print(t_trps);
