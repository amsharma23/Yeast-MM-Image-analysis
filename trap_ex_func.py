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
from all_funcs import crop_using_roi_tuple
from all_funcs import extract_roi_coordinates

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

from PIL import Image
from os import walk
from os import listdir
from os.path import isfile, join

from skimage.draw import ellipse
from skimage.measure import label, regionprops, regionprops_table
from skimage.transform import rotate

#from read_roi import read_roi_file
#from read_roi import read_roi_zip

import trap_ex_func as trapf
import cell_size_func as cellf

import shutil
#%%Functions for extracting traps

#Used
def remtraps(roi_od,seg_ims,path):
    [y_max,x_max] = np.shape(seg_ims[list(seg_ims.keys())[0]]);
    
    roi_seg_in_frames = getsinglecells(roi_od);#single cell segments deteced across in time frames
    roi_seg_in_frames = assignpositions(roi_seg_in_frames);#assigns each segmented single cell a position attribute

    
    roi_seg_in_mul_f = getmtimeframes(roi_seg_in_frames);#single cell in more than one frame
    roi_seg_in_all_f = getalltimeframes(roi_seg_in_mul_f); #single cell segments detected across all time
    roi_seg_in_ff = getcellsinff(roi_seg_in_all_f); #single cells in first frame
    [lftn,lftp] = getleftop(roi_seg_in_ff);# gets the top left most segmented cell and calculates the trap positions according to this
        
    crop_traps(lftp,roi_seg_in_ff,seg_ims,path,y_max,x_max);
    
     
           
def crop_traps(start_pos,roi_od,seg_imgs,path,ym,xm): #saves croped trap-removed images in a Crops/ folder in the path
    check =[0,0];
    nms = list(seg_imgs.keys());
    nms.sort();
    im = seg_imgs[nms[0]];
    chnls =0;
    x_strt = int(start_pos[0]);
    y_strt = int(start_pos[1]);
    tr_n=1;    
    while x_strt<xm: #start of channel should be under image max
        check[0] = x_strt;
        check[1] = y_strt;
        if((chnls%2) == 0): 
            while check[1]<ym: #trap position should under image max
                try:
                    val_yeast = im[check[1]+20,check[0]+35];#where we expect a yeast cell to be 
                    if(val_yeast!=0):
                        p_x = check[0];
                        p_y = check[1];
                        save_crop(seg_imgs,[p_x,p_y],nms,xm,ym,tr_n,chnls,path);
                        tr_n+=1;
                    
                    check[1] += 175;
                except IndexError:
                    check[1] +=175; 
                    continue;

            chnls+=1;
            x_strt += 230;
            
                
        else:
            while check[1]<ym:
                try:
                    val_yeast = im[check[1]+20,check[0]+35];#where we expect a yeast cell to be 
                    if(val_yeast!=0):
                        p_x = check[0];
                        p_y = check[1];
                        save_crop(seg_imgs,[p_x,p_y],nms,xm,ym,tr_n,chnls,path);
                        tr_n+=1;
                    
                    check[1] += 175;
                except IndexError:
                    check[1] +=175; 
                    continue;

            chnls+=1;
            x_strt += 720;
            
    if(tr_n==1): print("No filled traps were found for: "+path);


def save_crop(img,pos,nm,xmx,ymx,tn,chn,path):
    
    xl = pos[0]-55;
    if(xl<0): xl = 0;
    xr = pos[0]+140;
    if(xr>xmx): xr = xmx;
    yt = pos[1]-65;
    if(yt<0): yt = 0;
    yb = pos[1]+65;
    if(yb>ymx): yb = ymx;
    t=0;
    
    for imgs in nm:
        imm = img[imgs];
        val = imm[pos[1]-5:pos[1]+5,pos[0]-5:pos[0]+5];
        val = val[np.where(val!=0)];
        val_neigh = imm[pos[1]-10:pos[1]+10, pos[0]+75:pos[0]+95]; #neighbour trap is approx 85pixel length away
        val_neigh = val_neigh[np.where(val_neigh!=0)];
        if(len(val)!=0):
            imm[imm==val[0]] =0;
            #print(val[0]);
            
        if(len(val_neigh)!=0):
            imm[imm==val_neigh[0]]=0;
            #print(val_neigh[0]);
        
        
        
        #print(chn,tn,imgs);
        
        img_cr = imm[yt:yb,xl:xr];
        if(not os.path.exists(path+"/Crops")):
            os.mkdir(path+"/Crops");
        if(not os.path.exists(path+"/Crops/"+str(tn))):
            os.mkdir(path+"/Crops/"+str(tn));
        
        imwrite(path+"/Crops/"+str(tn)+"/"+str(t)+".tif", img_cr);#saves cropped image 
        
        t+=1;
    
        
        pics = [imread(pic) for pic in glob.glob(f"{path}/Crops/"+str(tn)+"/*.tif")];    
        immgs = [Image.fromarray(pic1) for pic1 in pics];
        immgs[0].save(path+"/Crops/"+str(tn)+"/kym.gif", save_all=True, append_images=immgs[1:], duration=50, loop=0);
        


def getsinglecells(roi_od):
    single_cells={};
    for el in roi_od.keys():
        if "single_cell" in el:
            single_cells[el] = (roi_od[el]);

    return(single_cells);


def assignpositions(roi_od):
  
    for el in roi_od.keys(): #assigns position to each segmented cell
        roi_od[el]['pos'] = [(max(roi_od[el]['x']) + min(roi_od[el]['x']))/2 , (max(roi_od[el]['y']) + min(roi_od[el]['y']))/2];
    
    return roi_od;

def getmtimeframes(roi_od): #gets segemented cells in multiple time frames; from Fiji plugin output
    
    unique_keys =[];
    tim_f = {};
    ml_t = [];
    roi_ret ={};
    for el in roi_od.keys():
        if (el.split(":")[0] not in unique_keys):
            unique_keys.append(el.split(":")[0]);
            tim_f[el.split(":")[0]] = 1;
        else:
            tim_f[el.split(":")[0]] +=1;
         
    
    for el in tim_f.keys():
        if(tim_f[el] > 1):
            ml_t.append(el);
            
                 
    for el in roi_od.keys():
        if(el.split(":")[0] in ml_t):
            roi_ret[el] = roi_od[el];
    return(roi_ret)


def getalltimeframes(roi_od):
    unique_keys =[];
    tim_f = {};
    al_t = [];
    roi_ret ={};
    for el in roi_od.keys():
        if (el.split(":")[0] not in unique_keys):
            unique_keys.append(el[0:3]);
            tim_f[el.split(":")[0]] = 1;
        else:
            tim_f[el.split(":")[0]] +=1;
            
    mx_t = max(tim_f.values());
    
    for el in tim_f.keys():
        if(tim_f[el] == mx_t):
            al_t.append(el);
            
                 
    for el in roi_od.keys():
        if(el.split(":")[0] in al_t):
            roi_ret[el] = roi_od[el];
            

            
    return(roi_ret)



def getcellsinff(roi_od):
    roi_ret ={};
    for el in roi_od.keys():
        if("-1" in el):
             roi_ret[el] = roi_od[el];

    return(roi_ret);




def left_top(roi_od): #finds the left top most segmented cell
    x_min = 0;
    cx=0;
    y_min = 0;
    cy=0;
    left_layer = [];
    left_tp = '';
    for el in roi_od.keys():
        xs = roi_od[el]['x'];
        if(cx == 0  or  min(xs)<=x_min): #if the first case or if finds cells much above 
            x_min = min(xs);
            cx+=1;
            
    for el in roi_od.keys():
        xs = roi_od[el]['x'];
        if((x_min-20)<min(xs)<(x_min+20)):
            left_layer.append(el);
    

            
    for el in left_layer:
        ys = roi_od[el]['y'];
        if(cy == 0  or  min(ys)<=y_min ): #if the first case or if finds cells much above 
            left_tp = el;
            y_min = min(ys);
            cy+=1;
    
    return(left_tp);

def getleftop(roi_od):         
    lft_name = left_top(roi_od);
    lftp_pos = [ (max(roi_od[lft_name]['x']) + min(roi_od[lft_name]['x']))/2 , (max(roi_od[lft_name]['y']) + min(roi_od[lft_name]['y']))/2];
    
    return[lft_name,lftp_pos]



#def get_fl_trap()









# def getchn(ltp,od,y_max): #to get all trap in channel
#     chn = {};
#     found_d = 1;
#     nxt_pos = [];
#     l = (len(od.keys()));
#     cur_y = ltp[1];
    
#     #print(ltp);
    
#     while (cur_y<y_max): 
#         ct = 0;
#         cn = 0;
#         found_d = 0;
#         nxt_pos =[];
#         for el in od.keys():
#             #print(ltp);
            
#             if ((ltp[0]+75<od[el]['pos'][0]<ltp[0]+95) and (ltp[1]-20<od[el]['pos'][1]<ltp[1]+20)):#for immediate neighbour
#                     chn[el] = od[el];
#                     ct += 1;
#                     if(ct==2): 
#                         ltp = nxt_pos; #moves to next trap
#                         if(len(ltp)==0): 
#                             print(el);
#                             break;
#                         break;
                       
#             elif((ltp[0]-20<od[el]['pos'][0]<ltp[0]+20) and (ltp[1]+160<od[el]['pos'][1]<ltp[1]+180)):#trap right below
#                     chn[el] = od[el];
#                     ct +=1;
#                     nxt_pos = od[el]['pos'];
#                     found_d =1;
#                     #print(nxt_pos);
#                     if((ct == 2)): 
#                         ltp = nxt_pos; #moves to next trap 
#                         if(len(ltp)==0): 
#                             print(el);
#                             break;
#                         break;       
#             elif(cn == l-1):
#                 if(found_d==1):
#                     ltp = nxt_pos;#moves to next trap
#                     if(len(ltp)==0): 
#                         print(el);
#                         break;
#                     break;
#                 else:
#                     ltp[1] += 170;#moves to position where we expect next trap
#                     cur_y = ltp[1];
#                     break;
#             cn+=1;           
           
            
#     return(chn)



# def jumpright(ltp,ltn,od):
#     pos =[];
#     name ='';
#     for el in od.keys():    
#         if ((ltp[0]+710<od[el]['pos'][0]<ltp[0]+730) and (ltp[1]-20<od[el]['pos'][1]<ltp[1]+20)):#for further channerl
#                 pos = od[el]['pos'];
#                 name = el;
#                 break;
  
#     return[pos,name];
    



# def moveright(ltp,ltn,od):
#     pos =[];
#     name='';
#     for el in od.keys():    
#         if ((ltp[0]+225<od[el]['pos'][0]<ltp[0]+245) and (ltp[1]-20<od[el]['pos'][1]<ltp[1]+20)):#for neighbour channel
#                 pos = od[el]['pos'];
#                 name = el;
#                 break;
  
#     return[pos,name];
  

    






# def gettrapcell(lft_name,lftp_pos,y_max,x_max,roi_od):
#     traps=[];
#     x_pos=[];
#     cnt =0;
    
  
    
#     for i in range(0,6): #remember this gives i's in [0..5]
        
#         if (i==0): 
#             #print(lftp_pos);
#             chnl = getchn(lftp_pos,roi_od,y_max);
#             chnl[lft_name] = roi_od[lft_name];
#             traps.append(chnl);
#         elif(i%2 == 0): 
#             [lftp_pos2,lft_name2] = jumpright(lftp_pos,lft_name,roi_od); 
#             if(lftp_pos2):
#                 lftp_pos = lftp_pos2; lft_name = lft_name2;
#                 chnl = getchn(lftp_pos,roi_od); 
#                 chnl[lft_name] = roi_od[lft_name];
#                 traps.append(chnl);
#         else: 
#             [lftp_pos1,lft_name1] = moveright(lftp_pos,lft_name,roi_od); 
#             if(lftp_pos1):

#                 lftp_pos = lftp_pos1; lft_name = lft_name1;
#                 chnl = getchn(lftp_pos,roi_od); 
#                 chnl[lft_name] = roi_od[lft_name];
#                 traps.append(chnl);
         
    
            
#     return(traps);



