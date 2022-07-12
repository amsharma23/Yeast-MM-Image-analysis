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

import shutil

#%%Functions for extracting traps


def gettraparr(roi_d):
    roi_seg_in_frames = getsinglecells(roi_d);#single cell segments deteced across in time frames

    roi_seg_in_frames = assignpositions(roi_seg_in_frames);#assigns each segmented single cell a position attribute

    roi_seg_in_mul_f = getmtimeframes(roi_seg_in_frames);#single cell in more than one frame

    roi_seg_in_all_f = getalltimeframes(roi_seg_in_mul_f); #single cell segments detected across all time

    roi_seg_in_ff = getcellsinff(roi_seg_in_all_f); #single cells in first frame

    [lftn,lftp] = getleftop(roi_seg_in_ff);# gets the top left most segmented cell and calculates the trap positions according to this


    roi_traps_arr = gettrapcell(lftn,lftp,roi_seg_in_ff);#segmented cells that are definitely traps

    return(roi_traps_arr)








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
    
    return(left_tp)



def getchn(ltp,od): #to get all trap in channel
    chn = {};
    found_d = 1;
    nxt_pos = [];
    l = (len(od.keys()));
    
    while found_d==1: 
        ct = 0;
        cn = 0;
        found_d = 0;
        nxt_pos =[];
        for el in od.keys():    
            if ((ltp[0]+75<od[el]['pos'][0]<ltp[0]+95) and (ltp[1]-20<od[el]['pos'][1]<ltp[1]+20)):#for immediate neighbour
                    chn[el] = od[el];
                    ct += 1;
                    if(ct==2): 
                        ltp = nxt_pos; #moves to next trap
                        
                        break;
                       
            elif((ltp[0]-20<od[el]['pos'][0]<ltp[0]+20) and (ltp[1]+160<od[el]['pos'][1]<ltp[1]+180)):#trap right below
                    chn[el] = od[el];
                    ct +=1;
                    nxt_pos = od[el]['pos'];
                    found_d =1;
                    #print(nxt_pos);
                    if((ct == 2)): 
                        ltp = nxt_pos; #moves to next trap 
                        break;       
            elif(cn == l-1):
                ltp = nxt_pos; #moves to next trap
                break;
            
            cn+=1;           
           
            
    return(chn)



def jumpright(ltp,ltn,od):
    pos =[];
    name ='';
    for el in od.keys():    
        if ((ltp[0]+710<od[el]['pos'][0]<ltp[0]+730) and (ltp[1]-20<od[el]['pos'][1]<ltp[1]+20)):#for further channerl
                pos = od[el]['pos'];
                name = el;
                break;
  
    return[pos,name];
    



def moveright(ltp,ltn,od):
    pos =[];
    name='';
    for el in od.keys():    
        if ((ltp[0]+225<od[el]['pos'][0]<ltp[0]+245) and (ltp[1]-20<od[el]['pos'][1]<ltp[1]+20)):#for neighbour channel
                pos = od[el]['pos'];
                name = el;
                break;
  
    return[pos,name];
  

def assignpositions(roi_od):
  
    for el in roi_od.keys(): #assigns position to each segmented cell
        roi_od[el]['pos'] = [(max(roi_od[el]['x']) + min(roi_od[el]['x']))/2 , (max(roi_od[el]['y']) + min(roi_od[el]['y']))/2];
    
    return roi_od;
    

def getleftop(roi_od):         
    lft_name = left_top(roi_od);
    lftp_pos = [ (max(roi_od[lft_name]['x']) + min(roi_od[lft_name]['x']))/2 , (max(roi_od[lft_name]['y']) + min(roi_od[lft_name]['y']))/2];
    
    return[lft_name,lftp_pos]





def gettrapcell(lft_name,lftp_pos,roi_od):
    traps=[];
    x_pos=[];
    cnt =0;
    
  
    
    for i in range(0,6): #remember this gives i's in [0..5]
        
        if (i==0): 
            chnl = getchn(lftp_pos,roi_od);
            chnl[lft_name] = roi_od[lft_name];
            traps.append(chnl);
        elif(i%2 == 0): 
            [lftp_pos2,lft_name2] = jumpright(lftp_pos,lft_name,roi_od); 
            if(lftp_pos2):
                lftp_pos = lftp_pos2; lft_name = lft_name2;
                chnl = getchn(lftp_pos,roi_od); 
                chnl[lft_name] = roi_od[lft_name];
                traps.append(chnl);
        else: 
            [lftp_pos1,lft_name1] = moveright(lftp_pos,lft_name,roi_od); 
            if(lftp_pos1):

                lftp_pos = lftp_pos1; lft_name = lft_name1;
                chnl = getchn(lftp_pos,roi_od); 
                chnl[lft_name] = roi_od[lft_name];
                traps.append(chnl);
         
    
            
    return(traps);


def getsinglecells(roi_od):
    single_cells={};
    for el in roi_od.keys():
        if "single_cell" in el:
            single_cells[el] = (roi_od[el]);

    return(single_cells);

def getalltimeframes(roi_od):
    max_t = 0;
    unique_keys =[];
    tim_f = {};
    al_t = [];
    roi_ret ={};
    for el in roi_od.keys():
        if (el[0:3] not in unique_keys):
            unique_keys.append(el[0:3]);
            tim_f[el[0:3]] = 1;
        else:
            tim_f[el[0:3]] +=1;
         
    mx_t = max(tim_f.values());
    
    for el in tim_f.keys():
        if(tim_f[el] == mx_t):
            al_t.append(el);
            
                 
    for el in roi_od.keys():
        if(el[0:3] in al_t):
            roi_ret[el] = roi_od[el];
            

            
    return(roi_ret)



def getmtimeframes(roi_od):
    roi_ret ={};
    for el in roi_od.keys():
        if("single_cell-" in el):
            roi_ret[el] = roi_od[el];
            
    return(roi_ret)




def getcellsinff(roi_od):
    roi_ret ={};
    l = len("single_cell-1");
    for el in roi_od.keys():
        if(el[(-1*l):] == "single_cell-1"):
             roi_ret[el] = roi_od[el];

    return(roi_ret);
