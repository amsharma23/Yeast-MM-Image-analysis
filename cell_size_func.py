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
import matplotlib.animation as animation
from matplotlib import image
import imageio
from PIL import Image
from scipy.spatial import ConvexHull


from os import walk
from os import listdir
from os.path import isfile, join


from read_roi import read_roi_file
from read_roi import read_roi_zip

import trap_ex_func as trapf
import cell_size_func as cellf

import collections
import shutil
#%% functions





    
def getycells(roi_check,roi_traps): #gets all the segmented yeast cells 
    roi_ret ={};
    roi_pos ={};
    for el in roi_check.keys():
        if(el not in roi_traps.keys()):
            roi_pos[el] = roi_check[el];
            
    for cel in roi_pos.keys():
        x = roi_pos[cel]['pos'][0];
        y = roi_pos[cel]['pos'][1];
        for el in roi_traps.keys():
            if ((roi_traps[el]['pos'][0]+25< x < roi_traps[el]['pos'][0]+45) and (roi_traps[el]['pos'][1]+10< y < roi_traps[el]['pos'][1]+30)):
                roi_ret[cel] = roi_pos[cel];
                break;
            elif((roi_traps[el]['pos'][0]-45< x < roi_traps[el]['pos'][0]-25) and (roi_traps[el]['pos'][1]+10< y < roi_traps[el]['pos'][1]+30)):
                roi_ret[cel] = roi_pos[cel];
                break;

    return(roi_ret);


def getarea(arr):
    ay = np.asarray(arr);
    area = ConvexHull(arr).area;
    
    return(area);
    
def getuniquecells(roi_od):
    uncel = {};
    names =[];
    for el in roi_od.keys():
        nm = el.split(':')[0];
        if (nm not in names):
            names.append(nm);     
            uncel[nm] = {'time_f' : 1};
            uncel[nm][uncel[nm]['time_f'] -1] = roi_od[el];
        else:
            uncel[nm]['time_f']+=1;
            uncel[nm][uncel[nm]['time_f'] -1] = roi_od[el];

    
    return(uncel)

def addarea(cells):
    for el in cells.keys():
        ts = cells[el]['time_f'];
        ims =[];
        for t in range(ts):
            x = cells[el][t]['x'];
            y = cells[el][t]['y'];
            cord = [];
            for ind in range(len(x)):
                cord.append([x[ind],-1*y[ind]]); #ConvexHull expects y coordinates to decrease as we go down(opposite to image array)
                   
            cells[el][t]['area'] = getarea(cord);
            
            
    return(cells)

def plot_AvT(cells,path):#saves area vs time plot
    fig, axs = plt.subplots();
    
    for el in cells.keys():
      ct=0;

        
      ts = cells[el]['time_f'];
      ar =[];
      tm =[];
      for t in range(ts):
          ar.append(cells[el][t]['area']);
          tm.append(15*t);
        
      if(len(ar)>=4): 
          if(not os.path.exists(path+"plots/"+el)):
              os.mkdir(path+"plots/"+el);
          axs.scatter(tm,ar);    
          plt.xlabel('minutes');
          plt.ylabel('pixel area');
          plt.savefig(path+"plots/"+el+"/"+str(el)+".png");
          axs.cla();
          ct+=1;
            
def savekym(cells,path):
    fig, axs = plt.subplots();
    for el in cells.keys():
            
        ts =cells[el]['time_f'];
        if(ts>=4):
            fig.clf();
            if(not os.path.exists(path+"plots/"+str(el)+"/kym")):
                os.mkdir(path+"plots/"+str(el)+"/kym");

            ims=[];
            ct=0;
            for t in range(ts):
                x = np.array(cells[el][t]['x']);
                y = np.array(cells[el][t]['y']);
                ax = plt.plot(x,y);
                plt.savefig(path+"plots/"+str(el)+"/kym/"+str(ct)+".png");   
                ims.append(imageio.imread(path+"plots/"+str(el)+"/kym/"+str(ct)+".png"));
                ct+=1;
            
            imageio.mimsave(path+"plots/"+str(el)+"/kym/"+"kym.gif", ims); 
            
            for fl in os.listdir(path+"plots/"+str(el)+"/kym/"):
               if fl.endswith(".png"): os.remove(path+"plots/"+str(el)+"/kym/"+fl);

            