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
from nd2reader import ND2Reader
import io
import math
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib import image as IM
import imageio
from PIL import Image

from scipy.spatial import ConvexHull
from skimage.draw import ellipse
from skimage.measure import label, regionprops, regionprops_table
from skimage.transform import rotate

import os
from os import listdir
from os.path import isfile, join

from read_roi import read_roi_file
from read_roi import read_roi_zip

import trap_ex_func as trapf
import cell_size_func as cellf

import collections
import shutil
#%% functions




def make_kym(pics,path): #makes a gif of volume est vs time
    vl = [];
    trj =[];
    immgs=[];
    fig = plt.figure();
    for pic in pics:
        pic_arr = imread(pic);
        label_img = label(pic_arr);
        regions = regionprops(label_img);
        vol_total = 0;#total volume across all regions
        fig, ax = plt.subplots(2);
        ax[0].imshow(pic_arr, cmap=plt.cm.gray);
        for props in regions:
            y0, x0 = props.centroid;
            orientation = props.orientation;
            x1 = x0 + math.cos(orientation) * 0.5 * props.axis_minor_length;
            y1 = y0 - math.sin(orientation) * 0.5 * props.axis_minor_length;
            x2 = x0 - math.sin(orientation) * 0.5 * props.axis_major_length;
            y2 = y0 - math.cos(orientation) * 0.5 * props.axis_major_length;
                
            ax[0].plot((x0, x1), (y0, y1), '-r', linewidth=2.5);
            ax[0].plot((x0, x2), (y0, y2), '-r', linewidth=2.5);
            ax[0].plot(x0, y0, '.g', markersize=15);

            minr, minc, maxr, maxc = props.bbox;
            bx = (minc, maxc, maxc, minc, minc);
            by = (minr, minr, maxr, maxr, minr);
            ax[0].plot(bx, by, '-b', linewidth=2.5);
                
            vol = (1/6)*math.pi*props.axis_major_length*(props.axis_minor_length**2)*(0.11**3);
            vol_total+=vol;
            
        vl.append(vol_total);
        ax[1].plot(np.array(range(len(vl)))*15,vl);
        ax[1].scatter(np.array(range(len(vl)))*15,vl,color='red');
        ax[1].set_ylabel('Vol μm^3');
        ax[1].set_xlabel('Time(mins)');
            
        img_buf = io.BytesIO();
        plt.savefig(img_buf, format='png');
        im =Image.open(img_buf);
        immgs.append(im);
        
    
    
    immgs_o = immgs[0];
    immgs_o.save(path+"/kym_vol.gif",format="GIF", append_images=immgs,save_all=True, duration=100, loop=0);
    img_buf.close();
    fig.clf();
    plt.close();
    return(vl);






def el_vol(path):
    traj={};
    t=1;
    path_crp = path+"Crops/"
    if(os.path.exists(path_crp)):#checks if the path provided actually has any traps with yeasts or not
        pics_folders = [join(path_crp,f) for f in listdir(path_crp) if not(isfile(join(path_crp,f)))];
        for fld in pics_folders:
            pics = [join(fld,f) for f in listdir(fld) if ('.tif' in f)];
            trj = make_kym(pics,fld);
            traj[t] = trj;
            t+=1;
        return(traj);
    else:
        print("No filled traps found in "+path);    
            
         

def plt_trj(trajs,path):
    fig = plt.figure();
    for el in trajs.keys():#goes over fovs
        for trjs in trajs[el]:#goes over trajectories 
            plt.plot(np.array(range(len(trjs)))*15,trjs);

    plt.savefig(path+"/trajs.tif",format="TIFF");
    print("Plotted trajectories for "+path);
    fig.clf();
    plt.close();


def plt_dist(trajs,path):
    fig = plt.figure();
    vols = [];
    for el in trajs.keys():#goes over fovs
        for trjs in trajs[el]:#goes over trajectories 
            for v in trjs:
                vols.append(v);
    print(len(vols));
    ser = pd.Series(np.array(vols));
    fig = ser.plot.hist(bins=20);
    plt.savefig(path+"/dist.tif",format='TIFF');

    fig.clf();
    plt.close();






















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

    
    return(uncel);

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
            
            
    return(cells);



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
          
          if(not os.path.exists(path+"plots/")):
              os.mkdir(path+"plots/");
          
          
          
          if(not os.path.exists(path+"plots/"+el)):
              os.mkdir(path+"plots/"+el);
          
            
          axs.scatter(tm,ar);    
          plt.xlabel('minutes');
          plt.ylabel('pixel area');
          plt.savefig(path+"plots/"+el+"/"+str(el)+".png");
          axs.cla();
          ct+=1;
            
def savekym(cells,path):
    fig, ax = plt.subplots();
    for el in cells.keys():
            
        ts =cells[el]['time_f'];
        if(ts>=4):
           
            
            if(not os.path.exists(path+"plots/"+el)):
                #print('here');
                os.mkdir(path+"plots/"+el);
            
            
            if(not os.path.exists(path+"plots/"+el+"/kym")):
                #print('here1');
                os.mkdir(path+"plots/"+el+"/kym");
            
            buf = io.BytesIO();
            ims=[];
            ct=0;
            for t in range(ts):
                fig.clf();
                x = np.array(cells[el][t]['x']);
                y = np.array(cells[el][t]['y']);
                ell = plt.plot(x,y);
                plt.savefig(buf,format='tif');
               
                ims.append(imageio.open(buf));
                buf.flush();
                ct+=1;
                
            buf.close();    
        
            imageio.mimsave(path+"plots/"+el+"/kym/kym.gif",ims);
            
            

            