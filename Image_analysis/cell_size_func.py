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
#from all_funcs import crop_using_roi_tuple
#from all_funcs import extract_roi_coordinates

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
import tqdm
#from read_roi import read_roi_file
#from read_roi import read_roi_zip

#import trap_ex_func as trapf
#import cell_size_func as cellf

import collections
import shutil
#%% functions

"""
Here the cell size functions can take in 2 types of segmented images-segmented images of FoV and segmented images of just the trap
The difference would be in the identification and verification of a segment being a trap or a cell.

I) In the case of the entire FoV image being segmented - the traps are located using their fixed relative positions
II) The trap images segmented - the trap volume is just eliminated from the volume (~180.12 μm^3)

"""

def make_kym(seg_pics,bf_pics,ar_m,ar_d,path,t_s): #makes a gif of volume est vs time
    


    trj =[];
    immgs=[];
    n = len(seg_pics.keys());
    t = np.arange(start=0,stop=n,step=1) * t_s/3600;

    fig = plt.figure(figsize=(8,4));
    bf_nm = list(bf_pics.keys());

    print('aaaaa');


    ax1 = fig.add_subplot(2,2,1);
    ax2 = fig.add_subplot(2,2,2);
    ax3 = fig.add_subplot(2,1,2);
    ax3.set_xlim([t[0],t[-1]]);


    for p,pkee in (enumerate(seg_pics.keys())):
        
        seg_pic_arr = seg_pics[pkee];
        seg_pic_arr[seg_pic_arr!=0] = 255;
        bf_pic_arr = bf_pics[bf_nm[p]];

        ax1.imshow(bf_pic_arr,cmap=plt.cm.gray);
        ax2.imshow(seg_pic_arr,cmap=plt.cm.gray);




        
        ax3.plot(t[0:p],ar_m[0:p],label='Mother',color='blue');
        ax3.scatter(t[0:p],ar_m[0:p],color='red',s=10);
        ax3.plot(t[0:p],ar_d[0:p],label='Daughter',color='orange');
        ax3.scatter(t[0:p],ar_d[0:p],color='blue',s=10);

        #ax3.legend();
        ax3.set_ylabel('Area $μm^2$');
        ax3.set_xlabel('Time(hrs)');
            
        img_buf = io.BytesIO();
        plt.savefig(img_buf, format='png');
        im =Image.open(img_buf);
        immgs.append(im);
    
    
    
    
    immgs_o = immgs[0];
    immgs_o.save(path+"/kym_vol.gif",format="GIF", append_images=immgs,save_all=True, duration=500, loop=0);
    img_buf.close();
    fig.clf();
    plt.close();







def ar_plot(images,n,t_s,path):

    ar_m =np.zeros((n,1));
    ar_d =np.zeros((n,1));

    for i,kee in  enumerate(images.keys()):

        im = images[kee];
        regions = regionprops(im);

        ars = [prop.area for prop in regions];
        ars = np.array(ars);
        try:
            ar_m[i] = np.amax(ars);
            ar_d[i] = ars[ars!=np.amax(ars)].sum();
        except:
            ar_m[i] = 0;
            ar_d[i] = 0;
        


    ar_m = ar_m * 0.325 * 0.325; #pixel size(andor)  = 6.5μm , 20x => 6.5/20 = 0.325μm
    ar_d = ar_d * 0.325 * 0.325;
    t = np.arange(start=0,stop=n,step=1) * t_s/3600;

    fig,ax = plt.subplots();

    ax.plot(t, ar_m,label='Mother');
    ax.scatter(t,ar_m,color='red',s=10);
    ax.plot(t, ar_d,label='Daughter');
    ax.scatter(t,ar_d,color='blue',s=10);
    ax.set_ylabel('Area $μm^2$');
    ax.set_xlabel('Time(hrs)');
    ax.legend();
    plt.savefig(path+"/area_count"+str(n)+".svg",format='svg');
    plt.savefig(path+"/area_count"+str(n)+".png",format='png')
    fig.clf();
    plt.close();



    return(ar_m,ar_d);


def el_vol(path,istrps):
    traj={};
    t=1;
    if(not istrps): path_crp = path+"Crops/"; #path of cropped segmented images
    else: 
        seg_im_pth  = [f for f in listdir(path) if ("seg" in f)];
        path_crp = path+seg_im_pth[0];

    
    if(os.path.exists(path_crp)):#checks if the path provided actually has any traps with yeasts or not
        pics_folders = [join(path_crp,f) for f in listdir(path_crp) if not(isfile(join(path_crp,f)))];
        for fld in pics_folders:
            pics = [join(fld,f) for f in listdir(fld) if ('.tif' in f)];
            #print(fld);
            trj = make_kym(pics,fld,istrps);
            traj[t] = trj;
            t+=1;
        return(traj);
    else:
        print("No filled traps found in "+path);    
            
def org_plot(path,n):
    orgn = np.zeros((18,1));
    
    flr_pth = path+"fluor/";
    for i in range(5,8):

        z_pth = flr_pth+str(i)+"/";
        images_pth = listdir(z_pth);
        images_pth.sort();
        
        j =0; 
        for files in images_pth:
            if(files!='.DS_Store') and ('.tif' in files):
                im = imread(z_pth+files);
                orgn[j] += np.sum(im);

                j+=1;
    x  = list(range(len(orgn)));
    x = [x[i]*15 for i in range(len(orgn))];
    fig,ax = plt.subplots();

    ax.plot(x, orgn);
    ax.scatter(x,orgn,color='red');
    ax.set_ylabel('AU');
    ax.set_xlabel('Time(mins)');
    plt.savefig(path+"/orgn_vol"+str(n)+".svg",format='svg');
    fig.clf();
    plt.close();

    orgnn =[orgn[i] for i in range(18)];
    return(orgnn);



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
                if(v!=0):
                    vols.append(v);
    print(len(vols));
    ser = pd.Series(np.array(vols));
    fig = ser.plot.hist(bins=50);
    plt.savefig(path+"/dist.tif",format='TIFF');

    fig.cla();
    plt.close();










