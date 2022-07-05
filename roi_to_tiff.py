
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:08:20 2022

@author: amansharma
REMEMBER TO HAVE YeaZ ACTIVATED IN CONDA ENV FOR USING read_roi
"""
#%% libraries
from tifffile import imwrite
#from all_funcs import crop_using_roi_tuple
#from all_funcs import extract_roi_coordinates
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
seg_zip_files = [f for f in listdir(zip_dir) if (isfile(join(zip_dir, f))  and ("seg" in f))];


zip_files = zip_files[1:];
zfs = zip_files.sort();


files =[];
for s in zip_files:
    files.append(s.replace('.zip', '')); 
    
st = "/Users/amansharma/Documents/Data/test/rois_im";
shutil.rmtree(st);
os.mkdir(st);

fld = [st+'/'+s for s in files];
for s in fld: os.mkdir(s);


#%% functions

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
        
        
        
        cropped = crop_using_roi_tuple(roi_coordinates,frames[0]);
        
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
        
    #print(np.shape(coordinates_list))
    return coordinates_list





def left_top(roi_od):
    x_min = 0;
    y_min = 0;
    left_layer = [];
    left_tp = '';
    for el in roi_od.keys():
        xs = roi_od[el]['x'];
        if(x_min == 0  or  min(xs)<x_min ): #if the first case or if finds cells much above 
            x_min = min(xs);
            
    for el in roi_od.keys():
        xs = roi_od[el]['x'];
        if((x_min-20)<min(xs)<(x_min+20)):
            left_layer.append(el);
    

            
    for el in left_layer:
        ys = roi_od[el]['y'];
        if(y_min == 0  or  min(ys)<y_min ): #if the first case or if finds cells much above 
            left_tp = el;
            y_min = min(ys);
    
    
    return(left_tp)



def getchn(ltp,od): #to get all trap in channel
    chn = [];
    nt_found = 0;
    nxt_pos = [];
    while nt_found==0: 
        ct = 0;
        for el in od.keys():    
            
            if ((ltp[0]+75<od[el]['pos'][0]<ltp[0]+95) and (ltp[1]-10<od[el]['pos'][1]<ltp[1]+10)):#for immediate neighbour
                    chn.append(el);
                    ct += 1;
                    if(ct==2): 
                        ltp = nxt_pos; #moves to next trap
                        break;
                       
            elif((ltp[0]-10<od[el]['pos'][0]<ltp[0]+10) and (ltp[1]+160<od[el]['pos'][1]<ltp[1]+180)):#trap right below
                    chn.append(el);
                    ct +=1;
                    nxt_pos = od[el]['pos'];
                    if(ct==2): 
                        ltp = nxt_pos; #moves to next trap 
                        break;
                       
        if(ct<2):
            nt_found =1;
            
            
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
  









def gettrapcell(roi_od):
    traps=[];
    x_pos=[];
    cnt =0;
    
    lft_name = left_top(roi_od);
    lftp_pos = [ (max(roi_od[lft_name]['x']) + min(roi_od[lft_name]['x']))/2 , (max(roi_od[lft_name]['y']) + min(roi_od[lft_name]['y']))/2];


    
    for el in roi_od.keys(): #assigns position to each "cell"
        roi_od[el]['pos'] = [(max(roi_od[el]['x']) + min(roi_od[el]['x']))/2 , (max(roi_od[el]['y']) + min(roi_od[el]['y']))/2];
        
    
    for i in range(0,4):
        
        if (i==0): 
            chnl = getchn(lftp_pos,roi_od);
            chnl.append(lft_name); 
            #print(len(chnl));
            traps.append(chnl);
        elif(i==2): 
            [lftp_pos2,lft_name2] = jumpright(lftp_pos,lft_name,roi_od); 
        
            if(lftp_pos2):
                lftp_pos = lftp_pos2; lft_name = lft_name2;
                chnl = getchn(lftp_pos,roi_od); 
                chnl.append(lft_name);
               # print(len(chnl));
                traps.append(chnl);
        else: 
            [lftp_pos1,lft_name1] = moveright(lftp_pos,lft_name,roi_od); 
            
            if(lftp_pos1):

                lftp_pos = lftp_pos1; lft_name = lft_name1;
                chnl = getchn(lftp_pos,roi_od); 
                chnl.append(lft_name); 
               # print(len(chnl));
                traps.append(chnl);
         
    
            
    return(traps)
                
 
    
    
    
    
    

#%% main function: 
cells =[]

for fls in files:

    if("seg" not in fls):    
        a=1;
        # roi_tuples = extract_roi_coordinates(zip_dir+fls+'.zip');    
        # fov = int(fls[-1]);
        # stf = fld[files.index(fls)];    
        # for j in range(0,len(roi_tuples)):
        #     roi_coordinates = roi_tuples[j][0];
        #     st = stf+"/roi_"+str(j+1);
        #     nd2_roi_extractor('/Volumes/GodardSSD/Microscope Images/Saransh/20220531_ylb128_gal2p.nd2',st,roi_coordinates,fov,j);
        
    elif("seg" in fls):
        
        roi_arr = read_roi_zip(zip_dir+fls+'.zip'); #this is an ordered dict with keys(get by roi_arr.keys()), and values are dict with name,x,y,n(no. of points) 
        trp_cel = gettrapcell(roi_arr);
        
        for el in roi_arr.keys():
            found =False;
            for ch in trp_cel:
                if el in ch: found = True;
                
            if not(found): cells.append(el)
                
                

print(len(cells))
                 
                
                
                
    
                
                
                
                
                
                
                
                