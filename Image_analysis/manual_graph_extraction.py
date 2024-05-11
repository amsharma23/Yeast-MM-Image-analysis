#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 15:19:55 2024

@author: amansharma
"""


#%%Library load
import re
import time
from nellie import logger
from nellie.im_info.im_info import ImInfo
from nellie.utils.general import get_reshaped_image
import numpy as np
import scipy.ndimage as ndi
from scipy.spatial import cKDTree
import pandas as pd
from magicgui import magicgui
import math
from tifffile import imwrite
from tifffile import imread
import os
from napari.settings import get_settings
import napari
from natsort import natsorted
import networkx as nx
import matplotlib.pyplot as plt
#%% Read all the image files
get_settings().application.ipy_interactive = True
bin_folder = "/Users/amansharma/Documents/Data/Nellie_test/GUI_test/time_series_test"
tif_files = os.listdir(bin_folder)
tif_files = [file.removesuffix('_nodewise.csv') for file in tif_files if file.endswith('_nodewise.csv')]
tif_files = natsorted(tif_files)

print(tif_files)
#%%
def get_float_pos(st):
    st = re.split(r'[ \[\]]', st)
    pos = [int(element) for element in st if element != '']

    return pos
def load_image_and_skel(bf,idx=0):
    print('File number is: '+str(tif_files[idx]))
    
    global node_path
    global nd_pdf #the nodewise dataframe
    raw_im_path = os.path.join(bf,'mito_richglu.ome.tiff')
    
    skel_im_path = os.path.join(bf,'nellie_output/mito_richglu.ome-ch0-im_pixel_class.ome.tif')
    
    ts_raw_imgs = imread(raw_im_path)
    ts_skel_imgs = imread(skel_im_path)
    color_ex = []
    tp = int(re.split(r'[ _ ]', tif_files[idx])[-1])
    print('Timpe Point is:'+str(tp))
    
    node_path = os.path.join(bf, tif_files[idx]+'_etxracted.csv')
    print(os.path.exists(node_path))
    raw_im = ts_raw_imgs[tp]
    skel_im = ts_skel_imgs[tp] 
    skel_im = np.transpose(np.nonzero(skel_im))
    
    if (os.path.exists(node_path)):
        nd_pdf = pd.read_csv(node_path)
        if(not(pd.isna(nd_pdf.index.max()))):
            
            
            pos_ex = nd_pdf['Position(ZXY)'].values
            deg_ex = nd_pdf['Degree of Node'].values.astype(int)
            pos_ex = [get_float_pos(el) for el in pos_ex]
            print(pos_ex)
            
            matching_indices = np.argwhere(np.all(skel_im[:, None, :] == pos_ex, axis=2))
            matching_indices = (matching_indices[:,0])
            face_color_arr = ['red' for i in range(len(skel_im))]
            print(len(face_color_arr))
            
            for ni,i in enumerate(matching_indices):
              
                if deg_ex[ni] == 1: color_ex.append('blue')
                else: color_ex.append('green')


        else:
            print('HERE')
            pos_ex = []
            deg_ex = []
            nd_pdf =pd.DataFrame(columns=['Degree of Node','Position(ZXY)'])
            node_path = os.path.join(bf, tif_files[idx]+'_etxracted.csv')
            face_color_arr = ['red' for i in range(len(skel_im))]
            print(len(face_color_arr))
            nd_pdf.to_csv(node_path,index=False)    
        
        
    else: 
        print('HERE')
        pos_ex = []
        deg_ex = []
        nd_pdf =pd.DataFrame(columns=['Degree of Node','Position(ZXY)'])
        node_path = os.path.join(bf, tif_files[idx]+'_etxracted.csv')
        face_color_arr = ['red' for i in range(len(skel_im))]
        print(len(face_color_arr))
        nd_pdf.to_csv(node_path,index=False)
    
    
    
    
    
    
    
            
    return raw_im, skel_im,face_color_arr,pos_ex,color_ex

idx = 0
curr_node_idx = 0 

raw_im, skel_im,fca,imp_l,imp_c = load_image_and_skel(bin_folder,idx)
viewer = napari.view_image(raw_im)
#viewer.add_image(skel_im,name='skeleton') 
points_layer = viewer.add_points(skel_im,size=3,face_color = fca) #
if(len(imp_l)!=0):
    p_l = viewer.add_points(imp_l,size=5,face_color=imp_c,name='imp_l')

@viewer.bind_key('b')
def load_junction(viewer):
    ind = list(viewer.layers[1].selected_data)[0]
    
    pos =(viewer.layers[1].data[ind])
    
    insert_loc = nd_pdf.index.max()
    

    if pd.isna(insert_loc):
        insert_loc = 0    
    else:
        insert_loc = insert_loc+1
    
    nd_pdf.loc[insert_loc,'Degree of Node'] = 3
    nd_pdf.loc[insert_loc,'Position(ZXY)'] = str(pos)
    if(len(viewer.layers)>2):
        pos_ex = list(viewer.layers[2].data)
        if (any(np.array_equal(pos, arr) for arr in pos_ex)):
            pos_ex[-1] = pos
            viewer.layers[2].data = pos_ex
            print(list(viewer.layers[2].data))
            color_ex = list(viewer.layers[2].face_color)
            
        else:
            pos_ex.append(pos)
            viewer.layers[2].data = pos_ex
            print(list(viewer.layers[2].data))
            color_ex = list(viewer.layers[2].face_color)
           
    else:        
        p_l = viewer.add_points(pos,size=5,face_color=[[0.,1.,0.,1.]],name='imp_l')
        nd_pdf.to_csv(node_path,index=False)
    
@viewer.bind_key('t')
def load_tip(viewer):
    ind = list(viewer.layers[1].selected_data)[0]
    
    pos =(viewer.layers[1].data[ind])
    
    insert_loc = nd_pdf.index.max()
    

    if pd.isna(insert_loc):
        insert_loc = 0    
    else:
        insert_loc = insert_loc+1
    
    nd_pdf.loc[insert_loc,'Degree of Node'] = 1
    nd_pdf.loc[insert_loc,'Position(ZXY)'] = str(pos)
    if(len(viewer.layers)>2):
        pos_ex = list(viewer.layers[2].data)
        print(type(pos_ex[0]))
        if (any(np.array_equal(pos, arr) for arr in pos_ex)):
            pos_ex[-1] = pos
            viewer.layers[2].data = pos_ex
            print(list(viewer.layers[2].data))
            color_ex = list(viewer.layers[2].face_color)
            
        else:
            pos_ex.append(pos)
            viewer.layers[2].data = pos_ex
            print(list(viewer.layers[2].data))
            color_ex = list(viewer.layers[2].face_color)
            
            
            color_ex[-1] = [0.,0.,1.,1.]
            viewer.layers[2].face_color = color_ex
            nd_pdf.to_csv(node_path,index=False)

    else:        
        p_l = viewer.add_points(pos,size=5,face_color=[[0.,0.,1.,1.]],name='imp_l')
        nd_pdf.to_csv(node_path,index=False)

@viewer.bind_key('r')
def remove_special_node(viewer):
    ind = list(viewer.layers[2].selected_data)
    nd_pdf0 = pd.read_csv(node_path)
    print(nd_pdf0)
    if(len(ind)==0):
        print('Please select points only from "imp_l"')
    else:
        pos =(viewer.layers[2].data[ind][0])
        ind = [get_float_pos(st) for st in list(nd_pdf0['Position(ZXY)'])]
        for ni,i  in enumerate(ind):
            if(all(x == y for x, y in zip(i, pos)) and len(pos) == len(i)):
                print(ni)
                nd_pdf0.drop((ni),inplace=True)
                nd_pdf0.to_csv(node_path,index=False)
                raw_im, skel_im,fca,imp_l,imp_c = load_image_and_skel(bin_folder,idx)
                viewer.layers.remove('imp_l')
                if(len(imp_l)!=0):
                    p_l = viewer.add_points(imp_l,size=5,face_color=imp_c,name='imp_l')

@viewer.bind_key('n')
def move_on(viewer):
    global idx

    idx = (idx +1)%len(tif_files)
    
    viewer.layers.clear()
    raw_im, skel_im, fca,imp_l,imp_c  = load_image_and_skel(bin_folder,idx)
    viewer.add_image(raw_im)
    points_layer = viewer.add_points(skel_im,size=3,face_color = fca) #
    if(len(imp_l)!=0):
        p_l = viewer.add_points(imp_l,size=5,face_color=imp_c,name='imp_l')
    else:
        viewer.layers.remove('imp_l')
