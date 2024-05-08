#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:57:19 2024

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
import networkx as nx
import matplotlib.pyplot as plt
#%% Read all the image files
get_settings().application.ipy_interactive = True
bin_folder = "/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/"
tif_files = os.listdir(bin_folder)
tif_files = [file.removesuffix('.tif') for file in tif_files if file.endswith(".tif")]
tif_files.sort()
print(tif_files)
#%%

def get_float_pos(st):
    st = re.split(r'[ \[\]]', st)
    pos = [float(element) for element in st if element != '']

    return pos

def add_point(viewer,point_pos_arr,feat,color_arr,new_point_pos):
    
    point_pos_arr.append(new_point_pos)
    color_arr.append('yellow')    
    viewer.add_points(point_pos_arr,features=feat,text = 'label' ,size=10,face_color = color_arr)

def load_image_and_skel(idx=0,curr_idx=0):
    print('IDX is: '+str(idx))
    
    global node_path
    global nd_pdf #the nodewise dataframe
    
    raw_im_path = '/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/'+tif_files[idx]+'.tif'
    skel_im_path = '/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/nellie_output/'+tif_files[idx]+'-ch0-im_pixel_class.ome.tif'
    node_path = '/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/'+tif_files[idx]+'_nodewise.csv'
    nd_pdf = pd.read_csv(node_path)
    
    imp_indxs =  nd_pdf['Degree of Node']!= (2 or 0) #only the tips and branch points are important to highlight
    indxs = np.arange(len(imp_indxs))[imp_indxs]
    np.append(indxs,curr_idx)
    
    imp_nds = nd_pdf.index.values[imp_indxs]
    imp_nds_pos_st = (nd_pdf.loc[imp_indxs,'Position(ZXY)'].values)
    imp_nds_deg = (nd_pdf.loc[imp_indxs,'Degree of Node'].values).astype(int)
    
      
            
      
    imp_nds_label = (nd_pdf.loc[imp_indxs,'Node#'].values).astype(str)
    imp_nds_cc = (nd_pdf.loc[imp_indxs,'CC(Island)#'].values).astype(str)
    
    node_label = [i+'-'+imp_nds_label[ni] for ni,i in enumerate(imp_nds_cc)]
    
    imp_nds_pos_fl = []
    curr_pos = get_float_pos(nd_pdf.loc[curr_idx,'Position(ZXY)'])
    for ni,i in enumerate(imp_nds_pos_st):
        pos = get_float_pos(i)
        imp_nds_pos_fl.append(pos)
    
    face_color_deg = []
    for deg in imp_nds_deg:
        print(type(deg))
        if(deg==1): face_color_deg.append('red')
        else: face_color_deg.append('green')

    face_color_deg.append('yellow')

     


    raw_im = imread(raw_im_path)
    skel_im = imread(skel_im_path)
    feat = {'degree': imp_nds_deg, 'label':node_label, 'index':indxs}
    
    
    return raw_im, skel_im, imp_nds_pos_fl, feat, face_color_deg



idx = 0
curr_node_idx = 0 

raw_im, skel_im, imp_nds_pos_fl, feat, face_color_deg = load_image_and_skel(idx,curr_node_idx)
viewer = napari.view_image(raw_im)
viewer.add_image(skel_im,name='skeleton') 
points_layer = viewer.add_points(imp_nds_pos_fl,features=feat,text = 'label' ,size=10,face_color = face_color_deg) 


@viewer.bind_key('\u2191')
def move_to_next_node(viewer):
    print('Move up')
    
        
@viewer.bind_key('r')
def remove_edge(viewer):
    print('Here')
    sel = list(viewer.layers[2].selected_data)
    indcs = list(viewer.layers[2].features['index'])
    if len(sel)==2:
        sel_indcs = [indcs[idx] for idx in sel]
        print('Selected nodes are:'+str(sel_indcs))
        neigh_0 = re.split(r'[ \[\]]',nd_pdf.loc[sel_indcs[0],'Neighbours'])
        neigh_0i =[int(element) for element in neigh_0 if element != '']
        neigh_1 = re.split(r'[ \[\]]', nd_pdf.loc[sel_indcs[1],'Neighbours'])
        neigh_1i =[int(element) for element in neigh_1 if element != '']
        
        if((sel_indcs[1] in neigh_0i) and (sel_indcs[0] in neigh_1i)):
            nd_pdf.loc[sel_indcs[0],'Degree of Node'] = (int(nd_pdf.loc[sel_indcs[0],'Degree of Node'])) - 1
            nneigh_0 = [element for element in neigh_0i if element != (sel_indcs[1])]
            nd_pdf.loc[sel_indcs[0],'Neighbours'] = str(nneigh_0)
            nd_pdf.loc[sel_indcs[1],'Degree of Node'] = (int(nd_pdf.loc[sel_indcs[1],'Degree of Node'])) - 1
            nneigh_1 = [element for element in neigh_1i if element != (sel_indcs[0])]
            nd_pdf.loc[sel_indcs[1],'Neighbours'] = str(nneigh_1)
            print(nd_pdf.loc[sel_indcs[0],'Neighbours'],nd_pdf.loc[sel_indcs[1],'Neighbours'])
            nd_pdf.to_csv(node_path, index=False)
            print("Modification saved")
                                                      
        else:
            print('Selected nodes are not neighbours')
        
             
    else:
        print('Please select only 2 nodes for edge removal')


@viewer.bind_key('n')
def move_on(viewer):
    global idx

    idx += 1
    
    viewer.layers.clear()
    raw_im, skel_im, imp_nds_pos_fl, feat, face_color_deg = load_image_and_skel(idx)
    viewer.add_image(raw_im)
    viewer.add_image(skel_im,name='skeleton') 
    points_layer = viewer.add_points(imp_nds_pos_fl,features=feat,text = 'label' ,size=10,face_color = face_color_deg) 
    