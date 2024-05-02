#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:57:19 2024

@author: amansharma
"""

#%%
import re
import time
from nellie import logger
from nellie.im_info.im_info import ImInfo
from nellie.utils.general import get_reshaped_image
import numpy as np
import scipy.ndimage as ndi
from scipy.spatial import cKDTree
import pandas as pd
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
 
def preview(idx=0):
    print('IDX is: '+str(idx))
    global node_path
    global nd_pdf    
    raw_im_path = '/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/'+tif_files[idx]+'.tif'
    skel_im_path = '/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/nellie_output/'+tif_files[idx]+'-ch0-im_pixel_class.ome.tif'
    node_path = '/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/'+tif_files[idx]+'_nodewise.csv'
    nd_pdf = pd.read_csv(node_path)
    
    imp_indxs =  nd_pdf['Degree of Node']!= (2 or 0)
    indxs = np.arange(len(imp_indxs))[imp_indxs]
    print(indxs)
    imp_nds = nd_pdf.index.values[imp_indxs]
    imp_nds_pos_st = (nd_pdf.loc[imp_indxs,'Position(ZXY)'].values)
    imp_nds_deg = (nd_pdf.loc[imp_indxs,'Degree of Node'].values).astype(int)
    imp_nds_label = (nd_pdf.loc[imp_indxs,'Node#'].values).astype(str)
    imp_nds_cc = (nd_pdf.loc[imp_indxs,'CC(Island)#'].values).astype(str)
    
    node_label = [i+'-'+imp_nds_label[ni] for ni,i in enumerate(imp_nds_cc)]
    
    imp_nds_pos_fl = []
    
    
    for ni,i in enumerate(imp_nds_pos_st):
        st_p = re.split(r'[ \[\]]', i)
        #deg = float(imp_nds_pos_deg[ni])
        pos = [float(element) for element in st_p if element != '']
        imp_nds_pos_fl.append(pos)
    

     


    raw_im = imread(raw_im_path)
    skel_im = imread(skel_im_path)
    feat = {'degree': imp_nds_deg, 'label':node_label, 'index':indxs}
    
    face_color_deg = ['red','green']
    return raw_im, skel_im, imp_nds_pos_fl, feat, face_color_deg



idx = 0

raw_im, skel_im, imp_nds_pos_fl, feat, face_color_deg = preview(idx)
viewer = napari.view_image(raw_im)
viewer.add_image(skel_im,name='skeleton') 
points_layer = viewer.add_points(imp_nds_pos_fl,features=feat,text = 'label' ,size=10, face_color = 'degree', face_color_cycle = face_color_deg) 


      
        
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
    raw_im, skel_im, imp_nds_pos_fl, feat, face_color_deg = preview(idx)
    viewer.add_image(raw_im)
    viewer.add_image(skel_im,name='skeleton') 
    points_layer = viewer.add_points(imp_nds_pos_fl,features=feat,text = 'label' ,size=10, face_color = 'degree', face_color_cycle = face_color_deg) 

            
           
    
    
#%% click handler
def on_click(layer, event):
    # Check if any points were clicked
    if event is not None and event.type == 'mouse_release':
        # Get the index of the clicked point
        index = event.data.view.index
        # Toggle the selected state of the clicked point
        points_layer.selected_data = [index]  # Select only the clicked point