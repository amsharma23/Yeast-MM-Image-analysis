#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified on Mon Apr 08 2024 : original code - https://github.com/aelefebv/nellie-supplemental/blob/main/case_studies/GNN/graph_frame.py ; https://github.com/aelefebv/nellie/tree/main/nellie

@author: amansharma
"""
#%% Library imports
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
import napari
import networkx as nx
import matplotlib.pyplot as plt
#%% Read all the image file and with the Nellie output folder
bin_folder = "/Users/amansharma/Documents/Data/Nellie_test/GUI_test/bin_30_40_manual/"
tif_files = os.listdir(bin_folder)
tif_files = [file.removesuffix('.tif') for file in tif_files if file.endswith(".tif")]
tif_files.sort()

print(tif_files)
#%% Functions used in the nodewise extraction of properties

def graph_gen(tree):#Generates a graph from the skeleton image file
    graph = nx.Graph()
    
    if (len(tree.neighbors)!=0):
        for nn,neigh in enumerate(tree.neighbors):
            existing_edges = set()
            for nh in neigh:
                ed = tuple(sorted([str(nn), str(nh)]))
                if ed not in existing_edges: 
                    graph.add_edge(str(nn), str(nh))
                    existing_edges.add(ed)
    
        return graph

    else:        
        return graph


    
def nodewise_props(tree, graph, sn,number_of_nodes,nodew_path,df_nodewise=None, visited=None):#Traverses the skeleton graph in a depth first manner
    
    if visited is None:
        visited = set()

    if df_nodewise is None:
        df_nodewise = pd.read_csv(nodew_path, encoding='utf-8')
        
    insert_loc = df_nodewise.index.max()
    

    if pd.isna(insert_loc):
        insert_loc = 0    
    else:
        insert_loc = insert_loc+1
    
    
    df_nodewise.loc[insert_loc, 'CC(Island)#'] = tree.label
    df_nodewise.loc[insert_loc, 'Node#'] = int(sn)
    df_nodewise.loc[insert_loc, 'Degree of Node'] = int(graph.degree(sn))
    df_nodewise.loc[insert_loc, 'Position(ZXY)'] = str(tree.voxel_idxs[int(sn)])
    print(type(list(graph.neighbors(sn))[0]))
    df_nodewise.loc[insert_loc, 'Neighbours'] = str(list(graph.neighbors(sn)))

    
    if((insert_loc+1)==number_of_nodes):
        print(df_nodewise)
        df_nodewise.to_csv(nodew_path,encoding='utf-8',index=False)
        
    visited.add(sn)

    for neighbor in graph[sn]:
        if neighbor not in visited:
            nodewise_props(tree, graph, neighbor, number_of_nodes,nodew_path, df_nodewise, visited)        

def save_graph_fig(tree,graph,bf,file):
    
    plt.figure()
    nx.draw_kamada_kawai(graph, with_labels=True, font_weight='bold')
    
    #pos = nx.spring_layout(graph)
    #nx.draw_networkx_nodes(graph, pos,label=graph.nodes, node_color = 'r', node_size = 50, alpha = 1)
    ax = plt.gca()
    
    
    plt.axis('off')
    file_path =  os.path.join(bf,file+'_'+str(tree.label)+'.png')
    plt.savefig(file_path,dpi=700)
    plt.title(str(file)+'_'+str(tree.label))
    #plt.show()
    ax.clear()
    plt.clf()
    
def collected_prop(tf,df_average,bf):

    graph_files = os.listdir(bf)
    graph_files = [file.removesuffix('.gml') for file in graph_files if file.endswith('.gml')]
    graph_files.sort()
    
    for nfil,fil in enumerate(graph_files):
        s_deg,nodes,tips, junc, loops = [0,0,0,0,0]
        graph = nx.read_gml(os.path.join(bf,fil+'.gml'))
        nodes = nodes + graph.number_of_nodes()
        
        deg_z_n = sum(1 for node, degree in graph.degree() if degree == 1)
        deg_th_n = sum(1 for node, degree in graph.degree() if degree >= 3)
        
        s_deg = s_deg+ sum(degree for node,degree in graph.degree)
        tips = tips + deg_z_n
        junc = junc + deg_th_n
        
        cyc_b = list(nx.cycle_basis(nx_graph))
        loops = loops + len(cyc_b)
        
        df_average.loc[nfil,'File#'] = fil 
        df_average.loc[nfil,'#Nodes'] = nodes  
        df_average.loc[nfil,'#Tips(Deg1)'] = tips    
        df_average.loc[nfil,'#Junctions(Deg3+)'] = junc    
        df_average.loc[nfil,'Avg Deg'] = s_deg/nodes
        df_average.loc[nfil,'#Loops'] = loops
    
    df_average.to_csv(os.path.join(bin_folder,'Collected_prop.csv'), encoding='utf-8')
    
#%%Part of the supplemental code of Nellie
class Tree:
    def __init__(self, label: int, voxel_idxs: np.ndarray, global_idxs):
        self.label = label
        self.voxel_idxs = voxel_idxs
        self.global_idxs = global_idxs
        self.neighbors = []
        self.start_node = None
        self.jump_distances = None
        self.nodelists = None
        self.multiscale_edge_list = None

    def get_neighbors(self):
        ckdtree = cKDTree(self.voxel_idxs)
        self.tree = ckdtree.tree
        self.neighbors = ckdtree.query_ball_point(self.voxel_idxs, r=1.733)  # a little over sqrt(3) this in term of voxel distance
        self.neighbors = [np.array(neighbor) for neighbor in self.neighbors] #list of neighbor of each node + node itself: node 0 nNeigh - [0,1](say)
        self.neighbors = [neighbor[neighbor != i] for i, neighbor in enumerate(self.neighbors)] #removes node itself: node 0 nNeigh - [1]

    def get_start_node(self):
        # pick the first node with only one neighbor. If none exists, pick the first node
        for i, neighbor in enumerate(self.neighbors):
            if len(neighbor) == 1:
                self.start_node = i
                return
        self.start_node = 0

#%%Obtain a ckdTree [nearest neighbour dataframe] and graph of the skeleton file - https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
for n_file,file in enumerate(tif_files):
    
    
    pix_calss_path = os.path.join(bin_folder,"nellie_output",file+"-ch0-im_pixel_class.ome.tif")
    pix_class = imread(pix_calss_path) #reads all the 3d coordinates and their values
    
    df_nodewise = pd.DataFrame(columns=['CC(Island)#','Node#','Degree of Node','Position(ZXY)','Neighbours'])
    df_nodewise.to_csv(os.path.join(bin_folder,file+'_nodewise.csv'),index=False)
    nodewise_path = os.path.join(bin_folder,file+'_nodewise.csv')
    
    tree_labels,_ = ndi.label(pix_class , structure=np.ones((3, 3, 3)))#finds and labels topological islands
    
    valid_coords = np.argwhere(tree_labels > 0) #coordinates that are non-zero
    valid_coord_labels = tree_labels[valid_coords[:, 0], valid_coords[:, 1], valid_coords[:, 2]]
    
    unique_labels = np.unique(tree_labels)
    
    
    trees=[]
    
    for label_num, label in enumerate(unique_labels):
    
                if label == 0:
                    continue
                global_idxs = np.argwhere(valid_coord_labels == label).flatten().tolist() #chooses one connected component
                tree = Tree(label, valid_coords[valid_coord_labels == label], global_idxs) #converst the 3D skeleton image to a ckdTree 
            
                tree.get_neighbors()
                tree.get_start_node()
                trees.append(tree)
    
    print('Number of trees are:'+str(len(trees)))
    ttl_nodes = 0
    for nt,tre in enumerate(trees):
        
        nx_graph = nx.Graph()
        start_node = str(tre.start_node)
        nx_graph = graph_gen(tre)
        print('Start node is: '+str(start_node))
        print(nx_graph)
        ttl_nodes = ttl_nodes + nx_graph.number_of_nodes()
        print('Total Number of nodes:'+str(ttl_nodes))
        
        if(nx_graph):
            print('Here for '+str(n_file)+' '+str(tre.label))
            nodewise_props(tre,nx_graph,start_node,ttl_nodes,nodewise_path)
            nx.write_gml(nx_graph, os.path.join(bin_folder,file+'_'+str(tre.label)+'.gml'))#saves graph as gml
            save_graph_fig(tre,nx_graph,bin_folder,file)
        
        #
        
    
 
#%%Collects the network level properties to a CSV
df_average = pd.DataFrame(columns=['File#','#CC(Islands)','#Nodes','#Tips(Deg1)','#Junctions(Deg3+)','#Loops','Avg Deg'])
collected_prop(tif_files,df_average, bin_folder)


