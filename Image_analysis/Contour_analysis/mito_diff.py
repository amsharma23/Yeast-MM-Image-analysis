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
from mpl_toolkits import mplot3d
import imageio
import numpy.ma as ma 
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

mit_cont_path = "/Users/amansharma/Documents/Data/2019_November_week5_liyana/Gluc/mito_cont";
img_nms_cont = os.listdir(mit_cont_path);
img_nms_cont.sort();


phto_conv_path= "/Users/amansharma/Documents/Data/2019_November_week5_liyana/Gluc/after_PC";
img_nms_pc = os.listdir(phto_conv_path);
img_nms_pc.sort();

imgs = [];
for nm in img_nms_cont:
	img = imread(os.path.join(mit_cont_path,nm));
	imgs.append(img);

gfp_contour_images = np.array(imgs);

imgs = []; 
for nm in img_nms_pc:
	img = imread(os.path.join(phto_conv_path,nm));
	imgs.append(img);
pc_images = np.array(imgs);

for i in range(len(img_nms_pc)):

	#print(np.shape(gfp_contour_images[i]));
	gfp_img = gfp_contour_images[i];

	

	inten_on_contour = ma.masked_array(pc_images[i],mask= gfp_img>27);
	
	mx =np.max(inten_on_contour);
	print(mx);
	mn =np.min(inten_on_contour);
	print(np.shape(inten_on_contour));

	x_arr = np.arange(0,len(inten_on_contour[0,:]),1);	
	y_arr = np.arange(0,len(inten_on_contour[:,0]),1)

	X,Y = np.meshgrid(x_arr,y_arr);


	zs = np.array(inten_on_contour[np.ravel(Y),np.ravel(X)]);
	Z = zs.reshape(X.shape);

	fig = plt.figure();
	ax = plt.axes();
	#ax.set_xlim(len(x_arr));
	#ax.set_ylim(len(y_arr));
	#ax.set_xlabel(r'X $\mu$m');
	#ax.set_ylabel(r'Y $\mu$m');
	#ax.set_zlabel('Intensity a.u.');
	plt.imshow(inten_on_contour,cmap='inferno');
	plt.colorbar();
	plt.savefig("/Users/amansharma/Documents/Data/2019_November_week5_liyana/Gluc/heat_map/"+str(i));
