#Import all packages

from nd2reader import ND2Reader
import matplotlib.pyplot as plt
from skimage import data
import napari
from skimage.data import astronaut
import cv2
import numpy as np
import tifffile as tiff
from nd2reader import ND2Reader
import numpy as np
import glob
import csv
import os

loc = "/Users/amansharma/Documents/Data/test/UNET/Training_data/BF/Trap1_aug/resize_and_rescaled_8bit";
loc1 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/MasksBF/trap1_aug/resize_and_rescaled_8bit";
dest = "/Users/amansharma/Documents/Data/test/UNET/Training_data/cellpose/train";

imgs = os.listdir(loc);
#print(imgs);
if('.DS_Store' in imgs): imgs.remove('.DS_Store');
imgs.sort();

msks = os.listdir(loc1);
if('.DS_Store' in msks): msks.remove('.DS_Store');
msks.sort();

nu=0;
for i in range(len(imgs)):
	img_n = imgs[i];
	img = tiff.imread(os.path.join(loc,img_n));

	tiff.imsave(os.path.join(dest,'trap100'+str(nu+i)+'_img.tif'),img);

for i in range(len(msks)):
	msk_n = msks[i];
	msk = tiff.imread(os.path.join(loc1,msk_n));

	tiff.imsave(os.path.join(dest,'trap100'+str(nu+i)+'_mask.tif'),msk);

