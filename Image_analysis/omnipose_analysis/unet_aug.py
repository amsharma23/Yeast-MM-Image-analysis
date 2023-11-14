#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sept  12 13:41:08 2022

@author: amansharma
For Data Augmentation to train UNET; reshape and resize for prediction
"""
import random
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
import tifffile as tiff
import os
import cv2 as cv
import albumentations as A
 
img_fld = "/Users/amansharma/Documents/Data/test/UNET/Prediction data/BF/Trap3-8bit";
#mask_fld = "/Users/amansharma/Documents/Data/test/UNET/Training data/MasksBF/trap2_aug/resize_and_rescaled_8bit";

#immg_fld = "/Users/amansharma/Documents/Data/test/UNET/Training data/BF/Trap2_aug/";
#massk_fld = "/Users/amansharma/Documents/Data/test/UNET/Training data/MasksBF/trap2_aug/";


img_fls = os.listdir(img_fld);
#img_fls.remove('.DS_Store');
img_fls.sort();
#mask_fls = os.listdir(mask_fld);
#mask_fls.remove('.DS_Store');
#mask_fls.sort();

img_mask_dict_loc = {};

for i in range(len(img_fls)):
#	img_mask_dict_loc[i] = [os.path.join(img_fld,img_fls[i]),os.path.join(mask_fld,mask_fls[i])];
	img_mask_dict_loc[i] = os.path.join(img_fld,img_fls[i]);

aug = A.Compose([
    A.VerticalFlip(p=0.5),              
    A.RandomRotate90(p=0.5),
    A.GridDistortion(p=0.5),             
    A.CLAHE(p=0.8),
    A.RandomBrightnessContrast(p=0.8),    
    A.RandomGamma(p=0.8)]);

rss = A.PadIfNeeded(min_height=128, min_width=128, p=1);


random.seed(11);

for i in img_mask_dict_loc.keys():
	#print(img_mask_dict_loc[i][0],img_mask_dict_loc[i][1]);
	im = tiff.imread(img_mask_dict_loc[i]);
	#msk = tiff.imread(img_mask_dict_loc[i][1]);

#	augmented = aug(image=im,mask=msk);
	rss_d = rss(image=im);
	rs_d = rss_d['image'];
	
	# if(not (os.path.exists(os.path.join(immg_fld,"augmented")))):
	# 	os.mkdir(os.path.join(immg_fld,"augmented"));
	# if(not (os.path.exists(os.path.join(massk_fld,"augmented")))):
	#  	os.mkdir(os.path.join(massk_fld,"augmented"));

#	aug_img = augmented['image'];
#	aug_msk = augmented['mask'];

	nm_im = img_mask_dict_loc[i].split('/');
	nm_im = nm_im[-1];
	
#	nm_msk = img_mask_dict_loc[i][1].split('/');
#	nm_msk = nm_msk[-1];


	tiff.imwrite(os.path.join("/Users/amansharma/Documents/Data/test/UNET/Prediction data/BF/Trap3-8bit_rr",nm_im),rs_d);
#	tiff.imsave(os.path.join(massk_fld,"augmented/"+nm_msk),aug_msk);





