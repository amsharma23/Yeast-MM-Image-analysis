#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:30:08 2022

@author: amansharma
For testing of UNET
Relying majorly on JunLab github training notebook and blog post by Vidushi Bhatia(https://medium.com/geekculture/u-net-implementation-from-scratch-using-tensorflow-b4342266e406)
Snippets used from blog post have been marked

"""

import random
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
#from new_unet_training import bce_dice_loss
from tensorflow.keras import layers
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Conv2DTranspose
from tensorflow.keras.layers import concatenate
from tensorflow.keras import models
from tensorflow.keras import losses
from tensorflow.keras.losses import binary_crossentropy
from sklearn.model_selection import train_test_split
import tifffile as tiff
import os
import cv2 as cv



model_path = "/Users/amansharma/Documents/Data/test/UNET/Model/2022_10_06_unet.hdf5";
unet = models.load_model(model_path,compile=False);
unet.summary();

imgs_fld = "/Users/amansharma/Documents/Data/test/UNET/Prediction data/BF/Trap3-8bit_rr/";
imgs_fls = os.listdir(imgs_fld);
imgs_flls = [os.path.join(imgs_fld,i) for i in imgs_fls if ('.tif') in i];
imgs_fls.sort();


imgs = [];
for i in imgs_flls:
	with tiff.TiffFile(i) as img_tiff:
			cur_img = img_tiff.asarray() / 256.0 ;
	cur_img = cur_img.reshape(128,128,1);
	cur_img = np.expand_dims(cur_img,axis=3);
	cur_img = cur_img[np.newaxis, ...];
	imgs.append(cur_img);

print(np.shape(imgs));

for i in range(len(imgs)): 
#	print(immgs.shape);	
	pred = unet.predict(imgs[i]);
	print(pred.shape);
	tiff.imwrite("/Users/amansharma/Documents/Data/test/UNET/Prediction data/BF/unet_pred/"+imgs_fls[i],pred[0,:,:,0]);
