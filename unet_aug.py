#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sept  12 13:41:08 2022

@author: amansharma
For Data Augmentation to train UNET
"""

import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
import tifffile
import os
import cv2 as cv
import albumentations as alb
 
img_fld = "/Users/amansharma/Documents/Data/test/UNET/Training data/BF/Trap1_aug/resize_and_rescaled";
mask_fld = "/Users/amansharma/Documents/Data/test/UNET/Training data/MasksBF/trap1_aug/resize_and_rescaled";

img_fls = os.listdir(img_fld);
img_fls.sort();
mask_fls = os.listdir(mask_fld);
mask_fls.sort();

img_mask_dict_loc = {};

for i in range(len(img_fls)):
	img_mask_dict_loc[i] = [os.path.join(img_fld,img_fls[i]),os.path.join(mask_fld,mask_fls[i])];

transf = alb.Compose([alb.augmentations.geometric.transforms.Flip,
					  alb.augmentations.geometric.transforms.OpticalDistortion 	])