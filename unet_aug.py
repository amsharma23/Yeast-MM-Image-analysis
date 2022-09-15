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

transf = alb.Compose([alb.])
