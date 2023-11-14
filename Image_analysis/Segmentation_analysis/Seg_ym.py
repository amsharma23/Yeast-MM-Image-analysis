#segementation with YeastMate python backend with GPU implementation

from tifffile import imwrite
from tifffile import imread
from PIL import Image
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import subprocess
from time import perf_counter_ns
from os import walk
from os import listdir
from os.path import isfile, join

from yeastmatedetector.inference import YeastMatePredictor

path = "/Users/amansharma/Documents/Data/test/Speed_test/FOV2/Cropped_ts/";
detections = [];
masks = [];
tstrt = perf_counter_ns();
for imgs in os.listdir(path):
	if(imgs!='.DS_Store'):

		img = imread("/Users/amansharma/Documents/Data/test/Speed_test/FOV2/FOV2_ts/1.tif")
		predictor = YeastMatePredictor('/Users/amansharma/yeastmate/models/yeastmate.yaml','/Users/amansharma/yeastmate/models/yeastmate_weights.pth');

		detection, mask = predictor.inference(img);
		detections.append(detection);
		masks.append(mask);

tstop = perf_counter_ns();
print((tstop-tstrt)/(10**9));
