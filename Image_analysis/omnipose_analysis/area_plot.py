#%% libraries
from tifffile import imwrite
from tifffile import imread
from PIL import Image
#from all_funcs import crop_using_roi_tuple
#from all_funcs import extract_roi_coordinates

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from os import walk
from os import listdir
from os.path import isfile, join
from pathlib import Path

#from read_roi import read_roi_file
#from read_roi import read_roi_zip

#import trap_ex_func as trapf
#import cell_size_func as cellf

import shutil

#%% code

loc_1 = '/Users/amansharma/Documents/Data/test/cellpose/results/masks_res_w_aug_corrected/trap3/tif/';
img_nms = [str(p) for p in Path(loc_1).rglob('*.tif')];

imgs = [imread(nm) for nm in img_nms];
time = range(len(imgs));
ar =[]
for t in time:
	img = np.array(imgs[t]);

	ar.append(np.count_nonzero(img==4)+np.count_nonzero(img==1)); #from lookin at the JPEG images can tell most of the time the mother segmentation number is 4


fog, ax = plt.subplots();

ax.scatter(time, ar);
ax.plot(time,ar);
plt.show();
