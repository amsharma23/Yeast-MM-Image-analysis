from pathlib import Path
import os
from cellpose_omni import io as IO 
from cellpose_omni import models
import matplotlib as mlp 
import matplotlib.pyplot as plt 
import tifffile as tif
import numpy as np

test_img_fld = '/Users/amansharma/Documents/Data/test/cellpose/test/trap3/';
img_nms = [str(p) for p in Path(test_img_fld).rglob("*.tif")];

imgs = [IO.imread(f) for f in img_nms];

n = range(len(imgs));

model_loc = '/Users/amansharma/Documents/Data/test/cellpose/train/models/cellpose_residual_on_style_on_concatenation_off_omni_nclasses_4_train_2023_01_20_18_46_49.387237_correct';
model = models.CellposeModel(gpu=True, model_type=model_loc);

masks, flows, styles = model.eval([imgs[i] for i in n],channels=[0,0],rescale=None,mask_threshold=-1,
                                  transparency=True,flow_threshold=0,omni=True,
                                  cluster=True, resample=True,verbose=True);

for idx,i in enumerate(n):
    maski = masks[idx];
    mlp.image.imsave('/Users/amansharma/Documents/Data/test/cellpose/results/masks_res_w_aug_corrected/trap3/jpeg/'+str(idx)+'.jpg',maski);
    tif.imwrite('/Users/amansharma/Documents/Data/test/cellpose/results/masks_res_w_aug_corrected/trap3/tif/'+str(idx)+'.tif',maski)

