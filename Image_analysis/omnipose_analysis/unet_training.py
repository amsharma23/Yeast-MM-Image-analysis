#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:30:08 2022

@author: amansharma
For training of UNET
Relying majorly on JunLab github training notebook and blog post by Vidushi Bhatia(https://medium.com/geekculture/u-net-implementation-from-scratch-using-tensorflow-b4342266e406)
Snippets used from blog post have been marked

"""

import random
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Conv2DTranspose
from tensorflow.keras.layers import concatenate
from tensorflow.keras.losses import binary_crossentropy
from sklearn.model_selection import train_test_split
import tifffile as tiff
import os
import cv2 as cv



#encoder function - blog snippet
def EncoderMiniBlock(inputs, n_filters=32, dropout_prob=0.0, max_pooling=True):
    """
    This block uses multiple convolution layers, max pool, relu activation to create an architecture for learning. 
    Dropout can be added for regularization to prevent overfitting. 
    The block returns the activation values for next layer along with a skip connection which will be used in the decoder
    """
    # Add 2 Conv Layers with relu activation and HeNormal initialization using TensorFlow 
    # Proper initialization prevents from the problem of exploding and vanishing gradients 
    # 'Same' padding will pad the input to conv layer such that the output has the same height and width (hence, is not reduced in size) 
    conv = Conv2D(n_filters, 
                  3,   # Kernel size   
                  activation='relu',
                  padding='same',
                  kernel_initializer='RandomUniform')(inputs)
    conv = Conv2D(n_filters, 
                  3,   # Kernel size
                  activation='relu',
                  padding='same',
                  kernel_initializer='RandomUniform')(conv)
    
    # Batch Normalization will normalize the output of the last layer based on the batch's mean and standard deviation
    conv = BatchNormalization()(conv, training=False)

    # In case of overfitting, dropout will regularize the loss and gradient computation to shrink the influence of weights on output
    if dropout_prob > 0:     
        conv = tf.keras.layers.Dropout(dropout_prob)(conv)

    # Pooling reduces the size of the image while keeping the number of channels same
    # Pooling has been kept as optional as the last encoder layer does not use pooling (hence, makes the encoder block flexible to use)
    # Below, Max pooling considers the maximum of the input slice for output computation and uses stride of 2 to traverse across input image
    if max_pooling:
        next_layer = tf.keras.layers.MaxPooling2D(pool_size = (2,2))(conv)    
    else:
        next_layer = conv

    # skip connection (without max pooling) will be input to the decoder layer to prevent information loss during transpose convolutions      
    skip_connection = conv
    
    return next_layer, skip_connection








#decoder function -- blog snippet
def DecoderMiniBlock(prev_layer_input, skip_layer_input, n_filters=32):
    """
    Decoder Block first uses transpose convolution to upscale the image to a bigger size and then,
    merges the result with skip layer results from encoder block
    Adding 2 convolutions with 'same' padding helps further increase the depth of the network for better predictions
    The function returns the decoded layer output
    """
    # Start with a transpose convolution layer to first increase the size of the image
    up = Conv2DTranspose(
                 n_filters,
                 (3,3),    # Kernel size
                 strides=(2,2),
                 padding='same')(prev_layer_input)

    # Merge the skip connection from previous block to prevent information loss
    merge = concatenate([up, skip_layer_input], axis=3)
    
    # Add 2 Conv Layers with relu activation and HeNormal initialization for further processing
    # The parameters for the function are similar to encoder
    conv = Conv2D(n_filters, 
                 3,     # Kernel size
                 activation='relu',
                 padding='same',
                 kernel_initializer='RandomUniform')(merge)
    conv = Conv2D(n_filters,
                 3,   # Kernel size
                 activation='relu',
                 padding='same',
                 kernel_initializer='RandomUniform')(conv)
    return conv


#combine the encoder-decoder blocks - blog snippet
def UNetCompiled(input_size=(128, 128,1), n_filters=32, n_classes=3):
   """
   Combine both encoder and decoder blocks according to the U-Net research paper
   Return the model as output 
   """
   # Input size represent the size of 1 image (the size used for pre-processing) 
   inputs = Input(input_size)
    
   # Encoder includes multiple convolutional mini blocks with different maxpooling, dropout and filter parameters
   # Observe that the filters are increasing as we go deeper into the network which will increasse the # channels of the image 
   cblock1 = EncoderMiniBlock(inputs, n_filters,dropout_prob=0.3, max_pooling=True)
   cblock2 = EncoderMiniBlock(cblock1[0],n_filters*2,dropout_prob=0.3, max_pooling=True)
   cblock3 = EncoderMiniBlock(cblock2[0], n_filters*4,dropout_prob=0.3, max_pooling=True)
   cblock4 = EncoderMiniBlock(cblock3[0], n_filters*8,dropout_prob=0.3, max_pooling=True)
   cblock5 = EncoderMiniBlock(cblock4[0], n_filters*16, dropout_prob=0.3, max_pooling=False) 
    
   # Decoder includes multiple mini blocks with decreasing number of filters
   # Observe the skip connections from the encoder are given as input to the decoder
   # Recall the 2nd output of encoder block was skip connection, hence cblockn[1] is used
   ublock6 = DecoderMiniBlock(cblock5[0], cblock4[1],  n_filters * 8)
   ublock7 = DecoderMiniBlock(ublock6, cblock3[1],  n_filters * 4)
   ublock8 = DecoderMiniBlock(ublock7, cblock2[1],  n_filters * 2)
   ublock9 = DecoderMiniBlock(ublock8, cblock1[1],  n_filters)

   # Complete the model with 1 3x3 convolution layer (Same as the prev Conv Layers)
   # Followed by a 1x1 Conv layer to get the image to the desired size. 
   # Observe the number of channels will be equal to number of output classes
   conv9 = Conv2D(n_filters,
                3,
                activation='relu',
                padding='same',
                kernel_initializer='RandomUniform')(ublock9)

   conv10 = Conv2D(n_classes, 1, padding='same')(conv9)
    
   # Define the model
   model = tf.keras.Model(inputs=inputs, outputs=conv10)

   return model




#making image mask dictionary - makes dictionary of images(noramlized and masks)
def make_img_msk_dic(img_path,mask_path):
	org_fls=[];
	msk_fls=[];
	ct=0;
	img_mask_dict = {};
	for st in os.listdir(img_path):
		if('.tif' in st): org_fls.append(st);
	org_fls.sort();

	for st in os.listdir(mask_path):
		if('.tif' in st): msk_fls.append(st);
	msk_fls.sort();

	for i in range(len(org_fls)):
		cur_img_path = os.path.join(org_fld1, org_fls[i])
		cur_mask_path = os.path.join(masks_fld1,msk_fls[i])
		with tiff.TiffFile(cur_img_path) as img_tiff:
			cur_img = img_tiff.asarray() / 256.
		with tiff.TiffFile(cur_mask_path) as mask_tiff:
			cur_mask = mask_tiff.asarray() / 256.
		img_mask_dict[ct] = [cur_img, cur_mask];
		img_mask_dict[ct] = [img_mask_dict[ct][0].reshape(128,128,1),img_mask_dict[ct][1].reshape(128,128,1)];
		ct=ct+1;

	return(img_mask_dict);


def VisualizeResults(index):
    img = X_valid[index]
    img = img[np.newaxis, ...]
    pred_y = unet.predict(img)
    pred_mask = tf.argmax(pred_y[0], axis=-1)
    pred_mask = pred_mask[..., tf.newaxis]
    fig, arr = plt.subplots(1, 3, figsize=(15, 15))
    print(pred_y.shape)
    #print(pred_mask[:,:,:]);
    arr[0].imshow(pred_y[0,:,:,0])
    arr[0].set_title('Processed Image')
    arr[1].imshow(pred_y[0,:,:,1])
    arr[1].set_title('Actual Masked Image ')
    arr[2].imshow(pred_y[0,:,:,2])
    arr[2].set_title('Predicted Masked Image ')
    fig.savefig("/Users/amansharma/Documents/Data/test/UNET/prediction.png");






#image loading with masks

org_fld1 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/BF/Trap1_aug/resize_and_rescaled_8bit";
masks_fld1 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/MasksBF/trap1_aug/resize_and_rescaled_8bit";

dict0 = make_img_msk_dic(org_fld1,masks_fld1);


aug_img_fld1 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/BF/Trap1_aug/augmented";
aug_masks_fld1 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/MasksBF/trap1_aug/augmented";
dict1 = make_img_msk_dic(aug_img_fld1,aug_masks_fld1);


org_fld2 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/BF/Trap2_aug/resize_and_rescaled_8bit";
masks_fld2 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/MasksBF/trap2_aug/resize_and_rescaled_8bit";

dict2 = make_img_msk_dic(org_fld2,masks_fld2);


aug_img_fld2 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/BF/Trap2_aug/augmented";
aug_masks_fld2 = "/Users/amansharma/Documents/Data/test/UNET/Training_data/MasksBF/trap2_aug/augmented";
dict3 = make_img_msk_dic(aug_img_fld2,aug_masks_fld2);


imgs=[];
masks=[];

for i in [dict0,dict1,dict2,dict3]:
	for j in i.keys():
		i[j][0] = i[j][0]
		imgs.append(i[j][0]);
		masks.append(i[j][1]);

imgs = np.array(imgs)
masks = np.array(masks)[...,0]

mask0 = masks == 0
mask1 = masks == 1
mask2 = masks == 2

masks = np.stack((mask0, mask1, mask2), axis = -1)

print(imgs.shape,masks.shape);


rd = random.randint(1,300);
X_train, X_valid, y_train, y_valid = train_test_split(imgs, masks, test_size=0.4,random_state=rd)
X_train = np.expand_dims(X_train, axis=3)
X_valid = np.expand_dims(X_valid, axis=3)

print(len(X_train),len(y_train));


unet = UNetCompiled(input_size=(128,128,1), n_filters=32, n_classes=3);

unet.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=.0001), 
             loss=tf.keras.losses.MeanSquaredError(),
              metrics=['accuracy']);
unet.summary()
print(y_train)
results = unet.fit(x=X_train, y=y_train, batch_size=8, epochs=50, validation_data=(X_valid, y_valid))

fig, axis = plt.subplots(1, 2, figsize=(20, 5))
axis[0].plot(results.history["loss"], color='r', label = 'train loss')
axis[0].plot(results.history["val_loss"], color='b', label = 'dev loss')
axis[0].set_title('Loss Comparison')
axis[0].legend()
axis[1].plot(results.history["accuracy"], color='r', label = 'train accuracy')
axis[1].plot(results.history["val_accuracy"], color='b', label = 'dev accuracy')
axis[1].set_title('Accuracy Comparison')
axis[1].legend()
fig.savefig("/Users/amansharma/Documents/Data/test/UNET/results.png");


print("Evaluate:");
unet.evaluate(X_valid,y_valid)

index = 0
VisualizeResults(0)
