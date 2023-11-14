#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:30:08 2022

@author: amansharma
For training of UNET
Relying majorly on JunLab github training notebook and blog post by Vidushi Bhatia(https://medium.com/geekculture/u-net-implementation-from-scratch-using-tensorflow-b4342266e406)

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
from tensorflow.keras import models
from tensorflow.keras import losses
from tensorflow.keras.losses import binary_crossentropy
from sklearn.model_selection import train_test_split
import tifffile as tiff
import os
import cv2 as cv

# define what happens at each layer
def conv_block(input_tensor, num_filters):
    encoder = layers.Conv2D(num_filters, (3, 3), padding='same')(input_tensor)
    encoder = layers.BatchNormalization()(encoder)
    encoder = layers.Activation('relu')(encoder)
    encoder = layers.Conv2D(num_filters, (3, 3), padding='same')(encoder)
    encoder = layers.BatchNormalization()(encoder)
    encoder = layers.Activation('relu')(encoder)
    return encoder

def encoder_block(input_tensor, num_filters):
    encoder = conv_block(input_tensor, num_filters)
    encoder_pool = layers.MaxPooling2D((2, 2), strides=(2, 2))(encoder)
    return encoder_pool, encoder

def decoder_block(input_tensor, concat_tensor, num_filters):
    decoder = layers.Conv2DTranspose(num_filters, (2, 2), strides=(2, 2), padding='same')(input_tensor)
    decoder = layers.concatenate([concat_tensor, decoder], axis=-1)
    decoder = layers.BatchNormalization()(decoder)
    decoder = layers.Activation('relu')(decoder)
    decoder = layers.Conv2D(num_filters, (3, 3), padding='same')(decoder)
    decoder = layers.BatchNormalization()(decoder)
    decoder = layers.Activation('relu')(decoder)
    decoder = layers.Conv2D(num_filters, (3, 3), padding='same')(decoder)
    decoder = layers.BatchNormalization()(decoder)
    decoder = layers.Activation('relu')(decoder)
    return decoder


#combine the encoder-decoder blocks - blog snippet
def UNetCompiled(target_size=(128, 128,1), n_filters=32, n_classes=3):
   """
   Combine both encoder and decoder blocks according to the U-Net research paper
   Return the model as output 
   """
   # make the layers
   inputs = layers.Input(shape=(target_size[0], target_size[1], 1))
   # 256
   encoder0_pool, encoder0 = encoder_block(inputs, n_filters)
   # 128
   encoder1_pool, encoder1 = encoder_block(encoder0_pool, n_filters*2)
   # 64
   encoder2_pool, encoder2 = encoder_block(encoder1_pool, n_filters*4)
   # 32
   encoder3_pool, encoder3 = encoder_block(encoder2_pool, n_filters*8)
   # 16
   center = conv_block(encoder3_pool, n_filters*16) # we were using 128 before
   # center
   # 32
   decoder3 = decoder_block(center, encoder3, n_filters*8)
   # 64
   decoder2 = decoder_block(decoder3, encoder2, n_filters*4)
   # 64
   decoder1 = decoder_block(decoder2, encoder1, n_filters*2)
   # 128
   decoder0 = decoder_block(decoder1, encoder0, n_filters)
   # 256
   outputs = layers.Conv2D(n_classes, (1, 1), activation='sigmoid')(decoder0)

   # make the model
   model = models.Model(inputs=[inputs], outputs=[outputs])

   return model

def dice_coeff(y_true, y_pred):
    smooth = 0.01 # originally 1. Make sure this is float. 
    score_factor = 2. # originally 2. Same, keep as float. 
    # Flatten
    y_true_f = tf.reshape(y_true, [-1]) # flattens tensor
    y_pred_f = tf.reshape(y_pred, [-1])
    intersection = tf.reduce_sum(y_true_f * y_pred_f) # sums the resulting product
    score = (score_factor * intersection + smooth) / (tf.reduce_sum(y_true_f) + tf.reduce_sum(y_pred_f) + smooth)
    return score

def dice_loss(y_true, y_pred):
    loss = 1 - dice_coeff(y_true, y_pred)
    return loss

def bce_dice_loss(y_true, y_pred):
    loss = losses.binary_crossentropy(y_true, y_pred) + dice_loss(y_true, y_pred)
    return loss


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
    print('Predictions:')

    img = X_valid[index]
    img = img[np.newaxis, ...]
    print(img.shape);
    pred_y = unet.predict(img)
    #print(pred_y.shape)
    pred_mask = tf.argmax(pred_y[0], axis=-1)
    pred_mask = pred_mask[..., tf.newaxis]
    fig, arr = plt.subplots(1, 3, figsize=(15, 15))
    print(pred_mask.shape)
    #print(pred_mask[:,:,:]);
    arr[0].imshow(X_valid[index][:,:,0,0])
    arr[0].set_title('Processed Image')
    arr[1].imshow(y_valid[index])
    arr[1].set_title('Actual Masked Image ')
    arr[2].imshow(pred_y[0,:,:])
    arr[2].set_title('Predicted Masked Image ')
    fig.savefig("/Users/amansharma/Documents/Data/test/UNET/prediction.png");
    tiff.imwrite("/Users/amansharma/Documents/Data/test/UNET/pred_img1.tif",pred_y[0,:,:]);        

    img = X_valid[index+6]
    img = img[np.newaxis, ...]
    pred_y = unet.predict(img)
    pred_mask = tf.argmax(pred_y[0], axis=-1)
    pred_mask = pred_mask[..., tf.newaxis]
    fig, arr = plt.subplots(1, 3, figsize=(15, 15))
    print(pred_mask.shape)
    #print(pred_mask[:,:,:]);
    arr[0].imshow(X_valid[index+6][:,:,0,0])
    arr[0].set_title('Processed Image')
    arr[1].imshow(y_valid[index+6])
    arr[1].set_title('Actual Masked Image ')
    arr[2].imshow(pred_y[0,:,:])
    arr[2].set_title('Predicted Masked Image ')
    fig.savefig("/Users/amansharma/Documents/Data/test/UNET/prediction1.png");
    tiff.imwrite("/Users/amansharma/Documents/Data/test/UNET/pred_img2.tif",pred_y[0,:,:]);






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
		imgs.append(i[j][0]);
		masks.append(i[j][1]);

imgs = np.array(imgs)
masks = np.array(masks)[...,0]

# mask0 = masks == 0
# mask1 = masks == 1
# mask2 = masks == 2

# masks = np.stack((mask0, mask1, mask2), axis = -1)

print(imgs.shape,masks.shape);


rd = random.randint(1,300);
X_train, X_valid, y_train, y_valid = train_test_split(imgs, masks, test_size=0.4,random_state=rd)
X_train = np.expand_dims(X_train, axis=3)
X_valid = np.expand_dims(X_valid, axis=3)

print(len(X_train),len(y_train));

save_model_path = "/Users/amansharma/Documents/Data/test/UNET/Model/2022_12_01_unet.hdf5";
#mdls_cb = tf.keras.callbacks.ModelCheckpoint(filepath=save_model_path,monitor='val_dice_loss',save_best_only=True,verbose=1);


unet = UNetCompiled(target_size=(128,128,1), n_filters=32, n_classes=1);

unet.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=.000 001), 
             loss=tf.keras.losses.MeanSquaredError(),
              metrics=['accuracy']);
unet.summary()

results = unet.fit(x=X_train, y=y_train, batch_size=8, epochs=50, validation_data=(X_valid, y_valid))

models.save_model(unet,save_model_path);

#dice = results.history['dice_loss']
#val_dice = results.history['val_dice_loss']

loss = results.history['loss']
val_loss = results.history['val_loss']

# epochs_range = range(epochs)
epochs_range = results.epoch

plt.figure(figsize=(8, 4))
plt.subplot(1, 2, 1)
plt.plot(epochs_range, dice, label='Training Dice Loss')
plt.plot(epochs_range, val_dice, label='Validation Dice Loss')
plt.legend(loc='upper right')
plt.title('Training and Validation Dice Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.xlim(0, None)
plt.ylim(0,1)

plt.subplot(1, 2, 2)
plt.plot(epochs_range, loss, label='Training Loss')
plt.plot(epochs_range, val_loss, label='Validation Loss')
plt.legend(loc='upper right')
plt.title('Training and Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.xlim(0, None)
plt.ylim(0,1)


plt.savefig("/Users/amansharma/Documents/Data/test/UNET/results.pdf", dpi=200)
#plt.show()
plt.close()
# fig.savefig("/Users/amansharma/Documents/Data/test/UNET/results.png");


print("Evaluate:");
unet.evaluate(X_valid,y_valid)

index = 0
VisualizeResults(0)





