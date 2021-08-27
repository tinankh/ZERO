#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import imageio
import numpy as np

mask1 = imageio.imread(sys.argv[1]) # forgery mask F
mask2 = imageio.imread(sys.argv[2]) # forgery mask M
img = imageio.imread(sys.argv[3]).astype('uint8') # luminance image
n,m = mask2.shape

if (np.sum(mask1) + np.sum(mask2) > 0):

    # The simulated version is smaller than the tested image
    diffn = mask1.shape[0] - mask2.shape[0]
    diffm = mask1.shape[1] - mask2.shape[1]

    # The simulated version is resized so the masks can be overlayed
    mask2_resized = np.append(np.zeros((n,diffm)), mask2, axis=1)
    mask2_resized = np.append(np.zeros((diffn,m+diffm)), mask2_resized, axis=0)
    mask2_resized[np.where(mask1 > 0)] = 0 # forgery mask F is chosen over forgery mask M

    # forgeries via foreign grid are red and via missing grid are blue
    colored = np.dstack((mask1,np.zeros((mask1.shape)),mask2_resized)).astype(np.uint8)

    # masks are overlayed over the luminance image
    img = (0.2*img).astype('uint8')
    result = np.dstack((img, img, img)).astype(np.uint8)
    result[np.where(colored>0)] = colored[np.where(colored>0)]

    imageio.imwrite('result_zero.png', result)

else:
    imageio.imwrite('result_zero.png', img)
