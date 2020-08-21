#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import imageio
import numpy as np

import matplotlib.pyplot as plt

color = 'hsv'
cmap = plt.get_cmap(color)


im = imageio.imread(sys.argv[1])
height, width = im.shape

# testim = np.zeros((height, width, 1))
# for n in range(1,64):
#     for i in range((n-1)*8,8*n):
#         for j in range(width):
#             testim[i,j] = n

# imageio.imwrite('testimage.tiff', testim)

output = np.zeros((height, width, 3))

for x in range(height):
    for y in range(width):
        if (im[x,y] == 255):
            output[x,y] = np.array((0, 0, 0))
        else:
            rgb = 255*np.array(cmap((im[x,y])/64)[:3])
            output[x,y] = rgb.astype(np.uint8)

imageio.imwrite('colored_votes_' + color + '.png', output)
