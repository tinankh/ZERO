#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import imageio
import numpy as np

import matplotlib.pyplot as plt
cmap = plt.get_cmap('Accent')


im = imageio.imread(sys.argv[1])
height, width = im.shape

output = np.zeros((height, width, 3))

for x in range(height):
    for y in range(width):
        if (im[x,y] == 255):
            output[x,y] = np.array((0, 0, 0))
        else:
            rgb = 255*np.array(cmap((im[x,y]+1)/64)[:3])
            output[x,y] = rgb.astype(np.uint8)
            # output[x,y] = (int(1 + 3*im[x,y]),
            #                int(1 + 2*im[x,y]),
            #                int(1+ im[x,y]))

imageio.imwrite('colored_votes.png', output)
