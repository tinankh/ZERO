#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import imageio
import numpy as np


im = imageio.imread(sys.argv[1])
height, width = im.shape

output = np.zeros((height, width, 3))

for x in range(width):
    for y in range(height):
        if (im[x,y] == 255):
            output[x,y] = (0, 0, 0)
        else:
            output[x,y] = (int(1 + 3*im[x,y]),
                           int(1 + 2*im[x,y]),
                           int(1+ im[x,y]))


imageio.imwrite('colored_votes.png', output)
