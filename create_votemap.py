#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import imageio
import numpy as np

import matplotlib.pyplot as plt

im = imageio.imread(sys.argv[1])

cmap1 = plt.get_cmap('tab20')
cmap2 = plt.get_cmap('tab20b')
cmap3 = plt.get_cmap('tab20c')
cmap4 = plt.get_cmap('Set3')

def colormap(v):
    v = v.copy()
    # swap 0 and 4 so that 0 is colored with green
    v0 = v == 0
    v4 = v == 4
    v[v0] = 4
    v[v4] = 0
    v2 = cmap1(v/20)
    v2[v >= 20] = cmap2((v[v >= 20] - 20)/20)
    v2[v >= 40] = cmap3((v[v >= 40] - 40)/20)
    v2[v >= 60] = cmap4((v[v >= 60] - 60)/20)
    return v2

output = 255 * colormap(im)[...,:3]
output[im == 255] = 0

output = output.astype(np.uint8)
imageio.imwrite('colored_'+ str(sys.argv[1]), output)
