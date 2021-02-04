#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:20:15 2021

@author: aaron
"""

import matplotlib.pyplot as plt
import numpy as np
import pyshtools as pysh

# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.75)

degrees = np.arange(101, dtype=float)
degrees[0] = np.inf
power = degrees**(-2)

clm = pysh.SHCoeffs.from_random(power, seed=12345)

# fig, ax = clm.plot_spectrum(show=False)

#fig, ax = clm.plot_spectrum2d(cmap_rlimits=(1.e-7, 0.1), show=False)

grid = clm.expand()
fig, ax = grid.plot(show=False)

