#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to plot output results in the interface. 

@author: maxime - maxime.bernard@uni-potsdam.de

"""

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib
import numpy as np


#-----------------------------------------------------------------------------
#-------------------------------- Functions ----------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def get_gradient_color(c1,c2,mix=0):
    """create a gradient color from color 1 to color 2. """
    c1 = np.array(matplotlib.colors.to_rgb(c1))
    c2 = np.array(matplotlib.colors.to_rgb(c2))
    return matplotlib.colors.to_hex((1-mix)*c1 + mix*c2)


def make_error_boxes(ax,xdata,ydata,xerror,yerror,facecolor='g',edgecolor='none',alpha=0.5,label='Portion extracted'):
    """ plot boxes for 4He/3He (see output.py)"""
   
    x_modif = [xdata.values[i-1]  if i > 0 else 0 for i in range(len(xdata.values))]
    Width = [xdata.values[i] - xdata.values[i-1] + xerror.values[i]  if i > 0 else xdata.values[i] + xerror.values[i] for i in range(len(xdata.values))]
    # Loop over data points, creat box from errors at each point
    errorboxes = [Rectangle((x, y - ye),W, 2*ye)
                  for x, y, W, ye in zip(x_modif, ydata.values, Width, yerror.values)]
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor,label=label)
    
    # Add collection to axes
    ax.add_collection(pc)