# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:51:28 2023

@author: Ryan Bosca
"""

import numpy

def _calc_piu(mins, maxs) -> float:
    """Calculate the percent integral uniformity

    Parameters:
        mins (float):minimum signal intensity
        maxs (float):maximum signal intensity

    Returns:
        float:Percent integral uniformity

    """

    return 100.*(1 - (maxs-mins)/(maxs+mins))


def _create_circular_mask(dim, c=None, r=None):
    """Createa a circular bitmask 

    Parameters:
      dim (list or tuple):2D image dimensions (height, width)
      c (list or tuple):circle center
      r (float):circle radius

    Returns:
        numpy.ndarry:Bitmask with circle of radius r at center c
    """

    h = dim[0]
    w = dim[1]

    # Calculate default center/radius based on input height/width
    if center is None:
        center = (int(w/2), int(h/2))
    if radius is None:
        radius = min(center[0], center[1], w-center[0], h-center[1])
    
    y, x = numpy.ogrid[:h, :w]
    dist = numpy.sqrt((x - center[0])**2 + (y - center[1])**2)

    # Set all values within the circle to 1 and those outside to 0
    dist[dist <= radius] = 1
    dist[dist > 1] = 0

    return bool(dist)