# -*- coding: utf-8 -*-
"""
Created on Mon May  1 14:06:36 2023

@author: 703355681
"""

import numpy
import SimpleITK

from dicomtools.acr.misc import _create_circular_mask, _calc_piu

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def _get_max_area(stats):

    area = [stats.GetNumberOfPixels(l) for l in stats.GetLabels()]
    if len(area) > 1:
        print(f"Found {len(area)} labels")
    idx = area.index(max(area))
    return idx+1, area[idx]



# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def calc_piu(img):

    #TODO: ultimately, this code will be part of a large collection of functions
    #      that calculate ACR phantom test results. A few issues need to be
    #      considered:
    #       (1) error checking
    #       (2) handling visualziations
    #       (3) managing dependencies
    #       (4) handling small and medium ACR phantoms
    #
    #      For now, the code assumes that the user has supplied an ITK image
    #      object from the series reader. This object should have a 11 axial
    #      MR images, with slice 7 representing the uniformity slice of the ACR's
    #      large MR phantom.

    # Get the image from slice 7
    sl7 = img[::, ::, 6]

    # Find the fill using IsoDataThreshold filter. Provide slice 7. Output mask
    # uses 0 for values outside (background) and 1 for values inside. This mask
    # is used to get the image foreground area with the larges area (by pixels)
    imgbin = SimpleITK.IsoDataThreshold(sl7, 0, 1, img.GetSize()[0])
    cc = SimpleITK.ConnectedComponent(imgbin)
    stats = SimpleITK.LabelIntensityStatisticsImageFilter()
    stats.Execute(cc, img)
    lbl, area = _get_max_area(stats)

    # Calculate the phantom center from the bounding box. This is a more 
    # reliable method of determining the phantom's center than using, for 
    # example the GetCentroid or GetCenterOfGravity methods
    bb = stats.GetBoundingBox(lbl)
    center = (bb[0] + bb[2] / 2,
              bb[1] + bb[3] / 2)

    # Calculate the "Large ROI" mask
    dimPix = img.GetSpacing()[:2]
    if (len(set(dimPix)) != 1):
        raise NotImplementedError("Non-isotropic resolution unsupported.")
    rPhysRoi = numpy.sqrt(200 / numpy.pi)  # 200 cm^2 ROI see large phantom guidance
    rPixRoi = int(rPhysRoi / dimPix[0] * 10)  # convert dimPix to mm
    arrLargeRoi = _create_circular_mask(img.GetSize(), center, rPixRoi)

    # Calculate the maxima ROI size
    #TODO: validate these calculations. I've noticed that using 1.1 gives an
    #      area more representative of a 1 cm^2 ROI
    rPhysRoi = numpy.sqrt(1.1 / numpy.pi)  # 1 cm^2 ROI
    rPixRoi = int(rPhysRoi / dimPix[0] * 10)  # convert dimPix to mm

    # Create a new version of the phantom array to be used in calculating the
    # min/max signal intensities using the Mean image filter
    #TODO: find a more robust way to determine the bit depth of the image
    #TODO: 2 arrays are likely not necessary. Use numpy.NaN for background
    #      instead of the existing logic
    arr = SimpleITK.GetArrayFromImage(img)
    arrMin = arr.copy()
    arrMax = arr.copy()
    arrMin[arrLargeRoi != 1] = 2**16 - 1  # assumed 16-bit image
    arrMax[arrLargeRoi != 1] = 0

    # Apply the mean filter with a circular kernel
    imgMin = SimpleITK.Mean(SimpleITK.GetImageFromArray(arrMin), radius=3*(rPixRoi, ))
    arrMin = SimpleITK.GetArrayFromImage(imgMin)
    imgMax = SimpleITK.Mean(SimpleITK.GetImageFromArray(arrMax), radius=3*(rPixRoi, ))
    arrMax = SimpleITK.GetArrayFromImage(imgMax)

    # Use the min/max to calculate the PIU
    return _calc_piu(arrMin.min(), arrMax.max())