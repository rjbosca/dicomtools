# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:23:12 2019

@author: Ryan
"""

from matplotlib import pyplot, lines
import numpy
from pathlib import Path
from scipy import ndimage  # for line profiles, until I can figure out SimpleITK
import SimpleITK

# Script testing variables
dAcr = Path(r'D:\images\ACR\BIS_HITACHI\20171010\00004DC9')
fDcm = [f for f in dAcr.glob("*")]
isPlot = True
isDebug = False


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def line_profile(img, idx1, idx2, n=1024):

    dim = img.GetDimension()

    # Determine the distance between the two points
    coor1 = img.TransformContinuousIndexToPhysicalPoint(idx1)
    coor2 = img.TransformContinuousIndexToPhysicalPoint(idx2)
    dist = numpy.sqrt(sum([(coor1[i] - coor2[i])**2 for i in range(dim)]))

    # Generate the coordinates (similar to MATLAB's meshgrid)
    coords = [numpy.linspace(idx1[ii], idx2[ii], n)
              for ii in range(img.GetDimension())]

    # Sample the image
    mapped = ndimage.map_coordinates(SimpleITK.GetArrayFromImage(img),
                                     numpy.vstack(coords[::-1]),
                                     order=1,
                                     mode='nearest')

    # Return the "x" and "y" vectors for plotting
    return numpy.linspace(0.0, dist, n), mapped


# -----------------------------------------------------------------------------
# Overlay setup
# -----------------------------------------------------------------------------
cmap = pyplot.get_cmap('jet_r')
cmap.set_bad(color='black', alpha=0)


#--------------
# Series reader
#---
reader = SimpleITK.ImageSeriesReader()
fNames = reader.GetGDCMSeriesFileNames(str(dAcr))
reader.SetFileNames(fNames)
image = reader.Execute()

# -----------------------------------------------------------------------------
# Sort the slices
# -----------------------------------------------------------------------------
seNum = []
seUid = []
slLoc = []
for f in fDcm:
    img = SimpleITK.ReadImage(str(f))
    slLoc.append(float(img.GetMetaData('0020|1041')))
    seNum.append(img.GetMetaData('0020|0011'))
    seUid.append(img.GetMetaData('0020|000e'))
if len(set(seUid)) > 1:
    raise NotImplementedError

# Create a second array with the sorted slice locations
slLocSort = slLoc
slLocSort.sort()

# Display the first slice
img = SimpleITK.ReadImage(str(fDcm[slLoc.index(min(slLoc))]))
imgArr = SimpleITK.GetArrayFromImage(img).squeeze()

# Show the binary image
binImg = SimpleITK.BinaryThreshold(img)
binImg = SimpleITK.BinaryNot(binImg)
binArr = SimpleITK.GetArrayFromImage(binImg).squeeze().astype(float)
binArr[binArr == 0] = numpy.NaN  # needed for overlay
#if isPlot and isDebug:
#    f1 = pyplot.figure()
#    ax = f1.add_axes([0, 0, 1, 1])
#    ax.imshow(imgArr, cmap='gray')
#    ax.set_axis_off()
#    im = ax.imshow(binArr,
#                   aspect='equal',
#                   cmap=cmap,
#                   alpha=0.5)
#    pyplot.show()

# Calculate the connected components of the binary image (also labels)
cc = SimpleITK.ConnectedComponent(binImg)
stats = SimpleITK.LabelIntensityStatisticsImageFilter()

# Apply those connected components to the DICOM image
stats.Execute(cc, img)

pixArea = [stats.GetNumberOfPixels(l) for l in stats.GetLabels()]
idxMax = stats.GetLabels()[pixArea.index(max(pixArea))]


# -----------------------------------------------------------------------------
# Geometric accuracy
# -----------------------------------------------------------------------------

# Calculate the background signal median
bkgThresh = stats.GetMedian(idxMax)/2
dim = img.GetSize()
if isPlot and isDebug:
    idxh = 87
    idxv = 80
    xh, yh = line_profile(img, (0, idxh, 0), (255, idxh, 0))
    th = numpy.where(yh >= bkgThresh)[0]
    xv, yv = line_profile(img, (idxv, 0, 0), (idxv, 255, 0))
    tv = numpy.where(yv >= bkgThresh)[0]

    f2 = pyplot.figure()
    ax1 = pyplot.subplot(2, 2, 1)
    ax1.imshow(imgArr, cmap='gray')
    ax1.add_line(lines.Line2D([0, dim[0]], [idxh, idxh], color='r'))
    ax1.set_axis_off()
    ax2 = pyplot.subplot(2, 2, 2)
    ax2.plot(xh, yh)
    ax2.add_line(lines.Line2D([xh[0], xh[-1]],
                              [bkgThresh, bkgThresh],
                              color='r'))
    ax2.add_line(lines.Line2D([xh[th[0]], xh[th[0]]],
                              [0, bkgThresh],
                              color='r'))
    ax2.add_line(lines.Line2D([xh[th[-1]], xh[th[-1]]],
                              [0, bkgThresh],
                              color='r'))
    ax3 = pyplot.subplot(2, 2, 3)
    ax3.imshow(imgArr, cmap='gray')
    ax3.add_line(lines.Line2D([idxv, idxv],
                              [0, dim[1]],
                              color='r'))
    ax3.set_axis_off()
    ax4 = pyplot.subplot(2, 2, 4)
    ax4.plot(xv, yv)
    ax4.add_line(lines.Line2D([xv[0], xv[-1]],
                              [bkgThresh, bkgThresh],
                              color='r'))
    ax4.add_line(lines.Line2D([xv[tv[0]], xv[tv[0]]],
                              [0, bkgThresh],
                              color='r'))
    ax4.add_line(lines.Line2D([xv[tv[-1]], xv[tv[-1]]],
                              [0, bkgThresh],
                              color='r'))
    pyplot.show()

# Find the maximum left-right
disth = numpy.zeros((dim[1], ))
maxh = 0
for idx in range(dim[1]):

    # Coordinates on the line
    c1 = (0, idx, 0)
    c2 = (dim[0]-1, idx, 0)
    xh, yh = line_profile(img, c1, c2, n=dim[0])
    th = numpy.where(yh >= bkgThresh)[0]
    if (len(th) > 1):
        disth[idx] = xh[th[-1]]-xh[th[0]]
        # Store the maximum width. Note that this only works for the horizontal
        # width
        if disth[idx] >= maxh:
            maxh = disth[idx]
        else:
            disth[idx] = 0
            break

maxh = numpy.where(disth == max(disth))[0]
idxh = maxh[0]
if (maxh.size > 1):
    idxh = maxh[int(maxh.size / 2)]


distv = numpy.zeros((dim[0], ))
for idx in range(dim[0]):

    xv, yv = line_profile(img, (idx, 0, 0), (idx, 255, 0), n=10000)
    tv = numpy.where(yv >= bkgThresh)[0]
    if (len(tv) > 1):
        distv[idx] = xv[tv[-1]]-xv[tv[0]]

if isPlot:
    f3 = pyplot.figure()
    ax1 = pyplot.subplot(1, 2, 1)
    ax1.imshow(imgArr, cmap='gray')
    ax1.plot(disth, range(dim[1]), color='red')
    ax2 = pyplot.subplot(1, 2, 2)
    ax2.imshow(imgArr, cmap='gray')
    ax2.plot(range(dim[0]), distv, color='red')


#TODO: slice 1 needs
#      (1) calculate the signal median - use this to determine FWHM
#      (2) slice thickness accuracy
#      (3) HCR
#      (4) slice position accuracy
#      (5) geometric accuracy


# -----------------------------------------------------------------------------
# PIU
# -----------------------------------------------------------------------------



if isPlot:
    f4 = pyplot.figure()
    ax1 = pyplot.subplot(1, 1, 1)
    ax1.imshow(imgArr, cmap='gray')