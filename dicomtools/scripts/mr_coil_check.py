# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:53:50 2019

@author: 703355681
"""

import matplotlib
from matplotlib import pyplot, lines
from dicomtools import vendormaps
from pydicom.errors import InvalidDicomError
import numpy
import pandas
import pydicom
import argparse
import SimpleITK

# Allows display from command prompt
try:
    matplotlib.use("TkAgg")
except ImportError:
    pass

# Script setup
isIdeDebug = True

coil_check = {'Series Number': [],
              'Coil Description': [],
              'Mean Signal': [],
              'Max Signal': [],
              'Min Signal': [],
              'Bkg Signal': [],
              'Bkg Std': [],
              'Freq FOV': [],
              'Ghost Sig': [],
              'SNR': [],
              'PIU': [],
              'PSG': [],
              'Mean Area': [],
              'Bkg Area': [],
              'Bkg Signal (L)': [],
              'Bkg Signal (R)': [],
              'Bkg Signal (T)': [],
              'Bkg Signal (B)': [],
              'Slice Loc': [],
              'Phase FOV': [],
              'Freq Enc Steps': [],
              'Phase Enc Steps': [],
              'Phase Enc Dir': [],
              'File': [],}


# -----------------------------------------------------------------------------
# Overlay setup
# -----------------------------------------------------------------------------
cmap = pyplot.get_cmap('jet_r')
cmap.set_bad(color='black', alpha=0)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def _calc_piu(minSig, maxSig):

    return 100.*(1 - (maxSig-minSig)/(maxSig+minSig))


def _finalize(args):

    # Write the file
    if args.out:
        df = pandas.DataFrame(coil_check)
        df.to_excel(str(args.out), index=False)
    print(coil_check)

    if args.show and not isIdeDebug:
        pyplot.show(block=True)
    elif args.show:
        pyplot.show()


def _filter_area_thresh(img, filt, cond, r=1, nIter=100, imgStats=None):

    # Set the filter radius and perform the first filter application
    filt.SetKernelRadius((r, )*3)
    imgFilt = filt.Execute(img)

    # Apply the filter at most nIter-1 times
    for idx in range(nIter-1):

        # Calc area
        cc = SimpleITK.ConnectedComponent(imgFilt)
        stats = SimpleITK.LabelIntensityStatisticsImageFilter()
        if (imgStats is not None):
            stats.Execute(cc, imgStats)
        else:
            stats.Execute(cc, imgFilt)
        lMax, area = _get_max_area(stats)

        # Apply stop condition
        if cond(area):
            break

        imgFilt = filt.Execute(imgFilt)

    print(f"Filter completed in n={idx} iterations")

    return imgFilt, stats, lMax


def _get_max_area(stats):

    area = [stats.GetNumberOfPixels(l) for l in stats.GetLabels()]
    if len(area) > 1:
        print(f"Found {len(area)} labels")
    idx = area.index(max(area))
    return idx+1, area[idx]


def calc_coil_metrics(hdr, img, args):

    # Get the image array
    arrImg = SimpleITK.GetArrayFromImage(img).squeeze().astype(float)

    # -------------------------------------------------------------------------
    # Find the background fill
    # -------------------------------------------------------------------------

    # Find the fill and calculate the area. Currently, using the image max and
    # 1/10th of the max as the upper and lower thresholds, respectively
    # TODO: this lower threshold works well for reasonable SNR (>100). How to
    #       adjust (check) if the SNR is low???
    # TODO: offer an iterative process to futher refine the fill signal
    #       detection
    imgBin = SimpleITK.IsoDataThreshold(img, 0, 1, img.GetSize()[0])
    ccOrig = SimpleITK.ConnectedComponent(imgBin)
    stats = SimpleITK.LabelIntensityStatisticsImageFilter()
    stats.Execute(ccOrig, img)
    lIdx, area = _get_max_area(stats)

    # Hueristic check for bad images. The idea: the detected area must have a
    # signal larger than a circumscribed square with half the width of the FOV.
    # Remove the file name from the list of files since no data will be stored.
    if area < (min(arrImg.shape)/4)**2:
        print("Unable to determine valid phantom signal...")
        idxPop = len(coil_check['Coil Description'])
        coil_check['File'].pop(idxPop)
        return

    # Initializes the binary label array and shows the image
    if args.show:
        arrLbl = SimpleITK.GetArrayFromImage(imgBin).squeeze()
        arrLbl = arrLbl.astype(float)
        pyplot.figure()
        axs = [pyplot.subplot(1, 2, 1)]
        axs[0].imshow(arrImg, cmap='gray')
        axs[0].set_axis_off()

    # -------------------------------------------------------------------------
    # Calculate signal mean
    # -------------------------------------------------------------------------

    # Get some information about the acquisition
    coil_check['Coil Description'].append(vendormaps.ReceiveCoil(h))
    coil_check['Slice Loc'].append(float(h.SliceLocation))
    coil_check['Freq FOV'].append(vendormaps.FieldOfView(h)[0])
    coil_check['Phase FOV'].append(vendormaps.FieldOfView(h)[1])
    coil_check['Phase Enc Dir'].append(vendormaps.PhaseEncodingDir(h))
    if (coil_check['Phase Enc Dir'][-1] == 'ROW'):
        coil_check['Freq Enc Steps'].append(h.AcquisitionMatrix[1])
        coil_check['Phase Enc Steps'].append(h.AcquisitionMatrix[2])
    else:
        coil_check['Freq Enc Steps'].append(h.AcquisitionMatrix[0])
        coil_check['Phase Enc Steps'].append(h.AcquisitionMatrix[3])

    # Erode the base binary image until 90% of the area is represented
    imgErd, statsErd, lbl = \
        _filter_area_thresh(imgBin,
                            SimpleITK.BinaryErodeImageFilter(),
                            lambda x: x/area <= 0.9,
                            imgStats=img)
    arrErd = SimpleITK.GetArrayFromImage(imgErd).squeeze().astype(bool)

    # Get the bounding box (this will be used to artifacts at the edge of the
    # FOV). bounding box = (x_min, y_min, z_min, x_width, y_width, z_width)
    bbErd = statsErd.GetBoundingBox(lbl)
    coorErd = [bbErd[0]-1, bbErd[0]+bbErd[3],
               bbErd[1]-1, bbErd[1]+bbErd[4]]

    # Update the binary label array
    if args.show:
        arrLbl[arrErd] += 1

        # Add the inner bounding box to the image
        axs[0].add_line(lines.Line2D(
            [0, img.GetSize()[0]], 2*[coorErd[2]], color='r'))  # top
        axs[0].add_line(lines.Line2D(
            [0, img.GetSize()[0]], 2*[coorErd[3]], color='r'))  # btm
        axs[0].add_line(lines.Line2D(
            2*[coorErd[0]], [0, img.GetSize()[1]], color='r'))  # left
        axs[0].add_line(lines.Line2D(
            2*[coorErd[1]], [0, img.GetSize()[1]], color='r'))  # right

    coil_check['Mean Signal'].append(statsErd.GetMean(lbl))
    coil_check['Mean Area'].append(stats.GetNumberOfPixels(lbl))

    # -------------------------------------------------------------------------
    # Calculate the background statistics
    # -------------------------------------------------------------------------

    # Dialate the base binary image until 120% of the area is represented
    imgDia, statsDia, lbl = \
        _filter_area_thresh(imgBin,
                            SimpleITK.BinaryDilateImageFilter(),
                            lambda x: x/area >= 1.2)
    arrDia = SimpleITK.GetArrayFromImage(imgDia).squeeze().astype(bool)

    # Update the binary label array with the dialated overlay
    if args.show:
        arrLbl[arrDia] += 1

    # Use the bounding box to determine the
    bbDia = statsDia.GetBoundingBox(lbl)
    coorDia = [bbDia[0]-1, bbDia[0]+bbDia[3],
               bbDia[1]-1, bbDia[1]+bbDia[4]]
    arr = numpy.zeros_like(arrDia, dtype=int)

    # Calculate the background signals (note that hte labels are important)
    if (coil_check['Phase Enc Dir'][-1] == 'ROW'):
        # Frequency portion
        arr[:coorDia[2], :] = 1  # top
        arr[coorDia[3]:, :] = 2  # bottom
        # Phase (i.e., ghosting) portion
        arr[coorErd[2]:coorErd[3], :coorDia[0]] = 3  # left
        arr[coorErd[2]:coorErd[3], coorDia[1]:] = 4  # right
        bkgSig = [arrImg[arr == 3].mean(),
                  arrImg[arr == 4].mean(),
                  arrImg[arr == 1].mean(),
                  arrImg[arr == 2].mean()]
    elif (coil_check['Phase Enc Dir'][-1] == 'COL'):
        # Frequency portion
        arr[:, :coorDia[0]] = 1  # left
        arr[:, coorDia[1]:] = 2  # right
        # Phase (i.e., ghosting) portion
        arr[:coorDia[2], coorErd[0]:coorErd[1]] = 3  # top
        arr[coorDia[3]:, coorErd[0]:coorErd[1]] = 4  # bottom
        bkgSig = [arrImg[arr == 1].mean(),
                  arrImg[arr == 2].mean(),
                  arrImg[arr == 3].mean(),
                  arrImg[arr == 4].mean()]
    else:
        raise NotImplementedError()
    coil_check['Bkg Signal (L)'].append(bkgSig[0])
    coil_check['Bkg Signal (R)'].append(bkgSig[1])
    coil_check['Bkg Signal (T)'].append(bkgSig[2])
    coil_check['Bkg Signal (B)'].append(bkgSig[3])

    # Update the binary label array with the backgound mask
    if args.show:

        # Update the background labels
        for i in range(1, 5):
            arrLbl[arr == i] = arrLbl.max() + 1

        # Add the bounding box lines on the image
        axs[0].add_line(lines.Line2D(
            [bbDia[0]-1, bbDia[0]+bbDia[3]], 2*[bbDia[1]-1]))  # top
        axs[0].add_line(lines.Line2D(
            [bbDia[0]-1, bbDia[0]+bbDia[3]], 2*[bbDia[1]+bbDia[4]]))  # btm
        axs[0].add_line(lines.Line2D(
            2*[bbDia[0]-1], [bbDia[1]-1, bbDia[1]+bbDia[4]]))  # left
        axs[0].add_line(lines.Line2D(
            2*[bbDia[0]+bbDia[3]], [bbDia[1]-1, bbDia[1]+bbDia[4]]))  # right

    # Use the label array to calculate some stats
    coil_check['Bkg Signal'].append(
                           arrImg[numpy.logical_or(arr == 1, arr == 2)].mean())
    coil_check['Bkg Std'].append(arrImg[arr != 0].std())
    coil_check['Bkg Area'].append(numpy.count_nonzero(arr))
    coil_check['Ghost Sig'].append(
                           arrImg[numpy.logical_or(arr == 3, arr == 4)].mean())

    # Calculate some metrics based on the background and mean signals
    coil_check['SNR'].append(coil_check['Mean Signal'][-1] /
                             coil_check['Bkg Std'][-1])
    coil_check['PSG'].append(100.*abs(sum(bkgSig[:2])-sum(bkgSig[-2:])) /
                             (2.*coil_check['Mean Signal'][-1]))

    # -------------------------------------------------------------------------
    # Calculate the ghosting signal
    # -------------------------------------------------------------------------

    # Get all background values outside of the dialated mask (this will help
    # avoid Gibbs ringing)
#    bbOrig = stats.GetBoundingBox(lIdx)
#    arr = SimpleITK.GetArrayFromImage(imgDia).squeeze().astype(bool)
#    if (coil_check['Phase Enc Dir'][-1] == "ROW"):
#        arr[:bbOrig[1]-1, :] = True
#        arr[sum(bbOrig[1:-1:3])+1:, :] = True
#    elif (coil_check['Phase Enc Dir'][-1] == "COL"):
#        arr[:, :bbOrig[0]-1] = True
#        arr[:, sum(bbOrig[0:-1:3])+1:] = True
#    else:
#        raise NotImplementedError()
#
    if args.show:
        arrLbl[arrLbl < 1] = numpy.NaN
        pyplot.imshow(arrLbl, cmap=cmap, alpha=0.5)
        # TODO: develop methods to detect/avoid RF leakage and
        #       truncation bands
        axs.append(pyplot.subplot(1, 2, 2))
        axs[1].hist(arrImg[numpy.invert(arrDia)], bins=150)

    # -------------------------------------------------------------------------
    # Calculate the min/max signal within the mean phantom signal
    # -------------------------------------------------------------------------

    # Calculate the ROI radius needed for the min/max signal calculation
    dimPix = img.GetSpacing()[:2]
    if (len(set(dimPix)) != 1):
        raise NotImplementedError("Non-isotrpoic resolution unsupported.")
    areaPhys = dimPix[0]**2 * numpy.array(img.GetSize()).prod()
    areaRoi = 0.15 * areaPhys / 100  # see 2015 ACR QC (pg 101)
    rPhysRoi = numpy.sqrt(areaRoi / numpy.pi)
    rPixRoi = int(rPhysRoi / dimPix[0])

    # Create a new version of the original image such that all areas outside
    # the eroded binary mask are a large negative number. This ensures that the
    # circular ROI used to calculate the average will not extend beyond the
    # mean signal ROI
    arrMax = arrImg.copy()
    arrMax[arrErd != 1] = -10 * arrImg.max()
    imgMax = SimpleITK.GetImageFromArray(arrMax)

    # Filter the image
    imgMax = SimpleITK.Mean(imgMax, radius=3*(rPixRoi, ))
    arrMax = SimpleITK.GetArrayFromImage(imgMax).squeeze()
    arrMax[arrMax < 0] = 0

    # Create a new version of the original image such that all areas outside
    # the eroded binary mask are a large negative number. This ensures that the
    # circular ROI used to calculate the average will not extend beyond the
    # mean signal ROI
    arrMin = arrImg.copy()
    arrMin[arrErd != 1] = 10**6 * arrImg.max()
    imgMin = SimpleITK.GetImageFromArray(arrMin)

    # Filter the image
    imgMin = SimpleITK.Mean(imgMin, radius=3*(rPixRoi, ))
    arrMin = SimpleITK.GetArrayFromImage(imgMin).squeeze()
    arrMin[arrMin > arrImg.max()] = 0

    coil_check['Max Signal'].append(arrMax.max())
    yMax, xMax = numpy.where(arrMax == coil_check['Max Signal'][-1])
    coil_check['Min Signal'].append(arrMin[arrMin != 0].min())
    yMin, xMin = numpy.where(arrMin == coil_check['Min Signal'][-1])

    # Calculate the PIU
    coil_check['PIU'].append(_calc_piu(coil_check['Min Signal'][-1],
                                       coil_check['Max Signal'][-1]))

    if args.show:
        # TODO: verify that the min/max ROIs are within the mean signal ROI

        # Show the ROIs on the eroded image
        axs[0].scatter(xMin, yMin, color='b')
        axs[0].scatter(xMax, yMax, color='r')


def _create_parser():

    parser = \
        argparse.ArgumentParser(prog=__file__,
                                description='Calculate MRI coil performance quality metrics',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', metavar='F', type=str, nargs='+',
                        help="Directory or file to analyze.")
    parser.add_argument('--out', help="Output Excel file")
    parser.add_argument('--show', action='store_true',
                        help=("When specified, verification images will "
                              "be shown."))

    return parser


def _parse_args(parser, args):

    from pathlib import Path

    # Validate the files
    files = []
    isWarned = False
    for f in args.files:

        f = Path(f)
        if f.is_dir():
            [files.append(ff) for ff in f.rglob('*') if ff.is_file()]
        elif f.is_file():
            files.append(f)
        elif (len(args.files) == 1):
            raise ImportError(f"Invalid file or directory: {f}")
        elif not isWarned:
            isWarned = True
            print("Ignoring invalid files or directories...")

    if not len(files):
        if (len(args.files) == 1):
            s = f'\t {args.files[0]}'
        else:
            s = ''.join([f"\t{ss}\n" for ss in s][0:-1])
        sErr = f"Unable to find any valid files or directories in:{s}"
        raise ImportError(sErr)

    args.files = files

    # Validate the ouptut file
    if args.out:
        fOut = Path(args.out)
        if '.xls' not in fOut.ext:
            raise ImportError()
        elif not fOut.parts()[-2].isdir():
            raise ImportError()

    return args


if __name__ == "__main__":

    parser = _create_parser()

    # For testing - comment otherwise
    if isIdeDebug:
        s = [r"C:\Users\703355681\Desktop\BIS_MR\Foot-Ankle",
             r"--show"]
        args = parser.parse_args(s)
    else:
        args = parser.parse_args()
    args = _parse_args(parser, args)

    # Read all  files so error checking can be performed along the way
    hdrs = []
    imgs = []
    isWarned = False
    for f in args.files:
        try:
            hdrs.append(pydicom.read_file(str(f)))
            imgs.append(SimpleITK.ReadImage(str(f)))
            coil_check['File'].append(str(f))
        except InvalidDicomError:
            if not isWarned:
                isWarned = True
                print("Ignoring invalid DICOM files...")

    for h, i in zip(hdrs, imgs):
        arrLbl = calc_coil_metrics(h, i, args)

    _finalize(args)
    # pyplot.figure()
    # pyplot.imshow(arrLbl)
