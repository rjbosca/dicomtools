# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 22:21:51 2017

@author: Ryan Bosca
"""


def ipp2plane(hdr, ang=10.) -> str:
    """Returns the slice/slab orientation string

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set
        ang (float): acceptance angle in degrees (default: 10)

    Returns:
        str:Orientation - one of Axial, Coronal, Sagittal, Oblique

    This function determines the closest possible plane based on the direction
    cosines. If the acquisition volume is within the acceptnace angle of any
    particular plane (Axial, Coronal, or Sagittal) that plane name is returned.
    """

    import numpy as np

    # Check the modality SOP class UID. Currently, the following storage
    # classes are supported (in the same order as if/elif)
    #   - MR Image Storage
    #   - Enhanced MR Image Storage
    sop_class_uid = hdr[0x0008, 0x0016].value
    if (sop_class_uid == "1.2.840.10008.5.1.4.1.1.4"):
        # Image orientation patient
        iop = hdr[0x0020, 0x0037].value
    elif (sop_class_uid == "1.2.840.10008.5.1.4.1.1.4.1"):
        # The Plane Orientation Sequence is stored in the Functional 
        # Groups Sequence.
        pfg = hdr[0x5200, 0x9230][0][0x0020, 0x9116][0]
        iop = np.array(pfg[0x0020, 0x0037].value, dtype=float)
    else:
        raise NotImplementedError(
            "Expected modality MR Iage Storage or Enhanced MR Image Storage")

    #TODO: consider making an "oblique-axial", "oblique-coronal", etc.
    #TODO: how large of an angle can be specified? In other words, at what 
    #      point does the angle correspond to a different plane?

    # The following dictionary is used to encode the image orientation
    dictOrient = {0: "Sagittal",
                  1: "Coronal",
                  2: "Axial",
                  3: "Oblique"}

    # Calculate the possible orientations. Calculating the quadrature of
    # the direction cosines provides a good estimate of the orientation
    d = np.vstack((iop - np.array([0,1,0,0,0,-1]),  # sagittal
                   iop - np.array([1,0,0,0,0,-1]),  # coronal
                   iop - np.array([1,0,0,0,1, 0])))  # axial
    d = np.sum(d*d, axis=1)

    # Set the output to oblique for angles greater than 45 deg.
    ind = np.argmin(d)
    if (d.min() > np.cos(ang)**2):
        ind = 3

    # Set the return value
    return dictOrient[ind]