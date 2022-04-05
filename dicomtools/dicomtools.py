# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 22:21:51 2017

@author: Ryan Bosca
"""

from typing import Union

from pathlib import Path

import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.misc import is_dicom


def get_dir_from_dicom(hdr: Union[str, Dataset, FileDataset]) -> str:
    """Converts DICOM meta-data to a valid directory

    Parameters
    ----------
    header : Path object, str, or dicomtools.header

    Returns
    -------
    dir : str
        Directory name generated from DICOM meta-data. Note that a
        directory name of "unknown" is generated when the requisite DICOM
        tags are not present.

    """

    # Determine what type of data the user provided, coverting when needed
    if (type(hdr) == str):
        hdr = Path(hdr)
        hdr = pydicom.dcmread(hdr)
    if (type(hdr) != Dataset) or (type(hdr) != FileDataset):
        raise NotImplementedError(f"Inavlid input of type: {type(hdr)}...")

    # Ensure that all fields necessary for converting the file information
    # to a valid directory name are present
    # TODO: in a future version of this function, there should be a way of
    #       specifying what data to use
    sNum = hdr[0x0020, 0x0011].value
    sDesc = hdr[0x0008, 0x0103e].value
    dOut = 'unknown'
    if sNum and sDesc:
        # Remove special characters
        sNum = ''.join([rc if rc.isalnum() else '_' for rc in str(sNum)])
        sDesc = ''.join([rc if rc.isalnum() else '_' for rc in sDesc])

        dOut = '--'.join([sNum, sDesc])

    return dOut


def mkdir_dicom(
    dicomDir: Union[Path, str], 
    useAllFiles: bool = False
    ) -> bool:
    """Rename a directory of DICOM files based on the meta-data

    Parameters
    ----------
    dicomDir : str or Path object
        Absolute to path of the directory to be renamed
    useAllFiles: bool
        When True, every file in the specified directory will be used to
        determine the target directory. The file will be moved to the
        target assuming a file of the same name does not exist at that
        location. By default, only the first DICOM file found is used to
        determine the target directory

    Returns
    -------
    success : bool
        True when the rename operation succeeds; otherwise False

    """

    # Force Path object from string
    dicomDir = Path(dicomDir)
    if not dicomDir.is_dir():
        raise ValueError("dicomDir must be a directory not a file...")

    # Check each file for DICOM data
    isSuccess = False
    for f in dicomDir.rglob("*"):

        # Only DICOM files are considered...
        if not f.is_file() or not is_dicom(f):
            continue

        try:

            # Get the new directory names
            d = get_dir_from_dicom(f)

            # Generate the new path name. Before making the directory via
            # the 'ensure_dir' method, the useAllFiles logic must be
            # evaluated. This is because of a potential error collison
            # with the 'rename' method.
            newPath = dicomDir / d
            newFile = newPath / f.name

            # Rename the directory or move the file depending on the user
            # option, returning the success
            if useAllFiles:
                newPath.mkdir(exist_ok=True)
                if newFile.exists():
                    continue
                else:
                    f.rename(newFile)
            else:
                dicomDir.rename(newPath)
                isSuccess = newPath.isdir()
                break

        except Exception as error:
            raise(error)

    # Remove the directory (if empty)
    lFiles = [f for f in dicomDir.iterdir()]
    if useAllFiles and not lFiles:
        dicomDir.remove()

    return isSuccess


def sort_dir(dicomDir: Union[Path, str]) -> None:
    """Sort a directory of DICOM images

    Parameters
    ----------
        dicomDir : str or pathlib object
            Full name of the directory containing DICOm files to be sorted

    """

    # Create Path object
    dicomDir = Path(dicomDir)
    if not dicomDir.is_dir():
        raise ValueError("dicomDir must be a directory not a file...")

    # Generate the files on the path
    lFiles = [f for f in dicomDir.rglob("*") if f.is_file()]

    for f in dicomDir.listdir():

        # Skip non-DICOM files
        if not is_dicom(f):
            continue

        # Read the file
        dcm = pydicom.dcmread(f)

        # Create the new directory
        newDir = dicomDir / get_dir_from_dicom(dcm)
        newDir.mkdir(exist_ok=True)
        f.rename(newDir / f.name)


# class dicom():
#     """Class for reading DICOM files

#     This class supports (via the gdcm package) reading DICOM meta-data and
#     image data. Currently, there is no support writing files

#     Parameters
#     ----------
#     file : str or py.path.local
#         Full file name of the DICOM file
#     autoLoad : bool
#         When True, the DICOM file data is imported on class instantiation.
#         Otherwise, the user is responsible for performing all read operations.

#     Attributes
#     ----------
#     file : py.path.local
#         DICOM file
#     header : dicomtools.header
#         Meta-data reference class

#     """

#     def __init__(self, file, **kwargs):

#         self._autoLoad = kwargs.get('autoLoad', True)
#         self.header = header(file, **kwargs)

#         # By setting the file, a large number of events will be performed
#         self.file = file

#     @property
#     def image(self):
#         """DICOM image array

#         """

#         if not hasattr(self, '_image'):
#             self._image = SimpleITK.ReadImage(self.file.strpath)

#         return SimpleITK.GetArrayFromImage(self._image)

#     def show(self, showSeries=False):
#         """Show the DICOM image(s)

#         Parameters
#         ----------
#         showSeries : bool
#             When True, any images in the same directory as the current DICOM
#             object will be displayed via a viewer that allows scrolling.
#             Default: False

#         """

#         if showSeries:

#             qApp = PyQt5.QtWidgets.QApplication([''])

#             v = viewer.DICOMviewer(self.file.dirpath().strpath)
#             v.setWindowTitle("DICOMviewer")
#             v.show()

#             sys.exit(qApp.exec_())

#         else:
#             # TODO: there should be some validation to ensure that the file
#             #       does, in fact, contain image data...
#             # TODO: the below code works if the DICOM image isn't a time series
#             img = SimpleITK.ReadImage(self.file.strpath)
#             imgArray = SimpleITK.GetArrayFromImage(img)[0, :, :]
#             fig = matplotlib.pyplot.Figure(frameon=False)
#             matplotlib.pyplot.imshow(imgArray, cmap='gray')
#             matplotlib.pyplot.show()
