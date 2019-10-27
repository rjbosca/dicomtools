# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 14:59:59 2017

@author: Ryan
"""

import os
import PyQt5
import numpy
import matplotlib
import SimpleITK

# Make sure that QT5 is being used. Note that the PyQt5 import must follow the
# matplotlib backend change.
matplotlib.use('Qt5Agg')
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow

# A few niceties from matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# from ipdb import set_trace
# PyQt5.QtCore.pyqtRemoveInputHook()
# set_trace()


# Define the GUI application
class DICOMviewer(QMainWindow):

    def __init__(self, dirDicom):

        super().__init__()
        self.init_app()
        self.init_canvas(dirDicom)

        # Update some window properties
        self.setWindowTitle("DICOM Viewer")
        self.show()

    def init_app(self):
        """Initialize the application window"""

        self.setAttribute(PyQt5.QtCore.Qt.WA_DeleteOnClose)
        self.setGeometry(300, 300, 500, 550)
        self.setWindowTitle("DICOM Viewer")

    def init_canvas(self, dirDicom):
        """Initialize the figure canvas within the application window"""

        # Generate a widget that will house the canvas
        self.canvasMain = QWidget(self)
        self.canvasMain.move(0, 0)
        self.canvasMain.resize(770, 770)

        # Initialize the file names
        if not os.path.isdir(dirDicom):
            raise NotADirectoryError(f'Invalid directory: {dirDicom}')
        self.canvasMain.dirDicom = dirDicom

        # Create a box that acts as the container for the canvas
        self.canvasMain = DICOMcanvas(self.canvasMain)

        # Change focus to the widget
        self.canvasMain.setFocus()

    def wheelEvent(self, event):

        # Get the mouse value
        val = event.angleDelta().y()
        val = round(val/abs(val))
        idx = val + self.canvasMain.idxFiles

        # Update the image index
        nFiles = self.canvasMain.image.GetSize()[2]
        if (idx < 0):
            idx = 0
        elif (idx > nFiles):
            idx = nFiles

        # Update the index and show the new image
        self.canvasMain.idxFiles = idx
        self.canvasMain.imshow()


# Define the figure canvas that will display the image
class DICOMcanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=5, dpi=100):

        # Initialize the figure that will be housed by the canvas
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_axes([0, 0.05, 1, 1])

        # Disable the axis ticks
        self.axes.get_xaxis().set_ticks([])
        self.axes.get_yaxis().set_ticks([])

        # Initialize the canvas
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   PyQt5.QtWidgets.QSizePolicy.Expanding,
                                   PyQt5.QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        # Initialize the image index and DICOM directory
        self.dirDicom = parent.dirDicom
        self.fileNames = \
            SimpleITK.ImageSeriesReader.GetGDCMSeriesFileNames(self.dirDicom)
        reader = SimpleITK.ImageSeriesReader()
        reader.SetFileNames(self.fileNames)
        self.image = reader.Execute()
        self.imageArray = SimpleITK.GetArrayFromImage(self.image)
        self.idxFiles = 0

        # Show the image
        self.imshow()

    def imshow(self):

        # Display the image
        fName = self.fileNames[self.idxFiles].replace(self.dirDicom, '')
        self.axes.imshow(self.imageArray[self.idxFiles, :, :], cmap='gray')
        self.axes.set_aspect('auto')
        self.axes.set_xlabel(f'Slice: {self.idxFiles+1}  File: ..{fName}')
        self.draw()
