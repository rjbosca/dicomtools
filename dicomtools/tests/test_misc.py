"""Test for misc.py"""

import pytest

from pydicom import Dataset

from dicomtools.misc import ipp2plane

class TestMisc:
    def test_ipp2plane(self):
        """Test the ipp2plane function"""

        # Create the MR Image Storage data set
        ds = Dataset()
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.4"

        #TODO: write tests to evalute user-specified angle other then 0

        # Test axial
        ds.ImageOrientationPatient = [1,0,0,0,1, 0]
        assert ipp2plane(ds) == "Axial"
        assert ipp2plane(ds, 0.) == "Axial"

        # Test coronal
        ds.ImageOrientationPatient = [1,0,0,0,0,-1]
        assert ipp2plane(ds) == "Coronal"
        assert ipp2plane(ds, 0.) == "Coronal"

        # Test sagittal
        ds.ImageOrientationPatient = [0,1,0,0,0,-1]
        assert ipp2plane(ds) == "Sagittal"
        assert ipp2plane(ds, 0.) == "Sagittal"

        # Test the value error for an invalide SOP Class
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.6.1"  # U/S image
        with pytest.raises(NotImplementedError):
            ipp2plane(ds)

        #TODO: write tests for Enhanced MR Image Storage