# -*- coding: utf-8 -*-


def BodyPart(hdr):
    """Returns the user-selected body part

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Body part

    """

    if (hdr[0x0008, 0x0016].value == "1.2.840.10008.5.1.4.1.1.4.1"):

        # Get the Frame Anatomic Sequence from the Shared Functional
        # Group. Note that the following code assumes a single item
        # in the sequence
        fra = hdr[0x5200, 0x9229][0][0x0020, 0x9071][0]
        if (len(fra[0x0008, 0x2218].value) == 1):
            for de in fra[0x0008, 0x2218][0]:
                if (de.description().lower() == "code meaning"):
                    return de.value
        else:
            raise NotImplementedError("Anatomic Region Sequence length not equal to 1")

    else:

        if hdr[0x0008, 0x2218].value:
            for de in hdr[0x0008, 0x2218].value:
                if (de.description().lower() == 'code meaning'):
                    return de.value

    # In the event that the Anatomic Region Sequence or code meaning within
    # that sequence is not found, simply return the Body Part Examined tag
    return hdr[0x0018, 0x0015].value


def CollimatorShape(hdr):
    """Retruns the x-ray equipment collimator shape

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (manufacturer == 'siemens') and (modality == 'DX'):
        return hdr[0x0018, 0x1147].value
    else:
        return hdr[0x0018, 0x1700].value


def FieldOfView(hdr):
    """Returns the MRI in-plane field of view as a tuple

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        tuple:(Frequency FOV, Phase FOV) - usually in mm

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        raise NotImplementedError("Unsupported modality")

    if (manufacturer == 'siemens'):

        # FIXME: validation should be performed using the same method as Philips
        freqFov, phaseFov = hdr[0x0051, 0x100c].value.split()[1].split(sep='*')

    elif (manufacturer == 'philips medical systems'):

        if (hdr[0x0008, 0x0016].value == "1.2.840.10008.5.1.4.1.1.4.1"):
            # The Pre-Fram Functional Group sequence contains the Pixel Measures
            # Sequence. Get the pixel spacing from this sequence and calculate
            # the FoV using the number of rows/columns
            pmsg = hdr[0x5200, 0x9230][0][0x0028, 0x9110][0]
            fovRow = float(pmsg[0x0028, 0x0030].value[0]) * \
                     float(hdr[0x0028, 0x0010].value)
            fovCol = float(pmsg[0x0028, 0x0030].value[1]) * \
                     float(hdr[0x0028, 0x0011].value)

            # The In-Plane Phase Encoding Direction is found in the MR FOV/
            # Geometry Sequence of the Shared Functional Group.
            peDir = PhaseEncodingDir(hdr)
            if (peDir.lower() == "row"):
                phaseFov = fovRow
                freqFov = fovCol
            else:
                phaseFov = fovCol
                freqFov = fovRow

        else:

            # Reconstruction Diameter (0018, 1100) - value is a copy of the
            # largest value of the Field of View. Confirmed in DICOM conformance
            # statements of the following platforms:
            #
            # - Intera: R2.6.3
            # - Achieva: R2.6.3, R3.2

            # With the above definition, one has to use the percent phase FOV to
            # determine which encoding direction is represented by the recon FOV
            pctPhFov = float(hdr[0x0018, 0x0094].value) / 100.

            if (pctPhFov >= 1):
                phaseFov = float(hdr[0x0018, 0x1100].value)
                freqFov = phaseFov / pctPhFov
            else:
                freqFov = float(hdr[0x0018, 0x1100].value)
                phaseFov = freqFov * pctPhFov

    elif (manufacturer == 'hitachi medical corporation'):

        # FIXME: I have no means by which to determine the FOV as there is
        # little to no documentation (including the conformance statement).
        # However, it appears that the Reconstruction Diameter (0018, 1100)
        # seems to be the frequency FOV and the "RectFOVRatio" (0029, 100d)
        # is the phase FOV ratio

        pctPhFov = float(hdr[0x0029, 0x100d].value) / 100.

        if (pctPhFov == 1):
            freqFov = float(hdr[0x0018, 0x1100].value)
            phaseFov = float(hdr[0x0018, 0x1100].value) * pctPhFov
        else:
            raise NotADirectoryError(f"Unknown FOV for FOV ratio: {pctPhFov}")

    else:
        raise NotImplementedError(f"Unknown manufacturer: {manufacturer}")

    return (float(freqFov), float(phaseFov))


def FilterMaterial(hdr):
    """Get x-ray producing equipment filter material

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        list:Filter material

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality == 'dx'):
        if (manufacturer == 'siemens'):
            # TODO: from Mobilett Mira DICOM conformance
            if 'cu' in hdr[0x0018, 0x1160].value.lower():
                return ['COPPER']
            elif 'cu' in hdr[0x0017, 0x100f].value.lower():
                return ['NONE']
            else:
                return []
        else:
            val = hdr[0x0018, 0x7050].value.upper()
            if '\\' in val:
                return val.split('\\')
            elif val:
                return [val]
            else:
                return []


def FilterThicknessMinimum(hdr):
    """Get x-ray producing equipment filter material thickness minimum

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        list:Filter thickness minimum

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality == 'dx'):
        fMat = FilterMaterial(hdr)
        if (manufacturer == 'siemens'):
            if ('NONE' in fMat):
                return [0]
            elif ('COPPER' in fMat):
                raise NotImplementedError()
            else:
                raise NotImplementedError()
        else:
            if fMat:
                fl = hdr[0x0018, 0x7054].value

                if not fl:
                    return ['']*len(fMat)

                if (type(fl) != list):
                    fl = [fl]

                if (len(fMat) == len(fl)):
                    return fl
                else:
                    raise NotImplementedError()

            else:
                return []


def FilterThicknessMaximum(hdr):
    """Get x-ray producing equipment filter material thickness maximum

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        list:Filter thickness maximum

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality == 'dx'):
        fMat = FilterMaterial(hdr)
        if (manufacturer == 'siemens'):
            if ('NONE' in fMat):
                return [0]
            elif ('COPPER' in fMat):
                raise NotImplementedError()
            else:
                raise NotImplementedError()
        else:
            if fMat:

                fl = hdr[0x0018, 0x7054].value

                if not fl:
                    return ['']*len(fMat)

                if (type(fl) != list):
                    fl = [fl]

                if (len(fMat) == len(fl)):
                    return fl
                else:
                    raise NotImplementedError()

            else:
                return []


def Grid(hdr):
    """Get x-ray producing equipment grid information

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Grid information (if available)

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality == 'dx'):
        s = hdr[0x0018, 0x1166].value.upper()
        if (manufacturer == 'ge healthcare'):
            if 'FOCUSED' in s:
                return s + f'\\{hdr[0x0018, 0x704c].value:.0f}'
            else:
                return s
        else:
            return s


def ImageProcessingDescription(hdr):
    """Get x-ray producing equipment image processing information

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        tuple:Image processing details

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality == 'dx'):
        if (manufacturer == 'siemens'):
            return (f'{hdr[0x0018, 0x1400].value}; Edge Enhancement Gain='
                    f'{hdr[0x0025, 0x100f].value}; Edge Enhancement Kernel='
                    f'{hdr[0x0025, 0x100e].value}; Harmonization Gain='
                    f'{hdr[0x0025, 0x100d].value}; Harmonization Kernel='
                    f'{hdr[0x0025, 0x100c].value}')
        else:
            return hdr[0x0018, 0x1400].value


def Modality(hdr):
    """Get modality

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Modality

    """

    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (manufacturer == 'imaging sciences international') and (
            hdr[0x0008, 0x1090].value in ['17-19DX']):
        return 'PX'
    else:
        return hdr[0x0008, 0x0060].value


def PhaseEncodingDir(hdr):
    """Get MRI phase encoding direction

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Phase encoding direction

    """

    modality = Modality(hdr).lower()

    if (modality != 'mr'):
        return

    if (hdr[0x0008, 0x0016].value == "1.2.840.10008.5.1.4.1.1.4.1"):
        # From the Shared Functional Group sequence, the In-Plane Phase
        # Encoding Direction can be found in the MR FOV/Geometry Sequence
        #FIXME: how to handle multi-slice data sets?
        geo = hdr[0x5200, 0x9229][0][0x0018, 0x9125]
        if (len(geo.value) != 1):
            raise NotImplementedError("Shared Functional Group sequence has "
                                      "too many itmes.")
        else:
            phEncDir = geo[0][0x0018, 0x1312].value
    else:
        # Try to find the phase encoding direction the easy way
        try:

            acqMat = hdr[0x0018, 0x1310].value  # acquisition matrix
            phEncDir = hdr[0x0018, 0x1312].value  # in-plane phase enc dir

            # Validate using the acquisition matrix
            if (phEncDir == 'ROW') and not (not acqMat[0] and acqMat[1]) or \
                    (phEncDir == 'COL') and not (acqMat[0] and not acqMat[1]):
                raise NotImplementedError("Phase encoding direction mismatch")

        except KeyError:

            # The "In-Plane Phase Encoding Direction" tag wasn't there.
            # Instead use the acquisition matrix to determine the direction.
            if acqMat[0] and not acqMat[1]:
                return 'COL'
            elif not acqMat[0] and acqMat[1]:
                return 'ROW'
            else:
                raise NotImplementedError("An unknown error occured.")

    return phEncDir


def Protocol(hdr):

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality == 'dx'):
        if (manufacturer == 'siemens'):
            # Per DICOM conformance
            val = hdr[0x0040, 0x0254].value
        elif (manufacturer == 'ge healthcare'):
            val = hdr[0x0018, 0x1030].value
        else:
            val = hdr[0x0008, 0x1032].value
            if (type(val) == list):
                val = val[0].value
        return val


def ReceiveCoil(hdr):
    """Get MR receive coil name

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Receive coil name (if available)

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        return

    if (manufacturer == 'philips medical systems'):
        if (hdr[0x0008, 0x0016].value == "1.2.840.10008.5.1.4.1.1.4.1"):
            # The Receive Coil Name can be found in the Shared Functional
            # Group sequence in the MR Receive Coil Sequence
            if len(hdr[0x5200, 0x9229].value) != 1:
                raise NotImplementedError("Shared Functional Group sequence has "
                                          "too many values.")
            elif (len(hdr[0x5200, 0x9229][0][0x0018, 0x9042].value) != 1):
                raise NotImplementedError("MR Receive Coil sequence has "
                                          "too many values.")
            else:
                sfg = hdr[0x5200, 0x9229][0][0x0018, 0x9042][0]
                rxCoilName = sfg[0x0018, 0x1250].value
                #FIXME: provide a means, as with older DICOM headers to
                # determine which coil elements were used.
        else:
            rxCoilName = hdr[0x0018, 0x1250].value
            # Check the Philips specific header for multi-coil selections
            de = hdr[0x2005, 0x140f][0]
            if (de[0x0018, 0x9048].value.lower() == "yes") and ((0x0018, 0x9047) in de):
                rxCoilName += f" ({de[0x0018, 0x9047].value})"
    elif (manufacturer == 'siemens'):
        # TODO: finish this...
        # The location depends on the Siemens CSA header version
        rxCoilName = 'Unknown'
        rxCoilStr = ''
        if (0x0018, 0x1250) in hdr:
            rxCoilName = hdr[0x0018, 0x1250].value
        if (hdr[0x0051, 0x1008].value.lower() == "image num 4"):
            rxCoilStr = hdr[0x0051, 0x100f].value

        # Check for a coil signature based on the coil string
        if (rxCoilName == 'Unknown') and rxCoilStr:

            # Body 18
            lBody18 = [f'BO{i}' for i in range(1, 4)]
            if any([cs in rxCoilStr for cs in lBody18]):
                rxCoilName = 'Body 18'

            # Body coil
            if ('BC' in rxCoilStr):
                rxCoilName = 'Body Coil'

            # Breast 2Ch
            lBreast2Ch = ['BL1', 'BR1']
            if any([cs in rxCoilStr for cs in lBreast2Ch]):
                rxCoilName = 'Breast 2Ch'

            # Breast 4Ch
            lBreast4Ch = ['BL2', 'BR2']
            if any([cs in rxCoilStr for cs in lBreast4Ch]):
                rxCoilName = 'Breast 4Ch'

            # Breast 16Ch
            lBreast16Ch = ['BL4', 'BR4']
            if any([cs in rxCoilStr for cs in lBreast16Ch]):
                rxCoilName = 'Breast 16Ch'

            # Flex large
            if ('FL' in rxCoilStr):
                rxCoilName = 'Flex Large'

            # Flex small
            if ('FS' in rxCoilStr):
                rxCoilName = 'Flex Small'

            # Foot/ankle 20
            lFootAnkle20 = ['FA', 'TO']
            if any([cs in rxCoilStr for cs in lFootAnkle20]):
                rxCoilName = 'Foot/Ankle 20'

            # Hand/wrist 16
            lHandWrist16 = [f'HW{i}' for i in range(1, 4)]
            if any([cs in rxCoilStr for cs in lHandWrist16]):
                rxCoilName = 'Hand/Wrist 16'

            # Head-neck 20 receive coil
            lHeadNeck20 = ['HE1', 'HE2', 'HE3', 'HE4', 'NE1', 'NE2']
            if any([cs in rxCoilStr for cs in lHeadNeck20]):
                rxCoilName = 'Head-Neck 20'

            # Peripheral angio
            lPeriphAngio = [f'PA{i}' for i in range(1, 7)]
            if any([cs in rxCoilStr for cs in lPeriphAngio]):
                rxCoilName = 'Periph Angio'

            # Shoulder large 16
            if ('SHL' in rxCoilStr):
                rxCoilName = 'Shoulder Large 16'

            # Shoulder small 16
            if ('SHS' in rxCoilStr):
                rxCoilName = 'Shoulder Small 16'

            # Spine 32
            lSpine32 = [f'SP{i}' for i in range(1, 9)]
            if any([cs in rxCoilStr for cs in lSpine32]):
                rxCoilName = 'Spine 32'

            # TxRx Knee 15
            if ('15K' in rxCoilStr):
                rxCoilName = 'TxRx Knee 15'

        # Update the coil string
        if rxCoilStr:
            rxCoilName += f" ({hdr[0x0051, 0x100f].value})"

    return rxCoilName


def SliceOrientation(hdr):
    """Get MR slice orientation

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Slice orientation

    """

    import numpy

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        return

    #TODO: verify all of this using DICOM conformance statements
    #TODO: implement the GE version
    if (manufacturer == 'philips medical systems'):
        if (hdr[0x0008, 0x0016].value == "1.2.840.10008.5.1.4.1.1.4.1"):

            # The following dictionary is used to determine the image
            # orientation
            dictOrient = {0: "Sagittal",
                          1: "Coronal",
                          2: "Axial",
                          3: "Oblique"}

            # The Plane Orientation Sequence is stored in the Functional 
            # Groups Sequence. Calculating the quadrature of the direction
            # cosines provides a good estimate of the orientation
            pfg = hdr[0x5200, 0x9230][0][0x0020, 0x9116][0]
            iop = numpy.array(pfg[0x0020, 0x0037].value, dtype=float)

            # Calculate the possible orientations
            d = numpy.vstack((iop - numpy.array([0,1,0,0,0,-1]),  # sagittal
                              iop - numpy.array([1,0,0,0,0,-1]),  # coronal
                              iop - numpy.array([1,0,0,0,1, 0])))  # axial
            d = numpy.sum(d*d, axis=1)

            # Set the output to oblique for angles greater than 45 deg.
            ind = numpy.argmin(d)
            if d.min() > numpy.cos(numpy.pi/4)**2:
                ind = 3

            # Set the return value
            val = dictOrient[ind]

        else:
            val = hdr[0x2001, 0x100b].value.lower()
    elif (manufacturer == 'siemens'):
        val = hdr[0x0051, 0x100e].value.lower()
    elif (manufacturer == 'hitachi medical corporation'):
        # Valid for V4.0 and V5.0
        val = hdr[0x0019, 0x1002].value.lower()
    else:
        raise NotImplementedError(f"Unknown manufacturer: {manufacturer}")

    # Convert the vendor values to a standard lexicon
    if ('tra' in val.lower()) or ('ax' in val).lower():
        val = 'Axial'
    elif 'cor' in val.lower():
        val = 'Coronal'
    elif 'sag' in val.lower():
        val = 'Sagittal'
    else:
        raise NotImplementedError(f"Unkown orientation: {val}")

    return val


def TransmitCoil(hdr):
    """Get MR transmit coil name

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Transmit coil name

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        return

    # TODO: see if the SAR fields might shed some light on the transmit
    if (manufacturer == 'philips medical systems'):
        if (hdr[0x0008, 0x0016].value == "1.2.840.10008.5.1.4.1.1.4.1"):
            # The Transmit Coil Name can be found in the Shared Functional
            # Group sequence in the MR Transmit Coil Sequence
            if len(hdr[0x5200, 0x9229].value) != 1:
                raise NotImplementedError("Shared Functional Group sequence has "
                                          "too many values.")
            elif (len(hdr[0x5200, 0x9229][0][0x0018, 0x9049].value) != 1):
                raise NotImplementedError("MR Receive Coil sequence has "
                                          "too many values.")
            else:
                sfg = hdr[0x5200, 0x9229][0][0x0018, 0x9049][0]
                txCoilName = sfg[0x0018, 0x9051].value
                #FIXME: provide a means, as with older DICOM headers to
                # determine which coil elements were used.
        else:
            # Check the Philips specific header for transmit coil name
            txCoilName = hdr[0x2005, 0x140f][0][0x0018, 0x9051].value
    elif (manufacturer == 'siemens'):
        txCoilName = hdr[0x0018, 0x1251].value
    else:
        raise NotImplementedError(
            f"Manufacturer '{manufacturer}' not supported...")

    return txCoilName
