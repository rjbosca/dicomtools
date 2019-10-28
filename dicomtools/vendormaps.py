# -*- coding: utf-8 -*-


def BodyPart(hdr):
    """Returns the user-selected body part

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Body part

    """

    if hdr[0x0008, 0x2218].value:
        for de in hdr[0x0008, 0x2218].value:
            if (de.description.lower() == 'code meaning'):
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
        tuple:Field of view (usually in mm)

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        raise NotImplementedError("Unsupported modality")

    if (manufacturer == 'siemens'):

        # FIXME: validation should be performed using the same method as Philips
        freqFov, phaseFov = hdr[0x0051, 0x100c].value.split()[1].split(sep='*')

    elif (manufacturer == 'philips medical systems'):

        freqFov = float(hdr[0x0018, 0x1100].value)  # recon diameter

        # Calculate the phase FOV
        pctPhFov = float(hdr[0x0018, 0x0094].value)
        phaseFov = freqFov * pctPhFov / 100

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

    manufacturer = hdr[0x0008, 0x0070].lower()

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

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        return

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
        rxCoilName = hdr[0x0018, 0x1250].value
        # Check the Philips specific header for multi-coil selections
        de = hdr[0x2005, 0x140f][0]
        if (de[0x0018, 0x9048].value.lower() == "yes") and ((0x0018, 0x9047) in de):
            rxCoilName += f" ({de[0x0018, 0x9047].value})"
    elif (manufacturer == 'siemens'):
        # TODO: finish this...
        # The location depends on the Siemens CSA header version
        rxCoilName = 'Unknown'
        if (0x0018, 0x1250) in hdr:
            rxCoilName = hdr[0x0018, 0x1250].value
        if (hdr[0x0051, 0x1008].value.lower() == "image num 4"):
            rxCoilName += f" ({hdr[0x0051, 0x100f].value})"

    return rxCoilName


def SliceOrientation(hdr):
    """Get MR slice orientation

    Parameters:
        hdr (pydicom.dataset.FileDataset):DICOM data set

    Returns:
        str:Slice orientation

    """

    modality = hdr[0x0008, 0x0060].value.lower()
    manufacturer = hdr[0x0008, 0x0070].value.lower()

    if (modality != 'mr'):
        return

    #TODO: verify all of this using DICOM conformance statements
    #TODO: implement the GE version
    if (manufacturer == 'philips medical systems'):
        val = hdr[0x2001, 0x100b].value.lower()
    elif (manufacturer == 'siemens'):
        val = hdr[0x0051, 0x100e].value.lower()
    else:
        raise NotImplementedError(f"Unknown manufacturer: {manufacturer}")

    # Convert the vendor values to a standard lexicon
    if 'tra' in val:
        val = 'Axial'
    elif 'cor' in val:
        val = 'Coronal'
    elif 'sag' in val:
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
        # Check the Philips specific header for transmit coil name
        txCoilName = hdr[0x2005, 0x140f][0][0x0018, 0x9051].value
    elif (manufacturer == 'siemens'):
        txCoilName = hdr[0x0018, 0x1251].value
    else:
        raise NotImplementedError(
            f"Manufacturer '{manufacturer}' not supported...")

    return txCoilName
