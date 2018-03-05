# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:04:54 2018

@author: 703355681
"""

import dicomtools


def CollimatorShape(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060].lower()
    manufacturer = hdr[0x0008, 0x0070].lower()
    if (manufacturer == 'siemens') and (modality == 'DX'):
        return hdr[0x0018, 0x1147]
    else:
        return hdr[0x0018, 0x1700]


def FilterMaterial(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060]
    manufacturer = hdr[0x0008, 0x0070]
    if modality.lower() == 'dx':
        if manufacturer.lower() == 'siemens':
            # TODO: from Mobilett Mira DICOM conformance
            if 'cu' in hdr[0x0018, 0x1160]:
                return ['COPPER']
            elif 'cu' in hdr[0x0017, 0x100f].lower():
                return ['NONE']
            else:
                raise NotImplementedError()
        else:
            val = hdr[0x0018, 0x7050].upper()
            if '\\' in val:
                return val.split('\\')
            elif val:
                return [val]
            else:
                return ['UNKNOWN']


def FilterThicknessMinimum(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060]
    manufacturer = hdr[0x0008, 0x0070]
    if modality.lower() == 'dx':
        fMat = FilterMaterial(hdr)
        if manufacturer.lower() == 'siemens':
            if ('NONE' in fMat):
                return [0]
            elif ('COPPER' in fMat):
                raise NotImplementedError()
            else:
                raise NotImplementedError()
        else:
            if fMat:
                fl = hdr[0x0018, 0x7054]
                if (type(fl) != list):
                    fl = [fl]
                return fl
            else:
                return []


def FilterThicknessMaximum(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060]
    manufacturer = hdr[0x0008, 0x0070]
    if modality.lower() == 'dx':
        fMat = FilterMaterial(hdr)
        if manufacturer.lower() == 'siemens':
            if ('None' in fMat):
                return [0]
            elif ('Copper' in fMat):
                raise NotImplementedError()
            else:
                raise NotImplementedError()
        else:
            if fMat:
                fl = hdr[0x0018, 0x7054]
                if (type(fl) != list):
                    fl = [fl]
                return fl
            else:
                return []


def Grid(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060]
    manufacturer = hdr[0x0008, 0x0070].lower()
    if (modality == 'dx'):
        s = hdr[0x0018, 0x1166].upper()
        if (manufacturer == 'ge healthcare'):
            if 'FOCUSED' in s:
                return s + f'\\{hdr[0x0018, 0x704c]:.0f}'
            else:
                return s
        else:
            return s


def ImageProcessingDescription(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060]
    manufacturer = hdr[0x0008, 0x0070]
    if modality.lower() == 'dx':
        if manufacturer.lower() == 'siemens':
            return (f'{hdr[0x0018, 0x1400]}; Edge Enhancement Gain='
                    f'{hdr[0x0025, 0x100f]}; Edge Enhancement Kernel='
                    f'{hdr[0x0025, 0x100e]}; Harmonization Gain='
                    f'{hdr[0x0025, 0x100d]}; Harmonization Kernel='
                    f'{hdr[0x0025, 0x100c]}')
        else:
            return hdr[0x0018, 0x1400]


def Protocol(hdr):
    assert(type(hdr) == dicomtools.header)
    modality = hdr[0x0008, 0x0060].lower()
    manufacturer = hdr[0x0008, 0x0070].lower()
    if (modality == 'dx'):
        if (manufacturer == 'siemens'):
            # Per DICOM conformance
            return hdr[0x0040, 0x0254]
        elif (manufacturer == 'ge healthcare'):
            return hdr[0x0018, 0x1030]
        else:
            try:
                return hdr[0x0008, 0x1032][0].value
            except Exception as e:
                raise e


def Modality(hdr):
    assert(type(hdr) == dicomtools.header)
    manufacturer = hdr[0x0008, 0x0070].lower()
    if (manufacturer == 'imaging sciences international') and (
            hdr[0x0008, 0x1090] in ['17-19DX']):
        return 'PX'
    else:
        return hdr[0x0008, 0x0060]
