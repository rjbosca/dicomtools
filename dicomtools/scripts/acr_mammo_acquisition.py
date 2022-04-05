# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:28:56 2018

@author: Ryan
"""

from pathlib import Path
from datetime import datetime

import pydicom

#TODO: build the ACR requirements into the excel spreadsheet

def store_or_pass(d, fld, val):
    if fld not in d:
        d[fld] = val


def _get_acq_params(hdr):

    is_vol = ('VOLUME' in hdr.ImageType)
    if is_vol:
        hdr = hdr.XRay3DAcquisitionSequence[0]
        
    # Use the DICOM header (or 3D acquisition sequence) to get the paramters
    kVp = _get_val(hdr, (0x0018, 0x0060))  # kVp
    if is_vol:
        mAs = _get_val(hdr, (0x0018, 0x9332))
        t = _get_val(hdr, (0x0018, 0x9328))  # Exposure time
    else:
        mAs = _get_val(hdr, (0x0018, 0x1153)) / 1000
        t = _get_val(hdr, (0x0018, 0x1150))  # Exposure time
    fs = _get_val(hdr, (0x0018, 0x1190))  # Focal spot(s)
    target = _get_val(hdr, (0x0018, 0x1191))  # Anode Target Material
    filt = _get_val(hdr, (0x0018, 0x7050))  # Filter material
    th = _get_val(hdr, (0x0018, 0x11a0))  # Body part thickness
    force = _get_val(hdr, (0x0018, 0x11a2))  # Compression force
    ang = _get_val(hdr, (0x0018, 0x1510))

    return kVp, t, mAs, fs, target, filt, th, force, ang


def _get_val(hdr, fld):

    try:
        if type(fld) == tuple:
            return hdr[fld].value
        elif type(fld) == str:
            return hdr.get(fld)
        else:
            raise NotImplementedError()
    except KeyError:
        return ""


def gen_protocol(args):

    import pandas

    # Initialize the protocol dictionary definition. Be sure to copy
    # this before updating.
    prot = {r'View': [],
            r'Type': [],
            r'Angle (deg)': [],
            r'Compression Force (N)': [],
            r'Compressed Breast Thickness (mm)': [],
            r'kVp': [],
            r'Time (ms)': [],
            r'mAs': [],
            r'Nominal Focal Spot (mm)': [],
            r'Tube Target': [],
            r'Filter': []
            }
    pat = {r'Accession': [],
           r'Date of exam': [],
           r'Age of patient': [],
           }

    exams = dict([])
    files = []
    sop = []  # stores SOP classes of import files

    for f in args.files:

        files.append(f.parent)

        # Attempt to load the DICOM file
        hdr = pydicom.dcmread(f)
        hdr.decode()

        # Verify the appropriate data content
        if (type(hdr) == pydicom.dicomdir.DicomDir) or (hdr.Modality == 'SR'):
            continue

        # Get some identifying info
        uid = hdr.StudyInstanceUID
        sop.append(hdr.SOPClassUID)

        # File successfully read - create a new container for this UID (or do
        # nothing if the UID exists)
        store_or_pass(exams, uid, [pat.copy(), prot.copy()])

        # Get the patient information
        dPat = exams[uid][0]
        sDate = datetime.strptime(hdr.StudyDate, "%Y%m%d")
        sDate = rf"{sDate.month}/{sDate.day}/{sDate.year}"
        if dPat[r'Date of exam'] and (sDate != dPat[r'Date of exam'][0]):
            raise NotImplementedError("Multiple exam dates...")
        elif not dPat[r'Date of exam']:
            pAge = (datetime.strptime(hdr.StudyDate, "%Y%m%d") -
                    datetime.strptime(hdr.PatientBirthDate, "%Y%m%d"))
            if (pAge.days > 365):
                pAge = round(pAge.days/365)
            else:
                raise NotImplementedError("Patient's age is less than 1 year.")
            dPat[r'Accession'].append("'" + str(hdr.AccessionNumber))
            dPat[r'Date of exam'].append(sDate)
            dPat[r'Age of patient'].append(pAge)

        # Get the dictionary and started adding some data
        dProt = exams[uid][1]

        # Get the series data
        view = (_get_val(hdr, (0x0020, 0x0062)) + " " + 
                _get_val(hdr, (0x0018, 0x5101)))
        dProt[r'View'].append(view)  # View position
        if ('VOLUME' in hdr.ImageType):
            dProt[r'Type'].append("3D")
        else:
            dProt[r'Type'].append("2D")
        kVp, t, mAs, fs, target, filt, th, force, ang = _get_acq_params(hdr)
        dProt[r'Angle (deg)'].append(ang)
        dProt[r'Compression Force (N)'].append(force)
        dProt[r'Compressed Breast Thickness (mm)'].append(th)
        dProt[r'kVp'].append(kVp)
        dProt[r'Time (ms)'].append(t)
        dProt[r'mAs'].append(mAs)
        dProt[r'Nominal Focal Spot (mm)'].append(fs)
        dProt[r'Tube Target'].append(target)
        dProt[r'Filter'].append(filt)

    if args.out:
        for ek in exams.keys():

            # Convert the patient info to a data frame
            dfPat = pandas.DataFrame(exams[ek][0])
            dfProt = pandas.DataFrame(exams[ek][1])

            # Create the Excel writer
            writer = pandas.ExcelWriter(
                                    args.dir / (dfPat['Accession'][0] + '.xlsx'),
                                    engine='xlsxwriter'
                                    )

            dfPat.to_excel(writer, 'Sheet1', index=False)
            dfProt.to_excel(writer, 'Sheet1', startrow=3, index=False)

            writer.save()


def _create_parser():

    import argparse

    parser = argparse.ArgumentParser(prog=__file__,
                                     description="Mammo DICOM protocol dump")
    parser.add_argument('path', metavar='F', type=str, nargs=1,
                        help=('Valid DICOM file (or directory) name '
                              'from which to generate a protocol'))
    parser.add_argument('--out', action='store_true', default=False,
                        help='Output mammo acquisitions to Excel file')

    return parser


def _process_args(args, parser):

    from pydicom.misc import is_dicom

    # Process the output flag first, as a target output file will needed to be
    # generated from the input files
    fl = Path(args.path[0])
    if fl.is_dir():  # add all DICOM files
        args.dir = fl
        dcms = [f for f in fl.rglob("*") if f.is_file() and is_dicom(f)]
    elif fl.is_file():  # evaluate single file
        args.dir = fl.parent
        dcms = [fl]
    else:
        raise NotImplementedError("An unknown error has occured while "
                                  f"parsing the file/dir input:\n{fl}")

    if not len(dcms):
        parser.error("At least one valid DICOM file or directory name must be"
                     " specified")
    else:
        args.files = dcms


if __name__ == "__main__":

    parser = _create_parser()
    args = parser.parse_args()
    _process_args(args, parser)
    gen_protocol(args)
