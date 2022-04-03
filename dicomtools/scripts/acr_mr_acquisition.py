# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:28:56 2018

@author: Ryan
"""

from pathlib import Path
import pydicom
import dicomtools
from datetime import datetime
from pydicom.errors import InvalidDicomError

#TODO: develop code for MV (Philips) to ensure that the correct
#      pixel pitch is determined. The Enhanced MR DICOM files
#      provides the wrong number of phase encodes
#TODO: build the ACR requirements into the excel spreadsheet
#TODO: add comments in, for example, the flip angle column 
#      heading to warn the user that the FA might mean either
#      excitation or refocusing...

def store_or_pass(d, fld, val):
    if fld not in d:
        d[fld] = val


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


def _is_mr_acq(hdr):

    # Check for MR or enhanced MR SOP class UIDs and a modality of MR
    is_mr_sop = any([hdr.SOPClassUID == '1.2.840.10008.5.1.4.1.1.4',
                     hdr.SOPClassUID == '1.2.840.10008.5.1.4.1.1.4.1'])
    is_mr = all([is_mr_sop, hdr.Modality.lower() == 'mr'])
    if not is_mr:
        return False

    ProtocolName = _get_val(hdr, (0x0018, 0x1030)).lower()
    ImageType = [val.lower() for val in hdr[0x0008, 0x0008].value]
    SeriesDescription = _get_val(hdr, (0x0008, 0x103E)).lower()

    # Check for common scout/localizer names in the series description or
    # selected protocol
    is_scout = any([("scout" in SeriesDescription),
                    ("scout" in ProtocolName),
                    ("survey" in SeriesDescription),
                    ("survey" in ProtocolName)])

    # Check for common non-acquisition data sets in the

    is_mip = any(['mip' in val for val in ImageType])
    is_bolus_track = any(['carebolus' in val for val in ImageType]) or \
                     ('carebolus' in ProtocolName)
    is_not_acq = is_scout or any([is_mip,
                                  is_bolus_track,
                                  'DERIVED' in hdr.ImageType,
                                  'SECONDARY' in hdr.ImageType])

    return all([not is_scout,not is_not_acq])


def _get_common_dir(files):

    return Path.home() / "Desktop"
#    files = set(files)
#    for f in files:
#        is_all = all([f in fl for fl in files])
#        if is_all:
#            # Found the commmon parent directory
#            return f
#
#    raise NotImplementedError("Unable to determine the "
#                              "common parent directory")


def _get_slice_info(hdr):

    slTh = 0.0
    slGap = 0.0

    if "SliceThickness" in hdr:
        slTh = float(hdr.SliceThickness)
        if "SpacingBetweenSlices" in hdr:
            slGap = float(hdr.SpacingBetweenSlices) - slTh
    else:
    
        for de in hdr[0x5200, 0x9230]:
            pms = de[0x0028, 0x9110][0]
            if pms[0x0018, 0x0050].value > pms[0x0018, 0x0088].value:
                slTh = float(pms[0x0018, 0x0088].value)
                slGap = 0.0
                break
            else:
                slTh = float(pms[0x0018, 0x0050].value)
                slGap =float(pms[0x0018, 0x0088].value) - slGap
                break

    return (slTh, slGap)


def _get_acq_matrix(hdr):

    if "AcquisitionMatrix" in hdr:
        acqMat = hdr.AcquisitionMatrix
    else:
        acqMat = hdr[0x5200, 0x9230][0][0x2005, 0x140f][0][0x0018, 0x1310].value

    return (max(acqMat[:2]), max(acqMat[-2:]))


def _get_acq_params(hdr):

    fa = 0.0
    te = 0.0
    tr = 0.0
    ti = ""
    nsa = 0.0

    if ("FlipAngle" in hdr) and ("RepetitionTime" in hdr) and \
        ("EchoTime" in hdr) and ("NumberOfAverages" in hdr):
        fa = hdr.FlipAngle
        tr = round(hdr.RepetitionTime, 4)
        te = round(hdr.EchoTime, 1)
        nsa = hdr.NumberOfAverages

    else:
        # MR Timing and Related Parameters sequence is found within the Shared
        # Functional Groups sequence. Info about the TR, FA, ETL, SAR, dB/dt,
        # etc. is stored here.
        mrt = hdr[0x5200, 0x9229][0][0x0018, 0x9112][0]
        fa = mrt.FlipAngle
        tr = round(mrt.RepetitionTime, 4)

        # MR Echo sequence is found within the Shared Functional Groups
        # sequence.
        mres = hdr[0x5200, 0x9230][0][0x0018, 0x9114][0]
        te = round(mres.EffectiveEchoTime, 1)

        # MR Averages sequence within the Shared Functional Groups sequence
        # appears to contain the number of signal averages
        mras = hdr[0x5200, 0x9229][0][0x0018, 0x9119][0]
        nsa = float(mras.NumberOfAverages)

    if ("InversionTime" in hdr):
        ti = round(hdr.InversionTime, 4)

    return (fa, te, tr, nsa, ti)


def gen_protocol(args):
    import pandas

    # Initialize the protocol dictionary definition. Be sure to copy
    # this before updating.
    prot = {r'Sequence name/type': [],
            r'Sequence #': [],
            r'Orientation': [],
            r'Dimension (2D/3D)': [],
            r'Slice thickness (mm)': [],
            r'Gap (mm)': [],
            r'FOVp': [],
            r'FOVf': [],
            r'Np': [],
            r'Nf': [],
            r'Area (mm^2)': [],
            r'# Acquisitions': [],
            r'TR (ms)': [],
            r'TE (ms)': [],
            r'Flip Angle': [],
            r'TI': [],
            r'Number of Images': [],
            }
    pat = {r'Accession': [],
           r'Date of exam': [],
           r'Age of patient': [],
           r'Weight of patient': [],
           }

    exams = dict([])
    files = []
    sop = []  # stores SOP classes of import files

    for f in args.files:

        files.append(f.parent)

        # Attempt to load the DICOM file
        try:
            hdr = pydicom.dcmread(f)
            hdr.decode()
            if 'DirectoryRecordSequence' in hdr:
                continue
            uid = hdr.StudyInstanceUID
            sop.append(hdr.SOPClassUID)
        except InvalidDicomError:
            print("Ignorning non-DICOM compliant files...")
            continue
        except Exception:
            raise ImportError("Unable to import DICOM file")

        # Verify the appropriate data content
        if not _is_mr_acq(hdr):
            continue

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
            dPat[r'Weight of patient'].append(hdr.PatientWeight)

        # Get the dictionary and started adding some data
        dProt = exams[uid][1]
        if hdr.SeriesNumber in dProt[r'Sequence #']:
            idx = dProt[r'Sequence #'].index(hdr.SeriesNumber)
            dProt[r'Number of Images'][idx] += 1
            continue

        # Get the series data
        dProt[r'Sequence name/type'].append(hdr.SeriesDescription)
        dProt[r'Sequence #'].append(hdr.SeriesNumber)
        dProt[r'Orientation'].append(dicomtools.SliceOrientation(hdr))
        dProt[r'Dimension (2D/3D)'].append(hdr.MRAcquisitionType)
        slTh, slGap = _get_slice_info(hdr)
        dProt[r'Slice thickness (mm)'].append(slTh)
        dProt[r'Gap (mm)'].append(slGap)
        fov = dicomtools.FieldOfView(hdr)
        dProt[r'FOVp'].append(round(fov[1], 1))
        dProt[r'FOVf'].append(round(fov[0], 1))
        nF, nP = _get_acq_matrix(hdr)
        dProt[r'Np'].append(nP)
        dProt[r'Nf'].append(nF)
        dProt[r'Area (mm^2)'].append(fov[0]*fov[1]/(nF*nP))
        fa, te, tr, nsa, ti = _get_acq_params(hdr)
        dProt[r'# Acquisitions'].append(nsa)
        dProt[r'TR (ms)'].append(tr)
        dProt[r'TE (ms)'].append(te)
        dProt[r'Flip Angle'].append(fa)
        dProt[r'TI'].append(ti)
        dProt[r'Number of Images'].append(1)

        # Store some general information about the exam
#        name = str(hdr.PatientName).replace('^',
#                                            ' ').strip().replace(' ', ', ')
#        eStr = '; '.join([f'Acc. #: {hdr.AccessionNumber}',
#                          f'MRN: {name}',
#                          f'date: {hdr.StudyDate}',
#                          f'name: {hdr.PatientName}',
#                          f'study desc.: {hdr.StudyDescription}'])

    if args.out:
        dOut = _get_common_dir(files)
        for ek in exams.keys():

            # Convert the patient info to a data frame
            dfPat = pandas.DataFrame(exams[ek][0])
            dfProt = pandas.DataFrame(exams[ek][1])

            # Create the Excel writer
            fXlsx = dfPat['Accession'][0] +  '.xlsx'
            writer = pandas.ExcelWriter(dOut / fXlsx.replace("'", ""),
                                        engine='xlsxwriter')

            dfPat.to_excel(writer, 'Sheet1', index=False)
            dfProt.to_excel(writer, 'Sheet1', startrow=3, index=False)

            writer.save()


def _create_parser():

    import argparse

    parser = argparse.ArgumentParser(prog=__file__,
                                     description="Mammo DICOM protocol dump")
    parser.add_argument('files', metavar='F', type=str, nargs='+',
                        help=('Valid DICOM file (or directory) names '
                              'from which to generate a protocol'))
    parser.add_argument('--out', action='store_true', default=False,
                        help='Output MR acquisitions to Excel file')

    return parser


def _process_args(args, parser):

    print("")

    # Process the output flag first, as a target output file will needed to be
    # generated from the input files
    validFiles = []
    for fl in args.files:
        fl = Path(fl)
        if not fl.exists():
            print(f"Ignoring invalid file or directory: {fl}")
            continue
        elif fl.is_dir():
            # All files, DICOM or not, are added to the list. DICOM checks
            # will occur later
            [validFiles.append(f) for f in fl.glob("**/*") if f.is_file()]
        elif fl.is_file():
            validFiles.append(fl)
        else:
            raise NotImplementedError("An unknown error has occured while "
                                      "parsing the file/dir input...")

    if not len(validFiles):
        print("")
        parser.error("At least one valid DICOM file or directory name must be"
                     " specified")
    else:
        args.files = validFiles


if __name__ == "__main__":

    parser = _create_parser()
    args = parser.parse_args()
    _process_args(args, parser)
    gen_protocol(args)
