# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:28:56 2018

@author: Ryan
"""

"""
The ACR requires the following exam protocol details:
    -Indication
    -Scanner acquisition settings: routine kV, mA/mAs/effective mAs, 
     collimation (N x T), pitch, rotation time, usage of radiation dose 
     reduction methods (AEC such as TCM, settigns for dose reudctions methods,
     etc)
    -Phase of respiration
    -Reconstruction settings: reconstructed image width (slice thickness),
     reconstruction interval, reconstruction kernel/filter, reconstructed FOV
    -Anatomical coverage (i.e. lung apices to lung bases, top of diaphragm to
     iliac crest, etc.)
    -IV contrast (with injectino rate and scan delay), if applicable
    -EKG gating (cardiac studies) policy
    
For the application materials, the following parameters are needed:
    -kV
    -mA
    -Time per rotation (s)
    -mAs (calculated by the System)
    -Effective mAs (or mAs per slice) as displayed by scanner
    -Scan FOV (cm)
    -Display FOV (cm)
    -Reconstructions Algorithm
    -Axial (A) or Helical (H)
    -# Data Channels used in a single rotation (N)
    -Z-axis collimation (T) in mmm
    -Table Increment (I): Axial Scans (mm) or Helical scans (mm/rotation)
    -Reconstructed image width (mm)
    -Reconsructed image interval (mm)
    -Dose Reduction Technique used
    -CTDIvol (recorded after scanning, not before, for each Scan Sequence)
    -DLP (Dose length product)
"""

from pathlib import Path
import copy
import pydicom
from datetime import datetime
from pydicom.errors import InvalidDicomError


def store_or_pass(d, fld, val):
    """
    Store or ignore input data into dictionary

    Parameters
    ----------
    d : dict
        Destination dictionary.
    fld : str
        Dictionary key, if not present in d, will store val in the new key.
    val
        Value to be stored at fld if fld is not present in d.

    Returns
    -------
    None.

    """
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


def _is_not_ct_acq(hdr):
    """
    Check that a DICOM header is acquired image data

    Parameters
    ----------
    hdr : pydicom.dataset.FileDataset
        DICOM header to be evaluated.

    Returns
    -------
    isAcq : bool
        True if DICOM header does not contain CT acquisition data.

    """

    # SeriesDescription = _get_val(hdr, (0x0008, 0x1030)).lower()
    # ProtocolName = _get_val(hdr, (0x0018, 0x1030)).lower()

    is_ct = all(['1.2.840.10008.5.1.4.1.1.2' in _get_val(hdr, (0x0008, 0x0016)),
                 hdr.Modality.lower() == 'ct',
                 'DirectoryRecordSequence' not in hdr])

    is_not_acq = any(['DERIVED' in _get_val(hdr, (0x0008, 0x0008)),
                      'SECONDARY' in _get_val(hdr, (0x0008, 0x0008)),
                      'VOLUME' in _get_val(hdr, (0x0008, 0x0008)),
                      'LOCALIZER' in _get_val(hdr, (0x0008, 0x0008))])

    return any([not is_ct,
                is_not_acq,
                ])


def gen_protocol(args):

    # Package imports
    import pandas

    # Initialize the protocol dictionary definition. **Be sure to copy this
    # before updating**
    prot = {r'Sequence #': [],
            r'Protocol name': [],
            r'Requested procedure': [],
            r'kVp': [],
            r'mA': [],
            r'Rotation Time': [],
            r'mAs': [],
            r'Effective mAs': [],
            r'Scan FOV': [],
            r'Display FOV': [],
            r'Reconstruction Algorithm': [],
            r'Scan Type': [],
            r'# of data channels (N)': [],
            r'Z-axis width (T) (mm)': [],
            r'Table increment (I)': [],
            r'Reconstructed image width (mm)': [],
            r'Reconstructed image interval (mm)': [],
            r'Slice Location': [],
            r'Dose reduction technique': [],
            r'CTDIvol': [],
            r'Pitch': []
            }
    pat = {r'Accession': [],
           r'Date of exam': [],
           r'Age of patient': [],
           r'Weight of patient': [],
           }

    exams = dict([])
    files = []
    sop = []  # stores SOP classes of imported files

    # Read all of the files
    for f in args.files:

        files.append(f.parent)

        # Attempt to load the DICOM file
        try:
            hdr = pydicom.dcmread(f)
            hdr.decode()
            if _is_not_ct_acq(hdr):  # ignore non-acquisition data
                continue
            uid = _get_val(hdr, (0x0020, 0x000d))
            suid = _get_val(hdr, (0x0020, 0x000e))
            sop.append(_get_val(hdr, (0x0008, 0x0016)))
        except InvalidDicomError:
            print("Ignorning non-DICOM compliant files...")
            continue
        except Exception:
            raise ImportError("Unable to import DICOM file")

        # File successfully read - create a new container for this series and
        # exam UID (or do nothing if the UID exists)
        store_or_pass(exams, uid, dict([]))
        store_or_pass(exams[uid], suid, [copy.deepcopy(pat), copy.deepcopy(prot), int(0)])

        # Some general info that will be used on multiple occassions
        dPat = exams[uid][suid][0]
        sDateDT = datetime.strptime(_get_val(hdr, (0x0008, 0x0020)), "%Y%m%d")
        sDate = rf"{sDateDT.month}/{sDateDT.day}/{sDateDT.year}"
        name = str(hdr.PatientName).replace('^',
                                            ' ').strip().replace(' ', ', ')
        mfr = _get_val(hdr, (0x0008, 0x0070)).lower()
        model = _get_val(hdr, (0x0008, 0x1090)).lower()
        ver = _get_val(hdr, (0x0018, 0x1020)).lower()
        col = float(_get_val(hdr, (0x0018, 0x9307)))
        sType = _get_val(hdr, (0x0018, 0x9302))

        # Get the patient information
        if dPat[r'Date of exam'] and (sDate != dPat[r'Date of exam'][0]):
            raise NotImplementedError("Multiple exam dates...")
        elif not dPat[r'Date of exam']:
            dob = _get_val(hdr, (0x0010, 0x0030))
            if dob:
                pAge = (sDateDT - datetime.strptime(dob, "%Y%m%d"))
                if (pAge.days > 365):
                    pAge = round(pAge.days/365)
                else:
                    pAge = round(pAge.days/365, 2)
                dPat[r'Age of patient'].append(pAge)
            else:
                dPat[r'Age of patient'].append(None)

            dPat[r'Accession'].append(_get_val(hdr, (0x0008, 0x0050)))
            dPat[r'Date of exam'].append(sDate)
            dPat[r'Weight of patient'].append(_get_val(hdr, "PatientWeight"))

        # Get the dictionary and start adding some data
        exams[uid][suid][2] += 1
        dProt = exams[uid][suid][1]

        # Get the series data
        if (mfr == 'toshiba') and (model == 'aquilion one'):
            dProt[r'Sequence #'].append(_get_val(hdr, (0x0020, 0x0011)))
            dProt[r'Protocol name'].append(_get_val(hdr, (0x0018, 0x1030)))
            dProt[r'Requested procedure'].append(_get_val(hdr, (0x0032, 0x1060)))
            dProt[r'Scan Type'].append(sType)
            dProt[r'kVp'].append(int(_get_val(hdr, (0x0018, 0x0060))))
            dProt[r'mA'].append(float(_get_val(hdr, (0x0018, 0x1151))))
            dProt[r'Rotation Time'].append(float(_get_val(hdr, (0x0018, 0x1150))))
            dProt[r'mAs'].append(float(_get_val(hdr, (0x0018, 0x1152))))
            dProt[r'Scan FOV'].append(float(_get_val(hdr, (0x0018, 0x0090))))
            dProt[r'Display FOV'].append(float(_get_val(hdr, (0x0018, 0x1100))))
            dProt[r'Reconstruction Algorithm'].append(_get_val(hdr, (0x0018, 0x1210)))
            dProt[r'Z-axis width (T) (mm)'].append(_get_val(hdr, (0x0018, 0x9306)))
            try:
                dProt[r'Table increment (I)'].append(abs(float(_get_val(hdr, (0x0018, 0x9310)))))
            except ValueError:
                dProt[r'Table increment (I)'].append(None)
            dProt[r'Reconstructed image width (mm)'].append(float(_get_val(hdr, (0x0018, 0x0050))))
            try:
                dProt[r'Reconstructed image interval (mm)'].append(float(_get_val(hdr, (0x7005, 0x1022))))
            except ValueError:
                dProt[r'Reconstructed image interval (mm)'].append(None)
            dProt[r'Dose reduction technique'].append(_get_val(hdr, (0x0018, 0x9323)))
            dProt[r'Slice Location'].append(float(_get_val(hdr, (0x0020, 0x01041))))
            try:
                dProt[r'CTDIvol'].append(float(_get_val(hdr, (0x0018, 0x9345))))
            except ValueError:
                dProt[r'CTDIvol'].append(None)
            dProt[r'Pitch'].append(_get_val(hdr, (0x0018, 0x9311)))

            # Calculated values
            dProt[r'# of data channels (N)'] = col/dProt[r'Z-axis width (T) (mm)'][-1]

        elif ((mfr == 'siemens') and (model == 'somatom definition as') and ('va48a' in ver)) or \
             ((mfr == 'siemens') and (model == 'somatom drive') and ('vb20a' in ver)) or \
             ((mfr == 'siemens') and (model == 'somatom edge plus') and ('vb20a' in ver)) or \
             ((mfr == 'siemens') and (model == 'somatom definition edge') and ('va48a' in ver)):
            dProt[r'Sequence #'].append(_get_val(hdr, (0x0020, 0x0011)))
            dProt[r'Protocol name'].append(_get_val(hdr, (0x0018, 0x1030)))
            dProt[r'Requested procedure'].append(_get_val(hdr, (0x0032, 0x1060)))
            dProt[r'Scan Type'].append(sType)
            dProt[r'kVp'].append(int(_get_val(hdr, (0x0018, 0x0060))))
            dProt[r'mA'].append(float(_get_val(hdr, (0x0018, 0x1151))))
            dProt[r'Rotation Time'].append(float(_get_val(hdr, (0x0018, 0x1150))))
            dProt[r'Effective mAs'].append(float(_get_val(hdr, (0x0018, 0x1152))))
            dProt[r'Scan FOV'].append(float(_get_val(hdr, (0x0018, 0x0090))))
            dProt[r'Display FOV'].append(float(_get_val(hdr, (0x0018, 0x1100))))
            dProt[r'Reconstruction Algorithm'].append(_get_val(hdr, (0x0018, 0x1210)))
            dProt[r'Z-axis width (T) (mm)'].append(float(_get_val(hdr, (0x0018, 0x9306))))
            try:
                dProt[r'Table increment (I)'].append(abs(float(_get_val(hdr, (0x0018, 0x9310)))))
            except ValueError:
                dProt[r'Table increment (I)'].append(None)
            dProt[r'Reconstructed image width (mm)'].append(float(_get_val(hdr, (0x0018, 0x0050))))
            dProt[r'Dose reduction technique'].append(_get_val(hdr, (0x0018, 0x9323)))
            dProt[r'Slice Location'].append(float(_get_val(hdr, (0x0020, 0x01041))))
            try:
                dProt[r'CTDIvol'].append(float(_get_val(hdr, (0x0018, 0x9345))))
            except ValueError:
                dProt[r'CTDIvol'].append(None)
            dProt[r'Pitch'].append(_get_val(hdr, (0x0018, 0x9311)))

            # Calculated values
            dProt[r'mAs'].append(dProt[r'mA'][-1]*dProt[r'Rotation Time'][-1])

        elif (mfr == 'ge medical systems') and (model == 'revolution evo'):
            dProt[r'Sequence #'].append(_get_val(hdr, (0x0020, 0x0011)))
            sProt = _get_val(hdr, (0x0018, 0x1030))
            sProt = sProt.replace("*", "")
            dProt[r'Protocol name'].append(sProt)
            #dProt[r'Requested procedure'].append(_get_val(hdr, (0x0032, 0x1060)))
            dProt.pop(r'Requested procedure', None)
            dProt[r'Scan Type'].append(_get_val(hdr, (0x0018, 0x0022)))
            dProt[r'kVp'].append(int(_get_val(hdr, (0x0018, 0x0060))))
            dProt[r'mA'].append(float(_get_val(hdr, (0x0018, 0x1151))))
            dProt[r'Rotation Time'].append(float(_get_val(hdr, (0x0018, 0x9305))))
            #dProt[r'Effective mAs'].append(float(_get_val(hdr, (0x0018, 0x1152))))
            dProt.pop(r'Effective mAs', None)
            dProt[r'Scan FOV'].append(float(_get_val(hdr, (0x0018, 0x0090))))
            dProt[r'Display FOV'].append(float(_get_val(hdr, (0x0018, 0x1100))))
            dProt[r'Reconstruction Algorithm'].append(_get_val(hdr, (0x0018, 0x1210)))
            dProt[r'Z-axis width (T) (mm)'].append(float(_get_val(hdr, (0x0018, 0x9306))))
            if (_get_val(hdr, (0x0018, 0x0022)) == 'CINE MODE'):
                dProt[r'Table increment (I)'].append('CINE')
            else:
                dProt[r'Table increment (I)'].append(abs(float(_get_val(hdr, (0x0018, 0x9310)))))
            dProt[r'Reconstructed image width (mm)'].append(float(_get_val(hdr, (0x0018, 0x0050))))
            dProt[r'Reconstructed image interval (mm)'].append(float(_get_val(hdr, (0x0018, 0x0088))))
            dProt[r'Dose reduction technique'].append(_get_val(hdr, (0x0053, 0x1040)))
            dProt[r'Slice Location'].append(float(_get_val(hdr, (0x0020, 0x01041))))
            try:
                dProt[r'CTDIvol'].append(float(_get_val(hdr, (0x0018, 0x9345))))
            except ValueError:
                dProt[r'CTDIvol'].append(None)
            dProt[r'Pitch'].append(_get_val(hdr, (0x0018, 0x9311)))

            # Calculated values
            dProt[r'mAs'].append(dProt[r'mA'][-1]*dProt[r'Rotation Time'][-1])
        else:
            raise NotImplementedError(f"Manufacturer or model not supported. Mfr: {mfr}; Model: {model}; Software version: {ver};")
        
        # Vendor independent calculations
        dProt[r'# of data channels (N)'].append(col/dProt[r'Z-axis width (T) (mm)'][-1])

        # Getting the DLP is probably best suited from the RDSR. The images
        # from an Aquilion ONE have a DLP tag (0x7005, 0x1040), but I can't
        # figure out how to turn the binary value into something useful.

        # Store some general information about the exam
#        eStr = '; '.join([f'Acc. #: {hdr.AccessionNumber}',
#                          f'MRN: {name}',
#                          f'date: {hdr.StudyDate}',
#                          f'name: {hdr.PatientName}',
#                          f'study desc.: {hdr.StudyDescription}'])

    if args.out:
        dOut = Path.home() / "Desktop"
        for ek in exams.keys():

            # Create the Excel writer and write the data
            dfPat = pandas.DataFrame(exams[ek][list(exams[ek].keys())[0]][0])
            writer = pandas.ExcelWriter(dOut / (dfPat['Accession'][0] + '.xlsx'),
                                        engine='xlsxwriter')
            sidx = 1

            for sk in exams[ek].keys():

                # Siemens scanners don't seem to provide a means by which to determine the
                # reconstructed image interval
                if ((mfr == 'siemens') and (model == 'somatom definition as') and ('va48a' in ver)) or \
                    ((mfr == 'siemens') and (model == 'somatom drive') and ('vb20a' in ver)) or \
                    ((mfr == 'siemens') and (model == 'somatom edge plus') and ('vb20a' in ver)) or \
                    ((mfr == 'siemens') and (model == 'somatom definition edge') and ('va48a' in ver)):
                    slLocs = exams[ek][sk][1]['Slice Location']
                    slLocs.sort()
                    if (len(slLocs) > 1):
                        slDiff = list(set([round(y-x, 1) for x, y in zip(slLocs[0::], slLocs[1::])]))
                        if (len(slDiff) > 1):
                            print(f"Unable to determine reconstructed slice interval - using max...\n",
                                  f"Found {len(slDiff)} different slice interavals.")
                            exams[ek][sk][1]['Reconstructed image interval (mm)'] = [max(slDiff)] * exams[ek][sk][2]
                        else:
                            exams[ek][sk][1]['Reconstructed image interval (mm)'] = slDiff * exams[ek][sk][2]
                    else:
                        exams[ek][sk][1]['Reconstructed image interval (mm)'] = exams[ek][sk][1]['Reconstructed image width (mm)']

                for k in exams[ek][sk][1].keys():
                    print(f'{k} length: {len(exams[ek][sk][1][k])}')

                # Convert the patient info to a data frame
                dfSer = pandas.DataFrame({'Series UID': [sk],
                                          'Number of Images': [exams[ek][sk][2]]})
                dfProt = pandas.DataFrame(exams[ek][sk][1])

                # Sort by slice location
                dfProt = dfProt.sort_values(by='Slice Location')
                shName = f'{dfProt["Sequence #"][0]}--{dfProt["Protocol name"][0]}'[:31]
                shName = shName.replace('/', '')

                dfPat.to_excel(writer, shName, index=False)
                dfSer.to_excel(writer, shName, startrow=3, index=False)
                dfProt.to_excel(writer, shName, startrow=6, index=False)

                # Increment
                sidx += 1

            writer.save()


def _create_parser():

    import argparse

    parser = argparse.ArgumentParser(prog=__file__,
                                     description="CT DICOM protocol dump")
    parser.add_argument('files', metavar='F', type=str, nargs='+',
                        help=('Valid DICOM file (or directory) names '
                              'from which to generate a protocol'))
    parser.add_argument('--out', action='store_true', default=False,
                        help='Output CT acquisitions to Excel file')

    return parser


def _process_args(args, parser):

    print("")

    # Process the output flag first, as a target output file will need to be
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

    # Parser setup and input parsing
    parser = _create_parser()
    args = parser.parse_args()
    _process_args(args, parser)

    # Generate the protocols
    gen_protocol(args)
