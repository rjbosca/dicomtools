6# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 08:39:36 2020

@author: 703355681
"""

from pathlib import Path
import pydicom


def get_slice(args):
    
    # Update the per-frame functional group to contain only the slice
    # of interest
    args.dcm[0x5200, 0x9230].value = \
        [args.dcm[0x5200, 0x9230][args.slice-1]]

    # Get the array for the slice of interest
    arrIm = args.dcm.pixel_array
    arrIm = arrIm[args.slice-1,]

    # Update the numnber of frames
    args.dcm.NumberOfFrames = 1

    # Write the array back to pixel data
    args.dcm.PixelData = arrIm.tobytes()

    # Before writing the image, verify that the labeling information is correct.
    # The "label" should contain City, State, and Zip
    print("Found the following institution information:")
    print(f"{args.dcm.InstitutionName}\n{args.dcm.InstitutionAddress}")

    # Write the new DICOM file
    args.dcm.save_as(args.file_out)


def _create_parser():

    import argparse

    parser = argparse.ArgumentParser(prog=__file__,
                                     description="Extract slice from DBT")
    parser.add_argument(
                    'file', type=str, nargs=1,
                    help='DICOM file from which to extract a slice'
                    )
    parser.add_argument('slice', help='Slice number to extract', type=int)

    return parser


def _process_args(args, parser):

    from pydicom.misc import is_dicom

    # Sanity check on input slice number
    if (args.slice < 1):
        raise ValueError(f"Invalid slice number: {args.slice}")

    # Validate the input file exists
    fl = Path(args.file[0])
    if not fl.is_file():  # evaluate single file
        raise ValueError("Valid file must be provided...")
    args.file = fl

    # Validate the output file does not exist
    f_name = fl.name.replace(fl.suffix, '') + f"_sl_{args.slice}" + fl.suffix
    args.file_out = fl.parent / f_name
    if args.file_out.is_file():
        raise ValueError(f"Output file already exists: {args.file_out}")

    # Load the DICOM and validate the SOP class UID and number of slices
    if not is_dicom(fl):
        raise ValueError(f"Invalid DICOM file: {fl}")
    args.dcm = pydicom.dcmread(fl)
    if (args.dcm.SOPClassUID != '1.2.840.10008.5.1.4.1.1.13.1.3'):
        raise ValueError(
            f"DICOM must be of class Breast Tomosynthesis Image Storage"
            )
    n = int(args.dcm.NumberOfFrames)
    if (args.slice > n):
        raise NotImplementedError(f"Invalid slice: {args.slice}\n"
                                  f"Slice must be between 1 and {n}")


if __name__ == "__main__":

    parser = _create_parser()
    args = parser.parse_args()
    _process_args(args, parser)
    get_slice(args)
