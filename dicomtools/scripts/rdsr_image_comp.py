# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:40:01 2022

@author: Ryan
"""


# Standard imports
import warnings
from pathlib import Path

import pandas  # used for creating output Excel
import pydicom   # DICOM interpreter


# Notes:
#--------
# Minimal error checking is performed. The following assumptions are made
#   - A single parent directory (variable: dRdsr)
#   - Parent directory contains only a single RDSR 
#   - All images are contained in the same parent directory
#   - "__ign__" directory will be ignored
#
#
#   Script Variables
# --------------------
# dRdsr: parent directory object containing the images/RDSR
# fOut: output Excel spreadsheet
# lMap: metadata-to-RDSR mapping
# tIgn: metadata elements to ignore


#--------------
# Script setup
#--------------

# Path definition as described above
d = Path(__file__).parent
dRdsr = d / "Sanford" / "20220301 - XRAY FOREAREM 2 VIEWS LT" / "DICOM"

# Output file definition
if (dRdsr.parts[-1] == "DICOM"):  # Sectra directory modifier
    fOut = dRdsr.parent / f"RDSR Comparison - {dRdsr.parts[-2]}.xlsx"
else:
    fOut = dRdsr.parent / f"RDSR Comparison - {dRdsr.parts[-1]}.xlsx"

# Metadata -> RDSR mappings. This is a list of tuples of the form:
#   ((group, element), RDSR code value)
# Any entry with a value of None is ignored in the lookup.
#TODO: include concept name in addition to concept code
lMap = [(None,             "113721"),  # Irradiation event type
        (None,             "123014"),  # Target region
        ((0x0018, 0x1030), "125203"),  # Acquisition protocol
        ((0x0011, 0x1044), "121106"),  # Table
        ((0x0018, 0x7060), None    ),  # Exposure control mode
        ((0x0018, 0x7062), None    ),  # Exposure control mode description
        ((0x0020, 0x0060), "111027"),  # Laterality
        ((0x0018, 0x5101), "111031"),  # Image view
        ((0x0018, 0x0060), "113733"),  # kVp
        ((0x0040, 0x0314), "111634"),  # Half value layer
        ((0x0018, 0x1151), "113734"),  # X-ray tube current
        ((0x0018, 0x1150), "113824"),  # Exposure time
        ((0x0018, 0x1153), None    ),  # Exposure
        ((0x0014, 0x4022), None    ),  # Pulse width
        ((0x0018, 0x1166), "111635"),  # Grid
        ((0x0018, 0x7040), None    ),  # Grid absorbing material
        ((0x0018, 0x7041), None    ),  # Grid spacing material
        ((0x0018, 0x7042), None    ),  # Grid thickness
        ((0x0018, 0x7046), None    ),  # Grid aspect ratio
        ((0x0018, 0x704c), None    ),  # Grid focal distance
        ((0x0018, 0x7050), "113757"),  # Added filtration material
        ((0x0018, 0x7052), "113758"),  # Filter thickness minimum
        ((0x0018, 0x7054), "113773"),  # Filter thickness maximum
        ((0x0018, 0x1110), "113750"),  # Distance source to detector
        ((0x0018, 0x1411), "113845"),  # Exposure index
        ((0x0018, 0x1412), "113846"),  # Target exposure indicator
        ((0x0018, 0x1413), "113847"),  # Deviation index
        ((0x0018, 0x115e), "122130"),  # Image area dose product
        ((0x0040, 0x8302), None    ),  # Entrance dose
        (None,             "113780"),  # Reference point definition
        (None,             "111636"),  # Entrance exposure at RP
        (None,             "113738"),  # Dose
        ((0x0018, 0x1149), "113788"),  # Field of view dimension(s)
        ((0x0018, 0x1149), "113789"),  # Field of view dimension(s)
        ((0x0018, 0x1147), None    ),  # Field of view shape
        ((0x0018, 0x7004), "113947"),  # Detector type
        ]

# DICOM elements to be ignored. Care should be taken in defining tags, as
# these elements are ignored when processing the DICOM metadata and RDSR
tIgn = [(0x0008, 0x0005),  # Specific character set
        (0x0008, 0x0016),  # SOP class UID
        (0x0008, 0x0018),  # SOP instance UID
        (0x0008, 0x0021),  # Series date
        (0x0008, 0x0023),  # Content date
        (0x0008, 0x0031),  # Series time
        (0x0008, 0x0033),  # Content time
        (0x0008, 0x0060),  # Modality
        (0x0008, 0x0090),  # Referring physician's name
        (0x0008, 0x103e),  # Series description
        (0x0010, 0x0010),  # Patient's name
        (0x0010, 0x0020),  # Patient ID
        (0x0020, 0x000d),  # Study instance UID
        (0x0020, 0x000e),  # Series instance UID
        (0x0020, 0x0010),  # Study ID
        (0x0020, 0x0011),  # Series number
        (0x0020, 0x0013),  # Instance number
        (0x0028, 0x3006),  # LUT data
        (0x6000, 0x3000),  # Overlay data
        (0x7fe0, 0x0010),  # Pixel data
        ]


#------------------
# Helper Functions
#------------------

def _decode_ds(ds, lde, ign=[]):
    """
    _decode_ds  Generate a list of data element info

    Parameters
    ----------
    ds : pydicom.dataset.FileDataset
        Dataset to convert into a list.
    lde : List
        Decoded data element list. The mutable list is not returned, but
        decoded data elements are appended to the list.
    ign : List
        List of tuples (group, element) to ignore

    Returns
    -------
    None. Appends formatted data element to lde.

    """


    def _decode_de(de, l):

        l.append([str(de.tag), str(de.name), str(de.VR), str(de.value)])


    for de in ds:

        if (de.tag in ign):
            continue

        _decode_de(de, lde)

        if (de.VR == "SQ"):
            lde[-1][-1] = f"{len(de.value)} item(s) ----"
            for d in de:
                _decode_ds(d, lde, ign)
            lde.append(["----"]*4)


def _find_container(ds, cv):
    """
    _find_cv  Find the dataset with a specific code value

    Parameters
    ----------
    ds : pydicom.dataset.Dataset
        Dataset in which to search for the specified code value.
    val : str
        Code value string.

    Returns
    -------
    None if code value not found. Otherwise, the dataset is returned.

    """

    # Return the input dataset if the code value exsits
    if (0x0008, 0x0100) in ds and (ds[0x0008, 0x0100].value == cv):
        return ds

    # Since the code value isn't in the provided dataset, search all data
    # elements that are sequences
    for de in ds:

        if (de.VR != "SQ"):
            continue

        # Search sequences. Any dataset containing the code value tag is
        # returned
        for dss in de:
            if _find_container(dss, cv) is not None:
                return dss


def _get_irrad_event(rdsr, uid):
    """
    _get_irrad_event  Gets an irradiation event dataset within an RDSR

    Parameters
    ----------
    rdsr : pydicom.dataset.FileDataset
        RDSR file imported using pydicom.
    uid : pydicom.uid.UID
        SOP instance UID that corresponds to a specific irradiation event.

    Returns
    -------
    None if irradiation event is not found within the RDSR. Otherwise, the
    dataset corresponding to that irradiation event is returned.

    """

    for ds in rdsr.ContentSequence:

        # Search all datasets in the content sequence. Those that have the
        # Irradiation Event X-ray Data code value will also have the referenced
        # SOP instance UID.
        dss = _find_container(ds, "113706")
        if dss is not None:  # Irradiation Event X-Ray Data

            # Find the Acquired Image SOP instance UID and compare
            dss = _find_container(ds, "113769")
            uidIrrad = dss.UID
            if (uidIrrad == uid):
                # Dataset found
                return ds


def _text_to_dicom(f):
    """
    _text_to_dicom  Decode a text file DICOM metadata dump

    Parameters
    ----------
    f : pathlib.Path
        Text file to interpret.

    Returns
    -------
    DICOM dataset.

    """

    def _parse_line(s):
    
        sDelim = [ls.strip() for ls in s.split("|")]
        if (len(sDelim) < 4):
            # Assume value output as string
            return (None, None, None, None, s)
        elif (len(sDelim) > 4):
            #TODO: handle case for value with delimiter
            pass
    
        # Complete the parse operation (VM is unused currently)
        (g, el, ln) = sDelim[0].split()
        vr = sDelim[2]
        val = sDelim[4]
    
        # Strip double quotes and check for truncation
        if (len(val) > 0):
            if (val[0] == '"'):
                val = val[1:]
            if (val[-1] == '"'):
                val = val[:-1]
            if (val[-3:] == "..."):
                val = None
                ln = None
    
        # Assume the group/element are in an appropriate format
        g = int(g, 16)
        el = int(el, 16)
    
        # Check for undefined length
        if (ln == "(undef.)"):
            ln = None
        elif (ln is not None):
            ln = int(ln)
    
        return g, el, ln, vr, val
    
    
    def _read_multiline_de(lFile, idx, val):
    
        for lf in lFile[idx:]:
    
            # Count delimiters
            nDelim =lf.count("|")
            if (nDelim == 4):
                if (val[-1] == '"'):
                    val = val[:-1]
                if (val[-3:] == "..."):
                    val = None
                return (idx, val)
            elif (nDelim > 4):
                raise NotImplementedError("BOO BAD!")
    
            # Combine the values
            if not len(lf):
                lf = "\n"
            val += lf
            idx += 1
    
    
    def _read_sq(lFile, idx):
    
        # Define sequence end keywords
        lSeqEnd = ["**** End sq", "---- End Seq"]
    
        # Generate a blank dataset
        dsSq = pydicom.Dataset()
    
        # Parse lines after current index
        lDs = []
        idxSq = 0
        for lf in lFile[idx + 1:]:
    
            # Increment the index
            idx += 1
    
            # Skip beginning line of a sequence or if a nested sequence was
            # read. For end item lines, store the current sequence dataset
            # , create a  blank dataset, and move to the next item
            if ("Item:" in lf) or (idxSq > idx):
                continue
            elif ("End Item" in lf):
                lDs.append(dsSq)
                dsSq = pydicom.Dataset()
                continue
            
            # Sequence termination
            if (lf in lSeqEnd):
                return (idx + 1, pydicom.sequence.Sequence(lDs))
    
            # Parse the line
            g, el, ln, vr, val = _parse_line(lf)
            if (vr == "SQ"):
                idxSq, val = _read_sq(lFile, idx + 1)
    
            #FIXME: error checking is needed...
    
            # Store the data element
            dsSq[g, el] = pydicom.DataElement((g, el), vr, val, ln)
    
    
    # Initialize the dataset
    ds = pydicom.Dataset()
    
    # Read the file
    l = []
    with open(f) as file:
        for line in file:
            l.append(line.strip())
    
    # Process the text
    val = None
    idxMl = -1  # multi-line read index
    idxSq = -1  # sequence read index
    g = None  # DICOM tag group value
    el = None  # DICOM element group value
    lDecode = l[4:]  # [4:] skips header
    for idx, line in enumerate(lDecode):
    
        # End of header reached
        isEoh = (len(line) > 13) and (len(set(line)) == 1)
        if isEoh:
            break
    
        # Continue for multi-line reads and sequence reads
        if (idxMl > idx):  # ignore multi-line reads
            continue
        else:  # reset index
            idxMl = 0
        if (idxSq > idx):  # ignore sequence lines
            continue
        else:  # reset index
            idxSq = 0
    
        # Count the number of delimiters in the line. If 0, the following code
        # reads a multi-line value. If more than 4, the delimiter is used in
        # the value
        nDelim = line.count("|")
        if (nDelim == 0):
            idxMl, val = _read_multiline_de(lDecode, idx, val)
            ds[g, el].value = val
            continue
    
        # Parse the DICOM group, DICOM element, length, value representation,
        # and value from the text line
        g, el, ln, vr, val = _parse_line(line)
    
        # Sequences require special processing
        if (vr == "SQ"):
            if (ln == 0):
                val = pydicom.sequence.Sequence([])
            else:
                idxSq, val = _read_sq(lDecode, idx)
    
        #TODO: is this needed???
        # Perform dictionary lookup
    
        # Create a data element
        #TODO: handle undefined length case...
        if ln is None:
            ds[g, el] = pydicom.DataElement((g, el), vr, val)
        else:
            ds[g, el] = pydicom.DataElement((g, el), vr, val, ln)

    return ds


#----------------------------
# File processing/validation
#----------------------------

# Validate the RDSR directory
if not dRdsr.is_dir():
    raise ImportError(f"Unable to locate the directory: '{dRdsr}'")

# Read all DICOM files in the exam directory
fDcms = [f for f in dRdsr.rglob("*")
         if f.is_file() and "__ign__" not in str(f)]
mds = []
for f in fDcms:
    ds = pydicom.dcmread(f, force=True)
    if not ds.file_meta:
        ds = _text_to_dicom(f)
    mds.append(ds)

# Ensure all files belong to the same study
uidStudy = [md.StudyInstanceUID for md in mds]
nSt = len(set(uidStudy))
if (nSt != 1):
    warnings.warn(f"Found {nSt} study UID(s). Expected only 1. Attempt to "
                  "process anyway...")

# Ensure exactly one RDSR is in the path
m = [md.Modality for md in mds]
nRdsr = m.count("SR")
if (nRdsr == 0):
    raise ImportError(f"No RDSR found in : '{dRdsr}'. Excepted exactly 1.")
elif (nRdsr > 1):
    sErr = [f"Found {nRdsr} RDSRs in: {dRdsr}",
            "Remove all except one of the following:"]
    _ = [sErr.append("\t" + str(fDcms[idx])) for idx, m in enumerate(m)
         if m == "SR"]
    raise ImportError("\n".join(sErr))

del nRdsr, nSt, m  # unused variables


#-----------------------------
# Process Image/RDSR Metadata
#-----------------------------

# List of dataframes and respective sheet names that will be written to the
# output Excel file
dfs = []

# Reformat the raw content
for idx, md in enumerate(mds):

    # RDSR index
    if (md.Modality == "SR"):
        iRdsr = idx

    # Initialize output list that will be converted to a data frame
    to_df = [["Tag", "Name", "VR", "Value"]]

    # Format the output to a pandas data frame
    _decode_ds(md, to_df, tIgn)

    dfs.append(pandas.DataFrame(to_df[1::], columns=to_df[0]))

# Generate the comparison
shName = len(dfs) * [None]
rdsr = mds[iRdsr]
for idx, md in enumerate(mds):

    to_df = [["DICOM Tag", "DICOM Data Element", "Value",
              "RDSR Tag/Code Value", "RDSR Data Element/Code Meaning", "Value"]]

    # Ignore the RDSR
    if (idx == iRdsr):
        shName[idx] = "RDSR"
        continue

    # Store the sheet name
    shName[idx] = f"{md.SeriesNumber}-{md.ViewPosition} {md.BodyPartExamined}"

    # Find the irradiation event (Content Sequence) that corresponds to the
    # current image
    ds = _get_irrad_event(rdsr, md.IrradiationEventUID)
    if (ds is None):
        warnings.warn("Unable to find irradiation event corresponding "
                      f"to image instance: {md.SOPInstanceUID}")
        continue

    # Append the comparison sheet
    shName.append(f"Comp-{shName[idx]}")

    # Get common data elements
    for de in md:
        if (de.tag in rdsr) and (de.tag not in tIgn):
            t = de.tag
            try:
                if (de.VR == "SQ"):
                    to_df.append([str(de.tag), 
                                  de.name, 
                                  de[0].CodeMeaning,
                                  str(rdsr[t].tag), 
                                  rdsr[t].name, 
                                  rdsr[t][0].CodeMeaning])
                else:
                    to_df.append([str(de.tag), 
                                  de.name, 
                                  str(de.value),
                                  str(rdsr[t].tag), 
                                  rdsr[t].name, 
                                  str(rdsr[t].value)])
            except Exception:
                continue

    # Map the DICOM tags and RDSR code values. If the tag or code value does
    # not exist in the metadata or RDSR, respectively,
    for lm in lMap:

        # Initialize output
        to_df.append(6 * [None])

        # DICOM image metadata mapping
        if (lm[0] is None):
            to_df[-1][:3] = 3 * ["No mapping"]
        elif (lm[0] in md):
            to_df[-1][:3] = [str(md[lm[0]].tag),
                             md[lm[0]].name,
                             str(md[lm[0]].value)]
        else:
            try:
                nm = pydicom.datadict.get_entry(lm[0])[2]
            except KeyError:
                nm = "???"
            to_df[-1][:3] = [f"({lm[0][0]:04x}, {lm[0][1]:04x})", nm, None]

        # RDSR code value mapping
        if (lm[1] is None):
            to_df[-1][3:] = 3 * ["No Mapping"]
        else:
            de = _find_container(ds, lm[1])
            if (de is not None) and (de.ValueType == "CONTAINER"):
                de = _find_container(de, lm[1])
            if (de is not None):
                if (de.ValueType == "TEXT"):
                    to_df[-1][3:] = [lm[1],
                                     de[0x0040, 0xa043][0].CodeMeaning,
                                     de.TextValue]
                elif (de.ValueType == "CODE"):
                    to_df[-1][3:] = [lm[1],
                                     de[0x0040, 0xa043][0].CodeMeaning,
                                     de[0x0040, 0xa168][0].CodeMeaning]
                elif (de.ValueType == "NUM"):
                    to_df[-1][3:] = [lm[1],
                                     de[0x0040, 0xa043][0].CodeMeaning,
                                     str(de[0x0040, 0xa300][0].NumericValue)]
                else:
                    raise ValueError(f"Unexpected value type: {de.ValueType}")
            else:
                to_df[-1][3:] = [lm[1], None, None]

    # Add the list to dataframe list
    dfs.append(pandas.DataFrame(to_df[1::], columns=to_df[0]))


#-----------------------------
# Write Image/RDSR Comparison
#-----------------------------

# Write to Excel
with pandas.ExcelWriter(fOut, mode="w") as writer:
    for idx, df in enumerate(dfs):
        df.to_excel(writer, sheet_name=shName[idx], index=False)