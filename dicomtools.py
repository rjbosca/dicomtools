# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 22:21:51 2017

@author: Ryan Bosca
"""

import py
import PyQt5
import gdcm
import SimpleITK
import sys
import matplotlib
import vendormaps
import viewer


class dicomMixin():
    """Mixin class to define common DICOM sub-class attributes

    Attributes
    ----------
    file : py.path.local
        DICOM file

    """

    def __init__(self):
        self._autoLoad = True

    @staticmethod
    def computehash(file, hashtype='md5', chunksize=524288):
        """Calculate the a file's hash value.

        This is a lazy wrapper for the class attribute computehash of
        py.path.local.

        Parameters
        ----------
            file : str or py.path.local
                Full file name
            hashtype : str
                The crytographic hash function to be used
            chunksize : int
                Chunk byte size used to determine the hash

        Returns
        -------
        str
            Hashed file hex digest

        Notes
        -----
            See also py.path.local:
                http://py.readthedocs.io/en/latest/path.html

        """

        file = fileStr2PyPath(file)
        return file.computehash(hashtype=hashtype, chunksize=chunksize)

    @property
    def file(self):
        """py.path.local: DICOM file"""
        return self.__file

    @file.setter
    def file(self, val):

        val = fileStr2PyPath(val)

        # Use GDCM to attempt reading the file
        self._gdcm_reader = gdcm.Reader()
        self._gdcm_reader.SetFileName(val.strpath)
        self.isDicom = self._gdcm_reader.Read()
        self._gdcm_file = self._gdcm_reader.GetFile()

        self.__file = val

    @staticmethod
    def isdicom(file):
        """Attempt to determine whether a given file is in the DICOM format.

        Parameters
        ----------
            file : str or py.path.local
                Full file name of the file to be checked

        Returns
        -------
        bool
            True if the file is a DICOM, False otherwise.

        """

        file = fileStr2PyPath(file)

        r = gdcm.Reader()
        r.SetFileName(file.strpath)
        return r.Read()

    @staticmethod
    def get_dir_from_dicom(hdr):
        """Converts DICOM meta-data to a valid directory

        Parameters
        ----------
        header : py.path.local, str, or dicomtools.header

        Returns
        -------
        dir : py.path.local
            Full directory name generated from DICOM meta-data

        """

        # Determine what type of data the user provided, coverting when needed
        if (type(hdr) == str) or (type(hdr) == py.path.local):
            hdr = header(hdr)
        assert (type(hdr) == header)

        # Ensure that all fields necessary for converting the file information
        # to a valid directory name are present
        # TODO: in a future version of this function, there should be a way of
        #       specifying what data to use
        sNum = hdr[0x0020, 0x0011]
        sDesc = hdr[0x0008, 0x0103e]
        if not(sNum) or not(sDesc):
            print(hdr)
            print("Unable to generate directory name!")
            return ''

        sNum = ''.join([rc if rc.isalnum() else '_' for rc in str(sNum)])
        sDesc = ''.join([rc if rc.isalnum() else '_' for rc in sDesc])

        # Remove special characters
        return '--'.join([sNum, sDesc])

    @staticmethod
    def mkdir_dicom(dicomDir):
        """Rename a directory of DICOM files based on the meta-data

        Parameters
        ----------
        dicomdir : str or py.path.local
            Absolute to path of the directory to be renamed

        Returns
        -------
        success : bool
            True when the rename operation succeeds; otherwise False

        """

        dicomDir = dirStr2PyPath(dicomDir)

        # Check each file for DICOM data
        for f in dicomMixin.seek_dicom(dicomDir, gen=True):
            try:
                # Get the new directory names
                d = dicomMixin.get_dir_from_dicom(f)

                # Rename the directory, returning the success
                newPath = py.path.local(dicomDir.dirname).join(d)
                dicomDir.rename(newPath)
                return newPath.isdir()

            except (py.error.EEXIST) as error:
                #TODO: this error checking is incomplete
                if (newPath == dicomDir):
                    return False
                else:
                    continue
            except Exception as error:
                raise(error)

    @staticmethod
    def seek_dicom(dicomDir, **kwargs):
        """Recursively seek all DICOM files in a specified directory

        Parameters
        ----------
            dicomDir : str or py.path.local
                Full name of the directory to be searched for DICOM files
            fil : str
                A filter (glob pattern or callable) that, if not matching the
                path/file will not yield that path/file. Default: None
            gen : bool
                When True (default), seek_dicom returns a generator that yields
                DICOM files. Otherwise, all DICOM files are returned as a list.

        Returns
        -------
            files : list or generator
                DICOM file(s) found in the specified directory and sub-
                directories. See also the optional input gen

        Notes
        -----
            For large directory trees, this function can require a significant
            amount of time to complete

        """

        dicomDir = dirStr2PyPath(dicomDir)

        if kwargs.get('gen', False):  # generator option
            return dicomMixin._seek_dicom_gen(dicomDir, **kwargs)
        else:  # get all files
            f = []
            for pp in dicomDir.visit(fil=kwargs.get('fil', None)):
                if pp.isfile() and dicomMixin.isdicom(pp):
                    f.append(pp)
            return f

    @staticmethod
    def _seek_dicom_gen(dicomDir, **kwargs):
        """Make a generator for recursively seeking DICOM files

        Notes
        -----
            See seek_dicom

        """

        for pp in dicomDir.visit(fil=kwargs.get('fil', None)):
            if pp.isfile() and dicomMixin.isdicom(pp):
                yield pp


class gdcmMixin():

    def _gdcm_gen(self, f, ds):
        """GDCM data element generator

        Parameters
        ----------
        fcn :
            a function object that is passed a GDCM data set object
        dataset :
            GDCM data set

        """

        iterDcm = ds.GetDES().begin()
        while not iterDcm.equal(ds.GetDES().end()):
            yield f(iterDcm.next())


class dicom(dicomMixin):
    """Class for reading DICOM files

    This class supports (via the gdcm package) reading DICOM meta-data and
    image data. Currently, there is no support writing files

    Parameters
    ----------
    file : str or py.path.local
        Full file name of the DICOM file
    autoLoad : bool
        When True, the DICOM file data is imported on class instantiation.
        Otherwise, the user is responsible for performing all read operations.

    Attributes
    ----------
    file : py.path.local
        DICOM file
    header: dicomtools.header
        Meta-data reference class

    """

    def __init__(self, file, **kwargs):

        self._autoLoad = kwargs.get('autoLoad', True)
        self.header = header(file, **kwargs)

        # By setting the file, a large number of events will be performed
        self.file = file

    @property
    def image(self):
        """DICOM image array

        """

        if not hasattr(self, '_image'):
            self._image = SimpleITK.ReadImage(self.file.strpath)

        return SimpleITK.GetArrayFromImage(self._image)

    def show(self, showSeries = False):
        """Show the DICOM image(s)

        Parameters
        ----------
        showSeries : bool
            When True, any images in the same directory as the current DICOM
            object will be displayed via a viewer that allows scrolling.
            Default: False

        """

        if showSeries:

            qApp = PyQt5.QtWidgets.QApplication([''])

            v = viewer.DICOMviewer(self.file.dirpath().strpath)
            v.setWindowTitle("DICOMviewer")
            v.show()

            sys.exit(qApp.exec_())

        else:
            # TODO: there should be some validation to ensure that the file
            #       does, in fact, contain image data...
            # TODO: the below code works if the DICOM image isn't a time series
            img = SimpleITK.ReadImage(self.file.strpath)
            imgArray = SimpleITK.GetArrayFromImage(img)[0, :, :]
            matplotlib.pyplot.imshow(imgArray)


class header(dicomMixin, gdcmMixin):
    """Class containing DICOM meta-data elements as attributes

    This class allows the user to control how DICOM meta-data are imported and
    is intended to be a helper class for the dicom class

    Parameters
    ----------
    autoLoad : bool
        When True, the DICOM file data is imported upon class instantiation.
        Otherwise, the user is responsible for performing all read operations.
        Note that when using the auto-load feature the read operation imports
        the entire DICOM header (i.e., the user cannot control how the import
        is performed).
    file : str or py.path.local
        Full file name of the DICOM file
    customLookup : bool
        When True, the requested DICOM meta-data will be filtered through any
        existing mapping function in vendormaps

    """

    def __init__(self, file, **kwargs):

        # Initialize the private group/element tuple property and _is_read
        self._tag_lookup = dict()
        self._is_read = False
        self._mangle = []

        # Set the user input
        self.file = file
        self._autoLoad = kwargs.get('autoLoad', True)
        self._customLookup = kwargs.get('customLookup', False)

        # Read the file now if auto-load is enabled
        if self._autoLoad:
            self.read()

    def __getitem__(self, key):

        # Validate the input and get the group/element
        if (len(key) == 2) or (type(key[0]) != int) or (type(key[1]) != int):

            # Attempt to get the attribute name and data element
            if key in self._tag_lookup:
                attName = self._tag_lookup[key]
                val = getattr(self, attName).value
            else:
                val = self.readsingle(key).value

        else:
            strErr = ('Key must be a 2-element tuple containing 2 hexadecimal '
                      f'str or int. Was {len(key)}-element {type(key)}.')
            raise ValueError(strErr)

        return val

    def __iter__(self):
        return (getattr(self, s) for s in self._tag_lookup.values())

    def __len__(self):
        return len(self._tag_lookup.keys())

    def __repr__(self):
        return ''

    def __str__(self):
        return '\n'.join([de.__str__() for de in tuple(self)])

    def _addelement(self, de):

        assert (type(de) == dataelement)
        if not de:
            return
        # TODO: verify that the description doesn't have other special
        #       characters

        # Add this data element to the tag look-up property
        att = de.desc2attr()
        tagNew = (de.group, de.element)

        # Store the new data element. There are a couple of cases to consider:
        #   (1) The requested attribute to be added already exists but with
        #       a different group/element from the old value. In this case, add
        #       the attribute name to the "_mangle" list, update the old tag
        #       with new corresponding description/group/element mangled name,
        #       mangle the new attribute name, and update the "_tag_lookup"
        #       dictionary for both new/old values
        #   (2) The requested attribute exists already (same group/element).
        #       Raise a NotImplementedError
        #   (3) The requested attribute does not already exists, but can be
        #       found in the "_mangle" list. Mangle the attribute with the
        #       corresponding group/element and store the data element
        if hasattr(self, att):
            tagOld = (getattr(self, att).group, getattr(self, att).element)
            if (tagNew != tagOld):
                self._mangle.append(att)
                oldAtt = f'{att}_{tagOld[0]:0{4}x}_{tagOld[1]:0{4}x}'
                self._tag_lookup[tagOld] = oldAtt
                setattr(self, oldAtt, getattr(self, att))
                delattr(self, att)
            else:
                strErr = (f'Attempted adding attribute "{att}", but this '
                          f'{type(self)} instance alredy has that attribute')
                raise NotImplementedError(strErr)
        if att in self._mangle:
            att = f'{att}_{tagNew[0]:0{4}x}_{tagNew[1]:0{4}x}'

        self._tag_lookup[tagNew] = att
        setattr(self, att, de)

    def readsingle(self, tag=()):
        """Reads only the requested tag

        Similar to 'read', except 'readtag' will read only the user-requested
        tag, appending that element, in place, to the header class instance. If
        the DICOM file has already been read, that tag will be

        Parameters
        ----------
        tag : tuple
            2-element tuples representing the DICOM group and element integer
            hex values

        Returns
        -------
        data : dicomtools.dataelement
            dicomtools data element representation of the requested DICOM data

        """

        assert (type(tag) == tuple) and (len(tag) == 2) and \
               (type(tag[0]) == int) and (type(tag[1]) == int)

        # Initialize the required gdcm string filter and get the data set/tag
        sf = gdcm.StringFilter()
        sf.SetFile(self._gdcm_file)

        # Generate and append the data element
        if not self._is_read:
            gdcmTag = gdcm.Tag(tag[0], tag[1])
            ds = self._gdcm_file.GetDataSet()
            if ds.FindDataElement(gdcmTag):
                de = ds.GetDataElement(gdcmTag)
                val = dataelement(de, sf)
                self._addelement(val)

        # Generate an empty data element if needed
        if 'val' not in locals():
            val = dataelement(gdcm.DataElement(), sf, tag=tag)  # default

        return val

    def read(self):
        """Reads the specified DICOM file

        Imports the entire DICOM header, modifying the header object in place.
        All DICOM elements found in the file will be added as attributes to
        the class instance.

        Returns
        -------
        success : bool
            True if the file was read successfully

        """

        if not self.isDicom:
            return self.isDicom

        if self._tag_lookup and not self._is_read:
            for k in self._tag_lookup.values():
                if hasattr(self, k):
                    delattr(self, k)
        elif self._is_read:
            return self._is_read

        # Initialize the required gdcm string filter and get the data set
        sf = gdcm.StringFilter()
        sf.SetFile(self._gdcm_file)
        ds = self._gdcm_file.GetDataSet()

        # A note on the following code... I attempted to use a simple WHILE
        # loop to read the DICOM header, but found little success. In fact,
        # using that control structure resulted in an infinite loop. Instead, a
        # generator is created (see the gdcmMixins method _gdcm_gen), which
        # does not cause the infinite loop) and allows finer iteration control.

        # Get all available header data
        for de in self._gdcm_gen(lambda d: dataelement(d, sf), ds):
            self._addelement(de)

        self._is_read = True
        return self._is_read


class dataelement(gdcmMixin):
    """Class for aggregating a DICOM meta-data elements

    Parameters
    ----------
    data : gdcm.DataElement
        GDCM data element
    filter : gdcm.StringFilter
        GDCM string filter. The file must be already be set


    Properties
    ----------
    element : int
        DICOM data element element number
    group : int
        DICOM data element group number
    items : list
        DICOM data sequence members, each is a dataelement instance
    parent : dicomtools.dataelement
        Parent data element. This is used only for sequences. An empty data
        element will be returned for non-sequences
    VR : str
        DICOM data element value representation
    value : str, int, float, or list
        DICOM data element value type casted according to the VR property. For
        a dictionary of the VR dependent data types, see the _VRdict property

    """

    # TODO: verify 'AE', 'AT', 'OB', 'OF', 'OW', 'SQ', 'UN', 'UT'
    _VRdict = {
            'AE': str,    # Application entity
            'AS': str,    # Age string
            'AT': str,    # Attribute tag
            'CS': str,    # Code string
            'DA': str,    # Date
            'DS': float,  # Decimal string
            'DT': str,    # Date/time
            'FD': float,  # Floating point double
            'FL': float,  # Floating point single
            'IS': int,    # Integer string
            'LO': str,    # Long string
            'LT': str,    # Long text
            'OB': int,    # Other byte
            'OF': float,  # Other float
            'OW': str,    # Other word
            'PN': str,    # Person's name
            'SH': str,    # Short string
            'SL': int,    # Signed long
            'SQ': list,   # Sequence
            'SS': int,    # Signed short
            'ST': str,    # Short text
            'TM': str,    # Time
            'UI': str,    # Unique identifier
            'UL': int,    # Unsigned long
            'UN': str,    # Unknown
            'US': int,    # Unsigned short
            'UT': str,    # Unlimited text
            '??': str,    # Unknown
               }

    def __init__(self, dEl, sFil, **kwargs):

        # Validate the inputs
        assert(type(dEl) == gdcm.DataElement)
        assert(type(sFil) == gdcm.StringFilter)

        self._converted_val = None
        self.items = []

        # The following is undocumented syntax for generating a null class
        # instance, which greatly simplifies the __getitem__ syntax of the
        # header class and the dataelement class retrieval of the "parent"
        # attribute.
        tag = kwargs.get('tag', ())
        if tag:
            assert (type(tag) == tuple) and (len(tag) == 2) and \
                (type(tag[0]) == int) and (type(tag[1]) == int)
            self.description = ''
            self.group = tag[0]
            self.element = tag[1]
            self._raw_val = ''
            self.VR = '??'
            return

        # Get the group/element and value representation
        hexGrEl = str(dEl.GetTag())
        self.group = int(hexGrEl[1:5], 16)
        self.element = int(hexGrEl[6:10], 16)
        self.VR = str(dEl.GetVR()).upper()

        # Get the raw data element (description and value) and initialize the
        # type cast data element storage
        val = sFil.ToStringPair(dEl)
        self.description = val[0].strip()
        self._raw_val = val[1].strip()

        # The following code will recursively find all data elements in a data
        # sequence. This is accomplished by initializing a loop through the
        # data sequence. Each item of the sequence is read and the nested GDCM
        # generator is initializes a new data sequence based the nested data
        # set. If the data set contains another sequence, all elements of that
        # sequence will be read accordingly before the parent sequence is
        # complete

        # Handle data sequences
        if (self.VR == 'SQ'):

            # Read the data elements of the sequence
            sq = dEl.GetValueAsSQ()
            for idx in range(sq.GetNumberOfItems()):
                item = sq.GetItem(idx+1)
                for de in self._gdcm_gen(lambda d: dataelement(d, sFil),
                                         item.GetNestedDataSet()):
                    if de:
                        self.items.append(de)

            # Update the '_parent' attribute for all items
            for de in self.items:
                de._parent = self

        # Initialize the level. Be setting the "_level" attribute at the end
        # of __init__, any sequences will correctly propagate the level to all
        # children data elements.
        self._level = 0

    def __bool__(self):
        return bool(self._raw_val) or bool(self.items)

    def __repr__(self):
        if self:
            return f'<{type(self).__name__}: {self.description}>'
        else:
            return f'<Empty {type(self).__name__}>'

    def __str__(self):
        if (self.VR != 'SQ'):
            sRep = (f'({self.group:#0{6}x}, {self.element:#0{6}x}):'
                    f'<{self.VR}>: {self.description} = {self.value}')
        else:
            sRep = (f'({self.group:#0{6}x}, {self.element:#0{6}x}):'
                    f'<{self.VR}>: {self.description}:\n')
            for i in self.items:
                sRep = sRep + '\t'*i._level + i.__str__() + '\n'
            sRep = sRep[:-1]
        return sRep

    @property
    def description(self):
        return self.__description

    @description.setter
    def description(self, val):
        if not val:
            self.__description = 'Internal Data'
        else:
            self.__description = val

    @property
    def _level(self):
        return self.___level

    @_level.setter
    def _level(self, val):
        self.___level = val
        if self.VR == 'SQ':
            for de in self.items:
                de._level = val + 1

    @property
    def parent(self):
        if hasattr(self, '_parent'):
            return self._parent
        else:
            self._parent = dataelement(tag=(0x0, 0x0))

    @property
    def value(self):
        """str, float, int, list : DICOM data element value interpreted based
        on the value representation (VR attribute)"""

        # If the data has already been converted, return the _converted_val
        # attribute
        if self._converted_val:
            return self._converted_val

        # Convert the data
        val = self._raw_val
        t = self._VRdict[self.VR]
        if (t == int) or (t == float):
            try:
                val = t(val)
            except ValueError:
                if '\\' in val:
                    val = val.split('\\')
                elif ',' in val:
                    val = val.split(',')
                else:
                    strErr = (f'Unknown separator in "{self.description}" with'
                              f' tag ({self.group},{self.element}): {val}')
                    raise NotImplementedError(strErr)
                val = [t(v) for v in val]

            # Store the (list of) numeric values to avoid processing the string
            # again. Note that this does not need to be done for strings, as
            # those values are simply returned
            self.__value = val

        elif (t == str):
            val = self._raw_val.encode('ascii', errors='ignore').decode()
        elif (t == list):
            val = self.items

        # Store and return the newly converted raw data element
        return val

    def desc2attr(self):
        """Convert a data element description into a proper attribute name.

        Returns
        -------
        str
            A string containing no special characters. For those elements
            containing the word 'Private' or character '?', will be converted
            to 'X_<group>_<element>' where X is the description mangled with
            the DICOM group and element hexadecimal strings.

        """

        s = self.description
        if (s.lower() == 'private creator') or \
                ('private' in s[0:7].lower()):
            s = (''.join([rc for rc in s if rc.isalnum()]) +
                 f'_{self.group:0{4}x}_{self.element:0{4}x}')
        elif not s or ('?' in s.lower()):
            s = f'Unknown_{self.group:0{4}x}_{self.element:0{4}x}'
        else:
            if '(' in s and ')' in s:
                s = s[:s.find('(')]
            s = ''.join([rc for rc in s if rc.isalnum()])

        return s


# -----------------------------------------------------------------------------
# --------------------------     Helper functions     -------------------------
# -----------------------------------------------------------------------------


def fileStr2PyPath(f):
    """Helper function to validate and convert file strings

    Parameters:
    -----------
        f (str or py.path.local): file to be validated/converted

    Returns:
    --------
        py.path.local
            File object

    """

    f = py.path.local(f)
    assert f.isfile()
    return f


def dirStr2PyPath(d):
    """Helper function to validate and convert directory strings

    Parameters:
    -----------
        d (str or py.path.local): file to be validated/converted

    Returns:
    --------
        py.path.local
            File object

    """

    d = py.path.local(d)
    assert d.isdir()
    return d
