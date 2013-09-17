# IDLSave - a python module to read IDL 'save' files
# Copyright (C) 2009 Thomas P. Robitaille
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import struct
import numpy as np

dtype_dict = {}
dtype_dict[1] = '>u1'
dtype_dict[2] = '>i2'
dtype_dict[3] = '>i4'
dtype_dict[4] = '>f4'
dtype_dict[5] = '>f8'
dtype_dict[7] = '|O'
dtype_dict[8] = '|O'
dtype_dict[12] = '>u2'
dtype_dict[13] = '>u4'
dtype_dict[14] = '>i8'
dtype_dict[15] = '>u8'


def align_32(f):
    pos = f.tell()
    if pos % 4 <> 0:
        f.seek(pos + 4 - pos % 4)


def skip_bytes(f, n):
    f.read(n)
    return


def read_bytes(f, n):
    return f.read(n)


def read_byte(f):
    return struct.unpack('>B', f.read(4)[0])[0]


def read_long(f):
    return struct.unpack('>l', f.read(4))[0]


def read_int16(f):
    return struct.unpack('>h', f.read(4)[2:4])[0]


def read_int32(f):
    return struct.unpack('>i', f.read(4))[0]


def read_int64(f):
    return struct.unpack('>q', f.read(8))[0]


def read_uint16(f):
    return struct.unpack('>H', f.read(4)[2:4])[0]


def read_uint32(f):
    return struct.unpack('>I', f.read(4))[0]


def read_uint64(f):
    return struct.unpack('>Q', f.read(8))[0]


def read_float32(f):
    return struct.unpack('>f', f.read(4))[0]


def read_float64(f):
    return struct.unpack('>d', f.read(8))[0]


class Pointer(object):

    def __init__(self, index):
        self.index = index
        return

rectype_dict = {}
rectype_dict[0] = "START_MARKER"
rectype_dict[1] = "COMMON_VARIABLE"
rectype_dict[2] = "VARIABLE"
rectype_dict[3] = "SYSTEM_VARIABLE"
rectype_dict[6] = "END_MARKER"
rectype_dict[10] = "TIMESTAMP"
rectype_dict[12] = "COMPILED"
rectype_dict[13] = "IDENTIFICATION"
rectype_dict[14] = "VERSION"
rectype_dict[15] = "HEAP_HEADER"
rectype_dict[16] = "HEAP_DATA"
rectype_dict[17] = "PROMOTE64"
rectype_dict[19] = "NOTICE"

struct_dict = {}


def read_string(f):
    length = read_long(f)
    if length > 0:
        chars = read_bytes(f, length)
        align_32(f)
    else:
        chars = None
    return chars


def read_string_data(f):
    '''Read a string'''
    length = read_long(f)
    if length > 0:
        length = read_long(f)
        string = read_bytes(f, length)
        align_32(f)
    else:
        string = None
    return string


def read_data(f, dtype):
    if dtype==1:
        if read_int32(f) <> 1:
            raise Exception("Error occured while reading byte variable")
        return read_byte(f)
    elif dtype==2:
        return read_int16(f)
    elif dtype==3:
        return read_int32(f)
    elif dtype==4:
        return read_float32(f)
    elif dtype==5:
        return read_float64(f)
    elif dtype==6:
        raise Exception("32-bit complex type not implemented")
    elif dtype==7:
        return read_string_data(f)
    elif dtype==8:
        raise Exception("Should not be here - please report this")
    elif dtype==9:
        raise Exception("32-bit complex type not implemented")
    elif dtype==10:
        return Pointer(read_int32(f))
    elif dtype==11:
        raise Exception("Object reference type not implemented")
    elif dtype==12:
        return read_uint16(f)
    elif dtype==13:
        return read_uint32(f)
    elif dtype==14:
        return read_int64(f)
    elif dtype==15:
        return read_uint64(f)
    else:
        raise Exception("Unknown IDL type: %i - please report this" % dtype)


def read_structure(f, array_desc, struct_desc):

    nrows = array_desc.nelements
    ncols = struct_desc.ntags
    columns = struct_desc.tagtable

    dtype = []
    for col in columns:
        if col.structure or col.array:
            dtype.append(((col.name.lower(), col.name), np.object_))
        else:
            if col.typecode in dtype_dict:
                dtype.append(((col.name.lower(), col.name),
                                    dtype_dict[col.typecode]))
            else:
                raise Exception("Variable type %i not implemented" %
                                                            col.typecode)

    structure = np.recarray((nrows, ), dtype=dtype)

    for i in range(nrows):
        for col in columns:
            dtype = col.typecode
            if col.structure:
                structure[col.name][i] = read_structure(f, \
                                            struct_desc.arrtable[col.name], \
                                            struct_desc.structtable[col.name])
            elif col.array:
                structure[col.name][i] = read_array(f, dtype, \
                                            struct_desc.arrtable[col.name])
            else:
                structure[col.name][i] = read_data(f, dtype)

    return structure


def read_array(f, typecode, array_desc):

    if typecode in [1, 3, 4, 5, 13, 14, 15]:

        if typecode == 1:
            nbytes = read_int32(f)
            if nbytes <> array_desc.nbytes:
                raise Exception("Error occured while reading byte array")

        # Read bytes as numpy array
        array = np.fromstring(f.read(array_desc.nbytes), \
                                dtype=dtype_dict[typecode])

    elif typecode in [2, 12]:

        # These are 2 byte types, need to skip every two as they are not packed

        array = np.fromstring(f.read(array_desc.nbytes*2), \
                                dtype=dtype_dict[typecode])[1::2]

    else:

        # Read bytes into list
        array = []
        for i in range(array_desc.nelements):
            dtype = typecode
            data = read_data(f, dtype)
            array.append(data)

        array = np.array(array, dtype=np.object_)

    # Reshape array if needed
    if array_desc.ndims > 1:
        dims = array_desc.dims[:array_desc.ndims]
        dims.reverse()
        array = array.reshape(dims)

    # Go to next alignment position
    align_32(f)

    return array


class Record(object):

    def __init__(self):
        self.end = False
        pass

    def read(self, f):

        self.recpos = f.tell()
        self.rectype = read_long(f)

        self.nextrec = read_uint32(f)

        if self.nextrec < self.recpos:
            self.nextrec = self.nextrec + read_uint32(f) * 2**32
        else:
            skip_bytes(f, 4)

        skip_bytes(f, 4)

        if not self.rectype in rectype_dict:
            raise Exception("Unknown RECTYPE: %i" % self.rectype)

        self.rectype = rectype_dict[self.rectype]

        if self.rectype in ["VARIABLE", "HEAP_DATA"]:

            if self.rectype == "VARIABLE":
                self.varname = read_string(f)
            else:
                self.heap_index = read_long(f)
                skip_bytes(f, 4)

            self.rectypedesc = TypeDesc().read(f)

            varstart = read_long(f)
            if varstart <> 7:
                raise Exception("VARSTART is not 7")

            if self.rectypedesc.structure:
                self.data = read_structure(f, self.rectypedesc.array_desc, \
                                            self.rectypedesc.struct_desc)
            elif self.rectypedesc.array:
                self.data = read_array(f, self.rectypedesc.typecode, \
                                            self.rectypedesc.array_desc)
            else:
                dtype = self.rectypedesc.typecode
                self.data = read_data(f, dtype)

        elif self.rectype == "TIMESTAMP":

            skip_bytes(f, 4*256)
            self.date = read_string(f)
            self.user = read_string(f)
            self.host = read_string(f)

        elif self.rectype == "VERSION":

            self.format = read_long(f)
            self.arch = read_string(f)
            self.os = read_string(f)
            self.release = read_string(f)

        elif self.rectype == "IDENTIFICATON":

            self.author = read_string(f)
            self.title = read_string(f)
            self.idcode = read_string(f)

        elif self.rectype == "NOTICE":

            self.notice = read_string(f)

        elif self.rectype == "HEAP_HEADER":

            self.nvalues = read_long(f)
            self.indices = []
            for i in range(self.nvalues):
                self.indices.append(read_long(f))

        elif self.rectype == "COMMONBLOCK":

            self.nvars = read_long(f)
            self.name = read_string(f)
            self.varnames = []
            for i in range(self.nvars):
                self.varnames.append(read_string(f))

        elif self.rectype == "END_MARKER":

            self.end = True

        elif self.rectype == "UNKNOWN":

            print "Skipping UNKNOWN record"

        elif self.rectype == "SYSTEM_VARIABLE":

            print "Skipping SYSTEM_VARIABLE record"

        else:

            raise Exception("RECTYPE=%s not implemented" % self.rectype)

        f.seek(self.nextrec)

        return self

    def __str__(self):

        string = ""

        if self.rectype == "VARIABLE":

            string += "Name: %s\n" % self.varname
            string += "Type: %i" % self.typedesc.typecode

        elif self.rectype == "TIMESTAMP":

            string += "Date: %s\n" % self.date
            string += "User: %s\n" % self.user
            string += "Host: %s" % self.host

        elif self.rectype == "VERSION":

            string += "Format: %s\n" % self.format
            string += "Architecture: %s\n" % self.arch
            string += "Operating System: %s\n" % self.os
            string += "IDL Version: %s" % self.release

        elif self.rectype == "IDENTIFICATON":

            string += "Author: %s\n" % self.author
            string += "Title: %s\n" % self.title
            string += "ID Code: %s" % self.idcode

        elif self.rectype == "NOTICE":

            string += self.notice

        elif self.rectype == "HEAP_HEADER":

            string += "todo"

        elif self.rectype == "COMMONBLOCK":

            string += "todo"

        return string

    def __repr__(self):
        return "<IDL " + self.rectype + " object>"


class TypeDesc(object):

    def __init__(self):
        pass

    def read(self, f):

        self.typecode = read_long(f)
        self.varflags = read_long(f)

        if self.varflags & 2 == 2:
            raise Exception("System variables not implemented")

        self.array = self.varflags & 4 == 4
        self.structure = self.varflags & 32 == 32

        # CHECK VARFLAGS HERE TO SEE IF ARRAY

        if self.structure:
            self.array_desc = ArrayDesc().read(f)
            self.struct_desc = StructDesc().read(f)
        elif self.array:
            self.array_desc = ArrayDesc().read(f)

        return self


class ArrayDesc(object):

    def __init__(self):
        pass

    def read(self, f):

        self.arrstart = read_long(f)

        skip_bytes(f, 4)

        self.nbytes = read_long(f)
        self.nelements = read_long(f)
        self.ndims = read_long(f)

        skip_bytes(f, 8)

        self.nmax = read_long(f)

        self.dims = []
        for d in range(self.nmax):
            self.dims.append(read_long(f))

        return self


class StructDesc(object):

    def __init__(self):
        pass

    def read(self, f):

        structstart = read_long(f)
        if structstart <> 9:
            raise Exception("STRUCTSTART should be 9")

        self.name = read_string(f)
        self.predef = read_long(f)
        self.ntags = read_long(f)
        self.nbytes = read_long(f)

        if self.predef & 1 == 0:

            self.tagtable = []
            for t in range(self.ntags):
                self.tagtable.append(TagDesc().read(f))

            for i in range(self.ntags):
                self.tagtable[i].name = read_string(f)

            self.arrtable = {}
            for t in self.tagtable:
                if t.array:
                    self.arrtable[t.name] = ArrayDesc().read(f)

            self.structtable = {}
            for t in self.tagtable:
                if t.structure:
                    self.structtable[t.name] = StructDesc().read(f)

            struct_dict[self.name] = (self.tagtable, \
                                        self.arrtable, self.structtable)

        else:

            if not self.name in struct_dict:
                raise Exception("PREDEF=1 but can't find definition")

            self.tagtable, self.arrtable, self.structtable = \
                                                        struct_dict[self.name]

        return self


class TagDesc(object):

    def __init__(self):
        pass

    def read(self, f):

        self.offset = read_long(f)
        self.typecode = read_long(f)
        tagflags = read_long(f)

        self.array = tagflags & 4 == 4
        self.structure = tagflags & 32 == 32
        self.scalar = self.typecode in dtype_dict
        # Assume '10'x is scalar

        return self


class IDLSaveFile(object):

    def __init__(self, filename):
        self.records = None
        self.variables = None
        self.filename = filename

    def parse(self, verbose=True):

        self.records = []
        self.variables = {}

        f = file(self.filename, 'rb')

        signature = read_bytes(f, 2)
        if signature <> 'SR':
            raise Exception("Invalid SIGNATURE: %s" % signature)

        recfmt = read_bytes(f, 2)
        if recfmt <> '\x00\x04':
            raise Exception("Invalid RECFMT: %s" % recfmt)

        rectypes = []
        while True:
            r = Record().read(f)
            self.records.append(r)
            rectypes.append(r.rectype)
            if r.end:
                break

        f.close()

        # Find heap data variables
        heap = {}
        for r in self.records:
            if r.rectype == "HEAP_DATA":
                heap[r.heap_index] = r.data

        for r in self.records:
            if r.rectype == "VARIABLE":
                while isinstance(r.data, Pointer):
                    r.data = heap[r.data.index]
                self.variables[r.varname.lower()] = r.data

        if verbose:

            for header in ['TIMESTAMP', 'VERSION', 'IDENTIFICATION']:
                if header in rectypes:
                    print "-"*50
                    pos = rectypes.index(header)
                    print self.records[pos]

            print "-"*50
            print "Successfully read %i records of which:" % \
                                                (len(self.records))
            for rt in set(rectypes):
                if rt <> 'END_MARKER':
                    print " - %i are of type %s" % (rectypes.count(rt), rt)
            print "-"*50

            if 'VARIABLE' in rectypes:
                print "Available variables:"
                for var in self.variables:
                    print " - %s [%s]" % (var, type(self.variables[var]))
                print "-"*50

    def __call__(self, key):
        if key.lower() in self.variables:
            return self.variables[key.lower()]
        else:
            raise Exception("No such variable: %s" % key)

    def __getitem__(self, key):
        if key.lower() in self.variables:
            return self.variables[key.lower()]
        else:
            raise Exception("No such variable: %s" % key)

    def __getattr__(self, key):
        if key.lower() in self.variables:
            return self.variables[key.lower()]
        else:
            raise Exception("No such variable: %s" % key)


def read(filename, verbose=True):
    s = IDLSaveFile(filename)
    s.parse(verbose=verbose)
    return s
