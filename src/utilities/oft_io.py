#---------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#---------------------------------------------------------------------------
#
# Class definition for manipulating Open FUSION Toolkit (OFT) structured binary files in Python
#
#---------------------------------------------------------------------------
from __future__ import print_function
import struct
import re
import argparse
eol_byte = '\n'.encode()


def decode_list(list):
    r'''Decode list for Python 3 compatibility'''
    return [val.decode("utf-8") for val in list]


class oft_histfile:
    r'''Class for loading and analyzing a structured OFT binary file'''
    def __init__(self,filename):
        r'''Load and parse binary file into Python representation'''
        self.filename = filename
        self.header = []
        self.nlines = 0
        self.data = {}
        #
        pad_size = 4
        stream_data = True
        with open(self.filename, 'rb') as handle:
            self.content = handle.read()
        data_reg = self.content.find("--- BEGIN DATA ---".encode())
        if data_reg == -1:
            self.setup_sizes_legacy()
            stream_data = False
        else:
            data_reg = data_reg + 19
            self.setup_sizes(data_reg)
        #
        for field in self.field_tags:
            if field is not None:
                self.data[field] = []
        #
        while (self.offset+1 < len(self.content)):
            # Look for end of line
            tmp1 = struct.unpack_from("=c", self.content, offset=self.offset-4)
            if tmp1[0] == eol_byte:
                self.offset += 1
            # Read line
            head_val = struct.unpack_from("i", self.content, offset=self.offset)
            self.offset += pad_size
            tmp = struct.unpack_from("=" + self.line_fmt, self.content, offset=self.offset)
            self.offset += self.line_length
            tail_val = struct.unpack_from("i", self.content, offset=self.offset)
            self.offset += pad_size
            self.nlines += 1
            if stream_data:
                if (head_val[0] != self.line_length) or (tail_val[0] != self.nfields):
                    print('Bad data line ({0})'.format(self.nlines))
                    raise IOError
            else:
                if (head_val[0] != self.line_length) or (tail_val[0] != self.line_length):
                    print('Bad data line ({0})'.format(self.nlines))
                    raise IOError
            if self.nfields > 1:
                k = 0
                for (j, field) in enumerate(self.field_tags):
                    if field is None:
                        k += 1
                        continue
                    nfields = self.field_repeats[j]
                    if nfields > 1:
                        self.data[field].append(tmp[k:k+nfields])
                    else:
                        self.data[field].append(tmp[k])
                    k += nfields
            else:
                self.data[self.field_tags[0]].append(zip(*[iter(tmp)]*self.dim[0]))

    def setup_sizes(self, data_reg):
        r'''Reader header information and setup binary reads'''
        base_dict = {}
        self.dim = None
        for line in self.content[:data_reg].split(eol_byte):
            line_tmp = line.decode()
            if (len(line_tmp) > 0) and (line_tmp[0] == '#'):
                self.header.append(line_tmp)
                continue
            fields = line_tmp.split(':')
            nfields = len(fields)
            if nfields == 1:
                continue
            elif nfields == 2:
                if len(fields[1]) == 0:
                    inner_dict = {}
                    base_dict[fields[0]] = inner_dict
                else:
                    if fields[0][2] == '-':
                        inner_dict[fields[0][4:]] = fields[1]
                    else:
                        base_dict[fields[0]] = fields[1]
        #
        self.nfields = int(base_dict['nfields'])
        self.offset = data_reg
        self.field_tags = []
        for field in base_dict['fields'].split():
            self.field_tags.append(field)
        self.field_descriptions = {}
        for field, description in base_dict['descriptions'].items():
            self.field_descriptions[field] = description.strip()
        #
        self.field_types = []
        self.field_sizes = []
        for field in base_dict['field_types'].split():
            if field == 'i4':
                self.field_types.append('i')
                self.field_sizes.append(4)
            elif field == 'i8':
                self.field_types.append('l')
                self.field_sizes.append(8)
            elif field == 'r4':
                self.field_types.append('f')
                self.field_sizes.append(4)
            elif field == 'r8':
                self.field_types.append('d')
                self.field_sizes.append(8)
            else:
                print('unknown type')
                raise IOError
        #
        self.field_repeats = []
        for field in base_dict['field_sizes'].split():
            self.field_repeats.append(int(field))
        #
        if self.nfields == 1:
            self.line_length = self.field_sizes[0]*self.dim[0]*self.dim[1]
            self.line_fmt = '{0}'.format(self.dim[0]*self.dim[1]) + self.field_types[0]
        else:
            self.line_length = 0
            self.line_fmt = ''
            for (ind, fsize) in enumerate(self.field_sizes):
                self.line_length += fsize*self.field_repeats[ind]
                self.line_fmt += self.field_types[ind]*self.field_repeats[ind]

    def setup_sizes_legacy(self):
        r'''Reader header information and setup binary reads using legacy format'''
        self.offset = 4
        tmp = struct.unpack_from("i",self.content,offset=self.offset)
        self.offset += 12
        self.nfields = int(tmp[0])
        if self.nfields == 1:
            tmp = struct.unpack_from("2i",self.content,offset=self.offset)
            self.offset += 16
            self.dim = tmp
        else:
            self.dim = None
        #
        if (self.dim != None) and (self.nfields > 1):
            print("Arrays not supported with multiple fields")
            raise ValueError
        #
        self.field_tags = []
        for _ in range(self.nfields):
            tmp = decode_list(struct.unpack_from("20s", self.content, offset=self.offset))
            self.offset += 20
            tag_name = tmp[0].strip()
            pattern = re.compile("[0-9]")
            if pattern.match(tag_name):
                tag_name = "d_" + tag_name
            self.field_tags.append(tag_name)
        self.offset += 8
        #
        self.field_types = []
        self.field_sizes = []
        self.field_repeats = []
        for _ in range(self.nfields):
            tmp_type = decode_list(struct.unpack_from("s",self.content,offset=self.offset))
            self.offset += 1
            #
            tmp = decode_list(struct.unpack_from("s",self.content,offset=self.offset))
            self.offset += 1
            self.field_sizes.append(int(tmp[0]))
            self.field_repeats.append(1)
            #
            if (tmp_type[0] + tmp[0]) == 'i4':
                self.field_types.append('i')
            elif (tmp_type[0] + tmp[0]) == 'i8':
                self.field_types.append('l')
            elif (tmp_type[0] + tmp[0]) == 'r4':
                self.field_types.append('f')
            elif (tmp_type[0] + tmp[0]) == 'r8':
                self.field_types.append('d')
            else:
                print('unknown type')
                raise IOError
        self.offset += 4
        #
        if self.nfields == 1:
            self.line_length = self.field_sizes[0]*self.dim[0]*self.dim[1]
            self.line_fmt = '{0}'.format(self.dim[0]*self.dim[1]) + self.field_types[0]
        else:
            self.line_length = 0
            self.line_fmt = ''
            for (ind, fsize) in enumerate(self.field_sizes):
                self.line_length += fsize
                self.line_fmt += self.field_types[ind]
        return

    def save_to_matlab(self,filename):
        r'''Convert data to MATLAB format'''
        try:
            import scipy.io as sio
        except:
            print("Unable to load SCIPY")
            raise
        else:
            print("Converting file to MATLAB format")
            sio.savemat(filename, self.data, oned_as='row')

    def save_to_hdf5(self,filename):
        r'''Convert data to HDF5 format'''
        try:
            import h5py
        except:
            print("Unable to load H5PY")
            raise
        else:
            print("Converting file to HDF5 format")
            with h5py.File(filename, 'w') as fid:
                for (ind, field) in enumerate(self.field_tags):
                    fid[field] = self.data[field]
                    fid[field].attrs['description'] = self.field_descriptions[field]

    def __repr__(self):
        r'''Print information about the file'''
        result = "\nOFT History file: {0}\n".format(self.filename)
        result += "  Number of fields = {0}\n".format(self.nfields)
        result += "  Number of entries = {0}\n".format(self.nlines)
        if self.dim != None:
            result += "  Field dimension = {0}\n".format(self.dim)
        result += "\n  Fields:\n"
        #
        for (ind, field) in enumerate(self.field_tags):
            if field is not None:
                result += '    {0}: {3} ({1}{2})\n'.format(field, self.field_types[ind], self.field_repeats[ind], self.field_descriptions[field])
        return result
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Pre-processing script for mesh files"
    parser.add_argument("--files", type=str, default=None, nargs='+', required=True, help="Files to view or convert")
    parser.add_argument("--convert_hdf5", action="store_true", default=False, help="Convert files to HDF5? (default: False)")
    parser.add_argument("--convert_matlab", action="store_true", default=False, help="Convert files to MATLAB? (default: False)")
    options = parser.parse_args()

    for file in options.files:
        hist_file = oft_histfile(file)
        file_prefix = file.split('.')[0]
        if options.convert_hdf5:
            hist_file.save_to_hdf5(file_prefix + ".h5")
        if options.convert_matlab:
            hist_file.save_to_matlab(file_prefix + ".mat")
        if not (options.convert_hdf5 or options.convert_matlab):
            print(hist_file)