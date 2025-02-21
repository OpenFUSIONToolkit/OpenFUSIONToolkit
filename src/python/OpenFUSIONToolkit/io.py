'''! General I/O functionality for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date July 2024
@ingroup doxy_oft_python
'''

## @file io.py
#
# General I/O functionality for Open FUSION Toolkit (OFT) Python interfaces
from __future__ import print_function
import os
import sys
import struct
import re
import subprocess
import numpy
import h5py
eol_byte = '\n'.encode()


class histfile:
    r'''! Class for loading and analyzing a structured OFT binary file'''
    def __init__(self,filename):
        r'''! Load and parse binary file into Python representation'''
        def decode_list(list):
            '''! Decode list for Python 3 compatibility'''
            return [val.decode("utf-8") for val in list]
        
        def setup_sizes(data_reg):
            r'''! Reader header information and setup binary reads'''
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

        def setup_sizes_legacy():
            r'''! Reader header information and setup binary reads using legacy format'''
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
        #
        self._filename = filename
        self.header = []
        self.nlines = 0
        self._data = {}
        #
        pad_size = 4
        stream_data = True
        with open(self._filename, 'rb') as handle:
            self.content = handle.read()
        data_reg = self.content.find("--- BEGIN DATA ---".encode())
        if data_reg == -1:
            setup_sizes_legacy()
            stream_data = False
        else:
            data_reg = data_reg + 19
            setup_sizes(data_reg)
        #
        for field in self.field_tags:
            if field is not None:
                self._data[field] = []
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
                        self._data[field].append(tmp[k:k+nfields])
                    else:
                        self._data[field].append(tmp[k])
                    k += nfields
            else:
                self._data[self.field_tags[0]].append(zip(*[iter(tmp)]*self.dim[0]))
        # Convert all data to NumPy arrays
        for key, data in self._data.items():
            self._data[key] = numpy.asarray(data)

    def save_to_matlab(self,filename):
        r'''! Convert data to MATLAB format'''
        try:
            import scipy.io as sio
        except:
            print("Unable to load SCIPY")
            raise
        else:
            print("Converting file to MATLAB format")
            sio.savemat(filename, self._data, oned_as='row')

    def save_to_hdf5(self,filename):
        r'''! Convert data to HDF5 format'''
        try:
            import h5py
        except:
            print("Unable to load H5PY")
            raise
        else:
            print("Converting file to HDF5 format")
            with h5py.File(filename, 'w') as fid:
                for (ind, field) in enumerate(self.field_tags):
                    fid[field] = self._data[field]
                    fid[field].attrs['description'] = self.field_descriptions[field]

    def get(self, keyname, value=None):
        return self._data.get(keyname,value)

    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def __getitem__(self, key):
        return self._data[key]
    
    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        r'''! Print information about the file'''
        result = "\nOFT History file: {0}\n".format(self._filename)
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


class XDMF_plot_mesh:
    '''! Class containing data for a single mesh from OFT plot files'''
    def __init__(self,mesh_obj):
        '''! Read in data and setup object

        @param mesh_obj HDF5 group corresponding to mesh (file needs to be open)
        '''
        self.type = mesh_obj['TYPE'][()]
        self.r = numpy.asarray(mesh_obj['R'])
        self.lc = numpy.asarray(mesh_obj['LC'])
        self.np = self.r.shape[0]
        self.nc = self.lc.shape[0]
        #
        self.static_fields = {}
        for field_name, field_obj in mesh_obj.get('0000',{}).items():
            self.static_fields[field_name] = numpy.asarray(field_obj)
        #
        self.times = []
        self.time_field_names = []
        self.time_fields = []
        for i in range(1,9998):
            timestep = mesh_obj.get('{0:04d}'.format(i+1),None)
            if timestep is None:
                break
            self.times.append(timestep['TIME'][0])
            step_dict = {}
            for field_name, field_obj in timestep.items():
                if field_name == 'TIME':
                    continue
                if field_name not in self.time_field_names:
                    self.time_field_names.append(field_name)
                step_dict[field_name] = numpy.asarray(field_obj)
            self.time_fields.append(step_dict)
        self.times = numpy.array(self.times)
        
    def get_field(self,name,time=None,timestep=None):
        '''! Get raw data associated with field at a given time or timestep

        @param name Name of field
        @param time Time value to evaluate field (linear interpolation)
        @param timestep Timestep index to retrieve
        '''
        if timestep is not None:
            if name not in self.time_field_names:
                raise KeyError('"{0}" is not one of the timestep fields'.format(name))
            if time is not None:
                raise ValueError('"time" and "timestep" cannot be specified simultaneously')
            if (timestep < 0) and (timestep > self.times.shape[0]-1):
                raise ValueError("Requested timestep outside available range [0,{0}]".format(self.times.shape[0]-1))
            if (name not in self.time_fields[timestep]):
                raise KeyError('"{0}" is not available at requested timestep {1}'.format(name,timestep))
            return self.time_fields[timestep][name]
        elif time is not None:
            if name not in self.time_field_names:
                raise KeyError('"{0}" is not one of the timestep fields'.format(name))
            if timestep is not None:
                raise ValueError('"time" and "timestep" cannot be specified simultaneously')
            for i in range(self.times.shape[0]-1):
                if (self.times[i] <= time) and (self.times[i+1] >= time):
                    if (name not in self.time_fields[i]):
                        raise KeyError('"{0}" is not available at required timestep {1}'.format(name,i))
                    if (name not in self.time_fields[i+1]):
                        raise KeyError('"{0}" is not available at required timestep {1}'.format(name,i+1))
                    return (self.time_fields[i][name]-self.time_fields[i+1][name])*(time-self.times[i])/(self.times[i+1]-self.times[i]) + self.time_fields[i][name]
            raise ValueError("Requested time outside available range [{0:.6E},{1:.6E}]".format(self.times[0],self.times[-1]))
        else:
            if name not in self.static_fields:
                raise KeyError('"{0}" is not one of the static fields'.format(name))
            return self.static_fields[name]
    
    def get_pyvista_grid(self):
        '''! Get pyvista representation of grid

        @returns `pyvista.UnstructuredGrid` object for grid
        '''
        try:
            import pyvista
        except ImportError:
            print('Failed to load "pyvista" package')
            raise
        if self.type == 31:
            celltype = pyvista.CellType.TETRA
            ncv = 4
        elif self.type == 32:
            celltype = pyvista.CellType.QUADRATIC_TETRA
            ncv = 10
        elif self.type == 33:
            celltype = pyvista.CellType.HEXAHEDRON
            ncv = 8
        elif self.type == 21:
            celltype = pyvista.CellType.TRIANGLE
            ncv = 3
        elif self.type == 22:
            celltype = pyvista.CellType.QUADRATIC_TRIANGLE
            ncv = 6
        elif self.type == 23:
            celltype = pyvista.CellType.QUAD
            ncv = 4
        elif self.type == 10:
            celltype = pyvista.CellType.LINE
            ncv = 2
        celltypes = numpy.array([celltype for _ in range(self.lc.shape[0])], dtype=numpy.int8)
        cells = numpy.insert(self.lc, [0,], ncv, axis=1)
        return pyvista.UnstructuredGrid(cells, celltypes, self.r)


class XDMF_plot_file:
    '''! Helper class for interacting with OFT plotting output files'''
    def __init__(self,filepath):
        '''! Load data from file and setup object

        @param filename Path to HDF5 plotting file
        '''
        self.filepath = filepath
        self._groups = {}
        with h5py.File(self.filepath,'r') as h5_file:
            for group_key, group_obj in h5_file.items():
                self._groups[group_key.lower()] = {}
                for mesh_key, mesh_obj in group_obj.items():
                    self._groups[group_key.lower()][mesh_key.lower()] = XDMF_plot_mesh(mesh_obj)
    
    def get(self, keyname, value=None):
        '''! Get plotting group (list of @ref XDMF_plot_mesh "meshes")

        @param keyname Name of group
        @param value Return value if not found
        '''
        return self._groups.get(keyname.lower(),value)

    def keys(self):
        '''! Get plotting groups in file'''
        return self._groups.keys()

    def items(self):
        '''! Return iterator over plotting group (name, value) pairs'''
        return self._groups.items()

    def __getitem__(self, key):
        return self._groups[key.lower()]
    
    def __iter__(self):
        return iter(self._groups)


def build_XDMF(path='.',repeat_static=False,pretty=False,legacy=False):
    '''! Build XDMF plot metadata files 

    @param path Folder to build XDMF files in (must include `oft_xdmf.XXXX.h5` or `dump.dat` files)
    @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
    @param pretty Use pretty printing (indentation) in XDMF files?
    @param legacy Use legacy XDMF script for processing `dump.dat` files?
    '''
    if legacy:
        cmd = [
            "{0}".format(sys.executable),
            "{0}".format(os.path.join(os.path.dirname(__file__),'..','build_xdmf-legacy.py'))
        ]
    else:
        cmd = [
            "{0}".format(sys.executable),
            "{0}".format(os.path.join(os.path.dirname(__file__),'..','build_xdmf.py'))
        ]
    if repeat_static:
        cmd.append("--repeat_static")
    if pretty:
        cmd.append("--pretty")
    subprocess.run(cmd,cwd=path)
    if legacy:
        return None
    else:
        return XDMF_plot_file(os.path.join(path,'oft_xdmf.0001.h5'))