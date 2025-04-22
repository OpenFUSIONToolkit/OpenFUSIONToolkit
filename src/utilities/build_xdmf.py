#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
#
# Driver script for creating XDMF descriptor files for Open FUSION Toolkit output
#
#------------------------------------------------------------------------------
from __future__ import print_function
import argparse
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import h5py


class xdmf_mesh:
    def __init__(self,id,filename,np,nc,ne=0,prefix=None,point_list="R",cell_list="LC"):
        self.id=id
        self.filebase=filename
        self.np=np
        self.nc=nc
        self.ne=ne
        self.prefix=prefix
        self.point_list=point_list
        self.cell_list=cell_list

    def set_type(self,type):
        if type==31:
            self.type="Tetrahedron"
            self.ncn=4
            self.nnodes=self.np
        elif type==32:
            self.type="Tet_10"
            self.ncn=10
            self.nnodes=self.np+self.ne
        elif type==21:
            self.type="Triangle"
            self.ncn=3
            self.nnodes=self.np
        elif type==22:
            self.type="Tri_6"
            self.ncn=6
            self.nnodes=self.np+self.ne
        elif type==33:
            self.type="Hexahedron"
            self.ncn=8
            self.nnodes=self.np
        elif type==23:
            self.type="Quadrilateral"
            self.ncn=4
            self.nnodes=self.np+self.ne
        elif type==10:
            self.type="Polyline"
            self.ncn=2
            self.nnodes=self.np

    def get_type(self):
        return self.type

    def get_ncn(self):
        return self.ncn

    def get_point_list(self):
        if self.prefix is None:
            return self.point_list
        else:
            return self.prefix + "/" + self.point_list

    def get_cell_list(self):
        if self.prefix is None:
            return self.cell_list
        else:
            return self.prefix + "/" + self.cell_list

    def count(self,type,centering):
        if centering == 'Node':
            output="{0}".format(self.nnodes)
        elif centering == 'Cell':
            output="{0}".format(self.nc)
        else:
            raise ValueError("Unknown field centering")
        #
        if type == 'Scalar':
            output += " 1"
        elif type == 'Vector':
            output += " 3"
        else:
            raise ValueError("Unknown field type")
        return output

    def filename(self):
        output=self.filebase+"."+self.id+".h5"
        return output


class xdmf_fields:
    def __init__(self,id,attributes,ts=0):
        self.id=id
        self.type=attributes['TYPE'][0].decode()
        self.loc=attributes['CENTERING'][0].decode()
        self.ts=ts


class xdmf_doc:
    def __init__(self,filename,padSize=4,prettyPrint=False):
        self.filename=filename
        self.meshes=[]
        self.fields=[]
        self.time=0.
        self.padSize=padSize
        self.prettyPrint=prettyPrint

    def get_nfields(self):
        return len(self.fields)

    def get_nmeshes(self):
        return len(self.meshes)

    def add_mesh(self,mesh):
        self.meshes.append(mesh)

    def add_field(self,field):
        self.fields.append(field)

    def set_time(self,time):
        self.time=time

    def start_mesh(self):
        # Begin combined mesh
        self.Domain=ET.SubElement(self.doc,"Domain")
        self.GridOuter=ET.SubElement(self.Domain,"Grid")
        self.GridOuter.set("Name", "Mesh")
        self.GridOuter.set("GridType", "Collection")
        self.GridOuter.set("CollectionType", "Spatial")
        # Set time value
        self.Time=ET.SubElement(self.GridOuter,"Time")
        self.Time.set("Value", "{0:.8E}".format(self.time))

    def insert_block(self,mesh):
        # Start block
        Block=ET.SubElement(self.GridOuter,"Grid")
        Block.set("Name", "block_{0}".format(mesh.id))
        # Start Cell list
        Top=ET.SubElement(Block,"Topology")
        Top.set("Type", mesh.get_type() )
        Top.set("NumberOfElements","{0}".format(mesh.nc))
        if mesh.get_type() == "Polyline":
            Top.set("NodesPerElement"," 2 ")
        # Add cell list reference
        LC=ET.SubElement(Top,"DataItem")
        LC.set("Dimensions", "{0} {1}".format(mesh.nc,mesh.get_ncn()))
        LC.set("NumberType", "Int")
        LC.set("Format", "HDF")
        LC.text = mesh.filename()+":/"+mesh.get_cell_list()
        # Start Point list
        Geom=ET.SubElement(Block,"Geometry")
        Geom.set("Type", "XYZ")
        # Add point list reference
        R=ET.SubElement(Geom,"DataItem")
        R.set("Dimensions", "{0} 3 ".format(mesh.nnodes))
        R.set("NumberType", "Float")
        R.set("Format", "HDF")
        R.text = mesh.filename()+":/"+mesh.get_point_list()
        # Return constructed block
        return Block

    def insert_field(self,mesh,field,block):
        Attr=ET.SubElement(block,"Attribute")
        Attr.set("Name", field.id)
        Attr.set("Type", field.type)
        Attr.set("Center", field.loc)
        #
        F=ET.SubElement(Attr,"DataItem")
        F.set("Dimensions",mesh.count(field.type,field.loc))
        F.set("NumberType", "Float")
        F.set("Format", "HDF")
        F.text = mesh.filename()+":/{0}/{1}/{2}".format(mesh.prefix,str(field.ts).zfill(self.padSize),field.id)

    def write_file(self):
        # Begin Document
        self.doc=ET.Element("Xdmf")
        self.doc.set('Version',"2.0")
        # Start mesh group
        self.start_mesh()
        # Write individual meshes
        for mesh in self.meshes:
            # Start block
            block = self.insert_block(mesh)
            for field in self.fields:
                self.insert_field(mesh, field, block)
        # Write document to file
        tree = ET.ElementTree(self.doc)
        if self.prettyPrint:
            ET.indent(tree, space="  ", level=0)
        tree.write(self.filename)

# ===============================================================
# Input parameters
# ===============================================================
parser = argparse.ArgumentParser()
parser.description = "Create Xdmf file for VisIt from Open FUSION Toolkit output."
parser.add_argument('-i', '--inprefix', type=str, default='oft_xdmf', help='Input file name prefix, default="oft_xdmf"')
parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Display debug information during Xdmf write.')
parser.add_argument('-s', '--size', type=int, default=20, help='Size of XML file cache before flush.')
parser.add_argument('-p', '--pretty', action="store_true", default=False, help='Print nicely formatted XML documents.')
parser.add_argument('-k', '--keep', action="store_true", default=False, help='Keep existing output files.')
parser.add_argument('--block_padding', type=int, default=4, help='Size block index padding.')
parser.add_argument('--repeat_static', action="store_true", default=False, help='Insert static fields into all timesteps.')
args = parser.parse_args()

inprefix=args.inprefix
verbose=args.verbose
cacheSize=args.size
repeatStatic=args.repeat_static
padSize=args.block_padding

if not args.keep:
    import os
    import glob
    print('Removing old Xdmf files')
    files = glob.glob('*.xmf')
    for filename in files:
        os.remove(filename)
    print('  Removed {0} files'.format(len(files)))

print()
print('Creating output files: {0}.{1}.h5'.format(inprefix,'X'*padSize))

xml_docs = []
with h5py.File("{0}.{1}.h5".format(inprefix,str(1).zfill(padSize)),'r') as h5_file:
    for group_key, group_obj in h5_file.items():
        print('  Found Group: {0}'.format(group_key))
        for mesh_key, mesh_obj in group_obj.items():
            print('    Found Mesh: {0}'.format(mesh_key))
            mesh_type = mesh_obj['TYPE'][()]
            # Create mesh description
            np = mesh_obj['R'].shape[0]
            nc = mesh_obj['LC'].shape[0]
            mesh=xdmf_mesh('0001',inprefix,np,nc,int(0),'{0}/{1}'.format(group_key,mesh_key),'R','LC')
            mesh.set_type(mesh_type)
            meshes = [mesh,]
            # Get other blocks if present
            nblocks_ds = mesh_obj.get('NBLOCKS')
            if nblocks_ds is not None:
                nblocks = nblocks_ds[0]
                print('      # of blocks: {0}'.format(nblocks))
                for i in range(nblocks):
                    if i == 0:
                        continue
                    with h5py.File("{0}.{1}.h5".format(inprefix,str(i+1).zfill(padSize)),'r') as h5_file2:
                        block_group = h5_file2.get(group_key,None)
                        if block_group is None:
                            raise ValueError("Group not found in block {0}".format(i+1))
                        block_mesh = block_group.get(mesh_key,None)
                        if block_mesh is None:
                            raise ValueError("Mesh not found in block {0}".format(i+1))
                        np = block_mesh['R'].shape[0]
                        nc = block_mesh['LC'].shape[0]
                        mesh=xdmf_mesh('{0:04d}'.format(i+1),inprefix,np,nc,int(0),'{0}/{1}'.format(group_key,mesh_key),'R','LC')
                        mesh.set_type(mesh_type)
                        meshes.append(mesh)
            # Add static fields
            static_fields = []
            xml_doc = xdmf_doc('{0}_{1}-static.xmf'.format(group_key,mesh_key),padSize,args.pretty)
            for mesh in meshes:
                xml_doc.add_mesh(mesh)
            for field_name, field_obj in mesh_obj.get('0000',{}).items():
                field=xdmf_fields(field_name,field_obj.attrs)
                xml_doc.add_field(field)
                static_fields.append(field)
            xml_docs.append(xml_doc)
            # Add timesteps
            for i in range(9998):
                xml_doc = xdmf_doc('{0}_{1}-{2:04d}.xmf'.format(group_key,mesh_key,i+1),padSize,args.pretty)
                for mesh in meshes:
                    xml_doc.add_mesh(mesh)
                timestep = mesh_obj.get('{0:04d}'.format(i+1),None)
                if timestep is None:
                    break
                xml_doc.set_time(timestep['TIME'][0])
                for field_name, field_obj in timestep.items():
                    if field_name == 'TIME':
                        continue
                    field=xdmf_fields(field_name,field_obj.attrs,i+1)
                    xml_doc.add_field(field)
                if repeatStatic:
                    for field in static_fields:
                        xml_doc.add_field(field)
                if len(xml_doc.fields) > 0:
                    xml_docs.append(xml_doc)

for xml_doc in xml_docs:
    xml_doc.write_file()
