#!/usr/bin/env python
#---------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#---------------------------------------------------------------------------
#
# Driver script for creating XDMF descriptor files for Open FUSION Toolkit output
#
#---------------------------------------------------------------------------
from __future__ import print_function
import argparse
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import h5py

# Indent ETREE document
def indent(elem, level=0):
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i
#
class xdmf_mesh:
    def __init__(self,id,filename,np,nc,ne=0,point_list="R",cell_list="LC"):
        self.id=id
        self.filebase=filename
        self.np=np
        self.nc=nc
        self.ne=ne
        self.point_list=point_list
        self.cell_list=cell_list
    #
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

    #
    def get_type(self):
        return self.type
    #
    def get_ncn(self):
        return self.ncn
    #
    def get_point_list(self):
        return self.point_list
    #
    def get_cell_list(self):
        return self.cell_list
    #
    def count(self,type):
        if type==11:
            output=" "+repr(self.nnodes)+" 1 "
            return output
        elif type==12:
            output=" "+repr(self.nc)+" 1 "
            return output
        elif type==31:
            output=" "+repr(self.nnodes)+" 1 "
            return output
        elif type==32:
            output=" "+repr(self.nc)+" 1 "
            return output
        elif type==21:
            output=" "+repr(self.nnodes)+" 3 "
            return output
        elif type==22:
            output=" "+repr(self.nc)+" 3 "
            return output
        elif type==41:
            output=" "+repr(self.nnodes)+" 3 "
            return output
        elif type==42:
            output=" "+repr(self.nc)+" 3 "
            return output
    #
    def filename(self):
        output=self.filebase+"."+self.id+".h5"
        return output
#
class xdmf_fields:
    def __init__(self,id,ts=0):
        self.id=id
        self.ts=ts
    #
    def set_type(self,type):
        self.tind=type
        if type==11:
            self.type="Scalar"
            self.loc="Node"
            self.fbase="scalar_dump"
        elif type==12:
            self.type="Scalar"
            self.loc="Cell"
            self.fbase="scalar_dump"
        elif type==31:
            self.type="Scalar"
            self.loc="Node"
            self.fbase="scalar_dump"
        elif type==32:
            self.type="Scalar"
            self.loc="Cell"
            self.fbase="scalar_dump"
        elif type==21:
            self.type="Vector"
            self.loc="Node"
            self.fbase="vector_dump"
        elif type==22:
            self.type="Vector"
            self.loc="Cell"
            self.fbase="vector_dump"
        elif type==41:
            self.type="Vector"
            self.loc="Node"
            self.fbase="vector_dump"
        elif type==42:
            self.type="Vector"
            self.loc="Cell"
            self.fbase="vector_dump"
#
class xdmf_doc:
    def __init__(self,filename,padSize=3,prettyPrint=False):
        self.filename=filename
        self.meshes=[]
        self.fields=[]
        self.time=0.
        self.padSize=padSize
        self.prettyPrint=prettyPrint
    #
    def get_nfields(self):
        return len(self.fields)
    #
    def get_nmeshes(self):
        return len(self.meshes)
    #
    def add_mesh(self,mesh):
        self.meshes.append(mesh)
    #
    def add_field(self,field):
        self.fields.append(field)
    #
    def set_time(self,time):
        self.time=time
    #
    def start_mesh(self):
        # Begin combined mesh
        self.Domain=ET.SubElement(self.doc,"Domain")
        self.GridOuter=ET.SubElement(self.Domain,"Grid")
        self.GridOuter.set("Name", "Mesh")
        self.GridOuter.set("GridType", "Collection")
        self.GridOuter.set("CollectionType", "Spatial")
        # Set time value
        self.Time=ET.SubElement(self.GridOuter,"Time")
        self.Time.set("Value", str(self.time))
    #
    def insert_block(self,mesh):
        # Start block
        Block=ET.SubElement(self.GridOuter,"Grid")
        Block.set("Name", "block_"+mesh.id)
        # Start Cell list
        Top=ET.SubElement(Block,"Topology")
        Top.set("Type", mesh.get_type() )
        Top.set("NumberOfElements"," "+str(mesh.nc)+" ")
        if mesh.get_type() == "Polyline":
            Top.set("NodesPerElement"," 2 ")
        # Add cell list reference
        LC=ET.SubElement(Top,"DataItem")
        LC.set("Dimensions", " "+str(mesh.nc)+" "+str(mesh.get_ncn())+" " )
        LC.set("NumberType", "Int")
        LC.set("Format", "HDF")
        LC.text = mesh.filename()+":/"+mesh.get_cell_list()
        # Start Point list
        Geom=ET.SubElement(Block,"Geometry")
        Geom.set("Type", "XYZ")
        # Add point list reference
        R=ET.SubElement(Geom,"DataItem")
        R.set("Dimensions", " "+str(mesh.nnodes)+" 3 ")
        R.set("NumberType", "Float")
        R.set("Format", "HDF")
        R.text = mesh.filename()+":/"+mesh.get_point_list()
        # Return constructed block
        return Block
    #
    def insert_field(self,mesh,field,block):
        Attr=ET.SubElement(block,"Attribute")
        Attr.set("Name", field.id)
        Attr.set("Type", field.type)
        Attr.set("Center", field.loc)
        #
        F=ET.SubElement(Attr,"DataItem")
        F.set("Dimensions",mesh.count(field.tind))
        F.set("NumberType", "Float")
        F.set("Format", "HDF")
        F.text = field.fbase+"."+mesh.id+".h5"+":/"+field.id+str(field.ts).zfill(self.padSize)
    # Write assembled XML document to file
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
        if self.prettyPrint:
            indent(self.doc)
        output=ET.tostring(self.doc).decode()
        with open(self.filename,'w') as fid:
            fid.write(output)

# ===============================================================
# Input parameters
# ===============================================================
parser = argparse.ArgumentParser()
parser.description = "Create Xdmf file for VisIt from Open FUSION Toolkit output."
parser.add_argument('-i', '--infile', type=str, default='oft_plot.0001.h5', help='Input file name, default="dump.dat"')
parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Display debug information during Xdmf write.')
parser.add_argument('-s', '--size', type=int, default=20, help='Size of XML file cache before flush.')
parser.add_argument('-p', '--pretty', action="store_true", default=False, help='Print nicely formatted XML documents.')
parser.add_argument('-k', '--keep', action="store_true", default=False, help='Keep existing output files.')
parser.add_argument('--repeat_static', action="store_true", default=False, help='Insert static fields into all timesteps.')
args = parser.parse_args()

infile=args.infile
verbose=args.verbose
cacheSize=args.size
repeatStatic=args.repeat_static
volFirst=True
surfFirst=True
padSize=4

if not args.keep:
    import os
    import glob
    print('Removing old Xdmf files')
    files = glob.glob('*.xmf')
    for filename in files:
        os.remove(filename)

print('Creating output files')

if verbose:
    print('Reading in output data')

xml_docs = []
with h5py.File(infile,'r') as h5_file:
    for key in h5_file:
        print(key)
        print(h5_file[key]['R'])
        print(h5_file[key]['LC'])
        mesh_type = h5_file[key]['TYPE'][()]
        # Create mesh description
        np = h5_file[key]['R'].shape[0]
        nc = h5_file[key]['LC'].shape[0]
        mesh=xdmf_mesh('0001','oft_plot',np,nc,int(0),'{0}/R'.format(key),'{0}/LC'.format(key))
        mesh.set_type(mesh_type)
        # Add static fields
        xml_doc = xdmf_doc('{0}_static.xmf'.format(key),padSize,args.pretty)
        xml_doc.add_mesh(mesh)
        for static_field in h5_file[key].get('static',{}):
            field=xdmf_fields(static_field)
            field.set_type(h5_file[key]['static'][static_field]['TYPE'][()])
            xml_doc.add_field(field)
        xml_docs.append(xml_doc)
        # Add timesteps
        for i in range(9999):
            xml_doc = xdmf_doc('{0}_{1:04d}.xmf'.format(key,i+1),padSize,args.pretty)
            xml_doc.add_mesh(mesh)
            timestep = h5_file[key].get('timestep{0:04d}'.format(i),None)
            if timestep is None:
                break
            xml_doc.set_time(timestep['TIME'][()])
            for field_name in timestep:
                if field_name == 'TIME':
                    continue
                field=xdmf_fields(field_name)
                field.set_type(timestep[field_name]['TYPE'][()])
                xml_doc.add_field(field)
            xml_docs.append(xml_doc)

for xml_doc in xml_docs:
    xml_doc.write_file()


