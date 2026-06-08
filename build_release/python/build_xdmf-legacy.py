#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
#
# Legacy driver script for creating XDMF descriptor files for Open FUSION Toolkit output
#
#------------------------------------------------------------------------------
from __future__ import print_function
import argparse
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

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
        if type==1:
            self.type="Tetrahedron"
            self.ncn=4
            self.nnodes=self.np
        elif type==2:
            self.type="Tet_10"
            self.ncn=10
            self.nnodes=self.np+self.ne
        elif type==11:
            self.type="Triangle"
            self.ncn=3
            self.nnodes=self.np
        elif type==12:
            self.type="Tri_6"
            self.ncn=6
            self.nnodes=self.np+self.ne
        elif type==3:
            self.type="Hexahedron"
            self.ncn=8
            self.nnodes=self.np
        elif type==13:
            self.type="Quadrilateral"
            self.ncn=4
            self.nnodes=self.np+self.ne
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
parser.add_argument('-i', '--infile', type=str, default='dump.dat', help='Input file name, default="dump.dat"')
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

def flush_io(vol_docs,surf_docs):
    global volFirst, surfFirst
    nVolFiles = len(vol_docs)
    # Check for non-empty volume file
    if volFirst:
        for doc in vol_docs:
            if doc.get_nfields() > 0:
                volFirst = False
                break
    # Check for non-empty boundary file
    if surfFirst:
        for doc in surf_docs:
            if doc.get_nfields() > 0:
                surfFirst = False
                break
    # Output volume files
    for doc in vol_docs:
        if doc.get_nfields() > 0 or volFirst:
            doc.write_file()
        volFirst = False
    # Output boundary files
    for doc in surf_docs:
        if doc.get_nfields() > 0 or surfFirst:
            doc.write_file()
        surfFirst = False

if not args.keep:
    import os
    import glob
    print('Removing old Xdmf files')
    files = glob.glob('*.xmf')
    for filename in files:
        os.remove(filename)

print('Creating output files')

f=open(infile,'r')
if verbose:
    print('Reading in output data')
line = f.readline()
line = f.readline()
vals=line.split()
order = int(vals[0])

if order==1:
    mesh_type=1
    bmesh_type=11
elif order==2:
    mesh_type=2
    bmesh_type=12
elif order==3:
    mesh_type=3
    bmesh_type=13

meshes=[]
bmeshes=[]
for line in f:
    vals=line.split()
    if vals==[]:
        if verbose:
            print('Mesh Data Complete')
        break
    mesh=xdmf_mesh(vals[0],"mesh",int(vals[1]),int(vals[2]),int(0),'R_vol','LC_vol')
    mesh.set_type(mesh_type)
    meshes.append(mesh)
    if int(vals[4])>0:
        bmesh=xdmf_mesh(vals[0],"mesh",int(vals[3]),int(vals[4]),int(0),'R_surf','LC_surf')
        bmesh.set_type(bmesh_type)
        bmeshes.append(bmesh)

vol_docs=[]
static_fields=[]
surf_docs=[]
surf_static_fields=[]
ind=-1
for line in f:
    vals=line.split()
    if vals==[]:
        break
    if vals[0]=='Time':
        ind=ind+1
        if ind==0:
            vol_doc=xdmf_doc('static.xmf',padSize,args.pretty)
            surf_doc=xdmf_doc('surf_static.xmf',padSize,args.pretty)
        else:
            vol_doc=xdmf_doc('out_'+str(ind).zfill(padSize)+'.xmf',padSize,args.pretty)
            surf_doc=xdmf_doc('surf_out_'+str(ind).zfill(padSize)+'.xmf',padSize,args.pretty)
        for mesh in meshes:
            vol_doc.add_mesh(mesh)
        for bmesh in bmeshes:
            surf_doc.add_mesh(bmesh)
        vol_doc.set_time(float(vals[2]))
        surf_doc.set_time(float(vals[2]))
        if verbose:
            print('Adding Time Step: {0}'.format(ind))
        for field in static_fields:
            vol_doc.add_field(field)
        for field in surf_static_fields:
            surf_doc.add_field(field)
        for line in f:
            vals=line.split()
            if vals==[]:
                break
            if vals[1]=='Data':
                continue
            if verbose:
                print(vals)
            field=xdmf_fields(vals[0],ind)
            field.set_type(int(vals[1]))
            if int(vals[1])<30:
                vol_doc.add_field(field)
                if repeatStatic:
                    if ind==0:
                        static_fields.append(field)
            else:
                surf_doc.add_field(field)
                if repeatStatic:
                    if ind==0:
                        surf_static_fields.append(field)
        vol_docs.append(vol_doc)
        surf_docs.append(surf_doc)
    #
    if len(vol_docs) == cacheSize:
        if verbose:
            print('Flushing')
        flush_io(vol_docs,surf_docs)
        vol_docs=[]
        surf_docs=[]

if vol_docs != []:
    if verbose:
        print('Flushing')
    flush_io(vol_docs,surf_docs)

if verbose:
    print('Field Data Complete')
f.close()
