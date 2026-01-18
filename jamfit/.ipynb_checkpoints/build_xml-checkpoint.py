import xml.etree.cElementTree as ET
import sys
import inspect

class Vcoil:
    def __init__(self, name, radius=None, Z = None, xyz = None, resistivity_per_length = None, sens_mask = None, npoints = None, iscustom = False):
        self.name = name
        self.resistivity_per_length = resistivity_per_length
        self.sens_mask = sens_mask
        self.npoints = npoints
        if iscustom: 
            if xyz is not None:
                self.xyz = xyz
                self.npoints = npoints
            else:
                raise ValueError("Custom coil must have a valid xyz = [x, y, z]")
            self.radius = None
        else: 
            self.xyz = None 
            self.radius = radius
            self.Z = Z

        
class Icoil: 
    def __init__(self, name, radius = None, Z=None, xyz = None, scale = None, sens_mask = None, npoints = None, iscustom = False):
        self.name = name
        self.sens_mask = sens_mask
        self.scale = scale
        self.custom = iscustom
        if iscustom: 
            if xyz is not None :
                self.xyz = xyz
                self.npoints = npoints 
            else:
                raise ValueError("Custom coil must have a valid xyz = [x, y, z]")
            self.radius = None
        else: 
            self.xyz = None 
            self.radius = radius
            self.Z = Z 



def create_oft(path, icoils, vcoils, resistivities): 
    '''
    Function solves filament and wall currents over time 
    
    Parameters:
    ----------
    path : string
        name and location of where file will be store i.e '/home/user/folder/oft.xml' 
    
    icoils : icoils
        should be a list of icoil objects [icoil1, icoil2, icoil3, ...]
        
    vcoils : vcoils
    should be a list of vcoil objects [vcoil1, vcoil2, vcoil3 , .... ]

    resistivites : float array
    should be list of resistivites by section [r1, r2, r3...] 
    
        
    Returns:
    ----------
    oft.xml file is created in specified directory. 
    '''
    oft_element = ET.Element("oft")
    xml_doc=ET.ElementTree(oft_element)
    thincurr_element = ET.SubElement(oft_element,"thincurr")
    eta_element = ET.SubElement(thincurr_element, "eta") 
    eta_element.text = ' '.join('{0:.4E}'.format(r) for r in resistivities) 
    if icoils: 
        icoil_element = ET.SubElement(thincurr_element, "icoils")
        for i in icoils: 
            if isinstance(i, Icoil):
                coil_set = ET.SubElement(icoil_element, "coil_set", attrib = {"name" : i.name})
                for xyz in i.xyz: 
                    coil_element = ET.SubElement(coil_set, "coil", attrib={"npoints": i.npoints})
                    if i.scale is None:                 
                        coil_element.set("scale","1.0")
                    else:     
                        coil_element.set("scale", i.scale)
                    coil_element.text = '{0:.6f}, {1:.6f}'.format(xyz[0], xyz[1], xyz[2])
            else: 
                raise ValueError("Non I-coil class in inputted I-coil list")
    if vcoils: 
        vcoil_element = ET.SubElement(thincurr_element, "vcoils")
        for v in vcoils: 
            if isinstance(v, Vcoil): 
                if v.custom: 
                    coil_set = ET.SubElement(vcoil_element, "coil_set", attrib = {"name" : v.name, "radius" : v.radius, "res_per_len" : v.res_per_len, "sens_mask" : v.sens_mask})
                    for xyz in v.xyz: 
                        coil_element = ET.SubElement(coil_set, "coil", attrib={"npoints": v.npoints})
                        if v.scale is None:                 
                            coil_element.set("scale","1.0")
                        else:     
                            coil_element.set("scale", i.scale)
                        coil_element.text = '{0:.6f}, {1:.6f}'.format(xyz[0], xyz[1], xyz[2])
                else:
                    coil_set = ET.SubElement(vcoil_element, "coil_set", attrib = {"name" : v.name})
                    for r_val, z_val in zip(v.radius, v.Z): 
                        coil_element = ET.SubElement(coil_set, "coil")
                        if v.scale is None:                 
                            coil_element.set("scale","1.0")
                        else:     
                            coil_element.set("scale", v.scale)
                        coil_element.text = '{0:.6f}, {1:.6f}'.format(r_val, z_val)
            else: 
                raise ValueError("Non V-coil class in inputted V-coil list")
    ET.indent(xml_doc, space="  ", level=0)
    xml_doc.write(path)
    print("oft file created") 
            

        
        