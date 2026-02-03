import xml.etree.ElementTree as ET


'''! xml editor for adding filament coil points
@authors Jamie Xia
@date Feb 2026
'''

# ===============================
# Top-level XML writer (base OFT_env class method)
# ===============================

def write_oft_xml(xml_blocks, path):
    """! Write OFT XML file from a list of XML block objects.
    @param xml_blocks: list of objects with build_XML methods
    @param path: string, output path for XML file
    """
    oft_element = ET.Element("oft")
    xml_doc = ET.ElementTree(oft_element)
    for xml_block in xml_blocks:
        xml_block.build_XML(oft_element)
    ET.indent(xml_doc, space="  ", level=0)
    xml_doc.write(path, encoding="utf-8", xml_declaration=True)
    print(f"XML file created at {path}")


# ===============================
# Icoil and Vcoil classes
# ===============================

class ThinCurr_Icoil:
    """! I-coil class for defining current-driven coils in OpenFUSIONToolkit"""
    
    def __init__(self, name, sens_mask=False):
        """! Initialize I-coil
        @param name: string, name of the coil set
        @param sens_mask: bool or array, sensitivity mask for coil points
        """
        self.name = name
    
        self.sens_mask = sens_mask
        self.subcoils = []  # List of (data, is_custom) tuples
    
    def add_subcoil(self, RZ=None, pts=None, scale=None, npoints=None):
        """! Add a subcoil to this I-coil set
        @param RZ: array of [R, Z] positions (for non-custom coils)
        @param pts: array of [x, y, z] positions (for custom coils)
        @param scale: float, scaling factor for coil current
        @param npoints: int, number of points (for custom coils)
        """
        if pts is not None:
            # Custom coil with 3D points
            self.subcoils.append({
                'is_custom': True,
                'pts': pts,
                'scale': scale if scale is not None else 1.0,
                'npoints': npoints
            })
        elif RZ is not None:
            # Standard RZ coil
            self.subcoils.append({
                'is_custom': False,
                'RZ': RZ,
                'scale': scale if scale is not None else 1.0
            })
        else:
            raise ValueError("Either RZ or pts must be provided")
    
    def build_XML(self, parent_tag):
        """! Build XML structure for this I-coil
        @param parent_tag: ET.Element, parent XML element to attach to
        """
        coil_set = ET.SubElement(parent_tag, "coil_set", attrib={"name": self.name})
        
        for subcoil in self.subcoils:
            if subcoil['is_custom']:
                # Custom 3D coil
                coil_element = ET.SubElement(
                    coil_set, 
                    "coil", 
                    attrib={"npoints": str(subcoil['npoints'])}
                )
                coil_element.set("scale", str(subcoil['scale']))
                # Format xyz points
                text_lines = []
                for xyz in subcoil['pts']:
                    text_lines.append('{0:.6f}, {1:.6f}, {2:.6f}'.format(xyz[0], xyz[1], xyz[2]))
                coil_element.text = "\n" + "\n".join(text_lines) + "\n"
            else:
                # Standard RZ coil
                for rz in subcoil['RZ']:
                    coil_element = ET.SubElement(coil_set, "coil")
                    coil_element.set("scale", str(subcoil['scale']))
                    coil_element.text = '{0:.6f}, {1:.6f}'.format(rz[0], rz[1])


class ThinCurr_Vcoil:
    """! V-coil class for defining voltage-driven coils in OpenFUSIONToolkit"""
    
    def __init__(self, name, resistivity_per_length=None, radius=None, sens_mask=False):
        """! Initialize V-coil
        @param name: string, name of the coil set
        @param resistivity_per_length: float, resistivity per unit length
        @param radius: float, coil radius
        @param sens_mask: bool or array, sensitivity mask for coil points
        """
        self.name = name
        self.resistivity_per_length = resistivity_per_length
        self.radius = radius
        self.sens_mask = sens_mask
        self.subcoils = []  # List of subcoil dictionaries
    
    def add_subcoil(self, RZ=None, pts=None, scale=None, npoints=None):
        """! Add a subcoil to this V-coil set
        @param RZ: array of [R, Z] positions (for non-custom coils)
        @param pts: array of [x, y, z] positions (for custom coils)
        @param scale: float, scaling factor for coil current
        @param npoints: int, number of points (for custom coils)
        """
        if pts is not None:
            # Custom coil with 3D points
            self.subcoils.append({
                'is_custom': True,
                'pts': pts,
                'scale': scale if scale is not None else 1.0,
                'npoints': npoints
            })
        elif RZ is not None:
            # Standard RZ coil
            self.subcoils.append({
                'is_custom': False,
                'RZ': RZ,
                'scale': scale if scale is not None else 1.0
            })
        else:
            raise ValueError("Either RZ or pts must be provided")
    
    def build_XML(self, parent_tag):
        """! Build XML structure for this V-coil
        @param parent_tag: ET.Element, parent XML element to attach to
        """
        # Build attributes for coil_set
        attrib = {"name": self.name}
        if self.radius is not None:
            attrib["radius"] = str(self.radius)
        if self.resistivity_per_length is not None:
            attrib["res_per_len"] = str(self.resistivity_per_length)
        if self. :
            attrib["sens_mask"] = str(self.sens_mask)
        
        coil_set = ET.SubElement(parent_tag, "coil_set", attrib=attrib)
        
        for subcoil in self.subcoils:
            if subcoil['is_custom']:
                # Custom 3D coil
                coil_element = ET.SubElement(
                    coil_set, 
                    "coil", 
                    attrib={"npoints": str(subcoil['npoints'])}
                )
                # Format xyz points
                text_lines = []
                for xyz in subcoil['pts']:
                    text_lines.append('{0:.6f}, {1:.6f}, {2:.6f}'.format(xyz[0], xyz[1], xyz[2]))
                coil_element.text = "\n" + "\n".join(text_lines) + "\n"
            else:
                # Standard RZ coil
                for rz in subcoil['RZ']:
                    coil_element = ET.SubElement(coil_set, "coil")
                    coil_element.set("scale", str(subcoil['scale']))
                    coil_element.text = '{0:.6f}, {1:.6f}'.format(rz[0], rz[1])


# ===============================
# Base ThinCurr XML class
# ===============================

class ThinCurr_XML:
    """! Container class for ThinCurr XML block"""
    
    name = "thincurr"
    
    def __init__(self):
        """! Initialize ThinCurr XML block"""
        self.icoils = []
        self.vcoils = []
        self.eta = None
    
    def add_Icoil(self, icoil):
        """! Add an I-coil to this ThinCurr block
        @param icoil: ThinCurr_Icoil object
        """
        if not isinstance(icoil, ThinCurr_Icoil):
            raise TypeError("Icoil must be of type ThinCurr_Icoil")
        self.icoils.append(icoil)
    
    def add_Vcoil(self, vcoil):
        """! Add a V-coil to this ThinCurr block
        @param vcoil: ThinCurr_Vcoil object
        """
        if not isinstance(vcoil, ThinCurr_Vcoil):
            raise TypeError("Vcoil must be of type ThinCurr_Vcoil")
        self.vcoils.append(vcoil)
    
    def set_eta(self, resistivities):
        """! Set resistivity values
        @param resistivities: array of float, resistivity values by section
        """
        self.eta = resistivities
    
    def build_XML(self, parent_tag):
        """! Build XML structure for ThinCurr block
        @param parent_tag: ET.Element, parent XML element to attach to
        """
        thincurr_element = ET.SubElement(parent_tag, "thincurr")
        
        # Add eta (resistivities) if present
        if self.eta is not None:
            eta_element = ET.SubElement(thincurr_element, "eta")
            eta_element.text = ' '.join('{0:.4E}'.format(r) for r in self.eta)
        
        # Add I-coils if present
        if self.icoils:
            icoil_element = ET.SubElement(thincurr_element, "icoils")
            for icoil in self.icoils:
                icoil.build_XML(icoil_element)
        
        # Add V-coils if present
        if self.vcoils:
            vcoil_element = ET.SubElement(thincurr_element, "vcoils")
            for vcoil in self.vcoils:
                vcoil.build_XML(vcoil_element)

