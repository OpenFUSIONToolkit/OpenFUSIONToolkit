#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! xml editor for adding filament coil points
@authors Jamie Xia
@date Feb 2026
'''
import xml.etree.ElementTree as ET
from .._core import bool_to_string
# ===============================
# Icoil and Vcoil classes
# ===============================

class ThinCurr_coil_set:
    def __init__(self, name, sens_mask=False):
        """! Initialize general coil set
        @param name Name of the coil set
        @param sens_mask Masking flag for coil (if `True`, coil is masked from sensors)
        """
        self.name = name
        self.sens_mask = sens_mask
        self.subcoils = []  # List of subcoil dictionaries
    
    def add_subcoil(self, RZ=None, pts=None, scale=None, npoints=None, hdf5_path=None):
        """! Add a subcoil to this V-coil set
        @param RZ [R, Z] position (for circular coils; centered on Z-axis)
        @param pts Array of [x, y, z] positions (for general coils)
        @param scale Scaling factor for coil current (`1.0` if not specified)
        @param npoints Number of points (for general coils)
        @param hdf5_path Path to HDF5 dataset containing coil points (alternative to `pts` and `RZ`)
        """
        if pts is not None: # General 3D coil
            if (RZ is not None) or (hdf5_path is not None):
                raise ValueError("Only one of `RZ`, `pts`, or `hdf5_path` should be provided")
            if npoints is None:
                npoints = len(pts)
            self.subcoils.append({'pts': pts,'scale': scale,'npts': npoints})
        elif RZ is not None: # Circular coil
            if (pts is not None) or (hdf5_path is not None):
                raise ValueError("Only one of `RZ`, `pts`, or `hdf5_path` should be provided")
            self.subcoils.append({'RZ': RZ,'scale': scale})
        elif hdf5_path is not None: # Coil points from HDF5 dataset
            if (RZ is not None) or (pts is not None):
                raise ValueError("Only one of `RZ`, `pts`, or `hdf5_path` should be provided")
            self.subcoils.append({'hdf5_path': hdf5_path,'scale': scale})
        else:
            raise ValueError("Either `RZ`, `pts`, or `hdf5_path` must be provided")
    
    def build_coil_XML(self, coil_set):
        """! Build XML structure for this V-coil
        @param coil_set Coil set XML element to attach to
        """
        for subcoil in self.subcoils:
            if 'pts' in subcoil:
                # Custom 3D coil - detected by presence of 'pts' key
                coil_element = ET.SubElement(coil_set, "coil", attrib={"npts": str(subcoil['npts'])})
                # Format xyz points
                text_lines = []
                for xyz in subcoil['pts']:
                    text_lines.append('{0:E}, {1:E}, {2:E}'.format(xyz[0], xyz[1], xyz[2]))
                coil_element.text = "\n" + "\n".join(text_lines) + "\n"
            elif 'RZ' in subcoil:
                coil_element = ET.SubElement(coil_set, "coil")
                coil_element.text = '{0:E}, {1:E}'.format(subcoil['RZ'][0], subcoil['RZ'][1])
            elif 'hdf5_path' in subcoil:
                coil_element = ET.SubElement(coil_set, "coil", attrib={"path": subcoil['hdf5_path']})
                if subcoil['scale'] is not None:
                    coil_element.set("scale", str(subcoil['scale']))
            # Handling common attributes
            if subcoil['scale'] is not None:
                coil_element.set("scale", str(subcoil['scale']))

class ThinCurr_Icoil(ThinCurr_coil_set):
    """! I-coil class for defining current-driven coils in OpenFUSIONToolkit"""
    
    def __init__(self, name, sens_mask=False):
        """! Initialize I-coil
        @param name Name of the coil set
        @param sens_mask Masking flag for coil (if `True`, coil is masked from sensors)
        """
        super().__init__(name, sens_mask)
    
    def build_XML(self, parent_tag):
        """! Build XML structure for this I-coil
        @param parent_tag Parent XML element to attach to
        """
        # Build attributes for coil_set
        attrib = {"name": self.name}
        if self.sens_mask:
            attrib["sens_mask"] = bool_to_string(self.sens_mask)
        coil_set = ET.SubElement(parent_tag, "coil_set", attrib=attrib)
        # Build subcoil elements
        self.build_coil_XML(coil_set)


class ThinCurr_Vcoil(ThinCurr_coil_set):
    """! V-coil class for defining voltage-driven coils in OpenFUSIONToolkit"""
    
    def __init__(self, name, resistivity_per_length, radius, sens_mask=False):
        """! Initialize V-coil
        @param name Name of the coil set
        @param resistivity_per_length Resistivity per unit length
        @param radius Coil radius 
        @param sens_mask Masking flag for coil (if `True`, coil is masked from sensors)
        """
        super().__init__(name, sens_mask)
        # Vcoil-specific fields
        self.resistivity_per_length = resistivity_per_length
        self.radius = radius
    
    def build_XML(self, parent_tag):
        """! Build XML structure for this V-coil
        @param parent_tag Parent XML element to attach to
        """
        # Build attributes for coil_set
        attrib = {"name": self.name,"radius": str(self.radius), "res_per_len": str(self.resistivity_per_length)}
        if self.sens_mask:
            attrib["sens_mask"] = bool_to_string(self.sens_mask)
        coil_set = ET.SubElement(parent_tag, "coil_set", attrib=attrib)
        # Build subcoil elements
        self.build_coil_XML(coil_set)


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
        self.eta_vol = None
        self.thickness = None
    
    def add_Icoil(self, icoil):
        """! Add an I-coil to this ThinCurr block
        @param icoil ThinCurr_Icoil object
        """
        if not isinstance(icoil, ThinCurr_Icoil):
            raise TypeError("Icoil must be of type ThinCurr_Icoil")
        self.icoils.append(icoil)
    
    def add_Vcoil(self, vcoil):
        """! Add a V-coil to this ThinCurr block
        @param vcoil ThinCurr_Vcoil object
        """
        if not isinstance(vcoil, ThinCurr_Vcoil):
            raise TypeError("Vcoil must be of type ThinCurr_Vcoil")
        self.vcoils.append(vcoil)
    
    def set_eta(self, resistivities):
        """! Set resistivity values
        @param resistivities Resistivity values by section
        """
        self.eta = resistivities

    def set_eta_vol(self, resistivities):
        """! Set volumetric resistivity values
        @param resistivities array of float, volumetric resistivity values by section
        """
        self.eta_vol = resistivities

    def set_thickness(self, thicknesses):
        """! Set thickness values
        @param thicknesses array of float, thickness values by section
        """
        self.thickness = thicknesses
    
    def build_XML(self, parent_tag):
        """! Build XML structure for ThinCurr block
        @param parent_tag Parent XML element to attach to
        """
        thincurr_element = ET.SubElement(parent_tag, "thincurr")
        
        # Add eta (resistivities) if present
        if self.eta is not None:
            eta_element = ET.SubElement(thincurr_element, "eta")
            eta_element.text = ' '.join('{0:.E}'.format(r) for r in self.eta)

        # Add eta_vol (volumetric resistivities) if present
        if self.eta_vol is not None:
            eta_vol_element = ET.SubElement(thincurr_element, "eta_vol")
            eta_vol_element.text = ' '.join('{0:.E}'.format(r) for r in self.eta_vol)

        # Add thickness if present
        if self.thickness is not None:
            thickness_element = ET.SubElement(thincurr_element, "thickness")
            thickness_element.text = ' '.join('{0:.E}'.format(t) for t in self.thickness)
        
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