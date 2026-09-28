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
import numpy
import h5py
from .._interface import bool_to_string, string_to_bool


def _split_delimited_float(string):
    if ',' in string:
        return [float(x.strip()) for x in string.split(',')]
    else:
        return [float(x.strip()) for x in string.split()]


class ThinCurr_coil_set:
    def __init__(self, name, sens_mask=False, resistivity_per_length=None, radius=None):
        """! Initialize general coil set
        @param name Name of the coil set
        @param sens_mask Masking flag for coil (if `True`, coil is masked from sensors)
        @param resistivity_per_length Resistivity per unit length (required if used as Vcoil)
        @param radius Coil radius (required if used as Vcoil)
        """
        self.name = name
        self.sens_mask = sens_mask
        self.subcoils = []  # List of subcoil dictionaries
        self.resistivity_per_length = resistivity_per_length
        self.radius = radius
        self.resistivity_per_length_list = []
        self.radius_list = []

    def add_circular_subcoil(self, RZ, scale=None, npoints=None, attr=None, resistivity_per_length=None, radius=None):
        """! Add a circular subcoil to this coil set
        @param RZ [R, Z] position (for circular coils; centered on Z-axis)
        @param scale Scaling factor for coil current (`1.0` if not specified)
        @param npoints Number of points for circular coil discretization (default: 181)
        @param resistivity_per_length Resistivity per unit length (required if used as Vcoil)
        @param radius Coil radius (required if used as Vcoil)
        """
        if npoints is None:
            npoints = 181  # Default number of points for circular coil
        theta = numpy.linspace(0, 2*numpy.pi, npoints)
        pts = numpy.column_stack([RZ[0] * numpy.cos(theta), RZ[0] * numpy.sin(theta), numpy.full(npoints, RZ[1])])
        self.add_subcoil(pts=pts, scale=scale, attr=attr, resistivity_per_length=resistivity_per_length, radius=radius)

    def add_subcoil(self, RZ=None, pts=None, scale=None, npoints=None, hdf5_path=None, attr=None, resistivity_per_length=None, radius=None):
        """! Add a subcoil to this V-coil set
        @param RZ [R, Z] position (for circular coils; centered on Z-axis)
        @param pts Array of [x, y, z] positions (for general coils)
        @param scale Scaling factor for coil current (`1.0` if not specified)
        @param npoints Number of points
        @param hdf5_path Path to HDF5 dataset containing coil points (alternative to `pts` and `RZ`)
        @param resistivity_per_length Resistivity per unit length (required if used as Vcoil)
        @param radius Coil radius (required if used as Vcoil)
        """
        if RZ is not None: # Circular coil
            if (pts is not None) or (hdf5_path is not None):
                raise ValueError("Only one of `RZ`, `pts`, or `hdf5_path` should be provided")
            self.add_circular_subcoil(RZ=RZ, scale=scale, npoints=npoints, attr=attr, resistivity_per_length=resistivity_per_length, radius=radius)
            return
        if attr is None:
            attr = {}
        if scale is not None:
            attr['scale'] = scale
        if resistivity_per_length is not None:
            attr['res_per_len'] = resistivity_per_length
        if radius is not None:
            attr['radius'] = radius
        if pts is not None: # General 3D coil
            if (RZ is not None) or (hdf5_path is not None):
                raise ValueError("Only one of `RZ`, `pts`, or `hdf5_path` should be provided")
            pts = numpy.asarray(pts)
            if npoints is None:
                npoints = pts.shape[0]
            self.subcoils.append({'pts': pts, 'npts': npoints, 'attr': attr})
        elif hdf5_path is not None: # Coil points from HDF5 dataset
            if (RZ is not None) or (pts is not None):
                raise ValueError("Only one of `RZ`, `pts`, or `hdf5_path` should be provided")
            self.subcoils.append({'hdf5_path': hdf5_path, 'attr': attr})
        else:
            raise ValueError("Either `RZ`, `pts`, or `hdf5_path` must be provided")
        self.resistivity_per_length_list.append(attr.get('res_per_len'))
        self.radius_list.append(attr.get('radius'))

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
            elif 'hdf5_path' in subcoil:
                coil_element = ET.SubElement(coil_set, "coil", attrib={"path": subcoil['hdf5_path']})
            # Handling attributes
            if subcoil.get('attr') is not None:
                for key, value in subcoil['attr'].items():
                    coil_element.set(key, str(value))

    def build_XML(self, parent_tag):
        """! Build XML structure for this V-coil
        @param parent_tag Parent XML element to attach to
        """
        if len(self.subcoils) == 0:
            raise ValueError("Coil set must contain at least one subcoil")
        # Build attributes for coil_set
        attrib = {}
        if self.name is not None:
            attrib["name"] = self.name
        if self.resistivity_per_length is not None:
            attrib["res_per_len"] = '{0:.6E}'.format(self.resistivity_per_length)
        if self.radius is not None:
            attrib["radius"] = '{0:.6E}'.format(self.radius)
        if self.sens_mask:
            attrib["sens_mask"] = bool_to_string(self.sens_mask)
        coil_set = ET.SubElement(parent_tag, "coil_set", attrib=attrib)
        # Build subcoil elements
        self.build_coil_XML(coil_set)

    @staticmethod
    def load_from_xml(coil_elem, fallback_name):
        """! Load coil from XML element
        @param coil_elem XML element representing the coil_set
        @param fallback_name Fallback name for the coil_set
        """
        name = coil_elem.attrib.get("name", fallback_name)
        sens_mask = string_to_bool(coil_elem.attrib.get("sens_mask", "false"))
        radius = float(coil_elem.attrib["radius"]) if "radius" in coil_elem.attrib else None
        resistivity_per_length = float(coil_elem.attrib["res_per_len"]) if "res_per_len" in coil_elem.attrib else None
        coil_set = ThinCurr_coil_set(name=name, sens_mask=sens_mask, resistivity_per_length=resistivity_per_length, radius=radius)
        for subcoil_elem in coil_elem.findall("coil"):
            subcoil_scale = float(subcoil_elem.attrib["scale"]) if "scale" in subcoil_elem.attrib else None
            subcoil_radius = float(subcoil_elem.attrib["radius"]) if "radius" in subcoil_elem.attrib else None
            subcoil_resistivity = float(subcoil_elem.attrib["res_per_len"]) if "res_per_len" in subcoil_elem.attrib else None
            if "path" in subcoil_elem.attrib:
                hdf5_path = subcoil_elem.attrib["path"]
                coil_set.add_subcoil(hdf5_path=hdf5_path, scale=subcoil_scale, resistivity_per_length=subcoil_resistivity, radius=subcoil_radius)
            else:
                npts = int(subcoil_elem.attrib.get("npts", 0))
                pts_text = subcoil_elem.text.strip()
                pts_lines = pts_text.splitlines()
                pts = numpy.array([_split_delimited_float(line) for line in pts_lines])
                if (npts == 0) or (pts.shape[0] == 1):
                    coil_set.add_circular_subcoil(RZ=pts[0], scale=subcoil_scale, resistivity_per_length=subcoil_resistivity, radius=subcoil_radius)
                else:
                    coil_set.add_subcoil(pts=pts, scale=subcoil_scale, npoints=npts, resistivity_per_length=subcoil_resistivity, radius=subcoil_radius)
        return coil_set

    def save_hdf5(self, h5_group, fallback_name):
        """! Save coil information to HDF5 group
        @param h5_group HDF5 group to save the coil information into
        @param fallback_name Fallback name for the coil_set
        """
        if len(self.subcoils) == 0:
            raise ValueError("Coil set must contain at least one subcoil")
        name = self.name if self.name is not None else fallback_name
        name_len = len(name)
        h5_group.create_dataset('NAME', data=numpy.array(name.encode('ascii'), dtype=f"S{name_len}"))
        h5_group.create_dataset('NCOILS', data=[len(self.subcoils),], dtype='i4')
        if self.sens_mask:
            h5_group.create_dataset('SENS_MASK', data=[1,], dtype='i4')
        # Save resistivity and radius information if available
        radius = [r if r is not None else self.radius for r in self.radius_list]
        resistivity_per_length = [r if r is not None else self.resistivity_per_length for r in self.resistivity_per_length_list]
        if any(r is None for r in radius) != any(r is None for r in resistivity_per_length):
            raise ValueError("Both radius and resistivity_per_length must be specified or neither must be specified.")
        if any(r is not None for r in radius):
            if all(r is not None for r in radius):
                h5_group.create_dataset('RADIUS', data=radius, dtype='f8')
            else:
                raise ValueError("Radius specified on only some subcoils, and no coilset value provided.")
        if any(r is not None for r in resistivity_per_length):
            if all(r is not None for r in resistivity_per_length):
                h5_group.create_dataset('RES_PER_LEN', data=resistivity_per_length, dtype='f8')
            else:
                raise ValueError("Resistivity per length specified on only some subcoils, and no coilset value provided.")
        # Save subcoil information
        scales = numpy.ones((len(self.subcoils),))
        for i, subcoil in enumerate(self.subcoils):
            subcoil_group = h5_group.create_group('coil{0:04d}'.format(i+1))
            if 'pts' in subcoil:
                subcoil_group.create_dataset('PTS', data=subcoil['pts'], dtype='f8')
            elif 'hdf5_path' in subcoil:
                source_filepath = subcoil['hdf5_path'].split(':')[0]
                source_dataset = subcoil['hdf5_path'].split(':')[1]
                with h5py.File(source_filepath, 'r') as source_file:
                    subcoil_group.create_dataset('PTS', data=source_file[source_dataset][()], dtype='f8')
            if subcoil.get('attr', {}).get('scale') is not None:
                scales[i] = subcoil['attr']['scale']
        h5_group.create_dataset('SCALES', data=scales, dtype='f8')


class ThinCurr_Icoil(ThinCurr_coil_set):
    """! I-coil class for defining current-driven coils in OpenFUSIONToolkit"""

    def __init__(self, name, sens_mask=False):
        """! Initialize I-coil
        @param name Name of the coil set
        @param sens_mask Masking flag for coil (if `True`, coil is masked from sensors)
        """
        super().__init__(name, sens_mask)


class ThinCurr_Vcoil(ThinCurr_coil_set):
    """! V-coil class for defining voltage-driven coils in OpenFUSIONToolkit"""

    def __init__(self, name, resistivity_per_length, radius, sens_mask=False):
        """! Initialize V-coil
        @param name Name of the coil set
        @param resistivity_per_length Resistivity per unit length
        @param radius Coil radius
        @param sens_mask Masking flag for coil (if `True`, coil is masked from sensors)
        """
        super().__init__(name, sens_mask=sens_mask, resistivity_per_length=resistivity_per_length, radius=radius)


def coil_from_hdf5(h5_group):
    """! Load coil information from HDF5 group
    @param h5_group HDF5 group containing the coil information
    """
    name = h5_group['NAME'][()]
    if isinstance(name, bytes):
        name = name.decode('ascii')
    sens_mask = bool(h5_group['SENS_MASK'][0]) if 'SENS_MASK' in h5_group else None
    if 'RES_PER_LEN' in h5_group:
        if 'RADIUS' not in h5_group:
            raise ValueError("RADIUS dataset must be present in HDF5 group if RES_PER_LEN is present")
        coil = ThinCurr_Vcoil(name=name,
                                resistivity_per_length=h5_group['RES_PER_LEN'][0],
                                radius=h5_group['RADIUS'][0],
                                sens_mask=sens_mask)
        resistivity_per_length_list = h5_group['RES_PER_LEN'][()]
        radius_list = h5_group['RADIUS'][()]
    else:
        coil = ThinCurr_Icoil(name=name,sens_mask=sens_mask)
        resistivity_per_length_list = None
        radius_list = None
    for i in range(h5_group['NCOILS'][0]):
        subcoil_group = h5_group['coil{0:04d}'.format(i+1)]
        scale = h5_group['SCALES'][i] if 'SCALES' in h5_group else None
        if resistivity_per_length_list is not None:
            coil.add_subcoil(pts=subcoil_group['PTS'][()],scale=scale,
                             resistivity_per_length=resistivity_per_length_list[i], radius=radius_list[i])
        else:
            coil.add_subcoil(pts=subcoil_group['PTS'][()],scale=scale)
    return coil


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

    @staticmethod
    def load(filename):
        """! Load ThinCurr information from XML file
        @param filename Path to the XML file
        """
        tree = ET.parse(filename)
        root = tree.getroot()
        thincurr = ThinCurr_XML()
        eta_vol_elem = root.find("./thincurr/eta_vol")
        if eta_vol_elem is not None:
            thincurr.set_eta_vol(_split_delimited_float(eta_vol_elem.text))
        eta_surf_elem = root.find("./thincurr/eta_surf")
        if eta_surf_elem is None:
            eta_surf_elem = root.find("./thincurr/eta")
        if eta_surf_elem is not None:
            thincurr.set_eta(_split_delimited_float(eta_surf_elem.text))
        thickness_elem = root.find("./thincurr/thickness")
        if thickness_elem is not None:
            thincurr.set_thickness(_split_delimited_float(thickness_elem.text))
        for icoil_elem in root.findall("./thincurr/icoils/coil_set"):
            thincurr.add_Icoil(ThinCurr_Icoil.load_from_xml(icoil_elem, fallback_name="Icoil{0}".format(len(thincurr.icoils)+1)))
        for vcoil_elem in root.findall("./thincurr/vcoils/coil_set"):
            thincurr.add_Vcoil(ThinCurr_Vcoil.load_from_xml(vcoil_elem, fallback_name="Vcoil{0}".format(len(thincurr.vcoils)+1)))
        return thincurr

    def add_Icoil(self, icoil):
        """! Add an I-coil to this ThinCurr block
        @param icoil ThinCurr_Icoil object
        """
        if not isinstance(icoil, ThinCurr_coil_set) or isinstance(icoil, ThinCurr_Vcoil):
            raise TypeError("Icoil must be of type ThinCurr_Icoil or ThinCurr_coil_set")
        self.icoils.append(icoil)

    def add_Vcoil(self, vcoil):
        """! Add a V-coil to this ThinCurr block
        @param vcoil ThinCurr_Vcoil object
        """
        if isinstance(vcoil, ThinCurr_coil_set):
            radius = [r if r is not None else vcoil.radius for r in vcoil.radius_list]
            resistivity_per_length = [r if r is not None else vcoil.resistivity_per_length for r in vcoil.resistivity_per_length_list]
            if any(r is None for r in radius) or any(r is None for r in resistivity_per_length):
                raise ValueError("Vcoil must have resistivity_per_length and radius defined for all subcoils or at the coilset level")
        else:
            raise TypeError("Vcoil must be of type ThinCurr_Vcoil or ThinCurr_coil_set (with resistivity_per_length and radius defined)")
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
            eta_element.text = ' '.join('{0:.6E}'.format(r) for r in self.eta)

        # Add eta_vol (volumetric resistivities) if present
        if self.eta_vol is not None:
            eta_vol_element = ET.SubElement(thincurr_element, "eta_vol")
            eta_vol_element.text = ' '.join('{0:.6E}'.format(r) for r in self.eta_vol)

        # Add thickness if present
        if self.thickness is not None:
            thickness_element = ET.SubElement(thincurr_element, "thickness")
            thickness_element.text = ' '.join('{0:.6E}'.format(t) for t in self.thickness)

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