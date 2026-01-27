import json
import numpy
import h5py
from matplotlib.path import Path
import matplotlib as plt
import xml.etree.ElementTree as ET
import sys
import inspect

'''! xml editor for adding filament coil points

@authors Jamie Xia
@date Jan 2026

'''
# ===============================
# XML utilities
# ===============================
def setup_filaments(resolution, fil_points_txt, xml_in, xml_out,  filename = None, tolerance = 1e-1, jsonfile = None, verbose = False, oft_env = None, tokamaker_meshfile = None): 
    ''' ! Sets up the filament points in an existing oft xml file.
    @param resolution : int, number of filaments in R and Z direction (resolution x
    resolution grid)
    @param fil_points_txt : string, path to temporary text file to store filament points
    @param xml_in : string, path to input oft xml file
    @param xml_out : string, path to output oft xml file with filaments added
    @param filename : string (optional), path to text file containing filament points           
    @param tolerance : float (optional), minimum distance from filament to boundary
    @param jsonfile : string (optional), path to json file containing cross section polygon
    @param verbose : bool (optional), if True, plot the filaments and limiter
    @param oft_env : string (optional), path to OpenFUSIONToolkit environment for plotting
    @param tokamaker_meshfile : string (optional), path to tokamaker mesh
    @returns points : array of (r, z) filament points
    @returns xml_out : string, path to output oft xml file with filaments added
    @returns R : array of R filament points
    @returns Z : array of Z filament points
    '''
    
    if jsonfile is None:
        raise ValueError("jsonfile must be provided to construct filaments")
    R, Z = construct_filaments_from_json(jsonfile, resolution, resolution, tol=tolerance)
    points = numpy.column_stack([R, Z])
    number_of_fils = len(points)
    numpy.savetxt(fil_points_txt, points, fmt='%.6f', delimiter= " ")
    edit_xml = ET.parse(xml_in)
    root = edit_xml.getroot()
    icoil_tag = root.find('./thincurr/icoils')
    if filename is None: # have the option to implement custom filament array 
        filename = fil_points_txt
    points = get_fil_points(filename)
    for r, z in points:
        coil_element = ET.SubElement(icoil_tag, 'coil_set')
        coil = ET.SubElement(coil_element, "coil")
        coil.set("scale", "1.0") #default scale value is 1.0
        coil.text = '{0:.6f}, {1:.6f}'.format(r, z)
    ET.indent(edit_xml, space="  ", level=0)
    edit_xml.write(xml_out, encoding="utf-8", xml_declaration=True) 
    print("XML file updated with filament points.") 
    if verbose: 
        from OpenFUSIONToolkit.TokaMaker import TokaMaker
        from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
        import matplotlib.pyplot as plt
        if oft_env is None or tokamaker_meshfile is None: 
            raise ValueError("oft_env and tokamaker_meshfile must be provided for verbose plotting")
        mygs = TokaMaker(oft_env)
        mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh(tokamaker_meshfile)
        mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
        mygs.setup(order = 2, F0=0.2518*1.23)
        limiter = mygs.lim_contour
        plt.plot(limiter[:,0],limiter[:,1],color='k')
        plt.scatter(R, Z, cmap='plasma', marker='o', edgecolors='k')
        plt.xlabel('r')
        plt.ylabel('z')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
    print(f"Constructed {number_of_fils} filaments and written to {xml_out}")
    return points, xml_out, R, Z
    
def get_fil_points(file):
    ''' ! Reads in a text file of filament points
    @param file : string, path to text file
    @returns points : list of (r, z) tuples
    '''
    points = []
    with open(file, "r") as f:
        for line in f:
            if line.strip():  # skip empty lines
                try:
                    r, z = map(float, line.split())
                    points.append((r, z))
                except ValueError:
                    print(f"Skipping invalid line: {line.strip()}")
    return points

# ===============================
# Geometry helpers
# ===============================
def get_cross_section(meshfile, jsonfile=None): 
    ''' ! Grabs the cross section from either a meshfile of json file
    @param meshfile .h5 file, a meshfile can be made in tokamaker of your machine (tokamak) 
    @param jsonfile (optional): a .json file, also can be made in tokamaker.  
    @returns R and Z arrays. When zipped together can be plotted as a cross section.  
    '''
    if jsonfile is None: 
        with h5py.File(meshfile, 'r') as f: 
            coords = f['mesh/R'][:]
            if coords.shape[1] == 3:  # Checks if it is in (X, Y, Z)
                R = numpy.sqrt(coords[:, 0]**2 + coords[:, 1]**2)  # Converts to R
                Z = coords[:, 2]  # Keeps the Z
            else:
                R, Z = coords[:, 0], coords[:, 1] 
    else: 
         with open(jsonfile,'r') as fid:
            geom = json.load(fid)
            R, Z = numpy.array(geom['limiter']).T  # Transpose to separate columns
    return R, Z #returns the cross section 

def point_to_segment_dist(px, py, x0, y0, x1, y1):
    """! Distance from point (px,py) to segment (x0,y0)-(x1,y1).
    @param px, py : float Point coordinates.
    @param x0, y0 : float Segment start coordinates
    @param x1, y1 : float Segment end coordinates
    @returns float Minimum distance from point to segment.
    """
    dx, dy = x1 - x0, y1 - y0
    if dx == dy == 0:
        # segment is a point
        return numpy.hypot(px - x0, py - y0)
    t = ((px - x0) * dx + (py - y0) * dy) / (dx*dx + dy*dy)
    t = numpy.clip(t, 0, 1)
    proj_x = x0 + t*dx
    proj_y = y0 + t*dy
    return numpy.hypot(px - proj_x, py - proj_y)

# ===============================
# Geometry / filament generation
# ===============================

def construct_filaments(meshfile, jsonfile, rows, columns): # make sure the resolution matches the g file resolution
    """ ! Constructs a grid of filaments inside a cross section polygon
    @param meshfile : str, path to tokamaker .h5 mesh file
    @param jsonfile : str, path to JSON file containing cross-section polygon under key 'limiter'. Can be made from tokamaker. 
    @param rows : int, number of filaments in R direction
    @param columns : int, number of filaments in Z direction
    @returns R_filaments, Z_filaments : arrays of float
    """
    R,Z = get_cross_section(meshfile, jsonfile) # grabs the cross section
    R_min, R_max = R.min(), R.max() # identifies the min and max r and z values (getting the size of the cross section) 
    Z_min, Z_max = Z.min(), Z.max() 
    R_grid = numpy.linspace(R_min, R_max, rows) 
    Z_grid = numpy.linspace(Z_min, Z_max, columns)
    R_mesh, Z_mesh = numpy.meshgrid(R_grid, Z_grid) # creates the grid of filaments based of resolution in R, Z space
    cross_section_boundary = numpy.column_stack((R, Z))
    centroid = numpy.mean(cross_section_boundary, axis=0) # the next few lines is organizing the filaments 
    angles = numpy.arctan2(cross_section_boundary[:, 1] - centroid[1], 
                        cross_section_boundary[:, 0] - centroid[0])                                                           
    sorted_indices = numpy.argsort(angles) # sorts points by angle around the centroid 
    sorted_boundary = cross_section_boundary[sorted_indices]                              
    boundary_path = Path(sorted_boundary) #creates a Path object from sorted points
    grid_points = numpy.column_stack((R_mesh.ravel(), Z_mesh.ravel())) #flattens the grid points 
    inside = boundary_path.contains_points(grid_points, None, -0.1) #checks if they are inside the cross section
    inside_mask = inside.reshape(R_mesh.shape)   # reshapes to match meshgrid shape
    return R_mesh[inside_mask], Z_mesh[inside_mask] # returns list of R and Z values that when put together make the grid of filaments


def construct_filaments_from_json(jsonfile, nR=10, nZ=10, tol=1e-1): 
    """ !   Path to JSON file containing cross-section polygon under key 'limiter'. Can be made from tokamaker. 
    Generate a grid of (R,Z) points inside a cross-section polygon and at least `tol` away from the boundary.
    @param jsonfile : str
    @returns R_points, Z_points : arrays of float
    """
    with open(jsonfile,'r') as fid:
        geom = json.load(fid)
        R_boundary, Z_boundary = numpy.array(geom['limiter']).T

    # Close polygon if needed
    if R_boundary[0] != R_boundary[-1] or Z_boundary[0] != Z_boundary[-1]:
        R_boundary = numpy.append(R_boundary, R_boundary[0])
        Z_boundary = numpy.append(Z_boundary, Z_boundary[0])

    polygon = Path(numpy.column_stack([R_boundary, Z_boundary]))

    # Bounding box
    Rmin, Rmax = numpy.min(R_boundary), numpy.max(R_boundary)
    Zmin, Zmax = numpy.min(Z_boundary), numpy.max(Z_boundary)

    # Candidate grid
    Rg = numpy.linspace(Rmin, Rmax, nR)
    Zg = numpy.linspace(Zmin, Zmax, nZ)
    RR, ZZ = numpy.meshgrid(Rg, Zg)
    points = numpy.column_stack([RR.ravel(), ZZ.ravel()])

    # Mask points inside polygon
    mask = polygon.contains_points(points)
    inside_points = points[mask]

    if tol > 0.0:
        keep = []
        for r, z in inside_points:
            # compute min distance to any polygon edge
            min_dist = numpy.min([
                point_to_segment_dist(r, z, R_boundary[i], Z_boundary[i],
                                      R_boundary[i+1], Z_boundary[i+1])
                for i in range(len(R_boundary)-1)
            ])
            if min_dist >= tol:
                keep.append([r, z])
        inside_points = numpy.array(keep)

    R_points = inside_points[:,0]
    Z_points = inside_points[:,1]

    return R_points, Z_points

# ===============================
# XML class definitions
# ===============================

class Vcoil:
    '''! Vcoil class for defining V-coils in OpenFUSIONToolkit'''
    def __init__(self, name, R=None, Z = None, xyz = None, scale = None, resistivity_per_length = None, sens_mask = None, npoints = None, iscustom = False):
        '''! Vcoil class for defining V-coils in OpenFUSIONToolkit
        @param name : string, name of the coil
        @param R : array of float, radial positions of coil points (for non-custom coils)
        @param Z : array of float, vertical positions of coil points (for non-custom coils)
        @param xyz : array of [x, y, z] positions of coil points (for custom coils)
        @param scale : float, scaling factor for coil current
        @param resistivity_per_length : float, resistivity per unit length of the coil (for custom coils)
        @param sens_mask : array of int, sensitivity mask for coil points 
        @param npoints : int, number of points in the coil (for custom coils)
        @param iscustom : bool, True if custom coil defined by xyz points, False if defined by R,Z
        '''
        self.name = name
        self.resistivity_per_length = resistivity_per_length
        self.sens_mask = sens_mask
        self.scale = scale
        self.npoints = npoints
        self.iscustom = iscustom
        if iscustom: 
            if xyz is not None:
                self.xyz = xyz
                self.R = R 
            else:
                raise ValueError("Custom coil must have a valid xyz = [x, y, z]")
        else: 
            self.xyz = None 
            self.R = R
            self.Z = Z

        
class Icoil: 
    '''! Icoil class for defining I-coils in OpenFUSIONToolkit'''
    def __init__(self, name, R = None, Z=None, xyz = None, scale = None, resistivity_per_length = None, sens_mask = None, npoints = None, iscustom = False):
        '''! Icoil class for defining I-coils in OpenFUSIONToolkit
        @param name : string, name of the coil
        @param R : array of float, radial positions of coil points (for non-custom coils)
        @param Z : array of float, vertical positions of coil points (for non-custom coils)
        @param xyz : array of [x, y, z] positions of coil points (for custom coils)
        @param scale : float, scaling factor for coil current
        @param resistivity_per_length : float, resistivity per unit length of the coil (for custom coils)
        @param sens_mask : array of int, sensitivity mask for coil points 
        @param npoints : int, number of points in the coil (for custom coils)
        @param iscustom : bool, True if custom coil defined by xyz points, False if defined by R,Z
        '''
        self.name = name
        self.sens_mask = sens_mask
        self.scale = scale
        self.npoints = npoints
        self.resistivity_per_length = resistivity_per_length
        self.iscustom = iscustom
        if iscustom: 
            if xyz is not None :
                self.xyz = xyz
                self.R = R 
            else:
                raise ValueError("Custom coil must have a valid xyz = [x, y, z]")
        else: 
            self.xyz = None 
            self.R = R
            self.Z = Z 


# ===============================
# Create XML function 
# ===============================

def create_xml(path, icoils, vcoils, resistivities): 
    '''! Function solves filament and wall currents over time 
    @param path : string, name and location of where file will be store i.e '/home/user/folder/oft.xml' 
    @param icoils : icoils, should be a list of icoil objects [icoil1, icoil2, icoil3, ...]       
    @param vcoils : vcoils, should be a list of vcoil objects [vcoil1, vcoil2, vcoil3 , .... ]
    @param resistivities : float array should be list of resistivites by section [r1, r2, r3...] 
    @returns oft.xml file is created in specified directory. 
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
                if i.iscustom: 
                    coil_set = ET.SubElement(icoil_element, "coil_set", attrib = {"name" : i.name})
                    for xyz in i.xyz: 
                        coil_element = ET.SubElement(coil_set, "coil", attrib={"npoints": str(i.npoints)})
                        if i.scale is None:                 
                            coil_element.set("scale","1.0")
                        else:     
                            coil_element.set("scale", str(i.scal))
                        coil_element.text = '{0:.6f}, {1:.6f}'.format(xyz[0], xyz[1], xyz[2])
                else: 
                    coil_set = ET.SubElement(icoil_element, "coil_set", attrib = {"name" : i.name})
                    for r_val, z_val in zip(i.R, i.Z): 
                        coil_element = ET.SubElement(coil_set, "coil")
                        if i.scale is None:                 
                            coil_element.set("scale","1.0")
                        else:     
                            coil_element.set("scale", i.scale)
                        coil_element.text = '{0:.6f}, {1:.6f}'.format(r_val, z_val)
            else: 
                raise ValueError("Non I-coil class in inputted I-coil list")
    if vcoils: 
        vcoil_element = ET.SubElement(thincurr_element, "vcoils")
        for v in vcoils: 
            if isinstance(v, Vcoil): 
                if v.iscustom: 
                    print(f'{v.name} and {v.R} and {v.xyz} and {v.resistivity_per_length} and {v.sens_mask} and {v.npoints} and {v.scale}')
                    coil_set = ET.SubElement(vcoil_element, "coil_set", attrib = {"name" : v.name, "R" : str(v.R), "res_per_len" : str(v.resistivity_per_length), "sens_mask" : str(v.sens_mask)})
                    coil_element = ET.SubElement(coil_set, "coil", attrib={"npoints": str(v.npoints)})
                    text_lines = []
                    for xyz in v.xyz: 
                        text_lines.append('{0:.6f}, {1:.6f}, {2:.6f}'.format(xyz[0], xyz[1], xyz[2]))
                    coil_element.text = "\n" + "\n".join(text_lines) + "\n"
                else:
                    coil_set = ET.SubElement(vcoil_element, "coil_set", attrib = {"name" : v.name})
                    for r_val, z_val in zip(v.R, v.Z): 
                        coil_element = ET.SubElement(coil_set, "coil")
                        if v.scale is None:                 
                            coil_element.set("scale","1.0")
                        else:     
                            coil_element.set("scale", str(v.scale))
                        coil_element.text = '{0:.6f}, {1:.6f}'.format(r_val, z_val)
            else: 
                raise ValueError("Non V-coil class in inputted V-coil list")
    ET.indent(xml_doc, space="  ", level=0)
    xml_doc.write(path)
    print("xml file created") 
            


# ===============================
# Plasma filament utilities
# =============================== 

def setup_synthetic_current(timepoints, ip_list, sigma_r, sigma_z, r0, z0, rmesh, zmesh):
    """! Sets up synthetic coil current data based on Gaussian distributions centered at (r0, z0) over time.
    @param timepoints: list of float Time points.
    @param ip_list: list of float Plasma current values at each time point.
    @param sigma_r: float Standard deviation in the radial direction.
    @param sigma_z: float Standard deviation in the vertical direction.
    @param r0: list of float Radial center positions at each time point.
    @param z0: list of float Vertical center positions at each time point.
    @param rmesh: np.ndarray Radial mesh grid.
    @param zmesh: np.ndarray Vertical mesh grid.
    @returns coil_curr: np.ndarray Shape (ntimes, nsensors+1) First column is time, remaining columns are sensor currents.
    """
    coil_curr = [] 
    for i in range(len(timepoints)): 
        gaussian_raw = numpy.exp(-((rmesh - r0[i])**2 / (2 * sigma_r**2) + (zmesh - z0[i])**2 / (2 * sigma_z**2)))
        gaussian_values = ip_list[i]*(gaussian_raw/numpy.sum(gaussian_raw))
        coil_curr.append(gaussian_values) 
    coil_curr = numpy.array(coil_curr)             
    time_column = numpy.array(timepoints).reshape(-1, 1)
    coil_curr = numpy.hstack((time_column, coil_curr))
    return coil_curr


def interpolate_total_current_test(coil_currs, nsteps, verbose=False): 
    """! Interpolates total current from coil sensor data to a higher-resolution time grid.
    @param coil_currs: np.ndarray, shape (ntimes, nsensors+1) irst column is time, remaining columns are sensor currents.
    @param nsteps: int Number of high-resolution steps between the first and last time.
    @returns high_res_time: np.ndarray
    @returns total_current_high_res: np.ndarray
    """
        
    times = coil_currs[:, 0]
    sensor_currents = coil_currs[:, 1:]
    total_current = numpy.sum(sensor_currents, axis=1)
    high_res_time = numpy.linspace(times[0], times[-1], nsteps + 1)
    total_current_high_res = numpy.interp(
        high_res_time,
        times,
        total_current
    )
    if verbose:
        plt.figure(figsize=(8, 5))
        plt.scatter(times, total_current, color='red', label='Original Data')
        plt.plot(high_res_time, total_current_high_res, color='blue',
                 label='Interpolated Total Current')
        plt.xlabel('Time [s]')
        plt.ylabel('Total Current [A]')
        plt.title('Total Current vs Time with Higher Resolution')
        plt.legend()
        plt.grid(True)
        plt.show()

    return high_res_time, total_current_high_res
