import json
import numpy as np
import h5py
from matplotlib.path import Path
import xml.etree.ElementTree as ET
import sys

sys.path.insert(0, '/Applications/OpenFUSIONToolkit/python')


# ===============================
# XML utilities
# ===============================
def setup_filaments(resolution, fil_points_txt, xml_in, xml_out,  filename = None, tolerance = 1e-1, jsonfile = None, verbose = False, oft_env = None, tokamaker_meshfile = None): 
    if jsonfile is None:
        raise ValueError("jsonfile must be provided to construct filaments")
    R, Z = construct_filaments_from_json(jsonfile, resolution, resolution, tol=tolerance)
    points = np.column_stack([R, Z])
    number_of_fils = len(points)
    np.savetxt(fil_points_txt, points, fmt='%.6f', delimiter= " ")
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
        vv = mygs.r
        plt.plot(limiter[:,0],limiter[:,1],color='k')
        r0, z0 = 0.0, 0.0
        sigma_r, sigma_z = 1, 1 
        plt.scatter(R, Z, cmap='plasma', marker='o', edgecolors='k')
        plt.xlabel('r')
        plt.ylabel('z')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
    print(f"Constructed {number_of_fils} filaments and written to {xml_out}")
    return points, xml_out, R, Z
    
def get_fil_points(file):
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
    '''
    grabs the cross section from either a meshfile of json file
    
    Parameters:
    ----------
    meshfile : .h5 file, a meshfile can be made in tokamaker of your machine (tokamak) 

    jsonfile (optional): a .json file, also can be made in tokamaker.  
        
    Returns:
    ----------
    R and Z arrays. When zipped together can be plotted as a cross section.  
    '''
    if jsonfile is None: 
        with h5py.File(meshfile, 'r') as f: 
            coords = f['mesh/R'][:]
            if coords.shape[1] == 3:  # Checks if it is in (X, Y, Z)
                R = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)  # Converts to R
                Z = coords[:, 2]  # Keeps the Z
            else:
                R, Z = coords[:, 0], coords[:, 1] 
    else: 
         with open(jsonfile,'r') as fid:
            geom = json.load(fid)
            R, Z = np.array(geom['limiter']).T  # Transpose to separate columns
    return R, Z #returns the cross section 

def point_to_segment_dist(px, py, x0, y0, x1, y1):
    """Distance from point (px,py) to segment (x0,y0)-(x1,y1)."""
    dx, dy = x1 - x0, y1 - y0
    if dx == dy == 0:
        # segment is a point
        return np.hypot(px - x0, py - y0)
    t = ((px - x0) * dx + (py - y0) * dy) / (dx*dx + dy*dy)
    t = np.clip(t, 0, 1)
    proj_x = x0 + t*dx
    proj_y = y0 + t*dy
    return np.hypot(px - proj_x, py - proj_y)

# ===============================
# Geometry / filament generation
# ===============================

def construct_filaments(meshfile, jsonfile, rows, columns): # make sure the resolution matches the g file resolution
    R,Z = get_cross_section(meshfile, jsonfile) # grabs the cross section
    R_min, R_max = R.min(), R.max() # identifies the min and max r and z values (getting the size of the cross section) 
    Z_min, Z_max = Z.min(), Z.max() 
    R_grid = np.linspace(R_min, R_max, rows) 
    Z_grid = np.linspace(Z_min, Z_max, columns)
    R_mesh, Z_mesh = np.meshgrid(R_grid, Z_grid) # creates the grid of filaments based of resolution in R, Z space
    cross_section_boundary = np.column_stack((R, Z))
    centroid = np.mean(cross_section_boundary, axis=0) # the next few lines is organizing the filaments 
    angles = np.arctan2(cross_section_boundary[:, 1] - centroid[1], 
                        cross_section_boundary[:, 0] - centroid[0])                                                           
    sorted_indices = np.argsort(angles) # sorts points by angle around the centroid 
    sorted_boundary = cross_section_boundary[sorted_indices]                              
    boundary_path = Path(sorted_boundary) #creates a Path object from sorted points
    grid_points = np.column_stack((R_mesh.ravel(), Z_mesh.ravel())) #flattens the grid points 
    inside = boundary_path.contains_points(grid_points, None, -0.1) #checks if they are inside the cross section
    inside_mask = inside.reshape(R_mesh.shape)   # reshapes to match meshgrid shape
    return R_mesh[inside_mask], Z_mesh[inside_mask] # returns list of R and Z values that when put together make the grid of filaments


def construct_filaments_from_json(jsonfile, nR=10, nZ=10, tol=1e-1): 
    """
    Input: jsonfile : str
        Path to JSON file containing cross-section polygon under key 'limiter'. Can be made from tokamaker. 
    Generate a grid of (R,Z) points inside a cross-section polygon
    and at least `tol` away from the boundary.
    """

    with open(jsonfile,'r') as fid:
        geom = json.load(fid)
        R_boundary, Z_boundary = np.array(geom['limiter']).T

    # Close polygon if needed
    if R_boundary[0] != R_boundary[-1] or Z_boundary[0] != Z_boundary[-1]:
        R_boundary = np.append(R_boundary, R_boundary[0])
        Z_boundary = np.append(Z_boundary, Z_boundary[0])

    polygon = Path(np.column_stack([R_boundary, Z_boundary]))

    # Bounding box
    Rmin, Rmax = np.min(R_boundary), np.max(R_boundary)
    Zmin, Zmax = np.min(Z_boundary), np.max(Z_boundary)

    # Candidate grid
    Rg = np.linspace(Rmin, Rmax, nR)
    Zg = np.linspace(Zmin, Zmax, nZ)
    RR, ZZ = np.meshgrid(Rg, Zg)
    points = np.column_stack([RR.ravel(), ZZ.ravel()])

    # Mask points inside polygon
    mask = polygon.contains_points(points)
    inside_points = points[mask]

    if tol > 0.0:
        keep = []
        for r, z in inside_points:
            # compute min distance to any polygon edge
            min_dist = np.min([
                point_to_segment_dist(r, z, R_boundary[i], Z_boundary[i],
                                      R_boundary[i+1], Z_boundary[i+1])
                for i in range(len(R_boundary)-1)
            ])
            if min_dist >= tol:
                keep.append([r, z])
        inside_points = np.array(keep)

    R_points = inside_points[:,0]
    Z_points = inside_points[:,1]

    return R_points, Z_points
