import json
import numpy as np
import h5py


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
    if jsonfile == None: 
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

def construct_filaments(jsonfile, rows, columns): # make sure the resolution matches the g file resolution
    R,Z = self.get_cross_section(jsonfile) # grabs the cross section
    R_min, R_max = R.min(), R.max() # identifies the min and max r and z values (getting the size of the cross section) 
    Z_min, Z_max = Z.min(), Z.max() 
    R_grid = np.linspace(R_min, R_max, row) 
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
