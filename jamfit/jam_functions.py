import numpy as np
import xml.etree.ElementTree as ET
import sys
import pandas as pd 
import os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from collections import OrderedDict

# home_dir = os.path.expanduser("~")
# oft_root_path = os.path.join(home_dir, "OpenFUSIONToolkit/install_release")
# os.environ["OFT_ROOTPATH"] = oft_root_path
# thincurr_python_path = os.getenv('OFT_ROOTPATH')
# if thincurr_python_path is not None:
#     sys.path.append(os.path.join(thincurr_python_path,'python'))


sys.path.insert(0, '/Applications/OpenFUSIONToolkit/python')
from OpenFUSIONToolkit.ThinCurr import ThinCurr
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr.sensor import Mirnov, save_sensors, circular_flux_loop
from OpenFUSIONToolkit.util import mu0
from OpenFUSIONToolkit.io import histfile
# input: input a list of arrays values that you want to find the index and value at that specified index within the inputted array at a specified time with a matching 
# time array
# output: gives the value and the index of that value within the array you are searching within
def getlist(insidelist, timearray, tofind, offset):
    indicies = [] 
    output = []
    for time in tofind:
        index = np.argmin(np.abs(timearray - (time + offset)))
        output.append(insidelist[index]) 
        indicies.append(index) 
    return output, indicies 

# goes through the xml file and grabs the names of each coil set while ignoring filaments (so fcoils and ecoils) 
def getcoilnames(xml): 
    tree = ET.parse(xml)
    root = tree.getroot()
    namelist = [] 
    for coil_set in root.findall(".//coil_set"):
        name = coil_set.get("name") 
        coil = [] 
        if name is not None: 
            namelist.append(name)
    return namelist


# applies a gassuian distribution based on an inputed ip list and spreads it with a sigmas distrubtion and centers it at r0 and z0 
# returns lisst of coil currents for filaments where the first index of every row is the time 
def setup_synthetic_current(timepoints, ip_list, sigma_r, sigma_z, r0, z0, rmesh, zmesh):
    coil_curr = [] 
    for i in range(len(timepoints)): 
        gaussian_raw = np.exp(-((rmesh - r0[i])**2 / (2 * sigma_r**2) + (zmesh - z0[i])**2 / (2 * sigma_z**2)))
        gaussian_values = ip_list[i]*(gaussian_raw/np.sum(gaussian_raw))
        coil_curr.append(gaussian_values) 
    coil_curr = np.array(coil_curr)             
    time_column = np.array(timepoints).reshape(-1, 1)
    coil_curr = np.hstack((time_column, coil_curr))
    return coil_curr
# applies a gaussian with a point beam runaway current at the edge of the plasma 
def setup_synthetic_current_runaway(timepoints, ip_list, sigma_r, sigma_z, r0, z0, sigma_r_small, sigma_z_small, r0_small, z0_small, rmesh, zmesh, A1=0.9, A2=0.2): 
    coil_curr = [] 
    for i in range(len(timepoints)):
        # Big central Gaussian
        g1 = np.exp(-(
            (rmesh - r0[i])**2 / (2 * sigma_r**2) +
            (zmesh - z0[i])**2 / (2 * sigma_z**2)
        ))
        # Smaller edge Gaussian (shifted toward edge)
        g2 = np.exp(-(
            (rmesh - r0_small[i])**2 / (2 * sigma_r_small**2) +
            (zmesh - z0_small[i])**2 / (2 * sigma_z_small**2)
        ))
        # Weighted sum of both
        gaussian_raw = A1 * g1 + A2 * g2 
        gaussian_values = ip_list[i]*(gaussian_raw/np.sum(gaussian_raw))
        coil_curr.append(gaussian_values) 
    coil_curr = np.array(coil_curr)             
    time_column = np.array(timepoints).reshape(-1, 1)
    coil_curr = np.hstack((time_column, coil_curr))
    return coil_curr
#i did not write this but i have tested it and it seems to work 
# returns the total sum of all filament currents given an inputed filament current array in time 
def interpolate_total_current(coil_currs, nsteps, verbose = False):
    """
    Interpolates total current from coil sensor data to a higher-resolution time grid.

    Parameters:
    - coil_currs: np.ndarray, shape (ntimes, nsensors+1)
        First column is time, remaining columns are sensor currents.
    - nsteps: int
        Number of high-resolution steps between the first and last time.

    Returns:
    - high_res_time: np.ndarray
    - total_current_high_res: np.ndarray
    """
    # Separate time and sensor currents
    times = coil_currs[:, 0]
    sensor_currents = coil_currs[:, 1:]

    # Generate high-resolution time grid
    high_res_time = np.linspace(times[0], times[-1], nsteps + 1)

    # Interpolate each sensor's current to high-res time
    interpolated_currents = np.array([
        interp1d(times, sensor, kind='linear')(high_res_time)
        for sensor in sensor_currents.T
    ]).T  # shape: (nsteps+1, nsensors)

    # Sum currents across sensors at each high-res time step
    total_current_high_res = np.sum(interpolated_currents, axis=1)

    # Plotting
    if verbose:
        plt.figure(figsize=(8, 5))
        plt.scatter(times, np.sum(sensor_currents, axis=1), color='red', label='Original Data')
        plt.plot(high_res_time, total_current_high_res, color='blue', label='Interpolated Total Current')
        plt.xlabel('Time [s]')
        plt.ylabel('Total Current [A]')
        plt.title('Total Current vs Time with Higher Resolution')
        plt.legend()
        plt.grid(True)
        plt.show()

    return high_res_time, total_current_high_res
# gives the r and z values of a flux loop location text file  
def read_points_flux(file_path,  sensor_names = None):
    points = []
    with open(file_path, 'r') as file:
        for line in file:
            if sensor_names is not None:
                name, r, z = line.strip().split()  # Split and ignore the first field (name)
                if name in sensor_names:
                    points.append([name,float(r), float(z)])
            else:
                _, r, z = line.strip().split()  # Split and ignore the first field (name)
                points.append([float(r), float(z)])  # Convert r and z to floats and append
    return points
# gives the name, r,z, polodial angle, and gamma angle from a mirnov magnetic location text file 
def read_points_mag(file_path, sensor_names = None):
    points = []
    with open(file_path, 'r') as file:
        for line in file:
            if sensor_names is not None: 
                name, r, z, pol, orient, _, _, _ = line.strip().split()  # Split and ignore the first field (name)
                if name in sensor_names:
                    points.append([name, float(r), float(z), float(orient)])
            else: 
                _, r, z, pol, orient, _, _, _ = line.strip().split()  # Split and ignore the first field (name)
                points.append([float(r), float(z), float(orient)])  # Convert r and z to floats and append
    return points
# adds only polodial sensors currently, returns a sensor object as well as sensor arrays for plotting 

def add_DIIID_sensors(magnetics_csv, all_sensors_dict): 
        """
        Reads a CSV of magnetic sensor locations and creates Mirnov and flux loop sensors.
        
        Parameters:
        - magnetics_csv: path to CSV containing sensor locations
        - all_sensors_dict: dict of all sensors to check against
        
        Returns:
        - sensors: list of Mirnov and flux loop sensor objects
        - plot_sensors: list of positions of sensors (Mirnov)
        - plot_sensors_flux: dict of flux loop plotting data
        """
        
        # Load CSV
        sensor_locs = pd.read_csv(magnetics_csv)
        sensor_locs["name"] = sensor_locs["name"].str.strip().str.upper()
        sensor_lookup = {row["name"]: row for _, row in sensor_locs.iterrows()} #if you need to lookup by name later
        df = pd.read_csv(magnetics_csv, encoding="utf-8-sig")
        sensor_dict = OrderedDict({col: df[col].to_numpy() for col in df.columns})
        
        # Initialize containers
        sensors = []
        plot_sensors = []
        plot_sensors_flux = {}
        orientations = []
        ordered_sensors = OrderedDict()
        
        # Order sensors based on all_sensors_dict
        for i in range(len(sensor_dict['name'])):
            name = sensor_dict['name'][i]
            if name in all_sensors_dict:
                ordered_sensors[name] = {
                    'R': sensor_dict['R'][i],
                    'Z': sensor_dict['Z'][i],
                    'type': sensor_dict['type'][i],
                    'gamma': sensor_dict['gamma'][i] if 'gamma' in sensor_dict else None,
                    'poloidal': sensor_dict['poloidal'][i] if 'poloidal' in sensor_dict else None
                }
        
        # Counters (optional)
        num_flux = 0
        num_mirnov = 0
        
        # Create sensor objects
        for s in all_sensors_dict.keys():
            s_clean = s.strip().upper()
            if s_clean in ordered_sensors:
                loc = ordered_sensors[s_clean]
                sensor_type = loc['type'].strip().lower()
                
                if sensor_type == 'mirnov':
                    tor_angle = loc['poloidal']
                    position = np.array([
                        loc['R'] * np.cos(tor_angle * np.pi / 180),
                        loc['R'] * np.sin(tor_angle * np.pi / 180),
                        loc['Z']
                    ])
                    
                    nx = np.cos(loc['gamma'] * np.pi / 180)
                    ny = 0
                    nz = np.sin(loc['gamma'] * np.pi / 180)
                    
                    cosphi = np.cos(tor_angle * np.pi / 180)
                    sinphi = np.sin(tor_angle * np.pi / 180)
                    
                    orientation = np.array([nx * cosphi - ny * sinphi,
                                            nx * sinphi + ny * cosphi,
                                            nz])
                    
                    sensors.append(Mirnov(position, orientation, name=s_clean))
                    num_mirnov += 1
                    plot_sensors.append(position)
                    orientations.append(orientation)
                    
                elif sensor_type == 'flux':
                    num_flux += 1
                    sensors.append(circular_flux_loop(loc['R'], loc['Z'], name=s_clean))
                    
                    theta = np.linspace(0.0, 2.0 * np.pi, 180)
                    x = loc['R'] * np.cos(theta)
                    y = loc['R'] * np.sin(theta)
                    z_plot = loc['Z'] * np.ones(len(y))
                    
                    plot_sensors_flux[s_clean] = {"x": x, "y": y, "z": z_plot}
        return sensors, plot_sensors, plot_sensors_flux, orientations


    
def load_data(ipfile, time_file, fcoil_file, fcoil_file_time, ecoil_file, ecoil_file_time, offset): 
    temp = np.load(ipfile)
    totalip = -temp["arr_0"] + offset

    temp = np.load(time_file)
    real_time = temp['arr_0']/1000
    temp = np.load(ecoil_file)
    temp1 = np.load(fcoil_file) 
    dict_temp = dict(temp)
    dict_temp1 = dict(temp1)
    coil_info = {**dict_temp1, **dict_temp}
    time_f = np.load(fcoil_file_time) 
    time_f = time_f['arr_0']/1E3
    time_e = np.load(ecoil_file_time) 
    time_e = time_e['arr_0']/1E3 
    return coil_info, totalip, real_time, time_f, time_e

def process_coils(coillist, coilcurr, coil_info, times, time_f, time_e): 
    numecoils = 0
    numfcoils = 0
    for coil in coillist: 
        if coil.startswith("F"):
            offset_f = time_f[0]
            indicies_f = [np.argmin(np.abs(time_f - (t + offset_f))) for t in times]
            temp_curr = coil_info[coil][:]
            indicies_f = [min(i, len(temp_curr) - 1) for i in indicies_f]  # Prevent out-of-bounds
            temp_curr = -temp_curr[indicies_f]
            newcol = np.array(temp_curr).reshape(-1, 1)
            coilcurr = np.hstack((coilcurr, newcol))
            numfcoils +=1 
        else: 
            offset_e = time_e[0]
            indicies_e = [np.argmin(np.abs(time_e - (t + offset_e))) for t in times]
            if coil.endswith("1") or coil.endswith("3") or coil.endswith("5") :
                temp_curr = coil_info['ECOILA'][:]
            else:
                temp_curr = coil_info['ECOILB'][:]
            temp_curr = -temp_curr[indicies_e] 
            newcol = np.array(temp_curr).reshape(-1, 1)
            coilcurr = np.hstack((coilcurr, newcol))
            numecoils +=1 

    return coilcurr, numecoils, numfcoils
