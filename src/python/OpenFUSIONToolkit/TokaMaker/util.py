#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! General utility and supporting functions for TokaMaker

@authors Chris Hansen
@date April 2024
@ingroup doxy_oft_python
'''
import struct
import numpy
from collections import OrderedDict
from .._interface import *
from ..util import read_fortran_namelist

## @cond
tokamaker_eval_green = ctypes_subroutine(oftpy_lib.tokamaker_eval_green,
    [c_int, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), c_double, c_double, ctypes_numpy_array(float64,1)])
## @endcond


def create_isoflux(npts, r0, z0, a, kappa, delta, kappaL=None, deltaL=None):
    r'''! Create isoflux points using simple analytic form

    @param npts Number of points to sample (evenly spaced in \f$\theta\f$)
    @param r0 Major radial position for magnetic axis
    @param z0 Vertical position for magnetic axis
    @param a Minor radius
    @param kappa Elongation (upper only if kappaL is set)
    @param delta Triangularity (upper only if deltaL is set)
    @param kappaL Lower elongation (default: kappa)
    @param deltaL Lower triangularity (default: delta)
    @result Point list [npts,2]
    '''
    kappaU = kappa
    deltaU = delta
    if kappaL is None:
        kappaL = kappaU
    if deltaL is None:
        deltaL = deltaU
    x0 = numpy.r_[r0, z0]
    isoflux_pts = numpy.zeros((npts,2))
    for i in range(npts):
        theta = i*2.0*numpy.pi/npts
        delta = ((deltaU + deltaL) + (deltaU - deltaL)*numpy.sin(theta))/2
        kappa = ((kappaU + kappaL) + (kappaU - kappaL)*numpy.sin(theta))/2
        isoflux_pts[i,:] = x0 + a*numpy.r_[numpy.cos(theta+numpy.arcsin(delta)*numpy.sin(theta)),kappa*numpy.sin(theta)]
    return isoflux_pts


def create_spline_flux_fun(npts,x,y,axis_bc=[1,0.0],edge_bc=[1,0.0],normalize=True):
    r'''! Build cubic spline flux function

    @param npts Number of points for definition
    @param x Location of spline "knots" in normalized flux
    @param y Value of flux function at spline "knots"
    @param axis_bc [SciPy BC specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html) on axis (\f$ \hat{\psi} = 0 \f$)
    @param edge_bc [SciPy BC specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html) on LCFS (\f$ \hat{\psi} = 1 \f$)
    @result Flux function definition dictionary
    '''
    try:
        from scipy.interpolate import CubicSpline
    except ImportError:
        print("Spline flux function requires SciPy")
        raise
    prof = CubicSpline(x,y,bc_type=[axis_bc,edge_bc])
    x_sample = numpy.linspace(0.0,1.0,npts)
    prof = {
        'type': 'linterp',
        'x': x_sample,
        'y': prof(x_sample)
    }
    if normalize:
        prof['y'] /= prof['y'][0]
    return prof



def create_power_flux_fun(npts,alpha,gamma):
    r'''! Build power law flux function of the form \f$ (1-\hat{\psi}^{\alpha})^{\gamma} \f$

    @param npts Number of points for definition
    @param alpha Inner exponent
    @param gamma Outer exponent
    @result Flux function definition dictionary
    '''
    psi_pts = numpy.linspace(0.0,1.0,npts)
    return {
        'type': 'linterp',
        'x': psi_pts,
        'y': numpy.power(1.0-numpy.power(psi_pts,alpha),gamma)
    }


def read_eqdsk(filename):
    '''! Read gEQDSK file

    @param filename Path to gEQDSK file
    @result Dictionary containing gEQDSK information
    '''
    def read_1d(fid, n):
        j = 0
        output = numpy.zeros((n,))
        for i in range(n):
            if j == 0:
                line = fid.readline()
            output[i] = line[j:j+16]
            j += 16
            if j == 16*5:
                j = 0
        return output

    def read_2d(fid, n, m):
        j = 0
        output = numpy.zeros((n, m))
        for k in range(n):
            for i in range(m):
                if j == 0:
                    line = fid.readline()
                output[k, i] = line[j:j+16]
                j += 16
                if j == 16*5:
                    j = 0
        return output
    # Read-in data
    eqdsk_obj = {}
    with open(filename, 'r') as fid:
        # Get sizes
        line = fid.readline()
        eqdsk_obj['case'] = line[:48]
        split_line = line[48:].split()
        eqdsk_obj['nr'] = int(split_line[-2])
        eqdsk_obj['nz'] = int(split_line[-1])
        # Read header content
        line_keys = [['rdim',  'zdim',  'rcentr',  'rleft',  'zmid'],
                     ['raxis', 'zaxis', 'psimag', 'psibry', 'bcentr'],
                     ['ip',    'skip',  'skip',   'skip',   'skip'],
                     ['skip',  'skip',  'skip',   'skip',   'skip']]
        for i in range(4):
            line = fid.readline()
            for j in range(5):
                if line_keys[i][j] == 'skip':
                    continue
                line_seg = line[j*16:(j+1)*16]
                eqdsk_obj[line_keys[i][j]] = float(line_seg)
        # Read flux profiles
        keys = ['fpol', 'pres', 'ffprim', 'pprime']
        for key in keys:
            eqdsk_obj[key] = read_1d(fid, eqdsk_obj['nr'])
        # Read PSI grid
        eqdsk_obj['psirz'] = read_2d(fid, eqdsk_obj['nz'],
                                        eqdsk_obj['nr'])
        # Read q-profile
        eqdsk_obj['qpsi'] = read_1d(fid, eqdsk_obj['nr'])
        # Read limiter count
        line = fid.readline()
        eqdsk_obj['nbbs'] = int(line.split()[0])
        eqdsk_obj['nlim'] = int(line.split()[1])
        # Read outer flux surface
        eqdsk_obj['rzout'] = read_2d(fid, eqdsk_obj['nbbs'], 2)
        # Read limiting corners
        eqdsk_obj['rzlim'] = read_2d(fid, eqdsk_obj['nlim'], 2)
    return eqdsk_obj


def read_ifile(filename):
    '''! Read i-file inverse equilibrium file

    @param filename Path to file
    @result Dictionary containing i-file information
    '''
    def read_array(content,offset,var_type,count):
        if var_type in ("i", "f"):
            var_size = 4
        elif var_type in ("l", "d"):
            var_size = 8
        else:
            raise ValueError("Invalid variable type")
        array_size = var_size*count
        #
        head_val = struct.unpack_from("i",content,offset=offset)
        if head_val[0] != array_size:
            raise ValueError("Dataframe size does not match array size")
        offset += 4
        body_val = struct.unpack_from("="+var_type*count,content,offset=offset)
        offset += array_size
        tail_val = struct.unpack_from("i",content,offset=offset)
        if head_val[0] != tail_val[0]:
            raise ValueError("Head and tail values disagree")
        offset += 4
        return numpy.array(body_val), offset
    #
    with open(filename, 'rb') as handle:
        content = handle.read()
    out_dict = {}
    offset = 0
    sizes, offset = read_array(content,offset,"i",2)
    out_dict["npsi"] = sizes[0]
    out_dict["ntheta"] = sizes[1]
    var_type = "d"
    try:
        out_dict["psi"], offset = read_array(content,offset,var_type,sizes[0])
    except ValueError:
        try:
            var_type = "f"
            out_dict["psi"], offset = read_array(content,offset,var_type,sizes[0])
        except ValueError:
            raise ValueError("Unable to determine float point datatype")
    out_dict["f"], offset = read_array(content,offset,var_type,sizes[0])
    out_dict["p"], offset = read_array(content,offset,var_type,sizes[0])
    out_dict["q"], offset = read_array(content,offset,var_type,sizes[0])
    R, offset = read_array(content,offset,var_type,sizes[0]*sizes[1])
    Z, offset = read_array(content,offset,var_type,sizes[0]*sizes[1])
    out_dict["R"] = R.reshape(sizes)
    out_dict["Z"] = Z.reshape(sizes)
    return out_dict


def eval_green(x,xc):
        r'''! Evaluate Green's function for a toroidal filament

        @param x Observation point [2]
        @param xc Coil location [:,2]
        @result \f$\psi(x)\f$ due to a coil with unit current [A] at xc
        '''
        n = x.shape[0]
        vals = numpy.zeros((n,),dtype=numpy.float64)
        r = x[:,0].copy()
        z = x[:,1].copy()
        tokamaker_eval_green(c_int(n),r,z,
            c_double(xc[0]),c_double(xc[1]),vals)
        return vals*mu0


def compute_forces_components(tMaker_obj,psi,cell_centered=False):
    r'''! Compute terms needed for evaluating forces in passively conducting regions

    @param tMaker_obj TokaMaker equilibrium object
    @param psi \f$ \psi \f$ corresponding to desired currents
    @param cell_centered Evaluate at cell centers instead of node points?
    @result J_cond, B_cond, mask, R
    '''
    # Get conductor currents at cell centers
    mask, J_cond = tMaker_obj.get_conductor_currents(psi,cell_centered=cell_centered)
    
    # Find points inside conducting regions
    pt_mask = numpy.zeros((tMaker_obj.r.shape[0],), dtype=numpy.int32)
    pt_mask[tMaker_obj.lc[mask,:]] = 1
    
    # Set psi and evaluate B-field in conducting regions
    psi_save = tMaker_obj.get_psi(normalized=False)
    tMaker_obj.set_psi(psi)
    field_eval = tMaker_obj.get_field_eval('B')
    B_cond = numpy.zeros((tMaker_obj.r.shape[0],3))
    for i in range(tMaker_obj.r.shape[0]):
        if pt_mask[i] == 0:
            continue
        B_cond[i,:] = field_eval.eval(tMaker_obj.r[i,:2])
    tMaker_obj.set_psi(psi_save) # Reset psi

    if cell_centered:
        # Convert to cell centered
        Bv_cond = (B_cond[tMaker_obj.lc[:,0],:] + B_cond[tMaker_obj.lc[:,1],:] + B_cond[tMaker_obj.lc[:,2],:])/3.0
        rcc = (tMaker_obj.r[tMaker_obj.lc[:,0],:] + tMaker_obj.r[tMaker_obj.lc[:,1],:] + tMaker_obj.r[tMaker_obj.lc[:,2],:])/3.0
        return J_cond, Bv_cond, mask, rcc
    else:
        return J_cond, B_cond, mask, tMaker_obj.r
def namelist_read(file0, silent=True, b_arr=False):
	#Return a dictionary with the parameters in the namelist file (file0)
	#b_arr refers to an array at the bottom of the file, if one exists

	f = open(str(file0), 'r')
	f_lines = f.readlines()

	datalines = [] #Initialize new dictionary of information

	if b_arr == True:
		brk_idx = [] #Find the last break and turn it into the bottom array
		equal_idx = [] #Find the last equals sign

		for i in range(len(f_lines)):
			brk_loc = f_lines[i].find('/')  #Find the lines with backslashes (line breaks)
			if brk_loc >= 0.0:
				brk_idx.append(i)

			equal_loc = f_lines[i].find('=')
			if equal_loc >= 0.0:
				equal_idx.append(i)

		f_lines[brk_idx[np.argmin(abs(brk_idx-np.max(equal_idx)))]] = 'b_arr = '

	#Prune everything after the comments and remove line breaks 
	end_idx = len(f_lines)

	for i in range(len(f_lines)):
		line = f_lines[i]

		if line.find('comment:') >= 0: #Find the comment at the end (if it exists and make this the end of the file)
			end_idx = i

		line_update = line.replace('\n','')  #Remove line breaks
		line_update = line_update.replace(',',' ')   #Remove commas
		line_update = line_update.replace('"', '')  #Remove comment/quotation marks
		comment_loc = line_update.find('!')  #Find the comment
		brk_loc = line_update.find('/')  #Find the lines with backslashes (line breaks)
		headr_loc = line_update.find('&') #and list headers (& symbols)

		if i < end_idx: #Only include lines before the end comment
			if headr_loc < 0: #Remove hardcoded linebreaks and list headers 
				if brk_loc < 0:
					if comment_loc < 0:  #If there are no comments in the line (-1)
						if line_update.strip() != '':
							datalines.append(line_update.strip())
					elif comment_loc >= 0:
						#ignore everything after the comment and get rid of white space (strip)
						line_update = line_update[:comment_loc].strip()
						if line_update != '':
							datalines.append(line_update)


	data_dict = OrderedDict()
	block_idx = []

	for i in range(len(datalines)):
		equal_loc = datalines[i].find('=')  #Find where the equal sign is in the array
		if equal_loc > 0:
			block_idx.append(i)
			key0 = datalines[i][:equal_loc].strip()
			data_dict[key0] = 0.0

	key_list = list(data_dict.keys())


	for i in range(len(key_list)):
		if key_list[i] != 'b_arr':
			#Concetenate each of the blocks (arrays that span multiple lines) of data
			if i < len(key_list)-1:
				mult_block = ' '.join(datalines[block_idx[i]:block_idx[i+1]])
				equal_loc = mult_block.find('=')
				mult_block = mult_block[equal_loc+1:].strip()
				#data_dict[key_list[i]] = np.array(mult_block.split(' '))#.astype('float')
				dupl_loc = mult_block.find('*')
				if dupl_loc >= 0:
					mult_block = mult_block.split()
					star_idx = [idx for idx, s in enumerate(mult_block) if '*' in s]

					for j in range(len(star_idx)):
						dupl_str = mult_block[star_idx[j]].split('*')
						try:
							dupl_arr = ' '.join(int(dupl_str[0])*[dupl_str[1]])
						except:
							continue
						mult_block[star_idx[j]] = dupl_arr

					mult_block = ' '.join(mult_block)

					try:
						data_dict[key_list[i]] = np.array(mult_block.split()).astype('float')

					except:
						data_dict[key_list[i]] = mult_block
				else:
					try:
						if len(mult_block.split()) > 1:  #If longer than one entry, make it an array
							data_dict[key_list[i]] = np.array(mult_block.split()).astype('float')
						else:    #if it is one entry, make it a float or an integer
							data_dict[key_list[i]] = eval(mult_block.split()[0])

					except:
						data_dict[key_list[i]] = mult_block

			#Try the same process for the last entry		
			else:
				mult_block = ' '.join(datalines[block_idx[i]:len(datalines)])  #On the last block, go to the end of the file
				equal_loc = mult_block.find('=')
				mult_block = mult_block[equal_loc+1:].strip()
				dupl_loc = mult_block.find('*')
				if dupl_loc >= 0:
					mult_block = mult_block.split()
					star_idx = [idx for idx, s in enumerate(mult_block) if '*' in s]

					for j in range(len(star_idx)):
						dupl_str = mult_block[star_idx[j]].split('*')
						try:
							dupl_arr = ' '.join(int(dupl_str[0])*[dupl_str[1]])
						except:
							continue
						mult_block[star_idx[j]] = dupl_arr

					mult_block = ' '.join(mult_block) 
					try:
						data_dict[key_list[i]] = np.array(mult_block.split()).astype('float')
					except:
						data_dict[key_list[i]] = mult_block
				else:
					try:
						data_dict[key_list[i]] = np.array(mult_block.split()).astype('float')
					except:
						data_dict[key_list[i]] = mult_block

		elif key_list[i] == 'b_arr':
			b_arr_len = len(datalines[block_idx[i]+1:len(datalines)])
			data_dict['b_arr'] = np.zeros((b_arr_len,6))

			for j in range(b_arr_len):
				try:
					b_line = np.array(re.split(r'(\s+)', datalines[block_idx[i]+1+j]))
					b_line_space_idx = []
					b_line_rmv_idx = []

					for k in range(len(b_line)):
						if b_line[k].isspace():
							b_line_space_idx.append(k)
					
					b_line_space_len = (lambda x:[len(l) for l in x])(b_line[b_line_space_idx])

					for l in range(len(b_line_space_idx)):
						if b_line_space_len[l] < 10: #== min(b_line_space_len):
							b_line_rmv_idx.append(b_line_space_idx[l])
						else:
							b_line[b_line_space_idx[l]] = 0.0

					b_line = np.delete(b_line, b_line_rmv_idx)

					#b_line = np.array(datalines[block_idx[i]+1+j].split()).astype('float')
					data_dict['b_arr'][j][0:len(b_line)] = b_line 


				except:
					data_dict['b_arr'][j][:] = np.nan

			if 'RF' not in key_list:
				try:
					data_dict['RF'] = data_dict['b_arr'][0:len(data_dict['TURNFC']),0]
					data_dict['ZF'] = data_dict['b_arr'][0:len(data_dict['TURNFC']),1]
					data_dict['WF'] = data_dict['b_arr'][0:len(data_dict['TURNFC']),2]
					data_dict['HF'] = data_dict['b_arr'][0:len(data_dict['TURNFC']),3]
					data_dict['AF'] = data_dict['b_arr'][0:len(data_dict['TURNFC']),4]
					data_dict['AF2'] = data_dict['b_arr'][0:len(data_dict['TURNFC']),5]
					Fcoil_in_barr = True

				except:
					print('Error in determining F-coil data')

			if 'RVS' not in key_list:
				try:
					VS_idx_start = len(data_dict['b_arr'])-len(data_dict['RSISVS']) #len(data_dict['TURNFC'])
					VS_idx_end = len(data_dict['b_arr'])#len(data_dict['RSISVS'])+VS_idx_start

					data_dict['RVS'] = data_dict['b_arr'][VS_idx_start:VS_idx_end,0]
					data_dict['ZVS'] = data_dict['b_arr'][VS_idx_start:VS_idx_end,1]
					data_dict['WVS'] = data_dict['b_arr'][VS_idx_start:VS_idx_end,2]
					data_dict['HVS'] = data_dict['b_arr'][VS_idx_start:VS_idx_end,3]
					data_dict['AVS'] = data_dict['b_arr'][VS_idx_start:VS_idx_end,4]
					data_dict['AVS2'] = data_dict['b_arr'][VS_idx_start:VS_idx_end,5]

					RVS_in_barr = True

				except:
					RVS_in_barr = False
					print('Error in determining vacuum vessel data')

			if 'RE' not in key_list:
				if data_dict['IECOIL'] == 1:
					try:
						if 'RF' in data_dict.keys():
							EC_idx_start = len(data_dict['TURNFC'])
						else:
							EC_idx_start = 0	

						if RVS_in_barr == True:
							EC_idx_end = VS_idx_start
						else:
							EC_idx_end = len(data_dict['b_arr'])

						data_dict['RE'] = data_dict['b_arr'][EC_idx_start:EC_idx_end,0]
						data_dict['ZE'] = data_dict['b_arr'][EC_idx_start:EC_idx_end,1]
						data_dict['WE'] = data_dict['b_arr'][EC_idx_start:EC_idx_end,2]
						data_dict['HE'] = data_dict['b_arr'][EC_idx_start:EC_idx_end,3]
						data_dict['ECID'] = data_dict['b_arr'][EC_idx_start:EC_idx_end,4]

					except:
						print('Error in determining E-coil data')


	if silent == False:
		for key in data_dict:
			print(key, data_dict[key])
			
	f.close()

	return data_dict

def read_mhdin(path, e_coil_names, f_coil_names):
	raw = namelist_read(path)
	machine_dict = OrderedDict()
	
	# Expand later
	keys = ['MPNAM2', 'LPNAME']
	for key in keys:
		names = raw[key].replace("'", " ")
		names = names.split()
		machine_dict[key] = names
		
	e_coil_dict = OrderedDict()
	for name in e_coil_names:
		e_coil_dict[name] = []

	raw['ECID'] = [x for x in raw['ECID'] if len(x.strip()) > 0]
	e_coil_vars = ['RE', 'ZE', 'WE', 'HE']
	for var in e_coil_vars:
		raw[var] = raw[var].split()

	for i in range(len(raw['ECID'])):
		if raw['ECID'][i].strip() == '':
			continue
		idx = int(raw['ECID'][i]) - 1
		e_coil_dict[e_coil_names[idx]].append([float(raw['RE'][i]), float(raw['ZE'][i]), float(raw['WE'][i]), float(raw['HE'][i])])
	machine_dict['ECOIL'] = e_coil_dict


	f_coil_vars = ['RF', 'ZF', 'WF', 'HF', 'TURNFC']
	for var in f_coil_vars:
		raw[var] = raw[var].split()

	f_coil_dict = OrderedDict()
	for i in range(len(f_coil_names)):
		f_coil_dict[f_coil_names[i]] = [float(raw['RF'][i]), float(raw['ZF'][i]), float(raw['WF'][i]), float(raw['HF'][i]), float(raw['TURNFC'][i])]
	machine_dict['FCOIL'] = f_coil_dict

	probe_angle_dict = OrderedDict()
	i = 0
	probe_angles = raw['AMP2'].split()
	for probe_name in machine_dict['MPNAM2']:
		probe_angle_dict[probe_name] = float(probe_angles[i])
		i = i + 1
	machine_dict['AMP2'] = probe_angle_dict

	return machine_dict

def read_kfile(path, e_coil_names, f_coil_names, machine_dict):	
	raw = namelist_read(path)
	
	probe_names = machine_dict['MPNAM2']
	probe_vals = raw['EXPMP2'].split()
	probe_errs = raw['BITMPI'].split()
	probe_selected = raw['FWTMP2'].split()
	probes_dict = OrderedDict()

	for i in range(len(probe_names)):
		if not float(probe_selected[i]):
			continue
		probes_dict[probe_names[i]] = [float(probe_vals[i]), float(probe_errs[i])]

	loop_names = machine_dict['LPNAME']
	loop_vals = raw['COILS'].split()
	loop_errs = raw['PSIBIT'].split()
	loop_selected = raw['FWTSI'].split()

	loops_dict = OrderedDict()
	for i in range(len(loop_names)):
		if not float(loop_selected[i]):
			continue
		loops_dict[loop_names[i]] = [float(loop_vals[i]) + float(raw['SIREF']), float(loop_errs[i])]
		
	e_coil_vals = raw['ECURRT'].split()
	e_coil_errs = raw['BITEC'].split()
	e_coil_dict = OrderedDict()
	for i in range(len(e_coil_names)):
		e_coil_dict[e_coil_names[i]] = [float(e_coil_vals[i]), float(e_coil_errs[i])]
	
	f_coil_vals = raw['BRSP'].split()
	f_coil_errs = raw['BITFC'].split()
	f_coil_dict = OrderedDict()

	for i in range(len(f_coil_names)):
		f_coil_dict[f_coil_names[i]] = [float(f_coil_vals[i]), float(f_coil_errs[i])]
	return probes_dict, loops_dict, e_coil_dict, f_coil_dict
