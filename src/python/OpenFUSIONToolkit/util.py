#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Helper functions for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
import os
import subprocess
import numpy
import h5py
from collections import OrderedDict
from ._interface import oftpy_dump_cov

# Common parameters
## Vacuum magnetic permeability
mu0 = numpy.pi*4.E-7
## Electron charge
eC = 1.60217663e-19


def run_shell_command(command, timeout=10, env_vars={}):
    '''! Run a shell command

    @param command Command and arguments to run
    @param timeout Timeout for command to complete [seconds]
    @param env_vars Modifications to runtime environment
    @returns `result` STDOUT, `errcode` Error code from command
    '''
    # Run shell command
    my_env = os.environ.copy()
    for key, val in env_vars.items():
        my_env[key] = val
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=my_env)
    # Wait for process to complete or timeout
    try:
        outs, _ = pid.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        pid.kill()
        outs, _ = pid.communicate()
        print("WARNING: Command timeout")
    errcode = pid.poll()
    result = outs.decode("utf-8")
    return result, errcode


def write_native_mesh(filename, r, lc, reg, nodesets=[], sidesets=[], ho_info=None, periodic_info=None):
    r'''Create a native HDF5 mesh file for OFT from the given mesh information

    @param filename Filename for mesh file
    @param r Points list [np,3]
    @param lc Cell list [nc,3] (1-based)
    @param reg Region list [nc]
    @param nodesets List of node sets
    @param sidesets List of side sets
    @param ho_info High-order grid information
    @param periodic_info Information for mesh periodicity
    '''
    print()
    print("Saving mesh: {0}".format(filename))
    with h5py.File(filename, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
        if len(nodesets) > 0:
            h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(nodesets),], dtype='i4')
            for i, node_set in enumerate(nodesets):
                h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set, dtype='i4')
        if len(sidesets) > 0:
            h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(sidesets),], dtype='i4')
            for i, side_set in enumerate(sidesets):
                h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set, dtype='i4')
        if ho_info is not None:
            h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
            h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1], dtype='i4')
            if ho_info[2] is not None:
                h5_file.create_dataset('mesh/ho_info/LF', data=ho_info[2], dtype='i4')
        if periodic_info is not None:
            h5_file.create_dataset('mesh/periodicity/nodes', data=periodic_info, dtype='i4')

def read_fortran_namelist(file0, silent=True, b_arr=False):
	r'''Return a dictionary with the parameters in the namelist file (file0)
	
	@param file0 File to read from
	@param silent If false, print dictionary entries after reading
	@param b_arr Refers to an array at the bottom of the file, if one exists
	@result Dictionary containing data from file0.
	'''

	f_lines = None
	with open(str(file0), 'r') as f:
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
			
	return data_dict
