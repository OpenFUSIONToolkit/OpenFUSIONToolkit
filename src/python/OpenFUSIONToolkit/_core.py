#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python interface for Open FUSION Toolkit common runtime functions

@authors Chris Hansen
@date Feb 2025
@ingroup doxy_oft_python
'''
import platform
import shutil
import tempfile
import numpy
from ._interface import *
from .util import run_shell_command


class OFT_env():
    '''! OpenFUSIONToolkit runtime environment class'''
    def __init__(self,debug_level=0,nthreads=2,unique_tempfiles='global',abort_callback=True):
        '''! Initialize OFT runtime object

        @param debug_level Level of debug printing (0-3)
        @param nthreads Number of threads for execution
        @param unique_tempfiles Method for temporary file creation ('global': Create unique folder
        in global temporary space, 'local_dir': Create unique folder in current working directory,
        'local_file': Use current working directory and append unique identifier to filenames,
        'none': Use non-unique names in local directory; can lead to conflict with multiple instances)
        @param abort_callback Use callback for "graceful" abort
        '''
        ## OS type
        self.os = platform.uname()[0]
        ## Number of physical CPUs (if multiple types are present only "performance" are counted)
        self.ncpus = None
        if self.os == 'Darwin':
            try:
                result, errcode = run_shell_command('sysctl -n hw.perflevel0.physicalcpu')
                if errcode == 0:
                    self.ncpus = int(result)
            except:
                pass
        elif self.os == 'Linux':
            try:
                result, errcode = run_shell_command('lscpu -p=CORE')
                if errcode == 0:
                    max_cpu = 0
                    for line in result.splitlines():
                        if line.startswith('#'):
                            continue
                        max_cpu = max(max_cpu,int(line))
                    self.ncpus = max_cpu+1
            except:
                pass
        if self.ncpus is None:
            print("WARNING: Could not detect number of physical cores.")
        else:
            if nthreads > self.ncpus:
                print('Warning: Request of {0} threads exceeds {1} physical cores detected by OFT (excluding "efficiency" cores)'.format(nthreads, self.ncpus))
                print("         If correct, this will significantly degrade performance.")
        ## ID of Python interpreter process
        self.pid = os.getpid()
        if unique_tempfiles == 'global':
            ## Directory for temporary files
            self.tempdir = os.path.join(tempfile.gettempdir(),'oft_{0}'.format(self.pid))
            try:
                os.mkdir(self.tempdir)
            except:
                print("Could not make temporary directory")
                raise
        elif unique_tempfiles == 'local_dir':
            self.tempdir = os.path.join(os.getcwd(),'oft_tmp-{0}'.format(self.pid))
            try:
                os.mkdir(self.tempdir)
            except:
                print("Could not make temporary directory")
                raise
        elif unique_tempfiles == 'local_file':
            self.tempdir = None
        elif unique_tempfiles == 'none':
            self.tempdir = None
            self.pid = None
            print("Warning: Using non-unique names/locations for temporary files can lead to conflicts if multiple python instances are used in the same directory")
        else:
            raise ValueError('Unknown value "{0}" for "unique_tempfiles"'.format(unique_tempfiles))
        ## Number of threads for execution
        self.nthreads = nthreads
        ## Debug level
        self.debug_level = debug_level
        ## Main input file groups
        self.oft_in_groups = {
            "runtime_options": {
                "debug": "{0}".format(debug_level),
            },
            "mesh_options": {
                "meshname": "'none'"
            }
        }
        ## Input filename for execution environment
        self.oft_ifile = self.unique_tmpfile('oftpyin')
        self.update_oft_in()
        # Initialize OFT
        slens = numpy.zeros((4,), dtype=numpy.int32)
        ifile_c = c_char_p(self.oft_ifile.encode())
        if abort_callback:
            oft_init(c_int(nthreads),ifile_c,slens,oft_python_abort)
        else:
            oft_init(c_int(nthreads),ifile_c,slens,c_void_p())
        ## General string size
        self.oft_slen = slens[1]
        ## Path string size
        self.oft_path_slen = slens[2]
        ## Error string size
        self.oft_error_slen = slens[3]
    
    def set_debug_level(self,debug_level):
        '''! Set debug verbosity level
        
        @param debug_level New value for debug level (must be in range [0,3])
        '''
        if (debug_level < 0) or (debug_level > 3):
            raise ValueError('Invalid value of "debug_level", must be in range [0,3]')
        oftpy_set_debug(c_int(debug_level))

    def set_num_threads(self,nthreads):
        '''! Set the number of OpenMP threads to use
        
        @param nthreads Number of threads to use for subsequent OpenMP parallel regions
        '''
        oftpy_set_nthreads(c_int(nthreads))
    
    def unique_tmpfile(self,filename):
        '''! Get unique temporary filename
        
        @param filename Base non-unique filename
        @result Unique filepath in suitable temporary location
        '''
        if self.tempdir is None:
            if self.pid is None:
                return '{0}'.format(filename)
            else:
                return '{0}-{1}'.format(filename,self.pid)
        else:
            return os.path.join(self.tempdir,filename)

    def path2c(self,path):
        '''! Convert general strings to C-compatible objects calls to OFT compiled API
        
        @param path Python path string
        @result `c_char_p` object containing path string value
        '''
        if len(path) > self.oft_path_slen:
            raise ValueError("Path length exceeds OFT library allowable lenght of {0}".format(self.oft_path_slen))
        return c_char_p(path.encode())

    def string2c(self,string):
        '''! Convert general strings to C-compatible objects calls to OFT compiled API
        
        @param string Python string
        @result `c_char_p` object containing string value
        '''
        if len(string) > self.oft_slen:
            raise ValueError("String length exceeds OFT library allowable lenght of {0}".format(self.oft_slen))
        return c_char_p(string.encode())

    def get_c_errorbuff(self):
        '''! Get properly-sized error string buffer for calls to OFT compiled API'''
        return create_string_buffer(b"",self.oft_error_slen)

    def update_oft_in(self):
        '''! Update input file with current settings (see @ref oft_in_groups)'''
        with open(self.oft_ifile, 'w+') as fid:
            for name, options in self.oft_in_groups.items():
                fid.write("&{0}\n".format(name))
                for option_name, option_value in options.items():
                    fid.write("  {0}={1}\n".format(option_name,option_value))
                fid.write("/\n\n")
    
    def __del__(self):
        '''! Destroy environment and cleanup known temporary files'''
        if self.tempdir is not None:
            try:
                shutil.rmtree(self.tempdir)
            except:
                print('Warning: unable to delete temporary directory "{0}"'.format(self.tempdir))
        else:
            try:
                os.remove(self.oft_ifile)
            except:
                print('Warning: unable to delete temporary file "{0}"'.format(self.oft_ifile))