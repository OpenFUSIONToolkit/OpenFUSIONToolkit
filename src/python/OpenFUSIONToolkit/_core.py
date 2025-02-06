import tempfile
import numpy
from ._interface import *


class OFT_env():
    '''! TokaMaker G-S solver class'''
    def __init__(self,debug_level=0,nthreads=2,use_tmpdir=True):
        '''! Initialize TokaMaker object

        @param debug_level Level of debug printing (0-3)
        @param nthreads Number of threads for execution
        '''
        ## ID of Python interpreter process
        self.pid = os.getpid()
        if use_tmpdir:
            ## Directory for temporary files
            self.tempdir = os.path.join(tempfile.gettempdir(),'oft_{0}'.format(self.pid))
            try:
                os.mkdir(self.tempdir)
            except:
                print("Could not make temporary directory")
                raise
        else:
            self.tempdir = None
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
        oft_init(c_int(nthreads),ifile_c,slens)
        ## General string size
        self.oft_slen = slens[1]
        ## Path string size
        self.oft_path_slen = slens[2]
        ## Error string size
        self.oft_error_slen = slens[3]
    
    def unique_tmpfile(self,filename):
        if self.tempdir is not None:
            return os.path.join(self.tempdir,filename)
        else:
            return '{0}-{1}'.format(filename,self.pid)

    def path2c(self,path):
        if len(path) > self.oft_path_slen:
            raise ValueError("Path length exceeds OFT library allowable lenght of {0}".format(self.oft_path_slen))
        return c_char_p(path.encode())

    def string2c(self,string):
        if len(string) > self.oft_slen:
            raise ValueError("String length exceeds OFT library allowable lenght of {0}".format(self.oft_slen))
        return c_char_p(string.encode())

    def get_c_errorbuff(self):
        return create_string_buffer(b"",self.oft_error_slen)

    def update_oft_in(self):
        '''! Update input file (`oftpyin`) with current settings'''
        with open(self.oft_ifile, 'w+') as fid:
            for name, options in self.oft_in_groups.items():
                fid.write("&{0}\n".format(name))
                for option_name, option_value in options.items():
                    fid.write("  {0}={1}\n".format(option_name,option_value))
                fid.write("/\n\n")
    
    def __del__(self):
        try:
            os.remove(self.oft_ifile)
        except:
            print('Warning: unable to delete temporary file "{0}"'.format(self.oft_ifile))