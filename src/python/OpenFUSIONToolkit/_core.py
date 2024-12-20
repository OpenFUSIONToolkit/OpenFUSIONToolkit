import numpy
from ._interface import *


class OFT_env():
    '''! TokaMaker G-S solver class'''
    def __init__(self,debug_level=0,nthreads=2):
        '''! Initialize TokaMaker object

        @param debug_level Level of debug printing (0-3)
        @param nthreads Number of threads for execution
        '''
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
        self.update_oft_in()
        # Initialize OFT
        slens = numpy.zeros((4,), dtype=numpy.int32)
        oft_init(c_int(nthreads),slens)
        ## General string size
        self.oft_slen = slens[1]
        ## Path string size
        self.oft_path_slen = slens[2]
        ## Error string size
        self.oft_error_slen = slens[3]
    
    def update_oft_in(self,filename='oftpyin'):
        '''! Update input file (`oftpyin`) with current settings'''
        with open(filename, 'w+') as fid:
            for name, options in self.oft_in_groups.items():
                fid.write("&{0}\n".format(name))
                for option_name, option_value in options.items():
                    fid.write("  {0}={1}\n".format(option_name,option_value))
                fid.write("/\n\n")