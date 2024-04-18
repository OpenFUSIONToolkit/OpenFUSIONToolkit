#!/usr/bin/python
from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt
import oft_io

dump = oft_io.oft_histfile('gem_flux.hist')
time = numpy.r_[dump.data['time']]*1.E3
flux = numpy.r_[dump.data['flux']]

fig, ax = plt.subplots(1,1,tight_layout=True)
ax.plot(time,flux[0]-flux)
ax.grid(True)
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Reconnected flux [Wb]')
plt.show()