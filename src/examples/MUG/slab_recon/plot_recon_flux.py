#!/usr/bin/python
from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt
from OpenFUSIONToolkit.io import histfile

dump = histfile('gem_flux.hist')
time = dump['time']*1.E3
flux = dump['flux']

fig, ax = plt.subplots(1,1,tight_layout=True)
ax.plot(time,flux[0]-flux)
ax.grid(True)
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Reconnected flux [Wb]')
plt.show()
