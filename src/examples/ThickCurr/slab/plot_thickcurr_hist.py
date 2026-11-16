#!/usr/bin/python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Plot magnetic energy from the ThickCurr slab example, compared with the decay of the slowest mode
import numpy
import matplotlib.pyplot as plt
from OpenFUSIONToolkit.io import histfile

tau = 42.05E-3  # Analytic decay time of slowest mode [s]

hist = histfile('thickcurr.hist')
time = hist['time']
energy = hist['energy']

fig, ax = plt.subplots(1,1,tight_layout=True)
ax.semilogy(time*1.E3,energy,label='ThickCurr')
ax.semilogy(time*1.E3,energy[-1]*numpy.exp(-2.0*(time-time[-1])/tau),'--',label=r'$e^{-2t/\tau}$')
ax.grid(True)
ax.legend()
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Magnetic energy [J]')
plt.show()
