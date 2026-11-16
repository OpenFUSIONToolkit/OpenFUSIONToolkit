#!/usr/bin/python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Plot flux through the hole of the torus and magnetic energy from the ThickCurr torus example
import matplotlib.pyplot as plt
from OpenFUSIONToolkit.io import histfile

flux_hist = histfile('torus_flux.hist')
energy_hist = histfile('thickcurr.hist')

fig, ax = plt.subplots(2,1,sharex=True,tight_layout=True)
ax[0].semilogy(flux_hist['time']*1.E3,flux_hist['flux'])
ax[0].set_ylabel('Flux through hole [Wb]')
ax[1].semilogy(energy_hist['time']*1.E3,energy_hist['energy'])
ax[1].set_ylabel('Magnetic energy [J]')
ax[1].set_xlabel('Time [ms]')
for ax_tmp in ax:
    ax_tmp.grid(True)
plt.show()
