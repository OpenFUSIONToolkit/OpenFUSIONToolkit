#!/usr/bin/python
from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt
from OpenFUSIONToolkit.io import histfile
import oft_mpl

mu0 = numpy.pi*(4.E-7)


def get_gr(time, total_energy, hist_file):
    # Get growth rate for linear run
    scal_fac = None
    is_linear = False
    for line in hist_file.header:
        if (line.count('linear') > 0) and (line.count('non-linear') == 0):
            is_linear = True
            print("  Linear run detected")
        if line.find('E0 =') > 0:
            scal_fac = float(line.split('=')[1])
            print("  Scale factor = {0}".format(scal_fac))
    if is_linear:
        dt = time[1] - time[0]
        if scal_fac is None:
            return time[1:], .5*numpy.log(total_energy[1:]/total_energy[:-1])/dt
        else:
            return time[1:], .5*numpy.log(total_energy[1:]/scal_fac)/dt
    else:
        return None


# Start script
nfiles = len(sys.argv) - 1
if nfiles > 0:
    files = sys.argv[1:]
else:
    print('No history files specified')
    sys.exit(1)
#
plot_id = 0
num_plots = []
for i in range(5):
    plot_id += 1
    num_plots.append(oft_mpl.plot(plot_id))
energy_plots = []
for i in range(3):
    plot_id += 1
    energy_plots.append(oft_mpl.plot(plot_id, subplots=2))
Tflux_plot = None
Tcurr_plot = None
Ti_plot = None
Ne_plot = None
Te_plot = None
derr_plot = None
gr_plot = None
# Loop over history files
for (i, file) in enumerate(files):
    # Read fields from MUG history file
    print('Reading file: "{0}"'.format(file))
    dump = histfile(file)
    ts = numpy.r_[dump['ts']]
    time = numpy.r_[dump['time']]*1.E3
    dt = time[1:-1] - time[0:-2]
    lits = numpy.r_[dump['lits']]
    nlits = numpy.r_[dump['nlits']]
    stime = numpy.r_[dump['stime']]
    mag_energy = numpy.r_[dump['men']]/mu0
    kin_energy = numpy.r_[dump['ven']]
    total_energy = mag_energy + kin_energy
    tflux = numpy.r_[dump['tflux']]
    tcurr = numpy.r_[dump['tcurr']]/mu0
    div_err = numpy.r_[dump['derr']]/(2.*mag_energy*mu0)
    gr_data = get_gr(time/1.E3, total_energy, dump)
    print()
    # Detect legacy field names
    if 'npart' in dump:
        ti_name = 'temp'
        ne_name = 'npart'
        te_name = 'tempe'
    else:
        ti_name = 'ti'
        ne_name = 'ne'
        te_name = 'te'
    # Plot solution time
    plot = num_plots[0]
    plot.add_plot(numpy.cumsum(stime), time, xlabel='Walltime [hr]', ylabel='Simulation Time [ms]')
    # Plot time step size
    plot = num_plots[1]
    plot.add_plot(time[:-2], dt*1.E3, xlabel='Time [ms]', ylabel='dt [us]', include_yzero=True)
    # Plot linear iteration count
    plot = num_plots[2]
    plot.add_plot(ts, lits, xlabel='Step index', ylabel='# of Linear iterations', include_yzero=True)
    # Plot non-linear iteration count
    plot = num_plots[3]
    plot.add_plot(ts, nlits, xlabel='Step index', ylabel='# of Non-Linear iterations', include_yzero=True)
    # Plot required solver time
    plot = num_plots[4]
    plot.add_plot(ts, stime, xlabel='Step index', ylabel='Solver time [s]', include_yzero=True)
    # Plot magnetic energy
    plot = energy_plots[0]
    plot.add_plot(time, mag_energy, ylabel='Magnetic energy [J]', subplot=1)
    plot.add_plot(time, mag_energy, xlabel='Time [ms]', ylabel='Magnetic energy [J]', subplot=2, type='semiy')
    # Plot kinetic energy
    plot = energy_plots[1]
    plot.add_plot(time, kin_energy, ylabel='Kinetic energy [J]', subplot=1)
    plot.add_plot(time, kin_energy, xlabel='Time [ms]', ylabel='Kinetic energy [J]', subplot=2, type='semiy')
    # Plot total energy
    plot = energy_plots[2]
    plot.add_plot(time, total_energy, ylabel='Total energy [J]', subplot=1)
    plot.add_plot(time, total_energy, xlabel='Time [ms]', ylabel='Total energy [J]', subplot=2, type='semiy')
    if abs(mag_energy).max() > 1.E-7:
        if Tflux_plot is None:
            plot_id += 1
            Tflux_plot = oft_mpl.plot(plot_id, subplots=2)
            plot_id += 1
            Tcurr_plot = oft_mpl.plot(plot_id, subplots=2)
        # Plot toroidal flux
        Tflux_plot.add_plot(time, tflux/1.E-3, ylabel='Toroidal Flux [mWb]', subplot=1)
        Tflux_plot.add_plot(time, abs(tflux)/1.E-3, xlabel='Time [ms]', ylabel='Toroidal Flux [mWb]', subplot=2, type='semiy')
        # Plot toroidal current
        Tcurr_plot.add_plot(time, tcurr/1.E3, ylabel='Toroidal Current [kA]', subplot=1)
        Tcurr_plot.add_plot(time, abs(tcurr)/1.E3, xlabel='Time [ms]', ylabel='Toroidal Current [kA]', subplot=2, type='semiy')
    if ti_name in dump:
        if Ti_plot is None:
            plot_id += 1
            Ti_plot = oft_mpl.plot(plot_id)
            plot_id += 1
            Ne_plot = oft_mpl.plot(plot_id)
        # Plot ion temperature
        Ti = numpy.r_[dump[ti_name]]
        Ti_plot.add_plot(time, Ti, xlabel='Time [ms]', ylabel='Average ion temperature [eV]', include_yzero=True)
        # Plot ion density
        Ne = numpy.r_[dump[ne_name]]
        Ne_plot.add_plot(time, Ne, xlabel='Time [ms]', ylabel='Average density [m^-3]', include_yzero=True)
    if te_name in dump:
        Te = numpy.r_[dump[te_name]]
        if Te.max() > 0.:
            if Te_plot is None:
                plot_id += 1
                Te_plot = oft_mpl.plot(plot_id)
            # Plot electron temperature
            Te_plot.add_plot(time, Te, xlabel='Time [ms]', ylabel='Average electron temperature [eV]', include_yzero=True)
    if div_err.max() > 0:
        if derr_plot is None:
            plot_id += 1
            derr_plot = oft_mpl.plot(plot_id)
        # Plot electron temperature
        derr_plot.add_plot(time, div_err, xlabel='Time [ms]', ylabel='Divergence error', type='semiy')
    if gr_data is not None:
        if gr_plot is None:
            plot_id += 1
            gr_plot = oft_mpl.plot(plot_id)
        # Plot growth rate
        gr_plot.add_plot(gr_data[0], gr_data[1], xlabel='Time [ms]', ylabel='Growth rate [1/s]', include_yzero=True)
# Display plots
plt.show()
