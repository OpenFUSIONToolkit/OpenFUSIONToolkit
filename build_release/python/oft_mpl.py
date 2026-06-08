from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def close_handler(event):
    if event.key == 'q':
        sys.exit(0)
    if event.key == 'w':
        plt.close(event.canvas.figure)

# def pick_handler(event):
#     if isinstance(event.artist, Line2D):
#         thisline = event.artist
#         xdata = thisline.get_xdata()
#         ydata = thisline.get_ydata()
#         minlist = numpy.sqrt((xdata - event.mouseevent.xdata)**2. + (ydata - event.mouseevent.ydata)**2.)
#         mind = numpy.unravel_index(minlist.argmin(), xdata.shape)
#         print('pick: x = {0}, y = {1}'.format(xdata[mind], ydata[mind]))

class plot:
    def __init__(self, figure=1, subplots=None):
        if subplots is None:
            self.fig, ax_tmp = plt.subplots(num=figure)
            self.ax = [ax_tmp]
        else:
            self.fig, self.ax = plt.subplots(subplots, 1, num=figure, sharex=True)
        self.cid = self.fig.canvas.mpl_connect('key_press_event', close_handler)
        # self.pid = self.fig.canvas.mpl_connect('pick_event', pick_handler)

    def add_plot(self, xdata=None, ydata=None, type='line', xlabel=None, ylabel=None,
                 label=None, subplot=1, grid=True, include_xzero=True, include_yzero=False):
        # Create axes for current plot
        curr_plot = self.ax[subplot-1]
        curr_plot.grid(grid)
        # Plot desired type
        if (xdata is not None) and (ydata is not None):
            if type == 'line':
                curr_plot.plot(xdata, ydata, picker=4, label=label)
            elif type == 'semix':
                curr_plot.semilogx(xdata, ydata, picker=4, label=label)
            elif type == 'semiy':
                curr_plot.semilogy(xdata, ydata, picker=4, label=label)
        # Force X=0 and Y=0 on plot
        if include_xzero or include_yzero:
            curr_plot.autoscale()
            if include_xzero:
                curr_plot.set_xlim(left=0)
            if include_yzero:
                curr_plot.set_ylim(bottom=0)
        # Add X-Axis label
        if (xlabel is not None) and (len(self.ax) == subplot):
            curr_plot.set_xlabel(xlabel)
        # Add Y-Axis label
        if ylabel is not None:
            curr_plot.set_ylabel(ylabel)

    def show_legend(self, subplot=-1):
        if subplot > 0:
            curr_plot = self.ax[subplot-1]
            curr_plot.legend()
        else:
            for curr_plot in self.ax:
                curr_plot.legend()
