import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
import os
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

class Figure:
    def __init__(self, xlabel=r"x", ylabel=r"y"):
        self.xmin = 0
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.nb_curves = 0
        self.xlabel = xlabel
        self.ylabel = ylabel
        ## Change RC default parameters >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        mpl.rcdefaults()
        mpl.rcParams['font.size'] = 24. # change the size of the font in every figure
        mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
        mpl.rcParams['axes.labelsize'] = 24.
        mpl.rcParams['xtick.labelsize'] = 24
        mpl.rcParams['ytick.labelsize'] = 24
        mpl.rcParams['xtick.direction'] = "in"
        mpl.rcParams['ytick.direction'] = "in"
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True
        mpl.rcParams['xtick.major.width'] = 0.6
        mpl.rcParams['ytick.major.width'] = 0.6
        mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
        mpl.rcParams['pdf.fonttype'] = 3
        ## Create figure and axes
        self.fig, self.axes = plt.subplots(1, 1, figsize = (9.2, 5.6)) # (1,1) means one plot, and figsize is w x h in inch of figure
        self.fig.subplots_adjust(left = 0.18, right = 0.82, bottom = 0.18, top = 0.95) # adjust the box of axes regarding the figure size
        ## set the distance between the ticks and their labels
        self.axes.tick_params(axis='x', which='major', pad=7, length=8)
        self.axes.tick_params(axis='y', which='major', pad=8, length=8)
        self.axes.tick_params(axis='x', which='minor', pad=7, length=4)
        self.axes.tick_params(axis='y', which='minor', pad=8, length=4)

    def add_plot(self, x, y, color="r", label="", ls="-", lw=2, marker="", ms=7, markevery=1):
        line = self.axes.plot(x, y, label = label, markevery=markevery)
        plt.setp(line, ls=ls, c=color, lw=lw, marker=marker, mfc=color, ms=ms, mec=color, mew=0)
        self.nb_curves +=1

    def add_label(self, label, xloc=0.22, yloc=0.86, ha="left"):
        self.fig.text(xloc, yloc, label, ha=ha)

    def add_legend(self, location=0, fontsize=16):
        self.axes.legend(loc = location, fontsize = fontsize, frameon = False,
                         numpoints=1, markerscale = 1, handletextpad=0.5)

    def add_hlines(self, y=0, color="k", ls="--", linewith=0.6):
        self.axes.axhline(y=y, ls=ls, c=color, linewidth=linewith)

    def add_vlines(self, xvalue=0, color="k", ls="--", linewith=0.6):
        self.axes.axvline(x=xvalue, ls=ls, c=color, linewidth=linewith)

    def set_xticks(self, xtics):
        """
        Sets xtics, which is the space between two ticks
        # mxtics is the space between two minor ticks
        """
        mxtics = xtics / 2.
        majorFormatter = FormatStrFormatter('%g') # put the format of the number of ticks
        self.axes.xaxis.set_major_locator(MultipleLocator(xtics))
        self.axes.xaxis.set_major_formatter(majorFormatter)
        self.axes.xaxis.set_minor_locator(MultipleLocator(mxtics))

    def set_yticks(self, ytics):
        """
        Sets ytics, which is the space between two ticks
        # mytics is the space between two minor ticks
        """
        mytics = ytics / 2.
        majorFormatter = FormatStrFormatter('%g') # put the format of the number of ticks
        self.axes.yaxis.set_major_locator(MultipleLocator(ytics))
        self.axes.yaxis.set_major_formatter(majorFormatter)
        self.axes.yaxis.set_minor_locator(MultipleLocator(mytics))

    def set_xlog(self):
        self.axes.set_xscale('log')

    def set_xlog(self):
        self.axes.set_yscale('log')

    def print_figure(self, show=True):
        ## Create x and y labels
        self.axes.set_xlabel(self.xlabel, labelpad = 8)
        self.axes.set_ylabel(self.ylabel, labelpad = 8)
        ## Create x and y limits
        self.axes.set_xlim(self.xmin, self.xmax)
        self.axes.set_ylim(self.ymin, self.ymax)
        if show:
            plt.show()

    def save_figure(self, filename, folder="../figures", extension=".pdf", bbox_inches="tight"):
        self.fig.savefig(folder + "/" + filename + extension, bbox_inches=bbox_inches)

