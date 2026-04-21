import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter
import os
from transport_analysis.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

## Info sample ////////////////////////////////////////////////////////////////#
SAMPLE_NAME = r"Al 2509B"
SAMPLE_LABEL = r"Al_2509B"


##############################################################################
### Plotting #################################################################
##############################################################################

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
mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                    # editing the text in illustrator


####################################################
## Plot Parameters #################################

##################################################################################

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

rhoxy=[]
RH=[]
FIELD= 14.0

Ts = [2,7,10,15,20,30,40,50,60,70,80,100,120,150,200,250,300]

for T in Ts:
    data = np.loadtxt(f"../data_analyzed/" + SAMPLE_LABEL+ f"_FSweep_rhoxy_RH_B=-14_to_14T_{T}K.txt", dtype="float", comments="#", skiprows=2,delimiter=',')

    filtered = data[ data[:,0] == FIELD ]

    if filtered.size !=0:
        rhoxy += [np.mean(filtered[:,1])]
        RH += [np.mean( filtered[:,2])]


    if filtered.size == 0:
        sorted_indices = np.argsort(data[:,0])
        x_sorted = data[sorted_indices,0]
        rhoxy_sorted = data[sorted_indices,1]
        RH_sorted = data[sorted_indices,2]
        rhoxy_interp = np.interp(FIELD, x_sorted,rhoxy_sorted)
        RH_interp = np.interp(FIELD, x_sorted,RH_sorted)

        rhoxy += [rhoxy_interp]   #rhoxy
        RH += [RH_interp]



rhoxy=np.array(rhoxy)
RH = np.array(RH)


"""
data_TS = np.loadtxt(f"../data_analyzed\LCCO_0p16_2505C_TSweep_rhoxy_RH_B=14T.txt", dtype="float", comments="#", skiprows=2,delimiter=',')
T_TS = data_TS[:,0]
rhoxy_TS = data_TS[:,1]
RH_TS = data_TS[:,2]
"""

fig_list=[]
figure = Figure(r"$T$ ( K )", r"$\rho_{xy}$ ( $\mu\Omega$ cm )")
figure.add_label(SAMPLE_NAME)
figure.add_label(r"$B$ = " + str(FIELD) + " T", yloc=0.79)
figure.add_hlines()
colors =  plt.cm.jet(np.linspace(0, 1, 3))
#figure.add_plot(T_TS, rhoxy_TS*1e8, color=colors[0],label=f"TS")
figure.add_plot(Ts, rhoxy*1e8, color="k",label=f"FS",marker=".",ms=10)
#figure.add_plot(rxx2["B"], rxx2["rhoxx"]*1e8, color="k",label="2.04K")
figure.add_legend(xloc=0.98,yloc=0.3)
#figure.xmin = 0
#figure.ymin = 0
#figure.ymax=0
figure.print_figure()
fig_list.append(figure.fig)


figure = Figure(r"$T$ ( K )", r"$R_{H}$ ( mm$^3$ / C )")
figure.add_label(SAMPLE_NAME)
figure.add_hlines()
figure.add_label(r"$B$ = " + str(FIELD) + " T", yloc=0.79)
colors =  plt.cm.jet(np.linspace(0, 1, 3))
#figure.add_plot(T_TS, RH_TS*1e9, color=colors[0],label=f"TS")
figure.add_plot(Ts, RH*1e9, color="k",label=f"FS",marker=".",ms=10)
#figure.add_plot(rxx2["B"], rxx2["rhoxx"]*1e8, color="k",label="2.04K")
figure.add_legend(xloc=0.98,yloc=0.3)
#figure.xmin = 0

#figure.ymin = 0
#figure.ymax = 0
figure.print_figure()
fig_list.append(figure.fig)




script_name = os.path.basename(__file__) # donne le nom du script avec l’extension .py
#figurename = "figures/" + SAMPLE_LABEL + script_name[0:-3] + ".pdf"
figurename = "../figures/" + SAMPLE_LABEL + f"_FSweep_Hall_B={FIELD}T.pdf"

data_to_txt = np.column_stack((Ts, RH,rhoxy))
np.savetxt("../data_analyzed/" + SAMPLE_LABEL+ f"_FSweep_rhoxy_Hall_B={FIELD}T.txt", data_to_txt, header="B (T), RH (m.Ohm/T), rhoxy (m.Ohm)", delimiter=",")



## Save in one PDF //////////////////////////////////////////////////////////////#
file_figures = PdfPages(figurename)

for fig in fig_list:
    file_figures.savefig(fig)

file_figures.close()