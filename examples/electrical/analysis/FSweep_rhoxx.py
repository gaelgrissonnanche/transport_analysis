from transport_analysis.electric import ResistivityB
from transport_analysis.electric import HallEffectB
from transport_analysis.figure import Figure
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import os
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Info sample ////////////////////////////////////////////////////////////////#
SAMPLE_NAME = r"Cu 2603B"
SAMPLE_LABEL = r"Cu_2603B"


## Which measurement? /////////////////////////////////////////////////////////#
DATE = "2026-03-20"
## Load data //////////////////////////////////////////////////////////////////#
Ts = [2,7,10,15,20,30,40,50,60,70,80,100,120,150,200,250,300]

rxx_dict = {}
rxy_dict = {}
RH_dict = {}

for TEMPERATURE in Ts:
    rxx = ResistivityB(SAMPLE_LABEL, DATE, TEMPERATURE)
    rxx.load_file()
    rxx.analysis_MR()
    rxx_dict[TEMPERATURE] = rxx

    rxy =  HallEffectB(SAMPLE_LABEL,DATE,TEMPERATURE)
    rxy.load_file()
    rxy.analysis_Hall()
    rxy_dict[TEMPERATURE] = rxy








## Figures //////////////////////////////////////////////////////////////////////#

fig_list = []


## Resistivity raw ------------------------------------------------------#
figure = Figure(r"$B$ ( T )", r"$\rho_{xx}$ ( $\mu\Omega$ cm )")
figure.add_label(SAMPLE_NAME)
colors =  plt.cm.turbo(np.linspace(0, 1, len(Ts)))
for i, TEMPERATURE in enumerate(Ts):
    figure.add_plot(rxx_dict[TEMPERATURE]["B_raw"], rxx_dict[TEMPERATURE]["rhoxx_raw"]*1e8, color=colors[i],label=str(TEMPERATURE) + "K")
figure.add_legend(yloc=1.0)
figure.xmin = -14
#figure.ymin = 25
#figure.ymax = 70
figure.print_figure()
fig_list.append(figure.fig)



## Resistivity symmetryzed ------------------------------------------------------#
figure = Figure(r"$B$ ( T )", r"$\rho_{xx}$ ( $\mu\Omega$ cm )")
figure.add_label(SAMPLE_NAME)

colors =  plt.cm.jet(np.linspace(0, 1, len(Ts)))
for i, TEMPERATURE in enumerate(Ts):
    figure.add_plot(rxx_dict[TEMPERATURE]["B"], rxx_dict[TEMPERATURE]["rhoxx"]*1e8, color=colors[i],label=str(TEMPERATURE) + "K")
figure.add_legend(yloc=1.0)
#figure.xmin = 0
#figure.ymin = 25
#figure.ymax = 70
figure.print_figure()
fig_list.append(figure.fig)




## Resistivity symmetryzed ------------------------------------------------------#
figure = Figure(r"$B$ ( T )", r"$\frac{\rho_{xx}(\rm H)-\rho_{xx}( \rm 0)}{\rho_{xx}(0)}$ ( % ) ")
#figure = Figure(r"$B$ ( T )", r"$\R_{xx}$ ( $\Omega$  )")
figure.add_label(SAMPLE_NAME)
#figure.add_label(r"$B$ = " + str(np.abs(FIELD)) + " T", yloc=0.79)
colors =  plt.cm.jet(np.linspace(0, 1, len(Ts)))
data = np.loadtxt(f"../data_analyzed/" + SAMPLE_LABEL+ f"_TSweep_rhoxx_B=0T.txt", dtype="float", comments="#", skiprows=2,delimiter=',')
sorted_indices = np.argsort(data[:,0])
T_sorted = data[sorted_indices,0]
rhoxx_sorted = data[sorted_indices,1]

for i, TEMPERATURE in enumerate(Ts):

    rxx0 = np.interp(TEMPERATURE, T_sorted, rhoxx_sorted)
    rxx0 = rxx_dict[TEMPERATURE]["rhoxx"][0]
    figure.add_plot(rxx_dict[TEMPERATURE]["B"], (rxx_dict[TEMPERATURE]["rhoxx"]-rxx0)/rxx0 *100, color=colors[i],label=str(TEMPERATURE) + "K")

figure.add_legend(yloc=1.0)
#figure.xmin = 0
#figure.ymin = 25
#figure.ymax = 70
figure.print_figure()
fig_list.append(figure.fig)



## rhoxy raw --------------------------------------------------#
figure = Figure(r"$B$ ( T )", r"$\rho_{xy}$ ( $\mu\Omega$ cm )")
figure.add_label(SAMPLE_NAME,yloc=0.22)
colors =  plt.cm.jet(np.linspace(0, 1, len(Ts)))
for i, TEMPERATURE in enumerate(Ts):
    #print(rxy_dict[TEMPERATURE]["B_raw"])
    figure.add_plot(rxy_dict[TEMPERATURE]["B_raw"], rxy_dict[TEMPERATURE]["rhoxy_raw"]*1e8, color=colors[i],label=str(TEMPERATURE) + "K")
figure.add_hlines()
figure.add_legend(yloc=1.0)
figure.xmin = -14
#figure.ymin = 0
figure.print_figure()
fig_list.append(figure.fig)




## rhoxy --------------------------------------------------#
figure = Figure(r"$B$ ( T )", r"$\rho_{xy}$ ( $\mu\Omega$ cm )")
#figure = Figure(r"$B$ ( T )", r"$R_{xy}$ ( $\Omega$  )")
figure.add_label(SAMPLE_NAME,yloc=0.22)
#figure.add_label(r"$B$ = " + str(np.abs(FIELD)) + " T", yloc=0.79)
colors =  plt.cm.jet(np.linspace(0, 1, len(Ts)))
for i, TEMPERATURE in enumerate(Ts):
    figure.add_plot(rxy_dict[TEMPERATURE]["B"], rxy_dict[TEMPERATURE]["rhoxy"]*1e8, color=colors[i],label=str(TEMPERATURE) + "K")
    #figure.add_plot(rxy_dict[TEMPERATURE]["B"], rxy_dict[TEMPERATURE]["rhoxy"], color=colors[i],label=str(TEMPERATURE) + "K")
figure.add_hlines()
#figure.add_plot(rxx2["B"], rxx2["rhoxx"]*1e8, color="k",label="2.04K")
figure.add_legend(yloc=1.0)
#figure.xmin = 0
#figure.ymin = 0
figure.print_figure()
fig_list.append(figure.fig)



## RH --------------------------------------------------#
figure = Figure(r"$B$ ( T )", r"$R_H$ ( mm$^3$ / C )")
#figure = Figure(r"$B$ ( T )", r"$R_{xy}$/B ")
figure.add_label(SAMPLE_NAME)
#figure.add_label(r"$B$ = " + str(np.abs(FIELD)) + " T", yloc=0.79)
colors =  plt.cm.jet(np.linspace(0, 1, len(Ts)))
for i, TEMPERATURE in enumerate(Ts):
    figure.add_plot(rxy_dict[TEMPERATURE]["B"], rxy_dict[TEMPERATURE]["RH"]*1e9, color=colors[i],label=str(TEMPERATURE) + "K")
    #figure.add_plot(rxy_dict[TEMPERATURE]["B"], rxy_dict[TEMPERATURE]["RH"], color=colors[i],label=str(TEMPERATURE) + "K")
figure.add_hlines()
#figure.add_plot(rxx2["B"], rxx2["rhoxx"]*1e8, color="k",label="2.04K")
figure.add_legend(yloc=1.0)
#figure.xmin = 0
figure.ymin = -0.05
figure.ymax = 0.01
figure.print_figure()
fig_list.append(figure.fig)



# Save to a text file with header
for i, TEMPERATURE in enumerate(Ts):
    data_to_txt_rhoxx = np.column_stack((rxx_dict[TEMPERATURE]["B"], rxx_dict[TEMPERATURE]["rhoxx"]))
    data_to_txt_rhoxy = np.column_stack((rxy_dict[TEMPERATURE]["B"], rxy_dict[TEMPERATURE]["rhoxy"],rxy_dict[TEMPERATURE]["RH"]))



    np.savetxt("../data_analyzed/" + SAMPLE_LABEL+ "_FSweep_rhoxx_B=-14_to_14T_"+ str(TEMPERATURE) + "K.txt", data_to_txt_rhoxx, header="B (T), rhoxx (m.Ohm)", delimiter=",")
    np.savetxt("../data_analyzed/" + SAMPLE_LABEL+ "_FSweep_rhoxy_RH_B=-14_to_14T_"+ str(TEMPERATURE) + "K.txt", data_to_txt_rhoxy, header="B (T), rhoxy (m.Ohm),RH (m.Ohm/T)", delimiter=",")




######################################################
script_name = os.path.basename(__file__) # donne le nom du script avec l’extension .py
figurename = script_name[0:-3] + ".pdf"

## Save in one PDF //////////////////////////////////////////////////////////////#
file_figures = PdfPages("../figures/" + SAMPLE_LABEL + "_" + figurename)

for fig in fig_list:
    file_figures.savefig(fig)

file_figures.close()
