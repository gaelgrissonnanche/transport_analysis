from transport_analysis.thermoelectric import SeebeckTemp
from transport_analysis.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Info sample //////////////////////////////////////////////////////////////////#
SAMPLE_NAME = r"SLIO $x$ = 0.158 G37"
SAMPLE_LABEL = r"SLIOx=0.158G37S14"

## Which measu rement? ///////////////////////////////////////////////////////////#
FIELD, DATE = [0, "2026-04-10"]
# FIELD, DATE = [-14, "2025-02-22"]

## Analysis /////////////////////////////////////////////////////////////////////#
seebeck_exp = SeebeckTemp(SAMPLE_LABEL, DATE, FIELD, measurement_type="Tsteps_Seebeck")
seebeck_exp.use_DC = True
seebeck_exp.load_file()
seebeck_exp.analysis()
seebeck_exp.save_file(["Tav", "S_AC", "S_DC", "dT_AC", "dT_phase", "Vs_phase", "dT_DC", "T0"])


## Figures //////////////////////////////////////////////////////////////////////#

fig_list = []

## S/T --------------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$S$ / $T$ ( $\mu$V / K$^2$ )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_hlines(y=0)
# figure.add_plot(seebeck_exp.data_raw[:,1], seebeck_exp.data_raw[:,3], color="b", label="raw")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["S_DC"] / seebeck_exp["Tav"] * 1e6, color="#000000", marker="s", ls="", label="DC")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["S_AC"] / seebeck_exp["Tav"] * 1e6, color="#ff0000", marker="o", ls="", label="AC")
figure.add_legend(location=4)
# figure.ymin = -0.1
# figure.ymax = 0.1
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)


## S/T vs logT ------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$S$ / $T$ ( $\mu$V / K$^2$ )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_hlines(y=0)
# figure.add_plot(seebeck_exp.data_raw[:,1], seebeck_exp.data_raw[:,3], color="b", label="raw")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["S_DC"] / seebeck_exp["Tav"] * 1e6, color="#000000", marker="s", ls="", label="DC")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["S_AC"] / seebeck_exp["Tav"] * 1e6, color="#ff0000", marker="o", ls="", label="AC")
figure.add_legend(location=4)
# figure.ymin = -0.1
# figure.ymax = 0.1
# figure.xmin = 5
figure.set_xlog()
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)

## S ----------------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$S$ ( $\mu$V / K )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_hlines(y=0)
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["S_DC"] * 1e6, color="#000000", label="DC", marker="s", ls="")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["S_AC"] * 1e6, color="#fc4c50", label="AC", marker="o", ls="")
figure.add_legend(location=3)
# figure.ymax = 3
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)

## Voltage S --------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$V_S$ ( $\mu$V )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_hlines(y=0)
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["Vs_R"] * 1e6, color="#fdac4d", marker="s", ls="")
# figure.ymax = 3
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)

## dT AC --------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$dT_{\rm AC}$ ( mK )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["dT_AC"] * 1e3, color="#4351fe", marker="s", ls="")
# figure.ymin = 0
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)

## dT_DC / Tav --------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$dT_{\rm DC} / T $ ( mK )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["dT_DC"] / seebeck_exp["Tav"]*100, color="#000000", marker="s", ls="")
# figure.ymin = 0
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)

## Phases -----------------------------------------------------------------------#
figure = Figure(r"$T$ ( K )", r"$\phi$ ( degrees )")
figure.add_label(SAMPLE_NAME, xloc=0.79, ha="right")
figure.add_hlines(y=0)
figure.add_hlines(y=-180)
figure.add_hlines(y=180)
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["dT_phase"], color="#4e5aff", label = "dT", marker="s", ls="")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["Vs_phase"], color="#44ff2d", label = "Vs", marker="s", ls="")
figure.add_plot(seebeck_exp["Tav"], seebeck_exp["dT_phase"]-seebeck_exp["Vs_phase"], color="#000000", label = "Diff", marker="s", ls="")
figure.add_legend(location=4)
# figure.ymax = 260
figure.add_label("$B$ = " + str(FIELD) + " T", xloc=0.79, yloc=0.79, ha="right")
figure.print_figure()
fig_list.append(figure.fig)


## Save in one PDF //////////////////////////////////////////////////////////////#
file_figures = PdfPages("../figures/" + seebeck_exp.file_analyzed + ".pdf")

for fig in fig_list:
    file_figures.savefig(fig)

file_figures.close()

