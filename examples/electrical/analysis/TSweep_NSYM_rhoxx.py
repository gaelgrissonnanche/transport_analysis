from transport_analysis.electric import ResistivityTemp
from transport_analysis.figure import Figure
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Info sample ////////////////////////////////////////////////////////////////#
SAMPLE_NAME = r"Cu 2603B"
SAMPLE_LABEL = r"Cu_2603B"

## Which measurement? /////////////////////////////////////////////////////////#
FIELD, DATE = [0, "2026-03-20"]
# FIELD, DATE = [14, "2026-03-21"]

## Load data //////////////////////////////////////////////////////////////////#
rxx = ResistivityTemp(SAMPLE_LABEL, DATE, FIELD)
rxx.load_file()
# rxx._filter_crazy("T", n=20)
rxx.analysis()
rxx.save_file(["T", "rhoxx", "Ixx"])

fig_rxx = Figure(r"$T$ ( K )", r"$\rho$ ( $\mu\Omega$ cm )")
fig_rxx.add_label(SAMPLE_NAME)
fig_rxx.add_label(r"$B$ = " + str(FIELD) + " T", yloc=0.79)
fig_rxx.add_plot(rxx["T"], rxx["rhoxx"]*1e8)
fig_rxx.xmin = 0
fig_rxx.ymin = 0
fig_rxx.print_figure()
## Use the same name for the figure file as for the analyzed data file
fig_rxx.save_figure(rxx.file_analyzed)

