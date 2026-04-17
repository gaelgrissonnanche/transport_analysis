import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample ////////////////////////////////////////////////////////////////#
SAMPLE_LABEL = r"SLIOx=0.158G37S14"

## Create dictionnaries of measurements info //////////////////////////////////#
info = {}     # keys are (field, date)

info[0, "2026-04-14"] = {
    "folder": "2026-04-14_Cooldown_0T_PEG=1mbar",
    "file": "2026-04-14-SLIO_x=0p158_G37_S14_T=300_to_2K_B=0T-Vout=2V_substrate-glass_PVTI=6mbar.txt",
    "comments": "",
    "columns":  {"T0": 1, "B":2,
                 "Vp_R": 6, "Vp_phase": 7, "Vm_R": 10, "Vm_phase": 11, "Vs_R": 14, "Vs_phase": 15,
                 "Vp_DC": 16, "Vm_DC": 17, "Vs_DC": 18,
                 },
    "signs": {"Tp": -1, "Tm": 1, "Vs": 1},
                      }

info[0, "2026-04-10"] = {
    "folder": "2026-04-10_Steps_PEG=120mbar",
    "file": "TS-0.00T-TOP-SLIOx=0.158G37S14-2026-04-10.dat",
    "comments": "",
    "columns":  {"T0": 0, "B":1,
                 "Vp_R": 9, "Vp_phase": 10, "Vm_R": 13, "Vm_phase": 14, "Vs_R": 17, "Vs_phase": 18,
                 "Vp_DC": 19, "Vm_DC": 20, "Vs_DC": 21,
                 "Vp_DC_OFF": 22, "Vm_DC_OFF": 23, "Vs_DC_OFF": 24,
                 },
    "signs": {"Tp": -1, "Tm": 1, "Vs": 1},
                      }


## Save info_dict with Pickle /////////////////////////////////////////////////#
with open("info_tsweep_" + SAMPLE_LABEL + ".pkl", "wb") as f:
    pickler_file = pickle.Pickler(f)
    pickler_file.dump(info)

