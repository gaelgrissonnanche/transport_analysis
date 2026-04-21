import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample ////////////////////////////////////////////////////////////////#
SAMPLE_LABEL = r"Cu_2603B"

## Create dictionnaries of measurements info //////////////////////////////////#
info = {}     # keys are (Temperature, date)


info[10, "2026-03-20"] = {
    "folder": "2026-03-20_FS_B=10Tto-10T",
    "file": "FS-10.00K-TOP-Cu_Modal-2026-2026-03-20.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 4, "Ixx":5, "Rxy": 7, "Ixy":8},
                        }

info[20, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-20.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[40, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-40.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[45, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-45.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[50, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-50.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[60, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-60.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[80, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-80.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[90, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-90.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[100, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-100.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[140, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-140.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[160, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-160.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[180, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-180.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[200, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-200.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[240, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-240.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[270, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-270.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})
info[300, "2026-03-20"] = info[10, "2026-03-20"].update({"file": "FS-300.00K-TOP-Cu_Modal-2026-2026-03-20.txt"})




## Save info_dict with Pickle ///////////////////////////
# //////////////////////#
with open("info_fsweep_" + SAMPLE_LABEL + ".pkl", "wb") as f:
    pickler_file = pickle.Pickler(f)
    pickler_file.dump(info)

