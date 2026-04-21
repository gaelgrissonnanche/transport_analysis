import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample ////////////////////////////////////////////////////////////////#
SAMPLE_LABEL = r"Cu_2603B"

## Create dictionnaries of measurements info //////////////////////////////////#
info = {}     # keys are (field, date)



info[0, "2026-03-20"] = {
    "folder": "2026-03-20- Resistivity cooldown B=0T",
    "file": "TS-300K-to-2K-Cu_2603B-Cu_2603A.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3}, # change geofactors
    "columns":  {"T": 1, "B":3, "Rxx": 4, "Ixx":5, "Rxy": 7, "Ixy":8},
                        }


#################################################################################################################
###################################### Magnetoresistance ########################################################
#################################################################################################################

#key [Temperature, date]
info[2, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=2K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[7, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=7K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }


info[10, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=10K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[15, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=15K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[20, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=20K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[30, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=30K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[40, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=40K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[50, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=50K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[60, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=60K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[70, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=70K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[80, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=80K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[100, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=100K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[120, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=120K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[150, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=150K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[200, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=200K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[250, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=250K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }

info[300, "2025-09-27"] = {
    "folder": "",
    "file": "2025-09-27 - Fieldsweep B=-14T_to_14T/2025-09-27-LSNO_0p12_2509B_Al_2509B_Cu_2509A2_Tsweep_B=14T_to_-14T_T=300K.txt",
    "comments": "",
    "geofactor": { "L" : 1.68*1e-3, "w" : 0.53*1e-3, "t" : 0.025*1e-3},
    "columns":  {"T": 1, "B":3, "Rxx": 10, "Ixx":11, "Rxy": 13, "Ixy":14},
                        }




## Save info_dict with Pickle ///////////////////////////
# //////////////////////#
with open("info_tsweep_" + SAMPLE_LABEL + ".pkl", "wb") as f:
    pickler_file = pickle.Pickler(f)
    pickler_file.dump(info)

