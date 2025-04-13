import numpy as np
from skimage.restoration import unwrap_phase
from scipy.optimize import brentq
import pickle
import scipy.ndimage as im # For removing crazy points
from .thermocouple import S_ther, E_ther
from numpy import sqrt, arctan2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

class SeebeckT:
    def __init__(self, samplelabel, date, field, measurement_type="Tsweep_Seebeck_NSYM"):
        self.samplelabel = samplelabel
        self.date        = date
        self.field       = field
        self.info        = None # contains all the info related to the measurement
        self.data_raw    = np.empty(1) # raw data array
        self.data        = {} # dict of all analyzed data
        self.L = None # length between contacts along x
        self.w = None # width of the sample along y
        self.t = None # thickness of the sample along z
        ## Files
        self.file_raw = None # store the name of the raw file
        self.folder_raw = None # store the name of the raw folder
        self.file_analyzed = None # store the name of the analyzed file
        self.measurement_type = measurement_type  # is used in the analyzed file name
        ## Do you want to use T_DC or T0
        self.use_DC = False
        ## Do you want to unwrap the phases
        self.unwrap = True
        ## Change sign of Seebeck
        self.sign_S = 1
        self.sign_Tp = 1
        self.sign_Tm = 1

    ## Methods >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def load_info(self):
        """
        Load the measurement pickle for the corresponding samplelabel. It contains
        the info regarding the measurements. The keys of info are (Field, Date), and
        return for each couple of (Field,Date) a particular dictionnary for each
        measurements. This dictionnary has for keys : \n

        "folder" : the name of the folder where data are stored in "data_r" \n
        "comments" : comments regarding this particular set of points \n
        "geofactor" : L the length, w the width and t the thickness (meters) \n
        "columns" : columns in the raw data file \n
        """
        with open("info_tsweep_" + self.samplelabel + ".pkl", "rb") as f: # 'w' means (over)write or create if doesn't exist
            unpickler_file = pickle.Unpickler(f)
            self.info = unpickler_file.load()[self.field, self.date] # keys of measurements are (H,date)

    def load_file(self, delimiter=","):
        """
        Loads the pickle and the data files
        """
        ## Load the pickle
        self.load_info()
        self.L = self.info["geofactor"]["L"] # length between contacts along x
        self.w = self.info["geofactor"]["w"] # width of the sample along y
        self.t = self.info["geofactor"]["t"] # thickness of the sample along z
        self.columns = self.info["columns"] # list of the columns in data file
        self.comments = self.info["comments"]
        ## Load the data file
        self.folder_raw = self.info["folder"]
        self.file_raw = self.info["file"]
        self.data_raw = np.loadtxt("../data_raw/" + self.folder_raw + "/"
                                   + self.file_raw,
                        dtype = "float", comments = "#", delimiter=delimiter)
        ## Extract raw data columns
        self.data = self.extract_raw_data(self.data_raw)

    def extract_raw_data(self, data_raw):
        data = {}
        data["T0"] = data_raw[:,self.columns["T0"]] # ref temperature
        # T+
        data["Vp_X"] = self.sign_Tp * data_raw[:,self.columns["Vp_X"]] * sqrt(2) # Real part of the voltage of the T+ thermocouple
        data["Vp_Y"] = self.sign_Tp * data_raw[:,self.columns["Vp_Y"]] * sqrt(2) # Im part of the voltage of the T+ thermocouple
        data["Vp_R"]   = np.absolute(data["Vp_X"] + 1j * data["Vp_Y"])
        data["Vp_phi"] = np.angle(data["Vp_X"] + 1j * data["Vp_Y"], deg=True) # degrees
        # T-
        data["Vm_X"] = self.sign_Tm * data_raw[:,self.columns["Vm_X"]] * sqrt(2) # Real part of the voltage of the T- thermocouple
        data["Vm_Y"] = self.sign_Tm * data_raw[:,self.columns["Vm_Y"]] * sqrt(2) # Im part of the voltage of the T- thermocouple
        data["Vm_R"]   = np.absolute(data["Vm_X"] + 1j * data["Vm_Y"])
        data["Vm_phi"] = np.angle(data["Vm_X"] + 1j * data["Vm_Y"], deg=True)  # degrees
        # Seebeck voltage
        data["Vs_X"] = self.sign_S * data_raw[:,self.columns["Vs_X"]] * sqrt(2) # Real part of the Seebeck voltage
        data["Vs_Y"] = self.sign_S * data_raw[:,self.columns["Vs_Y"]] * sqrt(2) # Im part of the Seebeck voltage
        data["Vs_R"]   = np.absolute(data["Vs_X"] + 1j * data["Vs_Y"])
        data["Vs_phi"] = np.angle(data["Vs_X"] + 1j * data["Vs_Y"], deg=True) # degrees
        ## DC
        data["Vp_DC_on"] = data_raw[:,self.columns["Vp_DC_on"]] # DC voltage of the T+ thermocouple, HEAT ON
        data["Vm_DC_on"] = data_raw[:,self.columns["Vm_DC_on"]] # DC voltage of the T- thermocouple, HEAT ON
        data["Vs_DC_on"] = data_raw[:,self.columns["Vs_DC_on"]] # DC Seebeck voltage, HEAT ON
        ## Extract current frequency
        data["freq"] = 2 * data_raw[:,self.columns["freq"]] # heat frequency
        data["I_amp"] = data_raw[:,self.columns["I_amp"]]
        data["I_offset"] = data_raw[:,self.columns["I_offset"]]
        data["R_heater"] = data_raw[:,self.columns["R_heater"]]
        ## Check for the magnetific field column
        try:
            data["B"] = data_raw[:,self.columns["B"]]  # try to see if it can load B
        except:
            data["B"] = np.ones_like(data["T0"]) * self.field
        try:
            data["Vp_DC_off"] = data_raw[:,self.columns["Vp_DC_off"]] # DC voltage of the T+ thermocouple, HEAT OFF
            data["Vm_DC_off"] = data_raw[:,self.columns["Vm_DC_off"]] # DC voltage of the T- thermocouple, HEAT OFF
            data["Vs_DC_off"] = data_raw[:,self.columns["Vs_DC_off"]] # DC Seebeck voltage, HEAT OFF
        except:
            data["Vp_DC_off"] = np.zeros_like(data["T0"])
            data["Vm_DC_off"] = np.zeros_like(data["T0"])
            data["Vs_DC_off"] = np.zeros_like(data["T0"])
        ## Remove background
        data["Vp_DC"]  = np.abs(data["Vp_DC_on"] - data["Vp_DC_off"])
        data["Vm_DC"]  = np.abs(data["Vm_DC_on"] - data["Vm_DC_off"])
        data["Vs_DC"]  = self.sign_S * (data["Vs_DC_on"] - data["Vs_DC_off"])
        return data

    def clean_up(self):
        """
        Clean the data by:
        - sorting the data in ascending order of temperature
        - removing the points with the same temperature values
        """
        T = self.data["T0"]
        ## Sort in respect to temperature
        index_sort = np.argsort(T)
        T = T[index_sort]
        for key in self.data.keys():
            self.data[key] = self.data[key][index_sort]
        ## Unique
        T, index_unique = np.unique(T, return_index=True)
        for key in self.data.keys():
            self.data[key] = self.data[key][index_unique]
        ## Save
        self.data["T0"] = T

    def calc_temp_ac(self):
        ## AC >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        self.data["Tp_R"]   = self._ac_temp_thermocouples(self.data["T0"], self.data["Vp_R"], self.data["Tp_DC"])
        self.data["Tm_R"]   = self._ac_temp_thermocouples(self.data["T0"], self.data["Vm_R"], self.data["Tm_DC"])
        ## Unwrap the phase
        if self.unwrap == True:
            self._unwrap_phase()
        else:
            self.data["Tp_phi"] = self.data["Vp_phi"] # degrees
            self.data["Tm_phi"] = self.data["Vm_phi"] # degrees

        Tp_X = self.data["Tp_R"] * np.cos(np.deg2rad(self.data["Tp_phi"]))
        Tp_Y = self.data["Tp_R"] * np.sin(np.deg2rad(self.data["Tp_phi"]))
        Tm_X = self.data["Tm_R"] * np.cos(np.deg2rad(self.data["Tm_phi"]))
        Tm_Y = self.data["Tm_R"] * np.sin(np.deg2rad(self.data["Tm_phi"]))
        ## AC thermal gradient
        self.data["dT_AC"] = sqrt((Tp_X - Tm_X)**2 + (Tp_Y - Tm_Y)**2)
        self.data["dT_phi"] = np.rad2deg(arctan2(self.data["Vp_Y"] - self.data["Vm_Y"],
                                                 self.data["Vp_X"] - self.data["Vm_X"]))

    def _ac_temp_thermocouples(self, T0, dV_AC, T_DC):
        # T_AC = np.zeros_like(T0)  # Initialize Tx array with zeros
        # # Vectorized brentq calculation for each element in the input arrays
        # if self.use_DC == False:
        #     T_DC = T0
        # for i in range(len(T0)):
        #     T_AC[i] = brentq(lambda x: x*S_ther(T_DC[i]) - dV_AC[i], 0, 1000)
        if self.use_DC == False:
            T_DC = T0
        T_AC   = dV_AC / S_ther(T_DC)
        return T_AC

    def calc_temp_dc(self):
        self.data["Tp_DC"]  = self._dc_temp_thermocouples(self.data["T0"], self.data["Vp_DC"])
        self.data["Tm_DC"]  = self._dc_temp_thermocouples(self.data["T0"], self.data["Vm_DC"])
        self.data["dT_DC"] = self.data["Tp_DC"] - self.data["Tm_DC"]
        self.data["Tav"] = (self.data["Tp_DC"] + self.data["Tm_DC"]) / 2
        ## Calculate Power DC
        self._power_dc()

    def _dc_temp_thermocouples(self, T0, dV_DC):
        # T_DC = np.zeros_like(T0)  # Initialize T_DC array with zeros
        # # Vectorized brentq calculation for each element in the input arrays
        # for i in range(len(T0)):
        #     T_DC[i] = brentq(lambda x: E_ther(x) - E_ther(T0[i]) - dV_DC[i], 0, 1000)
        T_DC = dV_DC / S_ther(T0) + T0
        return T_DC

    def _power_dc(self):
        """
        Calculate the amount of DC power based on the I_amp and I_offset
        """
        self.data["Q_DC"] = self.data["R_heater"] * ((self.data["I_amp"] / np.sqrt(2))**2 + self.data["I_offset"]**2)  # W
        self.data["j0"] = self.data["Q_DC"] / (self.w * self.t)

    def _unwrap_phase(self):
        self.data["Tp_phi"] = unwrap_phase(np.deg2rad(self.data["Vp_phi"] - 90))
        self.data["Tm_phi"] = unwrap_phase(np.deg2rad(self.data["Vm_phi"] - 90))
        self.data["Tp_phi"] = np.rad2deg(self.data["Tp_phi"]) + 90
        self.data["Tm_phi"] = np.rad2deg(self.data["Tm_phi"]) + 90
        ## To make sure that Tp_phi is always larger than Tm_phi and seperated by
        ## more than 180 degrees
        if self.data["Tp_phi"][0] < self.data["Tm_phi"][0]:
            self.data["Tp_phi"] += np.round((self.data["Tp_phi"][0] - self.data["Tm_phi"][0]) / 180) * 180

    def calc_seebeck(self):
        self.data["S_AC"] = self.data["Vs_R"] / self.data["dT_AC"] * np.cos(np.deg2rad(self.data["dT_phi"]-self.data["Vs_phi"]))
        self.data["S_DC"] = self.data["Vs_DC"] / self.data["dT_DC"]

    def create_file_name_analyzed(self):
        """
        Create the name of the analyzed file
        """
        self.file_analyzed = self.samplelabel + "_" + self.measurement_type + \
                    "_B=" + str(self.field) + 'T_' + str(self.date)

    def analysis(self, clean_up=True):
        """
        Perform the analysis
        """
        if clean_up:
            self.clean_up()
        self.calc_temp_dc()
        self.calc_temp_ac()
        self.calc_seebeck()
        self.create_file_name_analyzed()


    def save_file(self, to_save_list, folder='../data_analyzed', delimiter=',', extension=".txt"):
        """
        Save data in the analyzed folder after analysis
        - list_data: ["T0", "Tp", etc.] what needs to be saved in self.data
        - measurement: "TS_rxx_NSYM" for example
        """
        ## Create the array to save data
        data_s = np.empty((self.data["T0"].size, len(to_save_list)))
        for i, key in enumerate(to_save_list):
            data_s[:, i] = self.data[key]
        ## Save the data file
        np.savetxt(folder + '/' + self.file_analyzed + extension, data_s, delimiter=delimiter,
               fmt='%.7e', header = ",".join(to_save_list), comments = "#")

    ## Filters
    def _filter_crazy(self, key, n=50):
        ## Filter temperature points
        self.data[key] = im.filters.median_filter(self.data[key], size = n)
        # size correspond to the number of neighboors it takes to compare
        # the crazy points with the good ones. Then it corrects if it is crazy
        return self.data[key]

    ## Special Method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def __setitem__(self, key, value):
        """
        Set the ability to access the different data in the dictionnary self.data
        by just calling resistivity_obj["T"], resistivity_obj["rhoxx"], etc.
        """
        self.data[key] = value

    def __getitem__(self, key):
        """
        Get the ability to access the different data in the dictionnary self.data
        by just calling resistivity_obj["T"], resistivity_obj["rhoxx"], etc.
        """
        try:
            return self.data[key]
        except KeyError:
            print(f"{key} is not a defined")










    # def analysis(self):

    #     ## Thermal diffusivity >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #     # a = self.Vp_Y / self.Vp_X
    #     # b = self.Vm_Y / self.Vm_X
    #     # self.D_phi = unwrap_phase(np.arctan((a - b)/(1+a*b))*2)/2
    #     self.D_phi = np.deg2rad(self.Tp_phi - self.Tm_phi)
    #     self.l_phase = self.L / self.D_phi
    #     self.l_amp = self.L / np.log(self.Tp_R/self.Tm_R)
    #     self.alpha = (2*np.pi*self.freqs) * self.l_amp * self.l_phase / 2
    #     self.alpha_amp = (2*np.pi*self.freqs) * self.l_amp**2 / 2
    #     self.alpha_phase = (2*np.pi*self.freqs) * self.l_phase**2 / 2
    #     ## Kxx DC >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #     self.Kxx_DC = self.Q_DC / self.dTx_DC / (self.w * self.t / self.L)
    #     # print("x+ / l_amp")
    #     # print(self.xp / self.l_amp)
    #     # print("T+_R")
    #     # print(self.Tp_R)
    #     # print("x- / l_amp")
    #     # print(self.xm / self.l_amp)
    #     # print("T-_R")
    #     # print(self.Tm_R)
    #     # print("ratio")
    #     # print(1 / self.Tp_R * np.exp(- self.xp / self.l_amp))
    #     # print("ratio")
    #     # print(1 / self.Tm_R * np.exp(- self.xm / self.l_amp))
    #     self.Kxx_AC_p = self.j0 / self.Tp_R * np.exp(- self.xp / self.l_amp) / np.sqrt(1/self.l_amp**2 + 1/self.l_phase**2)
    #     self.Kxx_AC_m = self.j0 / self.Tm_R * np.exp(- self.xm / self.l_amp) / np.sqrt(1/self.l_amp**2 + 1/self.l_phase**2)
    #     # In theory Kxx_AC_p = Kxx_AC_m as a consequence of the analysis
    #     self.Kxx_AC   = np.sqrt(self.Kxx_AC_p * self.Kxx_AC_m)
    #     # self.Kxx_AC   = self.Kxx_AC_p
    #     # self.Kxx_AC   = self.Kxx_AC_m
    #     ## Specific heat >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #     self.Cp_DC = self.Kxx_DC / self.alpha / self.density # in J / K kg
    #     self.Cp_DC_P = self.Kxx_DC / self.alpha_amp / self.density # in J / K kg
    #     self.Cp_DC_Q = self.Kxx_DC / self.alpha_phase / self.density # in J / K kg
    #     self.Cp_AC = self.Kxx_AC_p / self.alpha / self.density # in J / K kg
    #     self.Cp_AC_P = self.Kxx_AC / self.alpha_amp / self.density # in J / K kg
    #     self.Cp_AC_Q = self.Kxx_AC / self.alpha_phase / self.density # in J / K kg



    # def save(self):
    #     data = np.vstack((self.Tav, self.Kxx_DC, self.Kxx_AC,
    #                       self.alpha, self.alpha_amp, self.alpha_phase,
    #                       self.Cp_DC, self.Cp_DC_P, self.Cp_DC_Q,
    #                       self.Cp_AC, self.Cp_AC_P, self.Cp_AC_Q,
    #                       self.Tp_R, self.Tp_phi, self.Tm_R, self.Tm_phi,
    #                       self.Tp_DC, self.Tm_DC,
    #                       self.freqs, self.I_amp, self.R_heater,
    #                       self.T0)).transpose()
    #     header = ",".join([
    #               "Tav[K]", "Kxx_DC[W/K*m]", "Kxx_AC[W/K*m]",
    #               "alpha[m^2/s]", "alpha_amp[m^2/s]", "alpha_phase[m^2/s]",
    #               "Cp_DC[J/K*kg]", "Cp_DC_amp[J/K*kg]", "Cp_DC_phase[J/K*kg]",
    #               "Cp_AC[J/K*kg]", "Cp_AC_amp[J/K*kg]", "Cp_AC_phase[J/K*kg]",
    #               "|T+|(K)", "Phi+(deg)", "|T-|(K)", "Phi-(deg)",
    #               "T+(K)", "T-(K)",
    #               "f(Hz)", "I_amp(A)", "R_heater(Ohm)",
    #               "T0(K)"
    #               ])
    #     np.savetxt('../data_analyzed/' + self.s_filename + ".dat",
    #                 data, fmt='%.7e',
    #                 header = header,
    #                 comments = "#", delimiter=',')

