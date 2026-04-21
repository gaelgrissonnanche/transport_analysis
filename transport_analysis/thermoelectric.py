import numpy as np
from skimage.restoration import unwrap_phase
from scipy.optimize import brentq
import pickle
import scipy.ndimage as im # For removing crazy points
from .thermocouple import Sther, Ether, Ether_inv
from numpy import sqrt, deg2rad, rad2deg
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

class SeebeckTemp:
    def __init__(self, samplelabel, date, field, measurement_type="Tsweep_Seebeck"):
        self.samplelabel = samplelabel
        self.date        = date
        self.field       = field
        self.info        = None # contains all the info related to the measurement
        self.data_raw    = np.empty(1) # raw data array
        self.data        = {} # dict of all analyzed data
        ## Files
        self.file_raw = None # store the name of the raw file
        self.folder_raw = None # store the name of the raw folder
        self.file_analyzed = None # store the name of the analyzed file
        self.measurement_type = measurement_type  # is used in the analyzed file name
        ## Do you want to use T_DC or T0
        self.use_DC = False

    ## Methods >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def load_info(self):
        """
        Load the measurement pickle for the corresponding samplelabel. It contains
        the info regarding the measurements. The keys of info are (Field, Date), and
        return for each couple of (Field,Date) a particular dictionnary for each
        measurements. This dictionnary has for keys : \n

        "folder" : the name of the folder where data are stored in "data_r" \n
        "comments" : comments regarding this particular set of points \n
        "columns" : columns in the raw data file \n
        """
        with open("info_tsweep_" + self.samplelabel + ".pkl", "rb") as f: # 'w' means (over)write or create if doesn't exist
            unpickler_file = pickle.Unpickler(f)
            self.info = unpickler_file.load()[self.field, self.date] # keys of measurements are (H,date)

    def load_file(self, delimiter=",", skiprows=0):
        """
        Loads the pickle and the data files
        """
        ## Load the pickle
        self.load_info()
        self.columns = self.info["columns"] # list of the columns in data file
        self.comments = self.info["comments"]
        ## Load the data file
        self.folder_raw = self.info["folder"]
        self.file_raw = self.info["file"]
        self.data_raw = np.loadtxt("../data_raw/" + self.folder_raw + "/"
                                   + self.file_raw,
                        dtype = "float", comments = "#", delimiter=delimiter, skiprows=skiprows)
        ## Change sign of signals
        self.sign_S = self.info["signs"]["Vs"]
        self.sign_Tp = self.info["signs"]["Tp"]
        self.sign_Tm = self.info["signs"]["Tm"]
        ## Extract raw data columns
        self.data = self.extract_raw_data(self.data_raw)

    def extract_raw_data(self, data_raw):
        data = {}
        data["T0"] = data_raw[:,self.columns["T0"]] # ref temperature
        # Voltage T+
        data["Vp_R"] = data_raw[:,self.columns["Vp_R"]] * sqrt(2) # modulus
        data["Vp_phase"] = data_raw[:,self.columns["Vp_phase"]] + rad2deg(np.angle(self.sign_Tp)) # phase
        data["Vp_phase"] = rad2deg(np.angle(np.exp(1j*deg2rad(data["Vp_phase"])))) # wrap back into -180 to 180
        # Voltage T-
        data["Vm_R"] =  data_raw[:,self.columns["Vm_R"]] * sqrt(2) # modulus
        data["Vm_phase"] = data_raw[:,self.columns["Vm_phase"]] + rad2deg(np.angle(self.sign_Tm)) # phase
        data["Vm_phase"] = rad2deg(np.angle(np.exp(1j*deg2rad(data["Vm_phase"])))) # wrap back into -180 to 180
        # Seebeck voltage
        data["Vs_R"] = data_raw[:,self.columns["Vs_R"]] * sqrt(2) # modulus
        data["Vs_phase"] = data_raw[:,self.columns["Vs_phase"]] + rad2deg(np.angle(self.sign_S))  # phase
        data["Vs_phase"] = rad2deg(np.angle(np.exp(1j*deg2rad(data["Vs_phase"])))) # wrap back into -180 to 180
        ## Check for the magnetific field column
        try:
            data["B"] = data_raw[:,self.columns["B"]]  # try to see if it can load B
        except:
            data["B"] = np.ones_like(data["T0"]) * self.field
        if self.use_DC == True:
            ## DC ON
            data["Vp_DC_ON"] = data_raw[:,self.columns["Vp_DC"]] # DC voltage of the T+ thermocouple, HEAT ON
            data["Vm_DC_ON"] = data_raw[:,self.columns["Vm_DC"]] # DC voltage of the T- thermocouple, HEAT ON
            data["Vs_DC_ON"] = data_raw[:,self.columns["Vs_DC"]] # DC Seebeck voltage, HEAT ON
            ## DC OFF
            data["Vp_DC_OFF"] = data_raw[:,self.columns["Vp_DC_OFF"]] # DC voltage of the T+ thermocouple, HEAT ON
            data["Vm_DC_OFF"] = data_raw[:,self.columns["Vm_DC_OFF"]] # DC voltage of the T- thermocouple, HEAT ON
            data["Vs_DC_OFF"] = data_raw[:,self.columns["Vs_DC_OFF"]] # DC Seebeck voltage, HEAT ON
            ## Remove background
            data["Vp_DC"]  = np.abs(data["Vp_DC_ON"] - data["Vp_DC_OFF"])
            data["Vm_DC"]  = np.abs(data["Vm_DC_ON"] - data["Vm_DC_OFF"])
            data["Vs_DC"]  = self.sign_S * (data["Vs_DC_ON"] - data["Vs_DC_OFF"])
        return data

    def analysis(self, clean_up=True, clean_against="T0"):
        """
        Perform the analysis
        """
        if clean_up:
            self.clean_up(clean_against)
        self.calc_temp_dc()
        self.calc_temp_ac()
        self.calc_seebeck()
        self.create_file_name_analyzed()

    def clean_up(self, clean_against="T0"):
        """
        Clean the data by:
        - sorting the data in ascending order of temperature
        - removing the points with the same temperature values
        """
        X_clean = self.data[clean_against]
        ## Sort in respect to temperature
        index_sort = np.argsort(X_clean)
        X_clean = X_clean[index_sort]
        for key in self.data.keys():
            self.data[key] = self.data[key][index_sort]
        ## Unique
        X_clean, index_unique = np.unique(X_clean, return_index=True)
        for key in self.data.keys():
            self.data[key] = self.data[key][index_unique]
        ## Save
        self.data[clean_against] = X_clean


    def calc_temp_ac(self):
        ## AC >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        self.data["Tp_R"]   = self._ac_temp_thermocouples(self.data["T0"], self.data["Vp_R"], self.data["Tp_DC"])
        self.data["Tm_R"]   = self._ac_temp_thermocouples(self.data["T0"], self.data["Vm_R"], self.data["Tm_DC"])
        self.data["Tp_phase"] = self.data["Vp_phase"] # degrees
        self.data["Tm_phase"] = self.data["Vm_phase"] # degrees

        Tp_j = self.data["Tp_R"] * np.exp(1j * deg2rad(self.data["Tp_phase"]))
        Tm_j = self.data["Tm_R"] * np.exp(1j * deg2rad(self.data["Tm_phase"]))

        ## AC thermal gradient
        self.data["dT_AC"] = np.absolute(Tp_j - Tm_j)
        self.data["dT_phase"] = rad2deg(np.angle(Tp_j - Tm_j))

    def _ac_temp_thermocouples(self, T0, dV_AC, T_DC):
        if self.use_DC == False:
            T_DC = T0
        T_AC = Ether_inv(Ether(T_DC) + dV_AC) - T_DC
        return T_AC

    def calc_temp_dc(self):
        if self.use_DC == True:
            self.data["Tp_DC"] = self._dc_temp_thermocouples(self.data["T0"], self.data["Vp_DC"])
            self.data["Tm_DC"] = self._dc_temp_thermocouples(self.data["T0"], self.data["Vm_DC"])
            self.data["dT_DC"] = self.data["Tp_DC"] - self.data["Tm_DC"]
            self.data["Tav"]   = (self.data["Tp_DC"] + self.data["Tm_DC"]) / 2
        else:
            self.data["Tp_DC"] = self.data["T0"]
            self.data["Tm_DC"] = self.data["T0"]
        self.data["dT_DC"] = self.data["Tp_DC"] - self.data["Tm_DC"]
        self.data["Tav"]   = (self.data["Tp_DC"] + self.data["Tm_DC"]) / 2

    def _dc_temp_thermocouples(self, T0, dV_DC):
        T_DC = Ether_inv(dV_DC + Ether(T0))
        return T_DC

    def calc_seebeck(self):
        self.data["S_AC"] = self.data["Vs_R"] / self.data["dT_AC"] * np.cos(np.deg2rad(self.data["dT_phase"]-self.data["Vs_phase"]))
        if self.use_DC == True:
            self.data["S_DC"] = self.data["Vs_DC"] / self.data["dT_DC"]
        else:
            self.data["S_DC"] = np.zeros_like(self.data["dT_DC"])

    def create_file_name_analyzed(self):
        """
        Create the name of the analyzed file
        """
        self.file_analyzed = self.measurement_type + \
                             "_" + "B=" + str(self.field) + 'T' + \
                             "_" + self.samplelabel + \
                             "_" + str(self.date)

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



class SeebeckField(SeebeckTemp):
    def __init__(self, samplelabel, date, temp, measurement_type="FS_Seebeck"):
        super().__init__(samplelabel=samplelabel,
                        date=date,
                        field=np.nan,
                        measurement_type=measurement_type)
        self.measurement_type = measurement_type
        self.samplelabel = samplelabel
        self.date        = date
        self.temp        = temp

    ## Methods >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def load_info(self):
        """
        Load the measurement pickle for the corresponding samplelabel. It contains
        the info regarding the measurements. The keys of info are (Field, Date), and
        return for each couple of (Field,Date) a particular dictionnary for each
        measurements. This dictionnary has for keys : \n

        "folder" : the name of the folder where data are stored in "data_r" \n
        "comments" : comments regarding this particular set of points \n
        "columns" : columns in the raw data file \n
        """
        with open("info_fsweep_" + self.samplelabel + ".pkl", "rb") as f: # 'w' means (over)write or create if doesn't exist
            unpickler_file = pickle.Unpickler(f)
            self.info = unpickler_file.load()[self.temp, self.date] # keys of measurements are (H,date)


    def create_file_name_analyzed(self):
        """
        Create the name of the analyzed file
        """
        self.file_analyzed = self.measurement_type + \
                            "_" + "T=" + str(self.temp) + 'K' + \
                            "_" + self.samplelabel + \
                            "_" + str(self.date)

