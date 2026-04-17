import numpy as np
import pickle
from sys import exit
import scipy.ndimage as im # For removing crazy points
from scipy.interpolate import interp1d
import os

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

class ResistivityT:
    def __init__(self, samplelabel, date, field, measurement_type="TS_rho_NSYM",subfolder = False):
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
        self.measurement_type = "TS_rho_NSYM" # is used in the analyzed file name
        self.subfolder = subfolder

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
        if self.subfolder == True:
            self.data_raw = np.loadtxt("../data_raw/" + self.folder_raw + "/"
                                   + self.file_raw,
                        dtype = "float", comments = "#", delimiter=delimiter)
        else:
            self.data_raw = np.loadtxt("../data_raw/" + self.folder_raw + "/"
                                   + self.file_raw,
                        dtype = "float", comments = "#", delimiter=delimiter)
        ## Extract raw data columns
        self.data = self.extract_raw_data(self.data_raw)


    def extract_raw_data(self, data_raw):
        data = {}
        data["T"] = data_raw[:,self.columns["T"]] # ref temperature
        data["Rxx"] = data_raw[:,self.columns["Rxx"]]  # resistance in Ohms
        data["Ixx"] = data_raw[:,self.columns["Ixx"]]  # current in Ohms
        try:
            data["B"] = data_raw[:,self.columns["B"]]  # try to see if it can load B
        except: pass
        return data

    def clean_up(self):
        """
        Clean the data by:
        - sorting the data in ascending order of temperature
        - removing the points with the same temperature values
        """
        T = self.data["T"]
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
        self.data["T"] = T

    def create_file_name_analyzed(self):
        """
        Create the name of the analyzed file
        """
        self.file_analyzed = self.samplelabel + "_" + self.measurement_type + \
                    "_B=" + str(self.field) + 'T_' + str(self.date)

    def calc_rhoxx(self):
        Rxx = self.data["Rxx"]
        rhoxx = Rxx * self.w * self.t / self.L # Ohm.m
        self.data["rhoxx"] = rhoxx

    def analysis(self, clean_up=True):
        """
        Perform the analysis
        """
        if clean_up:
            self.clean_up()
        self.calc_rhoxx()
        self.create_file_name_analyzed()
        return self.data["T"], self.data["rhoxx"], self.data["Rxx"]

    def save_file(self, to_save_list, folder='data_analyzed', delimiter=',', extension=".txt"):
        """
        Save data in the analyzed folder after analysis
        - list_data: ["T0", "Tp", etc.] what needs to be saved in self.data
        - measurement: "TS_rxx_NSYM" for example
        """
        ## Create the array to save data
        data_s = np.empty((self.data["T"].size, len(to_save_list)))
        for i, key in enumerate(to_save_list):
            data_s[:, i] = self.data[key]
        if self.subfolder == True:
            folder = '../' + folder
        ## Save the data file
        np.savetxt(folder + '/' + self.file_analyzed + extension, data_s, delimiter=delimiter,
               fmt='%.7e', header = ",".join(to_save_list), comments = "#")

    def _filter_crazy(self, key, n=50):
        ## Filter temperature points
        self.data[key] = im.filters.median_filter(self.data[key], size = n)
        # size correspond to
        # the number of neighboors it takes to compare the crazy points with
        # the good ones. Then it corrects if it is crazy

    ## Special Method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def __add__(self, obj_n):
        """
        When doing objxx_p + object_n = obj_sym.
        Create a new Resistivity instance with Rxx, Ixx and rhoxx symmetrized
        - obj_n is the Resistivity instance for negative fields
        """
        ## Creates a new object from the sum of the positive and negative field objects
        sym_obj = ResistivityT(self.samplelabel, self.date, np.abs(self.field))
        sym_obj.L = self.L
        sym_obj.w = self.w
        sym_obj.t = self.t
        ## Load the data
        T_p     = self.data["T"]
        Rxx_p   = self.data["Rxx"]
        Ixx_p   = self.data["Ixx"]
        T_n     = obj_n.data["T"]
        Rxx_n   = obj_n.data["Rxx"]
        Ixx_n   = obj_n.data["Ixx"]
        ## Make sure that the positive field is contained in the T range of the negative field
        index_range = (T_p > np.min(T_n)) * (T_p < np.max(T_n))
        T_p   = T_p[index_range]
        Rxx_p = Rxx_p[index_range]
        Ixx_p = Ixx_p[index_range]
        ## Interpolate the negative data on the positive field temperatures
        Rxx_n_i = np.interp(T_p, T_n, Rxx_n)
        Ixx_n_i = np.interp(T_p, T_n, Ixx_n)
        ## Symmetrize the data between positive and negative fields
        Rxx = (Rxx_p + Rxx_n_i) / 2
        Ixx = (Ixx_p + Ixx_n_i) / 2
        ## Fill up the symmetrized object
        sym_obj.data["T"]   = T_p
        sym_obj.data["Rxx"] = Rxx
        sym_obj.data["Ixx"] = Ixx
        ## Change the measurement type
        sym_obj.measurement_type = "TS_rho_SYM"
        return sym_obj

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



class HallEffectT(ResistivityT):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.measurement_type = "TS_RH_NSYM" # is used in the analyzed file name

        def extract_raw_data(self, data_raw):
            data = {}
            data["T"] = data_raw[:,self.columns["T"]] # ref temperature
            data["Rxy"] = data_raw[:,self.columns["Rxy"]]  # resistance in Ohms
            data["Ixy"] = data_raw[:,self.columns["Ixy"]]  # current in Ohms
            columns_phase = self.columns["Ixy"]+1
            data["phase_Rxy"] = data_raw[:, columns_phase]
            try:
                data["B"] = data_raw[:,self.columns["B"]]  # try to see if it can load B
            except: pass
            return data

        def calc_RH(self):
            Rxy = self.data["Rxy"] #* np.cos(self.data["phase_Rxy"])/abs(self.data["phase_Rxy"])
            rhoxy = Rxy * self.t # Ohm.m
            self.data["rhoxy"] = rhoxy
            self.data["RH"] = rhoxy / self.field

        def analysis(self, clean_up=True):
            """
            Perform the analysis
            """
            if clean_up:
                self.clean_up()
            self.calc_RH()
            self.create_file_name_analyzed()
            return self.data["T"], self.data["rhoxy"], self.data["RH"]


        ## Special Method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        def __sub__(self, obj_n):
            """
            When doing objxx_p - object_n = obj_sym.
            Create a new HalleffectT instance with Rxy, Ixy and rhoxy symmetrized
            - obj_n is the HalleffectT instance for negative fields
            """
            ## Creates a new object from the sum of the positive and negative field objects
            sym_obj = HallEffectT(self.samplelabel, self.date, np.abs(self.field))
            sym_obj.t = self.t
            ## Load the data
            T_p     = self.data["T"]
            Rxy_p   = self.data["Rxy"]
            Ixy_p   = self.data["Ixy"]
            T_n     = obj_n.data["T"]
            Rxy_n   = obj_n.data["Rxy"]
            Ixy_n   = obj_n.data["Ixy"]
            ## Make sure that the positive field is contained in the T range of the negative field
            index_range = (T_p > np.min(T_n)) * (T_p < np.max(T_n))
            T_p   = T_p[index_range]
            Rxy_p = Rxy_p[index_range]
            Ixy_p = Ixy_p[index_range]
            ## Interpolate the negative data on the positive field temperatures
            Rxy_n_i = np.interp(T_p, T_n, Rxy_n)
            Ixy_n_i = np.interp(T_p, T_n, Ixy_n)
            ## Symmetrize the data between positive and negative fields
            Rxx = (Rxy_p - Rxy_n_i) / 2
            Ixx = (Ixy_p - Ixy_n_i) / 2
            ## Fill up the symmetrized object
            sym_obj.data["T"]   = T_p
            sym_obj.data["Rxy"] = Rxx
            sym_obj.data["Ixy"] = Ixx
            ## Change the measurement type
            sym_obj.measurement_type = "TS_RH_SYM"
            return sym_obj


class ResistivityB(ResistivityT):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.measurement_type = "FS_MR_SYM" # is used in the analyzed file name

    def analysis_FS(self, clean_up=True):
        """
        Perform the analysis
        """
        if clean_up:
            self.clean_up()
        self.calc_rhoxx()
        self.create_file_name_analyzed()
        return self.data["B"], self.data["rhoxx"]

    def calc_MR(self, normalize=True):
        B = self.data["B"]
        Rxx = self.data["Rxx"]
        rhoxx = Rxx * self.w * self.t / self.L # Ohm.m
        B_mask = B>0
        B_plus = B[B_mask]
        B_minus = B[~B_mask]

        rhoxx_plus = rhoxx[B_mask]
        rhoxx_minus = rhoxx[~B_mask]

        min_val = -np.max(B_minus)
        max_val = -np.min(B_minus)
        mask = (B_plus >= min_val) & (B_plus <= max_val)
        B_plus_valid = B_plus[mask]
        rhoxx_plus = rhoxx_plus[mask]


        rhoxx_interp = interp1d(B_minus, rhoxx_minus)
        rhoxx_minus_interp = rhoxx_interp(-(B_plus_valid))


        rhoxx_sym = (rhoxx_plus+rhoxx_minus_interp)/2

        if normalize:
            rhoxx_sym = rhoxx_sym/rhoxx_sym[0]
        # B = np.concatenate((-B_plus_valid[::-1],B_plus_valid))
        # rhoxx_sym = np.concatenate((rhoxx_sym[::-1],rhoxx_sym))
        self.data["B"] = B_plus_valid
        self.data["rhoxx"] = rhoxx_sym
        self.data["rhoxx_raw"] = rhoxx
        self.data["B_raw"] = B



    def calc_Hall(self, normalize=True):
        B = self.data["B"]
        Rxy = self.data["Rxy"]
        rhoxy= Rxy *  self.t  # Ohm.m

        B_mask = B>0
        B_plus = B[B_mask]
        B_minus = B[~B_mask]

        rhoxy_plus = rhoxy[B_mask]
        rhoxy_minus = rhoxy[~B_mask]

        min_val = -np.max(B_minus)
        max_val = -np.min(B_minus)
        mask = (B_plus >= min_val) & (B_plus <= max_val)
        B_plus_valid = B_plus[mask]
        rhoxy_plus = rhoxy_plus[mask]


        rhoxy_interp = interp1d(B_minus, rhoxy_minus)
        rhoxy_minus_interp = rhoxy_interp(-(B_plus_valid))


        rhoxy_antisym = (rhoxy_plus-rhoxy_minus_interp)/2

        if normalize:
            rhoxy_sym = rhoxy_sym/rhoxy_sym[0]
        # B = np.concatenate((-B_plus_valid[::-1],B_plus_valid))
        # rhoxx_sym = np.concatenate((rhoxx_sym[::-1],rhoxx_sym))
        self.data["B"] = B_plus_valid
        self.data["rhoxy"] = rhoxy_sym
        self.data["RH"] = rhoxy_sym/B_plus_valid
        self.data["rhoxy_raw"] = rhoxy
        self.data["B_raw"]  = B





    def clean_up_MR(self):
        """
        Clean the data by:
        - sorting the data in ascending order of field
        - removing the points with the same field values
        """
        B = self.data["B"]
        ## Sort in respect to temperature
        index_sort = np.argsort(B)
        B = B[index_sort]
        for key in self.data.keys():
            self.data[key] = self.data[key][index_sort]
        ## Unique
        B, index_unique = np.unique(B, return_index=True)
        for key in self.data.keys():
            self.data[key] = self.data[key][index_unique]
        ## Save
        self.data["B"] = B



    def analysis_MR(self, clean_up = True, normalize=False):
        if clean_up:
            self.clean_up_MR()
        self.calc_MR(normalize=normalize)
        self.create_file_name_analyzed()
        return self.data["B"], self.data["rhoxx"]

    def analysis_Hall(self, clean_up = True, normalize = False):
        if clean_up:
            self.clean_up_MR()
        self.calc_Hall(normalize=normalize)
        self.create_file_name_analyzed()
        return self.data["B"], self.data["rhoxy"], self.data["RH"]


    def save_file(self, to_save_list, folder='../data_analyzed', delimiter=',', extension=".txt"):
        """
        Save data in the analyzed folder after analysis
        - list_data: ["T0", "Tp", etc.] what needs to be saved in self.data
        - measurement: "TS_rxx_NSYM" for example
        """
        ## Create the array to save data
        data_s = np.empty((self.data["B"].size, len(to_save_list)))
        for i, key in enumerate(to_save_list):
            data_s[:, i] = self.data[key]
        if self.subfolder == True:
            folder = '../' + folder
        ## Save the data file
        np.savetxt(folder + '/' + self.file_analyzed + extension, data_s, delimiter=delimiter,
               fmt='%.7e', header = ",".join(to_save_list), comments = "#")



class HallEffectB(ResistivityT):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.measurement_type = "TS_RH_NSYM" # is used in the analyzed file name

        def extract_raw_data(self, data_raw):
            data = {}
            data["T"] = data_raw[:,self.columns["T"]] # ref temperature
            data["Rxy"] = data_raw[:,self.columns["Rxy"]]  # resistance in Ohms
            data["Ixy"] = data_raw[:,self.columns["Ixy"]]  # current in Ohms
            column_phase = self.columns["Ixy"] + 1
            data["phase"] =  data_raw[:,column_phase]
            try:
                data["B"] = data_raw[:,self.columns["B"]]  # try to see if it can load B
            except: pass
            return data

        """
        def calc_RH(self):
            Rxy = self.data["Rxy"]
            rhoxy = Rxy * self.t # Ohm.m
            #print(rhoxy)

            self.data["rhoxy"] = rhoxy
            self.data["RH"] = rhoxy / self.field
        """

        def calc_RH(self, normalize=False):
            B = self.data["B"]
            Rxy = self.data["Rxy"] #*np.cos(self.data["phase"])/abs(np.cos(self.data["phase"]))
            rhoxy= Rxy  * self.t  # Ohm.m

            B_mask = B>0
            B_plus = B[B_mask]
            B_minus = B[~B_mask]

            rhoxy_plus = rhoxy[B_mask]
            rhoxy_minus = rhoxy[~B_mask]

            min_val = -np.max(B_minus)
            max_val = -np.min(B_minus)
            mask = (B_plus >= min_val) & (B_plus <= max_val)
            B_plus_valid = B_plus[mask]
            rhoxy_plus = rhoxy_plus[mask]


            rhoxy_interp = interp1d(B_minus, rhoxy_minus)
            rhoxy_minus_interp = rhoxy_interp(-(B_plus_valid))


            rhoxy_antisym = (rhoxy_plus-rhoxy_minus_interp)/2

            if normalize:
                rhoxy_sym = rhoxy_antisym/rhoxy_antisym[0]
            # B = np.concatenate((-B_plus_valid[::-1],B_plus_valid))
            # rhoxx_sym = np.concatenate((rhoxx_sym[::-1],rhoxx_sym))
            self.data["B"] = B_plus_valid
            self.data["rhoxy"] = rhoxy_antisym
            self.data["RH"] = rhoxy_antisym/B_plus_valid
            self.data["rhoxy_raw"] = rhoxy
            self.data["B_raw"]  = B


        def analysis_Hall(self, clean_up=True):
            """
            Perform the analysis
            """
            #print("Test")
            #if clean_up:
            #    self.clean_up()
            self.calc_RH()
            self.create_file_name_analyzed()
            return self.data["B"], self.data["rhoxy"], self.data["RH"]


        ## Special Method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        def __sub__(self, obj_n):
            """
            When doing objxx_p - object_n = obj_sym.
            Create a new HalleffectT instance with Rxy, Ixy and rhoxy symmetrized
            - obj_n is the HalleffectT instance for negative fields
            """
            ## Creates a new object from the sum of the positive and negative field objects
            sym_obj = HallEffectT(self.samplelabel, self.date, np.abs(self.field))
            sym_obj.t = self.t
            ## Load the data
            T_p     = self.data["T"]
            Rxy_p   = self.data["Rxy"]
            Ixy_p   = self.data["Ixy"]
            T_n     = obj_n.data["T"]
            Rxy_n   = obj_n.data["Rxy"]
            Ixy_n   = obj_n.data["Ixy"]
            ## Make sure that the positive field is contained in the T range of the negative field
            index_range = (T_p > np.min(T_n)) * (T_p < np.max(T_n))
            T_p   = T_p[index_range]
            Rxy_p = Rxy_p[index_range]
            Ixy_p = Ixy_p[index_range]
            ## Interpolate the negative data on the positive field temperatures
            Rxy_n_i = np.interp(T_p, T_n, Rxy_n)
            Ixy_n_i = np.interp(T_p, T_n, Ixy_n)
            ## Symmetrize the data between positive and negative fields
            Rxx = (Rxy_p - Rxy_n_i) / 2
            Ixx = (Ixy_p - Ixy_n_i) / 2
            ## Fill up the symmetrized object
            sym_obj.data["T"]   = T_p
            sym_obj.data["Rxy"] = Rxx
            sym_obj.data["Ixy"] = Ixx
            ## Change the measurement type
            sym_obj.measurement_type = "TS_RH_SYM"
            return sym_obj

