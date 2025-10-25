import os
import sys

class options:
    # constructor
    def __init__(self, fname):
        # option defaults
        self.cellsNumber = 64
        self.cellThickness = 10000
        self.photoBinsNumber = 100
        self.network = ""
        self.energyMin = None
        self.energyMax = None
        self.wavelengMin = None
        self.wavelengMax = None
        self.zenith_angle = 60.0
        self.TOA_para = 0.5
        self.plotRates = True
        self.useEntropyProduction = False
        self.species = []
        self.usePhotochemistry = True
        self.useHescape = False
        self.useWaterRemoval = False
        self.useAerosolformation = False
        self.gravity_species = {}
        self.constant_species = []
        self.drydep_species = {}
        self.emission_species = {}

        # required casting for options
        integerType = ["cellsNumber", "cellThickness", "photoBinsNumber"]
        floatType = ["energyMin", "energyMax", "wavelengMin", "wavelengMax","zenith_angle", "TOA_para"]
        boolType = [
            "plotRates",
            "useEntropyProduction",
            "usePhotochemistry",
            "useHescape",
            "useWaterRemoval",
            "useAerosolformation",
        ]
        listType = ["species", "constant_species"]
        dictionaryType = ["emission_species", "drydep_species","gravity_species"]

        # check if file exists
        if not os.path.isfile(fname):
            print("ERROR: option file not found:", fname)
            sys.exit()

        print("reading option file", fname)

        # open and read option file
        with open(fname, "r") as fh:
            for row in fh:
                srow = row.strip()
                if not srow or srow.startswith("#"):
                    continue

                # Split only at the first '=' to avoid breaking values with '=' inside
                parts = [x.strip() for x in srow.split("=", 1) if x.strip()]
                if len(parts) != 2:
                    continue  # skip malformed lines

                option, value = parts

                if option in integerType:
                    try:
                        value = int(value)
                    except ValueError:
                        value = 0

                elif option in floatType:
                    try:
                        value = float(value)
                    except ValueError:
                        value = 0.0

                elif option in boolType:
                    value = value.upper() == "T"

                elif option in listType:
                    # Handle empty or missing values safely
                    if not value:
                        value = []
                    else:
                        # Filter out empty strings to avoid ['']
                        value = [x.strip() for x in value.split(",") if x.strip()]

                elif option in dictionaryType:
                    # Safely parse key:value pairs
                    if not value:
                        value = {}
                    else:
                        value_pairs = [x.strip() for x in value.split(",") if x.strip()]
                        parsed_dict = {}
                        for pair in value_pairs:
                            if ":" in pair:
                                k, v = pair.split(":", 1)
                                parsed_dict[k.strip()] = v.strip()
                        value = parsed_dict

                # Assign to class attribute
                setattr(self, option, value)
