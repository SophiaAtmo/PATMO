import os
import sys


class options:

    # *******************
    # constructor
    def __init__(self, fname):
        # option defaults
        self.cellsNumber = 64
        self.photoBinsNumber = 100
        self.network = ""
        self.energyMin = None
        self.energyMax = None
        self.plotRates = True
        self.doTopology = True
        self.plotThermochemistry = True
        self.useEntropyProduction = False
        self.species = []
        self.usePhotochemistry = True

        # required casting for options
        integer_type = ["cellsNumber", "photoBinsNumber"]
        float_type = ["energyMin", "energyMax"]
        bool_type = ["plotRates", "useEntropyProduction",
                     "usePhotochemistry", "useReverse", "doTopology", "plotThermochemistry"]
        list_type = ["species"]

        # check if file exists
        if not os.path.isfile(fname):
            print "ERROR: option file not found: " + fname
            sys.exit()

        print "reading option file "+fname
        # open option file
        fh = open(fname, "rb")
        for row in fh:
            srow = row.strip()
            if srow == "":
                continue
            if srow.startswith("#"):
                continue
            (option, value) = [x.strip() for x in srow.split("=") if x != ""]
            if option in integer_type:
                value = int(value)
            if option in float_type:
                value = float(value)
            if option in bool_type:
                value = (value.upper() == "T")
            if option in list_type:
                value = [x.strip() for x in value.split(",")]
            # convert option string into variable
            setattr(self, option, value)
        fh.close()
