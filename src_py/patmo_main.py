import patmo_string


# *****************
def buildMain(network, options):
    # mass per nuclei functions
    allMassFunctions = ""
    for atomName in network.getAtoms():
        if ("(" in atomName) or ("_" in atomName):
            continue
        functionName = "patmo_getTotalMassNuclei_" + atomName
        functionNameUtils = "getTotalMassNuclei_" + atomName
        allMassFunctions += "!***************************\n"
        allMassFunctions += "function " + functionName + "()\n"
        allMassFunctions += " use patmo_utils\n"
        allMassFunctions += " implicit none\n"
        allMassFunctions += " real*8::" + functionName + "\n\n"
        allMassFunctions += functionName + " = " + functionNameUtils + "() \n\n"
        allMassFunctions += "end function\n\n"

    # replace pragma
    pragmaList = ["#PATMO_massNucleiFunctions"]
    replaceList = [allMassFunctions]

    # condition pragmas
    ifPragmas = ["#IFPATMO_use_opacity", "#IFPATMO_usePhotochemistry"]
    ifConditions = [options.usePhotochemistry, options.usePhotochemistry]

    patmo_string.fileReplaceBuild("src_f90/patmo.f90", "build/patmo.f90",
                                  pragmaList, replaceList, ifPragmas, ifConditions)
