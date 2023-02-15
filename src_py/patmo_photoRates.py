import patmo_string


def buildPhotoRates(network):
    allRates = ""
    for reaction in network.photoReactions:
        allRates += "!" + reaction.getVerbatim() + "\n"
        allRates += "krate(:," + str(reaction.index) + ") = " + reaction.rate + "\n\n"

    # replace commons pragma
    pragmaList = ["#PATMO_photoRates"]
    replaceList = [allRates]

    patmo_string.fileReplaceBuild("src_f90/patmo_photoRates.f90", "build/patmo_photoRates.f90",
                                  pragmaList, replaceList)
