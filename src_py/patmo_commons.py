import patmo_string


def buildCommons(network, options):
    allReactions = network.reactions + network.photoReactions + network.reverseReactions
    # dictionary for commons
    common_dict = []
    common_dict.append(["reactionsNumber", len(allReactions)])
    common_dict.append(["chemReactionsNumber", len(network.reactions)])
    common_dict.append(["photoReactionsNumber", len(network.photoReactions)])
    common_dict.append(["reverseReactionsNumber", len(network.reverseReactions)])
    common_dict.append(["chemSpeciesNumber", len(network.getSpecies())])
    common_dict.append(["speciesNumber", len(network.getSpecies()) + 2])
    common_dict.append(["positionTgas", len(network.getSpecies()) + 1])
    common_dict.append(["positionDummy", len(network.getSpecies()) + 2])
    common_dict.append(["cellsNumber", options.cellsNumber])
    common_dict.append(["photoBinsNumber", options.photoBinsNumber])

    # convert commons dictionary into f90 declarations
    all_commons = ""
    for variableName, value in common_dict:
        all_commons += "integer,parameter::" + variableName + " = " + str(value) + "\n"

    # commons for species index
    for species in network.getSpecies():
        all_commons += "integer,parameter::" + species.label + " = " + str(species.index) + "\n"

    photoPartners = []
    for reaction in network.photoReactions:
        photoPartners.append(reaction.reactants[0].label)
    photochemPartners = ""
    if len(photoPartners):
        photochemPartners = "integer,dimension(photoReactionsNumber)::photoPartnerIndex = (/"
        photochemPartners += (",".join(photoPartners))
        photochemPartners += "/)"

    indexArrays = ""
    indexDict = dict()
    maxReactants = network.getMaxReactants()
    maxProducts = network.getMaxProducts()
    for i in range(maxReactants):
        indexDict["indexReactants" + str(i + 1)] = []
    for i in range(maxProducts):
        indexDict["indexProducts" + str(i + 1)] = []
    for reaction in allReactions:
        for i in range(maxReactants):
            if i < len(reaction.reactants):
                indexDict["indexReactants" + str(i + 1)].append(reaction.reactants[i].label)
            else:
                indexDict["indexReactants" + str(i + 1)].append("positionDummy")
        for i in range(maxProducts):
            if i < len(reaction.products):
                indexDict["indexProducts" + str(i + 1)].append(reaction.products[i].label)
            else:
                indexDict["indexProducts" + str(i + 1)].append("positionDummy")

    indexArrays = ""
    for varName, idxs in indexDict.iteritems():
        indexArrays += "integer,parameter,dimension(reactionsNumber)::" + varName + " = (/" + (
            ",&\n".join(idxs)) + "/)\n"

    # replace commons pragma
    pragmaList = ["#PATMO_commons", "#PATMO_photochemPartners", "#PATMO_reaction_arrays"]
    replaceList = [all_commons, photochemPartners, indexArrays]

    patmo_string.fileReplaceBuild("src_f90/patmo_commons.f90", "build/patmo_commons.f90",
                                  pragmaList, replaceList)
