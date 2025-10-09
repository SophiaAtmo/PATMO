import patmo_string

def buildCommons(network,options):

	allReactions = network.reactions + network.photoReactions + network.reverseReactions
	#dictionary for commons
	commonDict = []
	commonDict.append(["reactionsNumber", len(allReactions)])
	commonDict.append(["chemReactionsNumber", len(network.reactions)])
	commonDict.append(["photoReactionsNumber", len(network.photoReactions)])
	commonDict.append(["reverseReactionsNumber", len(network.reverseReactions)])
	commonDict.append(["chemSpeciesNumber", len(network.getSpecies())])
	commonDict.append(["speciesNumber", len(network.getSpecies())+2])
	commonDict.append(["positionTgas", len(network.getSpecies())+1])
	commonDict.append(["positionDummy", len(network.getSpecies())+2])
	commonDict.append(["cellsNumber", options.cellsNumber])
	commonDict.append(["photoBinsNumber", options.photoBinsNumber])
	
	#convert commons dictionary into f90 declarations
	allCommons = ""
	for variableName,value in commonDict:
		allCommons += "integer,parameter::"+variableName+" = "+str(value)+"\n"

	#commons for species index
	for species in network.getSpecies():
		allCommons += "integer,parameter::"+species.label+" = "+str(species.index)+"\n"

	photoPartners = []
	for reaction in network.photoReactions:
		photoPartners.append(reaction.reactants[0].label)
	photochemPartners = ""
	if(len(photoPartners)):
		photochemPartners = "integer,dimension(photoReactionsNumber)::photoPartnerIndex = (/"
		photochemPartners += (",".join(photoPartners))
		photochemPartners += "/)"

	indexArrays = ""
	indexDict = dict()
	maxReactants = network.getMaxReactants()
	maxProducts = network.getMaxProducts()
	for i in range(maxReactants):
		indexDict["indexReactants"+str(i+1)] = []
	for i in range(maxProducts):
		indexDict["indexProducts"+str(i+1)] = []
	for reaction in (allReactions):
		for i in range(maxReactants):
			if(i<len(reaction.reactants)):
				indexDict["indexReactants"+str(i+1)].append(reaction.reactants[i].label)
			else:
				indexDict["indexReactants"+str(i+1)].append("positionDummy")
		for i in range(maxProducts):
			if(i<len(reaction.products)):
				indexDict["indexProducts"+str(i+1)].append(reaction.products[i].label)
			else:
				indexDict["indexProducts"+str(i+1)].append("positionDummy")

	indexArrays = ""
	for varName,idxs in indexDict.items():
		indexArrays += "integer,parameter,dimension(reactionsNumber)::"+varName+" = (/"+(",&\n".join(idxs))+"/)\n"

	#replace commons pragma
	pragmaList = ["#PATMO_commons","#PATMO_photochemPartners","#PATMO_reaction_arrays"]
	replaceList = [allCommons,photochemPartners,indexArrays]

	patmo_string.fileReplaceBuild("src_f90/patmo_commons.f90", "build/patmo_commons.f90", \
		pragmaList, replaceList)


