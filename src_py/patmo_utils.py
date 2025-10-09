import shutil
import patmo_string

#*****************
def buildUtils(network):
	allMass = allNames = ""
	#species mass and names
	for species in network.getSpecies():
		allMass += "getSpeciesMass("+species.label+") = " \
			+ patmo_string.f90DoubleCompact(species.mass) + "\n"
		allNames += "getSpeciesNames("+species.label+") = \"" \
			+ species.name + "\"\n"

	reactants = []
	weightedDegree = ""
	#add reactants array in flux
	for i in range(network.getMaxReactants()):
		reactants.append("n(indexReactants"+str(i+1)+"(i))")
		degreeVar = "degree(indexReactants"+str(i+1)+"(i))"
		weightedDegree += degreeVar +" = "+degreeVar+" + 1\n"
	for i in range(network.getMaxProducts()):
		degreeVar = "degree(indexProducts"+str(i+1)+"(i))"
		weightedDegree += degreeVar +" = "+degreeVar+" + 1\n"
	fluxReactants = (" &\n* ".join(reactants))

	#mass per nuclei functions
	allMassFunctions = ""
	for atomName in network.getAtoms():
		if(("(" in atomName) or ("_" in atomName)): continue
		functionName = "getTotalMassNuclei_"+atomName
		allMassFunctions += "!***************************\n"
		allMassFunctions += "function "+functionName+"()\n"
		allMassFunctions += " use patmo_commons\n"
		allMassFunctions += " use patmo_parameters\n"
		allMassFunctions += " implicit none\n"
		allMassFunctions += " integer::icell\n"
		allMassFunctions += " real*8::"+functionName+"\n"
		allMassFunctions += " real*8::m(speciesNumber)\n\n"
		allMassFunctions += " m(:) = getSpeciesMass()\n\n"
		allMassFunctions += " "+functionName+" = 0d0\n\n"
		allMassFunctions += "  do icell=1,cellsNumber\n"
		for species in network.getSpecies():
			if(not(atomName in species.atomDict)): continue
			mult = ""
			if(species.atomDict[atomName]>1): mult = " * "+patmo_string.f90DoubleCompact(species.atomDict[atomName])
			allMassFunctions += functionName + " = " + functionName + " + m(" +species.label+") * nall(icell," +species.label+")"+mult+"\n"
		allMassFunctions += "  end do\n\n"
		allMassFunctions += "end function\n\n"

	pragmaList = ["#PATMO_mass","#PATMO_speciesNames","#PATMO_flux","#PATMO_weightedDegree","#PATMO_massNucleiFunctions"]
	replaceList = [allMass,allNames,fluxReactants,weightedDegree,allMassFunctions]
	patmo_string.fileReplaceBuild("src_f90/patmo_utils.f90", "build/patmo_utils.f90", \
		pragmaList, replaceList)

#*****************
#copy a list of files from source to destination folder
def copyListTo(sourceFolder,destinationFolder,fileList):
	for fname in fileList:
		shutil.copyfile(patmo_string.pathFormat(sourceFolder)+fname, \
			patmo_string.pathFormat(destinationFolder)+fname)

#*****************
#copy a single file
def copyTo(sourceFile,destinationFile):
	shutil.copyfile(sourceFile, destinationFile)


#*****************
def isNumber(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
