import shutil
import patmo_string

#*****************
def buildMain(network,options):

	#mass per nuclei functions
	allMassFunctions = ""
	for atomName in network.getAtoms():
		if(("(" in atomName) or ("_" in atomName)): continue
		functionName = "patmo_getTotalMassNuclei_"+atomName
		functionNameUtils = "getTotalMassNuclei_"+atomName
		allMassFunctions += "!***************************\n"
		allMassFunctions += "function "+functionName+"()\n"
		allMassFunctions += " use patmo_utils\n"
		allMassFunctions += " implicit none\n"
		allMassFunctions += " real*8::"+functionName+"\n\n"
		allMassFunctions += functionName + " = "+functionNameUtils+"() \n\n"
		allMassFunctions += "end function\n\n"

	#:::2017/9/11 addition:::
	allPhotoReactions = ""
	allJValues = "write(22,*) i"
	for reaction in network.photoReactions:
		allPhotoReactions += ", " + reaction.getVerbatim().replace(" ", "")
		allJValues += ", krate(i,"+str(reaction.index)+")"

	allReactions = ""
	allReactionRates = "write(22,*) i"
	for reaction in (network.reactions + network.photoReactions + network.reverseReactions):
		allReactions += ", " + reaction.getVerbatim().replace(" ", "")
		allReactionRates += ", &\n        " + reaction.getRHS().replace("n(", "nall(").replace(":", "i")

	#replace pragma
	pragmaList = ["#PATMO_massNucleiFunctions", "#PATMO_JValueReactions", "#PATMO_JValues", "#PATMO_DumpReactions", "#PATMO_DumpAllReactionRates"]
	replaceList = [allMassFunctions, allPhotoReactions, allJValues, allReactions, allReactionRates]

	#condition pragmas
	ifPragmas = ["#IFPATMO_use_opacity","#IFPATMO_usePhotochemistry", "#IFPATMO_useHescape", "#IFPATMO_useHescape_dump"]
	ifConditions = [options.usePhotochemistry, options.usePhotochemistry, options.useHescape, options.useHescape]

	patmo_string.fileReplaceBuild("src_f90/patmo.f90", "build/patmo.f90", \
		pragmaList, replaceList, ifPragmas, ifConditions)
