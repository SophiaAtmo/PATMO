import patmo_string
import patmo_options
import patmo_commons
import math

def buildPhotoRates(network,options):

	allRates = ""
	for reaction in network.photoReactions:
		allRates +=  "!"+reaction.getVerbatim()+"\n"
		allRates +=  "krate(:,"+str(reaction.index)+") = "+reaction.rate+"\n\n"
	
	dE = 0.0
	dE += (options.wavelengMax-options.wavelengMin)/options.photoBinsNumber
	res =""
	res += "dE = " +str(dE) 
	
	zenith = 0.0
	zenith += math.cos(math.radians(float(options.zenith_angle)))
	mu = f"mu = {zenith:.6f}"
	
	coef = 0.0
	coef += (float(options.TOA_para))
	coef =f"coef = {coef:.6f}"
	#replace commons pragma
	pragmaList = ["#PATMO_photoRates","#PATMO_resolution", "#PATMO_zenith_angle","#PATMO_TOA_para"]
	replaceList = [allRates,res,mu,coef]

	patmo_string.fileReplaceBuild("src_f90/patmo_photoRates.f90", "build/patmo_photoRates.f90", \
		pragmaList, replaceList)
