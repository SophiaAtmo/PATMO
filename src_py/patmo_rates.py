import patmo_string

def buildRates(network):

	networkHas3body = False
	allRates = ""
	for reaction in network.reactions:
		allRates +=  "!"+reaction.getVerbatim()+"\n"
		rate = "krate(icell,"+str(reaction.index)+") = " + reaction.rate
		if("#VARIABLE" in reaction.rate):
			rate = reaction.rate.replace("#VARIABLE","krate(icell,"+str(reaction.index)+")")
			networkHas3body = True
		allRates +=  rate + "\n\n"

	#replace commons pragma
	pragmaList = ["#PATMO_rates"]
	replaceList = [allRates]

	#condition pragmas
	ifPragmas = ["#IFPATMO_has3body"]
	ifConditions = [networkHas3body]

	patmo_string.fileReplaceBuild("src_f90/patmo_rates.f90", "build/patmo_rates.f90", \
		pragmaList, replaceList, ifPragmas, ifConditions)
