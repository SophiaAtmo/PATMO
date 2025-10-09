import sys
import patmo_string

def buildReverseRates(network):

	Tmin = 0e0
	Tmax = 1e8
	reverseRates = ""
	for reverse in network.reverseReactions:
		rate = reverse.rate.replace("#VARIABLE","krate(i,"+str(reverse.index)+")")
		reverseRates += "!"+reverse.getVerbatim()+"\n"
		reverseRates += rate+"\n"
		Tmin = max(Tmin,reverse.Tmin)
		Tmax = min(Tmax,reverse.Tmax)

	Tmin = patmo_string.f90DoubleCompact(Tmin)
	Tmax = patmo_string.f90DoubleCompact(Tmax)

	pragmaList = ["#PATMO_reverseRates","#PATMO_Tmin","#PATMO_Tmax"]
	replaceList = [reverseRates,Tmin,Tmax]
	patmo_string.fileReplaceBuild("src_f90/patmo_reverseRates.f90", "build/patmo_reverseRates.f90", \
		pragmaList, replaceList)

