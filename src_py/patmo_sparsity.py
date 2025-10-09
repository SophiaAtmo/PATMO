import patmo_string

def buildSparsity(network):

	maxReactants = network.getMaxReactants()
	maxProducts = network.getMaxProducts()

	sparsity = ""
	for i in range(maxReactants):
		idxReactant = "(indexReactants"+str(i+1)+"(j)-1)*cellsNumber+i"
		for j in range(maxReactants):
			idxR = "(indexReactants"+str(j+1)+"(j)-1)*cellsNumber+i"
			sparsity += "Ms("+idxR+", &\n"+idxReactant+") = 1\n"
		sparsity += "\n"
		for j in range(maxProducts):
			idxP = "(indexProducts"+str(j+1)+"(j)-1)*cellsNumber+1"
			sparsity += "Ms("+idxP+", &\n"+idxReactant+") = 1\n"
		sparsity += "\n"

	#replace pragma
	pragmaList = ["#PATMO_sparsity"]
	replaceList = [sparsity]

	patmo_string.fileReplaceBuild("src_f90/patmo_sparsity.f90", "build/patmo_sparsity.f90", \
		pragmaList, replaceList)
