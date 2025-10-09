import patmo_string

def buildJacobian(network):

	#get number max number of reactants and products
	maxReactants = network.getMaxReactants()
	maxProducts = network.getMaxProducts()

	jacobian = ""
	#write the flux for chemistry, e.g. ddotn/dn2 = k*n1
	for i in range(maxReactants):
		jacobian += "flux"+str(i+1)+"(j) = krate(j,i)"
		for j in range(maxReactants):
			if(i==j): continue
			jacobian += " &\n * n(j,indexReactants"+str(j+1)+"(i))"
		jacobian += "\n"

	#write d(Rdot)/dR and d(Pdot)/dR
	#Jacobian for chemistry (loop on reactants, i.e. d/dR)
	for i in range(maxReactants):
		iRdot = "(indexReactants"+str(i+1)+"(i)-1)*cellsNumber+j"
		#loop on reactants i.e. Rdot
		for j in range(maxReactants):
			iR = "(indexReactants"+str(j+1)+"(i)-1)*cellsNumber+j"
			RRjac = "pd_vec("+iR+", &\n"+iRdot+")"
			jacobian += RRjac + " = &\n"+ RRjac + " - flux"+str(i+1)+"(j)\n"
		jacobian += "\n"
		#loop products i.e. Pdot
		for j in range(maxProducts):
			iP = "(indexProducts"+str(j+1)+"(i)-1)*cellsNumber+j"
			RPjac = "pd_vec("+iP+", &\n"+iRdot+")"
			jacobian += RPjac + " = &\n"+ RPjac + " + flux"+str(i+1)+"(j)\n"
		jacobian += "\n"

	#replace pragma
	pragmaList = ["#PATMO_jacobian"]
	replaceList = [jacobian]

	patmo_string.fileReplaceBuild("src_f90/patmo_jacobian.f90", "build/patmo_jacobian.f90", \
		pragmaList, replaceList)


