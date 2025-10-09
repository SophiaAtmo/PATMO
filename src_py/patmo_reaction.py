import sys,inspect
import patmo_species
import patmo_string
import patmo_error
#import patmo_photoFlux

class reaction:

	#******************
	def __init__(self):
		self.reactants = []
		self.products = []
		self.rate = None
		self.verbatim = None
		self.fileSafeVerbatim = None
		self.RHS = None
		self.index = None
		self.photoIndex = None
		self.xsec = None
		self.xsecEnergy = None
		self.xsecEnergySpan = None
		self.hash = None
		self.hashReverse = None
		self.reverseOriginalIndex = None
		self.htmlTableRow = None
		self.htmlVerbatim = None
		self.catalyzer = None
		self.is3body = False
		self.isPseudo3body = False #3body written as a 2body
		self.TgasMin = None
		self.TgasMax = None
		self.evaluationKIDA = None

	#*****************
	#make an html page for this reaction
	def makeHtmlPage(self):
		fout = open("htmlDocs/rate"+str(int(1e6)+self.index)+".html","w")
		fout.write(patmo_string.readFile("htmlDocs/header.src")+"\n")
		fout.write(self.getHtmlVerbatim()+"<br><br>\n")
		fout.write("<img src=\"ratePNGs/rate"+str(int(1e6)+self.index)+".png\" width=\"500px\">\n")
		fout.write("<img src=\"rateReversePNGs/rate"+str(int(1e6)+self.index)+".png\" width=\"500px\">\n")
		fout.write(patmo_string.readFile("htmlDocs/footer.src")+"\n")
		fout.close()

	#******************
	#return true if species name is in the reaction
	def isContainingSpeciesName(self,speciesName):
		if(self.isContainingSpeciesNameReactants(speciesName)): return True
		if(self.isContainingSpeciesNameProducts(speciesName)): return True
		#if nothing found in the above functions returns false
		return False

	#******************
	#return true if species name is contained in reactants
	def isContainingSpeciesNameReactants(self,speciesName):
		#loop on reactants
		for react in self.reactants:
			if(speciesName==react.name): return True
		#if nothing found in the loop above returns false
		return False

	#******************
	#return true if species name is contained in products
	def isContainingSpeciesNameProducts(self,speciesName):
		#loop on products
		for react in self.products:
			if(speciesName==react.name): return True
		#if nothing found in the loop above returns false
		return False

	#******************
	#prepare an html table raw for this reaction
	def getHtmlTableRow(self,highlightSpecies=None):
		#set empty spaces
		rowReactants = ["" for i in range(3)]
		rowProducts = ["" for i in range(6)]
		#replace spaces with html names
		for i in range(len(self.reactants)):
			boldLeft = boldRight = ""
			if(highlightSpecies==self.reactants[i].name):
				boldLeft = "<b>"
				boldRight = "</b>"
			rowReactants[i] = boldLeft+self.reactants[i].getHtmlName()+boldRight
		for i in range(len(self.products)):
			boldLeft = boldRight = ""
			if(highlightSpecies==self.products[i].name):
				boldLeft = "<b>"
				boldRight = "</b>"
			rowProducts[i] = boldLeft+self.products[i].getHtmlName()+boldRight
		#if photochemistry second reactant is h\nu
		if(self.photoIndex!=None): rowReactants[1] = "<i>h&nu;</i>"
		#put the pieces together
		rowString = "<tr><td>"+str(self.index)+"<td>"+("<td>+<td>".join(rowReactants))+"<td>&#8594;<td>" \
			+ ("<td>+<td>".join(rowProducts)) \
			+ "<td><a href=\"rate"+str(int(1e6)+self.index)+".html\">link</a>"
		#loop while string has empty table cells with plus signs
		while("<td>+<td><td>" in rowString):
			rowString = rowString.replace("<td>+<td><td>","<td><td><td>")
		return rowString

	#********************
	#verbatim in html format
	def getHtmlVerbatim(self):
		if(self.htmlVerbatim!=None): return self.htmlVerbatim
		self.htmlVerbatim = (" + ".join([x.getHtmlName() for x in self.reactants])) \
			+" &#8594; "+(" + ".join([x.getHtmlName() for x in self.products]))
		return self.htmlVerbatim

	#******************
	#check mass conservation
	#def checkMassConservation(self):
#		massReactants = sum([x.mass for x in self.reactants])
#		massProducts = sum([x.mass for x in self.products])
#		electronMass = 9.10938356e-28 #g
#		if(abs(massReactants-massProducts)>electronMass/2e0):
#			print "ERROR: mass conservation problem!"
#			for x in (self.reactants+self.products):
#				print x.name, x.mass
#			print self.getVerbatim()
#			sys.exit()

	#******************
	#check charge conservation
	def checkChargeConservation(self):
		chargeReactants = sum([x.charge for x in self.reactants])
		chargeProducts = sum([x.charge for x in self.products])
		if(chargeReactants+chargeProducts!=0):
			print ("ERROR: charge conservation problem!")
			print (self.getVerbatim())
			sys.exit()

	#******************
	#get reaction hash (e.g. R|R||P|P|P format)
	def getHash(self):
		if(self.hash!=None): return self.hash
		reactants = sorted([x.name for x in self.reactants])
		products = sorted([x.name for x in self.products])
		self.hash = ("|".join(reactants))+"||"+("|".join(products))

		return self.hash

	#******************
	#get reverse reaction hash (e.g. P|P|P||R|R format)
	def getHashReverse(self):
		if(self.hashReverse!=None): return self.hashReverse
		ahash = self.getHash().split("||")
		self.hashReverse = ahash[1]+"||"+ahash[0]
		return self.hashReverse

	#******************
	#get the verbarim version of the reaction
	def getVerbatim(self):
		if(self.verbatim!=None): return self.verbatim
		self.verbatim = (" + ".join([x.name for x in self.reactants])) + " -> " \
			+ (" + ".join([x.name for x in self.products]))
		return self.verbatim

	#******************
	#get a file name safe version of the verbatim reaction (underscores)
	def getFileSafeVerbatim(self):
		if(self.fileSafeVerbatim!=None): return self.fileSafeVerbatim
		self.fileSafeVerbatim = ("_".join([x.name for x in self.reactants])) + "__" \
			+ ("_".join([x.name for x in self.products]))
		return self.fileSafeVerbatim

	#*****************
	#get reverse rate assuming loaded thermochemistry.
	# returns the f90 function with loop on cells and T limits
	def getReverse(self,thermochemistry):
		reverse = reaction()
		reverse.reactants = self.products
		reverse.products = self.reactants
		reverse.reverseOriginalIndex = self.index

		#check if species are in thermochemistry
		for species in (reverse.reactants+reverse.products):
			if(not(species.name in thermochemistry)):
				print ("ERROR: species "+species.name+" missing in thermochemistry data!")
				sys.exit()

		nCoef = 7 #number of coefficients expected
		upCoef = [0e0 for i in range(nCoef)]
		lowCoef = [0e0 for i in range(nCoef)]
		Tmin = 0e0 #minimum temperature
		Tmax = 1e8 #maximum temperature
		Tmid = 0e0 #temperature between lower-T coef an upper-T
		#note: reactants here are forward's products
		for species in reverse.reactants:
			#higher Tmin and lower Tmax are limits for all polynomials
			Tmin = max(Tmin,thermochemistry[species.name]["Tmin"])
			Tmax = min(Tmax,thermochemistry[species.name]["Tmax"])
			#first round store Tmid and Tmin
			if(Tmid==0e0): Tmid = thermochemistry[species.name]["Tmid"]
			if(Tmin==0e0): Tmin = thermochemistry[species.name]["Tmin"]
			if(Tmid!=thermochemistry[species.name]["Tmid"]):
				#check if mid value is different from the others
				print ("ERROR: Tmid of "+species.name+" is different from the other polynomials")
				print (self.getVerbatim())
				sys.exit()
			if(Tmin!=thermochemistry[species.name]["Tmin"]):
				#check if min value is different from the others
				print ("ERROR: Tmin of "+species.name+" is different from the other polynomials")
				print (self.getVerbatim())
				sys.exit()
			for i in range(nCoef):
				lowCoef[i] += thermochemistry[species.name]["lowCoef"][i]
				upCoef[i] += thermochemistry[species.name]["upCoef"][i]
		#note: products here are forward reactants
		for species in reverse.products:
			#higher Tmin and lower Tmax are limits for all polynomials
			Tmin = max(Tmin,thermochemistry[species.name]["Tmin"])
			Tmax = min(Tmax,thermochemistry[species.name]["Tmax"])
			#check if mid value is different from the others
			if(Tmid!=thermochemistry[species.name]["Tmid"]):
				print ("ERROR: Tmid of "+species.name+" is different from the other polynomials")
				print (self.getVerbatim())
				sys.exit()
			#check if min value is different from the others
			if(Tmin!=thermochemistry[species.name]["Tmin"]):
				print ("ERROR: Tmin of "+species.name+" is different from the other polynomials")
				print (self.getVerbatim())
				sys.exit()
			for i in range(nCoef):
				lowCoef[i] -= thermochemistry[species.name]["lowCoef"][i]
				upCoef[i] -= thermochemistry[species.name]["upCoef"][i]

		#store temperature ranges
		reverse.Tmin = Tmin
		reverse.Tmax = Tmax
		reverse.Tmid = Tmid

		#coefficient polynomial multiplier with sign (eqn.6 in Burcat+2005)
		multPoly = [1e0, 1e0/2e0, 1e0/6e0, 1e0/12e0, 1e0/20e0, -1e0, 1e0]
		#change sign because kR = kF/kC
		multPoly = [-x for x in multPoly]
		for i in range(nCoef):
			lowCoef[i] *= multPoly[i]
			upCoef[i] *= multPoly[i]

		#add Tgas dependence to polynomial
		polyTmult = ["*(lnTgas(i)-1d0)","*Tgas(i)","*Tgas2(i)","*Tgas3(i)","*Tgas4(i)","*invTgas(i)",""]
		polyLow = []
		polyUp = []
		#put Tgas dependence and coefficients together
		for i in range(nCoef):
			polyLow.append(patmo_string.f90DoubleCompact(lowCoef[i])+polyTmult[i])
			polyUp.append(patmo_string.f90DoubleCompact(upCoef[i])+polyTmult[i])
		#join polyomials, replace +- signs, add line returns for readibility
		polyLowJoint = (" + ".join(polyLow)).replace(" + -", " - ").replace(" + ", " &\n+ ").replace(" - ", " &\n- ")
		polyUpJoint = (" + ".join(polyUp)).replace(" + -", " - ").replace(" + ", " &\n+ ").replace(" - ", " &\n- ")

		diffFactor = ""
		#note here is (kb*T/P0)**(P-R) where P and R are referred to forward
		# so for the reverse is the opposite (R-P)
		if(len(reverse.reactants)!=len(reverse.products)):
			diffFactor = "*(1.3806488d-22*Tgas(i))**(" \
				+str(len(reverse.reactants)-len(reverse.products))+")"

		reverse.rateLow = "exp("+polyLowJoint+")"+diffFactor
		reverse.rateUp = "exp("+polyUpJoint+")"+diffFactor


		#reverse rate consists ina a loop on cells to find Tgas limits
		reverse.rate = "do i=1,cellsNumber\n"
		reverse.rate += " if(Tgas(i)<"+patmo_string.f90DoubleCompact(Tmid)\
			+".and.Tgas(i).ge."+patmo_string.f90DoubleCompact(Tmin)+") then\n"
		reverse.rate += "  #VARIABLE = krate(i,"+str(self.index)+")*"+reverse.rateLow+"\n"
		reverse.rate += " elseif(Tgas(i).ge."+patmo_string.f90DoubleCompact(Tmid)\
			+".and.Tgas(i).le."+patmo_string.f90DoubleCompact(Tmax)+") then\n"
		reverse.rate += "  #VARIABLE = krate(i,"+str(self.index)+")*"+reverse.rateUp+"\n"
		reverse.rate += " else\n"
		reverse.rate += "  #VARIABLE = 0d0\n"
		reverse.rate += " end if\n"

		#print if pseudo 3body reaction is found (i.e. 3body written as 3body)
		# and divide by ntot
		if(self.isPseudo3body):
			print ("pseudo 3body found: "+self.getVerbatim())
			reverse.rate += " !divided because pseudo-3body\n"
			reverse.rate += " #VARIABLE = #VARIABLE / ntot(i)\n"

		reverse.rate += "end do\n"

		return reverse

	#*******************
	def plotPhotorate(self):
		rate = 0e0
		planck_eV = 4.135667516e-15 #eV*s
		#flux = getBB()
		for ienergy in range(len(self.xsecEnergySpan)):
			rate += flux[ienergy]*self.xsecEnergySpan[ienergy]*self.xsec[ienergy]/self.xsecEnergy[ienergy]/planck_eV
			fout.write(str(self.xsecEnergy[ienergy])+" "+str(self.xsec[ienergy])+"\n")


	#******************
	#dump the cross section to a fileName conained in folder, with specific postpone string
	# if fileName not present generate from file name safe verbatim string
	def dumpXsec(self,folder,fileName=None,postpone=""):
		if(fileName==None):
			fileName = self.getFileSafeVerbatim()

		fout = open(patmo_string.pathFormat(folder)+fileName+postpone,"w")
		fout.write("#xsec (eV,cm2) for "+self.getVerbatim()+"\n")
		for ienergy in range(len(self.xsecEnergy)):
			fout.write(str(self.xsecEnergy[ienergy])+" "+str(self.xsec[ienergy])+"\n")
		fout.close()

	#******************
	#get the RHS of the reaction (F90)
	def getRHS(self):
		if(self.RHS!=None): return self.RHS
		self.RHS = "krate(:,"+str(self.index)+")"
		for reactant in self.reactants:
			self.RHS += "*n(:,"+reactant.label+")"
		return self.RHS

	#******************
	#parse split file line using split format
	def parse(self,rowList,formatList):
		#check if split line and format have the same length
		if(len(rowList)!=len(formatList)):
			print (formatList)
			print (rowList)
			sys.exit("ERROR: format list has different size than row list")
		#loop on split lines
		for i in range(len(rowList)):
			fmt = formatList[i] #format element
			part = rowList[i] #line element
			#skip empty parts
			if(part==""): continue
			#reactant
			if(fmt=="R"):
				mySpecies = patmo_species.species()
				mySpecies.parse(part)
				self.reactants.append(mySpecies)
			#product
			elif(fmt=="P"):
				mySpecies = patmo_species.species()
				mySpecies.parse(part)
				self.products.append(mySpecies)
			#skip index
			elif(fmt=="idx"):
				continue
			#rate in F90 format
			elif(fmt=="rate"):
				self.rate = part.replace(";",",")
				#guess if this is a 3body written as a 2body
				if("ntot(icell)" in self.rate): self.isPseudo3body = True
			else:
				sys.exit("ERROR: unknown format element "+fmt)

	#**********************
	#parse a KIDA reaction, given reactants, products, and other data
	def parseKIDA(self,reactants,products,data):
		#parse reactants
		for reactant in reactants:
			mySpecies = patmo_species.species()
			mySpecies.parse(reactant)
			self.reactants.append(mySpecies)

		#parse products
		for product in products:
			mySpecies = patmo_species.species()
			mySpecies.parse(product)
			self.products.append(mySpecies)

		#store formula type
		formulaType = int(data[9])
		#store temperature
		(TgasMin,TgasMax) = [float(x) for x in data[7:9]]
		#change limits if necessary
		if(TgasMin<0e0): TgasMin = 0e0
		if(TgasMax>9999.): TgasMax = 1e8
		#store limits
		self.TgasMin = TgasMin
		self.TgasMax = TgasMax
		#read coefficients
		(ka,kb,kc) = [float(x) for x in data[:3]]

		#evaluation from KIDA
		self.evaluationKIDA = int(data[12])

		#defulat rate is none to rise error later
		krate = None
		#different types of reactions according to
		# http://kida.obs.u-bordeaux1.fr/help
		if(formulaType==1):
			#CR
			pass
		elif(formulaType==2):
			#photo draine
			pass
		elif(formulaType==3):
			krate = patmo_string.f90DoubleCompact(ka)
			if(kb!=0): krate += "*(Tgas/3d2)**("+str(float(kb))+")"
			if(kc!=0): krate += "*exp(-"+patmo_string.f90DoubleCompact(kc)+"/Tgas)"
		#ionpol1, a*b*(0.62+0.4767*c*sqrt(300/T))
		elif(formulaType==4):
			k1 = patmo_string.f90DoubleCompact(0.62*ka*kb)
			k2 = patmo_string.f90DoubleCompact(0.4767*ka*kb*kc)
			krate = k1+" + "+k2+"*sqrt(3d2/Tgas)"
		#ionpol2, a*b*(1+c*0.0967*sqrt(300/T)+c*c*300/10.526/T)
		elif(formulaType==5):
			k1 = patmo_string.f90DoubleCompact(ka*kb)
			k2 = patmo_string.f90DoubleCompact(0.0967*ka*kb*kc)
			k3 = patmo_string.f90DoubleCompact(ka*kb*kc**2*300./10.526)
			krate = k1+" + "+k2+"*sqrt(3d2/Tgas) + "+k3+"/Tgas"

		#rise error if nothing found
		if(krate==None):
			print ("ERROR: unkown reaction type "+str(formulaType))
			print (reactants)
			print (products)
			print (data)
			sys.exit()

		#copy reaction rate found to object attribute, replace also double signs
		self.rate = krate.replace("--","+").replace("+-","-").replace("-+","-")

	#*******************
	def parseKIDA3b(self,reactants,products,data):

		#set 3body attribute
		self.is3body = True

		#parse reactants
		for reactant in reactants:
			mySpecies = patmo_species.species()
			mySpecies.parse(reactant)
			self.reactants.append(mySpecies)
			#look for catalyzer (reactant present in products)
			if(reactant in products):
				#if catalyzer is not empty when found rise error (multiple catalyzer!)
				if(self.catalyzer!=None):
					print ("ERROR: too many catalyzers in 3body reaction!")
					print (reactants,products)
					patmo_error.trigError(__file__,inspect.currentframe())
				#store catalyzer
				self.catalyzer = mySpecies

		#check if catalyzer found
		if(self.catalyzer==None):
			print ("ERROR: catalyzer not found in 3body reaction!")
			print (reactants,products)
			patmo_error.trigError(__file__,inspect.currentframe())

		#parse products
		for product in products:
			mySpecies = patmo_species.species()
			mySpecies.parse(product)
			self.products.append(mySpecies)

		fa = patmo_string.f90DoubleCompact(data["fa"])
		fb = patmo_string.f90DoubleCompact(data["fb"])
		fc = patmo_string.f90DoubleCompact(data["fc"])
		fd = patmo_string.f90DoubleCompact(data["fd"])
		ka_low = patmo_string.f90DoubleCompact(data["ka_low"])
		kb_low = str(data["kb_low"])
		kc_low = patmo_string.f90DoubleCompact(data["kc_low"])
		ka_inf = patmo_string.f90DoubleCompact(data["ka_inf"])
		kb_inf = str(data["kb_inf"])
		kc_inf = patmo_string.f90DoubleCompact(data["kc_inf"])
		D = "0.14d0"
		if(data["fa"]==0e0):
			Fc = "Fc(icell) = 1d0"
			C = "-0.4d0"
			N = "0.75d0"
		else:
			Fc = "Fc(icell) = (1d0-"+fa+")"
			if(data["fb"]<1e99): Fc += "*exp(-Tgas/"+fb+")"
			if(data["fc"]>1e-99): Fc += "*"+fa+"*exp(-Tgas/"+fc+")"
			if(data["fd"]<1e99): Fc += "*exp(-"+fd+"/Tgas)"
			C = "-0.4d0 - 0.67*log10(Fc(icell))"
			N = "0.75d0 - 1.27*log10(Fc(icell))"
		y = "y(icell) = log10(Pr(icell)) &\n + ("+C+")/(("+N+") &\n - ("+D+")*log10(Pr(icell)) &\n + ("+C+"))"
		klow = "klow(icell) = "+ka_low
		if(data["kb_low"]!=0): klow += "*(Tgas/3d2)**("+kb_low+")"
		if(data["kc_low"]!=0): klow += "*exp(-"+kc_low+"/Tgas)"
		kinf = "kinf(icell) = "+ka_inf
		if(data["kb_inf"]!=0): kinf += "*(Tgas/3d2)**("+kb_inf+")"
		if(data["kc_inf"]!=0): kinf += "*exp(-"+kc_inf+"/Tgas)"
		#pr = "Pr(icell) = klow(icell)*nAll(icell,"+self.catalyzer.label+")/kinf(icell)"
		#written assuming that [M] is in the ODE
		pr = "Pr(icell) = klow(icell)/kinf(icell)"

		if(data["Pmin"]==-9999 and data["Pmax"]==9999):
			limitsPressureBegin = limitsPressureEnd = ""
		elif(data["Pmin"]==-9999 and  data["Pmax"]!=9999):
			limitsPressureBegin = "if(Pr(icell).le."+patmo_string.f90DoubleCompact(data["Pmax"])+") then\n"
			limitsPressureEnd = "end if"
		elif(data["Pmin"]!=-9999 and  data["Pmax"]==9999):
			limitsPressureBegin = "if(Pr(icell).ge."+patmo_string.f90DoubleCompact(data["Pmin"])+") then\n"
			limitsPressureEnd = "end if"
		else:
			limitsPressureBegin = "if(Pr(icell).ge."+patmo_string.f90DoubleCompact(data["Pmin"]) \
				+ " .and. Pr(icell).le."+patmo_string.f90DoubleCompact(data["Pmax"])+") then\n"
			limitsPressureEnd = "end if"

		self.rate = klow+"\n"
		self.rate += kinf+"\n"
		self.rate += pr+"\n"
		self.rate += limitsPressureBegin
		self.rate += Fc+"\n"
		self.rate += y+"\n"
		self.rate += "F(icell) = Fc(icell)**(1d0/(1d0+y(icell)**2))\n"
		self.rate += "#VARIABLE = kinf(icell)*(Pr(icell)/(1d0+Pr(icell)))*F(icell)\n"
		self.rate += limitsPressureEnd


