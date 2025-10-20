from math import log10,exp,log,sqrt
from subprocess import Popen, PIPE
import os,sys,shutil,inspect
import matplotlib.pyplot as plt
import patmo_reaction
import patmo_string
import patmo_species
import patmo_error

class network:


	#****************
	def __init__(self,options):
		self.reactions = []
		self.photoReactions = []
		self.reverseReactions = []
		self.species = None
		self.speciesNames = None
		self.maxReactants = None
		self.maxProducts = None
		self.ghosts = []
		self.thermochemistry = dict()
		self.fileThermochemistry = "data/thermochemistry/thermochemistry.dat"
		self.loadThermochemistry()
		self.options = options

		self.loadKIDA(options.species)

	#****************
	#load a network from file passed via options object
	def loadNetwork(self):
		fname = self.options.network
		#if network is not set skip loading network from file
		if(fname==""):
			#if also empty species list rise error
			if(self.options.species==[]):
				print ("ERROR: in option file")
				sys.exit()
			return
		fh = open(fname,"r")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			if(srow.startswith("@format:")):
				aformat = [x.strip() for x in srow.replace("@format:","").split(",") if x!=""]
				continue
			if(srow.startswith("@ghost:")):
				ghosts = [x.strip() for x in srow.replace("@ghost:","").split(",") if x!=""]
				for ghost in ghosts:
					mySpecies = patmo_species.species()
					mySpecies.parse(ghost)
					self.ghosts.append(mySpecies)
				continue

			arow = [x.strip() for x in srow.split(",")]
			myReaction = patmo_reaction.reaction()
			myReaction.parse(arow,aformat)
			myReaction.index = len(self.reactions) + 1
			self.reactions.append(myReaction)

	#*******************
	#load KIDA database from a list of species (e.g. H,H2,CO,...)
	def loadKIDA(self,speciesReferenceList):
		self.loadKIDA2body(speciesReferenceList)
		self.loadKIDA3body(speciesReferenceList)

	#**********************
	def loadKIDA3body(self,speciesReferenceList):
		kidaFile3body = "data/kida/kida_reac_tb_2015-03-19_1.csv"
		fh = open(kidaFile3body,"r")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			#read KIDA file row (even if F90-like format can be split by spaces)
			arow = [x.strip() for x in srow.split(";")]
			reactants = arow[:3]
			products = [x for x in arow[3:8] if x!=""]

			#check if reaction has to be skipped (non-required species are present)
			skipThisReaction = False
			#loop on found species
			for species in (reactants+products):
				#if species is not in reaction discard this
				if(not(species in speciesReferenceList)):
					skipThisReaction = True
					break
			#skip reaction if required by condition above
			if(skipThisReaction): continue


			data = dict()

			(data["ka_low"], data["kb_low"], data["kc_low"]) = \
				[float(x) for x in arow[9:12]]
			data["formula_klow"] = int(arow[15])
			(data["ka_inf"], data["kb_inf"], data["kc_inf"]) = \
				[float(x) for x in arow[16:19]]
			data["formula_kinf"] = int(arow[22])
			(data["fa"], data["fb"], data["fc"], data["fd"]) = \
				[float(x) for x in arow[23:27]]
			(data["Tmin"],data["Tmax"]) = [float(x) for x in arow[30:32]]
			(data["Pmin"],data["Pmax"]) = [float(x) for x in arow[32:34]]
			recommend = int(arow[37])

			#skip not recommended
			if(recommend==0): continue

			#prepare reaction object by parsing data
			myReaction = patmo_reaction.reaction()
			myReaction.parseKIDA3b(reactants,products,data)
			myReaction.index = len(self.reactions) + 1
			self.reactions.append(myReaction)
			print ("found in KIDA (3body): " + myReaction.getVerbatim())

		fh.close()

	#**********************
	def loadKIDA2body(self,speciesReferenceList):
		#kida file
		kidaFile = "data/kida/kida.uva.2014.dat"
		fh = open(kidaFile,"r")

		reactionDictionary = dict()
		#loop on kida database to extract reactions mathing criteria
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			#read KIDA file row (even if F90-like format can be split by spaces)
			arow = [x.strip() for x in srow.split(" ") if x.strip()!=""]
			#first two parts are reactants (no 3-body in this file)
			reactants = arow[:2]
			#loop assuming that species are until the first number found in the row
			for i in range(len(arow)):
				iProdMax = i #store index
				#break when number found
				if(patmo_string.isNumber(arow[i])): break
			#products are from the 3rd part to the last non-numeric
			products = arow[2:iProdMax]

			#check if reaction has to be skipped (non-required species are present)
			skipThisReaction = False
			#loop on found species
			for species in (reactants+products):
				#if species is not in reaction discard this
				if(not(species in speciesReferenceList)):
					skipThisReaction = True
					break
			#skip reaction if required by condition above
			if(skipThisReaction): continue
			#after reactants and products everything is "data"
			data = arow[iProdMax:]
			#prepare reaction object by parsing data
			myReaction = patmo_reaction.reaction()
			myReaction.parseKIDA(reactants,products,data)
			myReaction.index = len(self.reactions) + 1

			myHash = myReaction.getHash()
			#check if the reaction is already present in the dictionary
			if(myHash in reactionDictionary):
				#check for contiguity (note that in the
				contiguousA = (abs(reactionDictionary[myHash].TgasMin-myReaction.TgasMax)<=1e0)
				contiguousB = (abs(reactionDictionary[myHash].TgasMax-myReaction.TgasMin)<=1e0)
				#if contiguity found write rate with if statements
				if(contiguousA or contiguousB):
					if(contiguousA):
						rateIf = "if(Tgas<"+patmo_string.f90DoubleCompact(myReaction.TgasMax)+") then\n"
					if(contiguousB):
						rateIf = "if(Tgas>"+patmo_string.f90DoubleCompact(myReaction.TgasMin)+") then\n"
					rateIf += "  #VARIABLE = " + myReaction.rate+"\n"
					rateIf += "else\n"
					rateIf += "  #VARIABLE = " + reactionDictionary[myHash].rate+"\n"
					rateIf += "end if\n"
					myReaction.rate = rateIf
				else:
					#error if not contiguos
					print ("ERROR: multiple reactions with non-contiguous temperature limits!") 
					print (myReaction.getVerbatim())
					patmo_error.trigError(__file__,inspect.currentframe())


			reactionDictionary[myHash] = myReaction
			print ("found in KIDA: " + myReaction.getVerbatim())

		#add found reactions to the attribute list
		for (rhash,myReaction) in reactionDictionary.items():
			myReaction.index = len(self.reactions)+1
			self.reactions.append(myReaction)

		fh.close()


	#*******************
	#check charge, mass, sources, and sinks for the current network
	def checkAll(self):

		#check if at least one reaction has been loaded
		if(len(self.reactions+self.photoReactions+self.reverseReactions)==0):
			print ("ERROR: it seems you didn't load any reaction!")
			patmo_error.trigError(__file__,inspect.currentframe())

		#loop on reaction to check mass and charge conservation
		for reaction in (self.reactions+self.photoReactions+self.reverseReactions):
			reaction.checkChargeConservation()
			#reaction.checkMassConservation()

		#loop on reaction to store hashes to check if reactions are unique
		# Note: no photochemical and reverse reactions included
		hashes = []
		for reaction in (self.reactions):
			hashes.append(reaction.getHash())
		hashes = sorted(hashes)

		rhashOld = None
		multipleReactionHash = []
		#check for multiple reactions (fast method: contiguity in sorted list)
		for rhash in hashes:
			if(rhash==rhashOld):
				multipleReactionHash.append(rhash)
			rhashOld = rhash

		#if multiple reactions found trig error
		if(len(multipleReactionHash)>0):
			print ("ERROR: multiple reaction(s) found:")
			#loop on uniqe list of reaction hashes
			for rhash in list(set(multipleReactionHash)):
				#get verbatim from hash
				verbatim = patmo_string.hash2verbatim(rhash)
				print (" " + verbatim +"\t x "+ str(multipleReactionHash.count(rhash)+1))
			patmo_error.trigError(__file__,inspect.currentframe())


		#get dictionary for sinks and sources
		sList = self.getSinkSourceList()

		#if sink found write it with corresponding reactions
		if(len(sList["sink"])>0):
			print ("WARNING: sinks found:")
			#loop on sinks
			for sink in sList["sink"]:
				print (" "+sink)
				#loop on reactions
				for reaction in (self.reactions+self.photoReactions):
					#check if sink in reactants or products
					if(sink in [x.name for x in (reaction.reactants+reaction.products)]):
						print ("  "+reaction.getVerbatim())

		#if source found write it with corresponding reactions
		if(len(sList["source"])>0):
			print ("WARNING: sources found:")
			#loop on sources
			for source in sList["source"]:
				print (" "+source)
				#loop on reactions
				for reaction in (self.reactions+self.photoReactions):
					#check if source in reactants or products
					if(source in [x.name for x in (reaction.reactants+reaction.products)]):
						print ("  "+reaction.getVerbatim())

	#*****************
	#get the list of sink and sources
	def getSinkSourceList(self):
		allReactants = []
		allProducts = []
		#loop on reactions (note: exclude reverse!) to store all
		# reactants and products
		for reaction in (self.reactions+self.photoReactions):
			allReactants += [x.name for x in reaction.reactants]
			allProducts += [x.name for x in reaction.products]
		#uniqe list of reactants and products
		allReactants = list(set(allReactants))
		allProducts = list(set(allProducts))

		#prepare a dictionary to store sinks and sources
		sList = {"source":[], "sink":[]}
		#loop on reactants
		for reactant in allReactants:
			#if reactant not present in products store it
			if(not(reactant in allProducts)): sList["source"].append(reactant)
		#loop on products
		for product in allProducts:
			#if product not present in reactants store it
			if(not(product in allReactants)): sList["sink"].append(product)
		return sList

	#*******************
	#load thermochemistry from a Burcat's-like table
	# http://garfield.chem.elte.hu/Burcat/burcat.html
	def loadThermochemistry(self):

		#file format is from Burcat's paper: "The 7 term polynomials actually
		# include 15 constants. The first set of 7 constants belongs to the
		# Tmid-Tmax K polynomial, the second set of 7 constants belongs to
		# the Tmin-Tmid K polynomial, and the fifteenth constant is
		# H298/R = DeltafH298/R. The latter value (and the corresponding
		# position within the polynomial format) is not used by
		# most other programs."
		fh = open(self.fileThermochemistry,"r")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			#split lines at spaces (for header)
			arow = [x.strip() for x in srow.split(" ") if x!=""]
			#split lines with fix format (for data)
			frow = [row[i*15:(i+1)*15] for i in range(5)]
			#read header
			if(arow[-1]=="1"):
				data = dict()
				speciesName = arow[0]
				data["Tmid"] = float(arow[-2])
				data["Tmax"] = float(arow[-3])
				data["Tmin"] = float(arow[-4])
			#read other lines
			if(arow[-1]=="2"):
				upCoef = [float(x) for x in frow]
			if(arow[-1]=="3"):
				coef = [float(x) for x in frow]
				upCoef += coef[:2]
				lowCoef = coef[2:]
			if(arow[-1]=="4"):
				coef = [float(x) for x in frow]
				lowCoef += coef[:-1]
				#store data
				data["lowCoef"] = lowCoef
				data["upCoef"] = upCoef
				#store data in dictionary
				self.thermochemistry[speciesName] = data
		fh.close()

	#*********************
	#get the maximum number of reactants
	def getMaxReactants(self):
		if(self.maxReactants!=None): return self.maxReactants
		self.maxReactants = 0
		for reaction in (self.reactions+self.photoReactions):
			self.maxReactants = max(self.maxReactants,len(reaction.reactants))
		return self.maxReactants

	#*********************
	#get the maximum number of products
	def getMaxProducts(self):
		if(self.maxProducts!=None): return self.maxProducts
		self.maxProducts = 0
		for reaction in (self.reactions+self.photoReactions):
			self.maxProducts = max(self.maxProducts,len(reaction.products))
		return self.maxProducts

	#*********************
	#get all the species (as objects)
	def getSpecies(self):
		if(self.species!=None): return self.species
		specDict = dict()
		for reaction in self.reactions:
			for reactant in reaction.reactants:
				specDict[reactant.name] = reactant
			for product in reaction.products:
				specDict[product.name] = product
		self.species = []
		for name,species in specDict.items():
			species.index = len(self.species) + 1
			self.species.append(species)
		#add ghosts to species
		for ghost in self.ghosts:
			ghost.index = len(self.species) + 1
			self.species.append(ghost)
		return self.species

	#*********************
	#get all the species (as names)
	def getSpeciesNames(self):
		if(self.speciesNames!=None): return self.speciesNames
		self.speciesNames = []
		for species in self.getSpecies():
			self.speciesNames.append(species.name)
		return self.speciesNames

	#************************
	def getAtoms(self):
		allAtoms = []
		for species in self.getSpecies():
			allAtoms += species.atomDict.keys()
		return list(set(allAtoms))

	#*********************
	#find needed reverse (e.g. exclude phtochemistry and rev already in the ntw)
	def reverseCheck(self):
		reactionsData = []
		hashReference = []
		#loop on reactions to find hash and reverse hash
		for reaction in self.reactions:
			r = dict()
			r["hash"] = reaction.getHash()
			r["hashRev"] = reaction.getHashReverse()
			r["reaction"] = reaction
			reactionsData.append(r)
			#store hash as reference
			hashReference.append(reaction.getHash())

		reverseNeeded = []
		#loop on reactions to find if reverse is already in the forward list
		for rdata in reactionsData:
			#skip 3body
			#if(len(r["reaction"].products)>3): continue
			#check if reaction is present in the reference (if not do reverse)
			if(not(rdata["hashRev"] in hashReference) or not(self.options.useEntropyProduction)):
				reverseNeeded.append(rdata["reaction"])
			else:
				print ("ERROR: reverse already present:",rdata["reaction"].getVerbatim())
				print (" At this stage of the developement you need to remove it from the network!")
				patmo_error.trigError(__file__,inspect.currentframe())
		return reverseNeeded

	#***********************
	def doReverse(self):
		reverseNeeded = self.reverseCheck()
		for reaction in reverseNeeded:
			reverse = reaction.getReverse(self.thermochemistry)
			reverse.index = len(self.reactions+self.photoReactions+self.reverseReactions)+1
			self.reverseReactions.append(reverse)

	#*********************
	#prepare html docs
	def makeHtmlDocs(self):
		shutil.copyfile("htmlDocs/index.src","htmlDocs/index.html")
		self.makeHtmlSpeciesList()
		self.makeHtmlReactionsList()
		self.createTopology()
		self.plotRates()
		self.plotReverseRates()

	#*********************
	#prepare html species list
	def makeHtmlSpeciesList(self):
		fout = open("htmlDocs/speciesMenu.html","w")
		fout.write(patmo_string.readFile("htmlDocs/header.src")+"\n")
		fout.write("<table width=\"40%\">\n")
		speciesList = self.getSpecies()
		speciesList = sorted(speciesList, key=lambda x:x.name)
		for species in speciesList:
			#write table entry
			fout.write("<tr><td>"+species.getHtmlName()+"<td><a href=\"species_"+species.name+".html\">details</a><br>\n")
			#create reaction page
			species.makeHtmlPage(self)
		fout.write("</table>\n")
		fout.write(patmo_string.readFile("htmlDocs/footer.src")+"\n")
		fout.close()

	#*********************
	#prepare html reactions list
	def makeHtmlReactionsList(self):
		fout = open("htmlDocs/reactionsMenu.html","w")
		fout.write(patmo_string.readFile("htmlDocs/header.src")+"\n")
		fout.write("<table width=\"40%\">\n")
		for reaction in (self.reactions+self.photoReactions):
			#write table entry
			fout.write(reaction.getHtmlTableRow()+"\n")
			#create reaction page
			reaction.makeHtmlPage()
		fout.write("</table>\n")
		fout.write(patmo_string.readFile("htmlDocs/footer.src")+"\n")
		fout.close()

	#**********************
	#plot rates to PNGs
	def plotRates(self):

		#ouput folder
		outFolder = "htmlDocs/ratePNGs/"
		#create output folder if not present
		if(not(os.path.exists(outFolder))):
			os.makedirs(outFolder)

		#maximum span of a plot (to avoid small values)
		maxOrder = 6
		#number of grid points in Tgas (log)
		npoints = 30
		#min temperature
		Tmin = 3e0
		#max temperature
		Tmax = 1e4

		logTmin = log10(Tmin)
		logTmax = log10(Tmax)
		#compute grid points in log space
		aTgas = [1e1**(i*(logTmax-logTmin)/(npoints-1)+logTmin) for i in range(npoints)]


		#loop on reactions
		for reaction in self.reactions:
			plt.clf()

			#replace shortcuts
			rate = patmo_string.shortcutReplace(reaction.rate)

			ydata = []
			#evaluation OK flag
			evalOK = False
			#loop on Tgas to evaluate rate
			for Tgas in aTgas:
				evalRate = rate.replace("tgas",str(Tgas))
				#try to evaluate, if error skip
				try:
					ydata.append(eval(evalRate))
					evalOK = True
				except:
					ydata.append(0e0)
			#skip plot if no evaluation
			if(not(evalOK)): continue
			print ("plotting: "+reaction.getVerbatim())
			#avoid plot that spans many orders of magnitude
			if(max(ydata)/(min(ydata)+1e-40)>1e1**maxOrder):
				plt.ylim(max(ydata)/1e1**maxOrder,max(ydata)*1e1)
			#increase y span for constant plots
			if(max(ydata)==min(ydata)):
				plt.ylim(min(ydata)/1e1,max(ydata)*1e1)
			#set title and others
			plt.title(reaction.getVerbatim())
			plt.grid(True)
			plt.xlabel('Tgas/K')
			plt.ylabel('rate')
			#plot rate
			plt.loglog(aTgas,ydata)
			#save to PNG file
			fnamePNG = outFolder+"rate"+str(int(1e6+reaction.index))+".png"
			plt.savefig(fnamePNG)

	#**********************
	#plot rates to PNGs
	def plotReverseRates(self):

		#ouput folder
		outFolder = "htmlDocs/rateReversePNGs/"
		#create output folder if not present
		if(not(os.path.exists(outFolder))):
			os.makedirs(outFolder)

		#maximum span of a plot (to avoid small values)
		maxOrder = 6
		#number of grid points in Tgas (log)
		npoints = 30

		#compute grid points in log space

		refRates = dict()
		for reaction in self.reactions:
			refRates[reaction.index] = reaction.rate
		for reaction in self.reverseReactions:
			logTmin = log10(reaction.Tmin)
			logTmax = log10(reaction.Tmax)
			aTgas = [1e1**(i*(logTmax-logTmin)/(npoints-1)+logTmin) for i in range(npoints)]
			aTgas += [reaction.Tmin]
			aTgas += [reaction.Tmid]
			aTgas += [reaction.Tmax]
			aTgas = sorted(aTgas)
			plt.clf()

			ydata = []
			ydata2 = []
			#evaluation OK flag
			evalOK = False
			#loop on Tgas to evaluate rate
			for Tgas in aTgas:
				if(Tgas<reaction.Tmid):
					rate = "("+refRates[reaction.reverseOriginalIndex]+")*"+reaction.rateLow
				else:
					rate = "("+refRates[reaction.reverseOriginalIndex]+")*"+reaction.rateUp
				rate = rate.replace("&","").replace("\n","")
				rate = patmo_string.shortcutReplace(rate)

				rate2 = refRates[reaction.reverseOriginalIndex]
				rate2 = rate2.replace("&","").replace("\n","")
				rate2 = patmo_string.shortcutReplace(rate2)

				evalRate = rate.replace("tgas",str(Tgas))
				evalRate2 = rate2.replace("tgas",str(Tgas))
				#try to evaluate, if error skip
				try:
					ydata.append(eval(evalRate))
					ydata2.append(eval(evalRate2))
					evalOK = True
				except:
					ydata.append(0e0)
					ydata2.append(0e0)
			#skip plot if no evaluation
			if(not(evalOK)): continue
			print ("ploting reverse: "+reaction.getVerbatim())
			#avoid plot that spans many orders of magnitude
			#if(max(ydata)/(min(ydata)+1e-40)>1e1**maxOrder and max(ydata)<1e0):
			#	plt.ylim(max(ydata)/1e1**maxOrder,max(ydata)*1e1)
			#increase y span for constant plots
			#if(max(ydata)==min(ydata)):
			#	plt.ylim(min(ydata)/1e1,max(ydata)*1e1)
			plt.ylim(min(ydata+ydata2)+1e-30, max(ydata+ydata2)*1e1)
			#set title and others
			plt.title(reaction.getVerbatim() + " (reverse)")
			plt.grid(True)
			plt.xlabel('Tgas/K')
			plt.ylabel('rate')
			#plot rate
			plt.loglog(aTgas,ydata)
			plt.loglog(aTgas,ydata2,"r--")

			rateMin = "("+refRates[reaction.reverseOriginalIndex]+")*"+reaction.rateLow
			rateMin = rateMin.replace("&","").replace("\n","")
			rateMin = patmo_string.shortcutReplace(rateMin)
			rateMin = eval(rateMin.replace("tgas",str(reaction.Tmin)))

			rateMid = "("+refRates[reaction.reverseOriginalIndex]+")*"+reaction.rateLow
			rateMid = rateMid.replace("&","").replace("\n","")
			rateMid = patmo_string.shortcutReplace(rateMid)
			rateMid = eval(rateMid.replace("tgas",str(reaction.Tmid)))

			rateMax = "("+refRates[reaction.reverseOriginalIndex]+")*"+reaction.rateUp
			rateMax = rateMax.replace("&","").replace("\n","")
			rateMax = patmo_string.shortcutReplace(rateMax)
			rateMax = eval(rateMax.replace("tgas",str(reaction.Tmax)))

			plt.loglog([reaction.Tmin,reaction.Tmid,reaction.Tmax], \
				[rateMin,rateMid,rateMax],"ro")

			#save to PNG file
			fnamePNG = outFolder+"rate"+str(int(1e6+reaction.reverseOriginalIndex))+".png"
			plt.savefig(fnamePNG)

	#*********************
	def buildODE(self,options):

		ODEdict = dict()
		for reaction in (self.reactions+self.photoReactions+self.reverseReactions):
			for reactant in reaction.reactants:
				key = "dn(:,"+reactant.label+")"
				if(not(key in ODEdict)): ODEdict[key] = []
				ODEdict[key].append(" - " + reaction.getRHS())
			for product in reaction.products:
				key = "dn(:,"+product.label+")"
				if(not(key in ODEdict)): ODEdict[key] = []
				ODEdict[key].append(" + " + reaction.getRHS())
		
		fullODE = ""
		for ode, RHS in ODEdict.items():
			species = ode.split(":")[1].split(")")[0].split("_")[2].strip()  # Extract species from the key
			if species not in options.constant_species:
				fullODE += ode + " = &\n" + (" &\n".join(RHS)) + "\n\n"
			
		const_spec = ""
		for species in options.constant_species:
			const_spec +="dn(:,patmo_idx_"+species+") = 0d0"+"\n"

		#if no species 
		if(self.options.constant_species=={} or self.options.drydep_species=={} or self.options.emission_species=={}):
				print ("Warning: No constant_species or dry deposition or emission data in option file")
				
		
		drydep = ""
		for idx, val in options.drydep_species.items():
			drydep += (
				f"if (n(1,patmo_idx_{idx}) > {val}/{options.cellThickness}) then\n"
				f"  dn(1,patmo_idx_{idx}) = dn(1,patmo_idx_{idx}) - ({val}/{options.cellThickness}) * n(1,patmo_idx_{idx})\n"
				f"end if\n"
			)
			if idx in ["CH4", "O2"]:
				drydep += (
					f"{idx}Flux = -{options.cellThickness} * dn(1, patmo_idx_{idx})\n"
					f"dn(1,patmo_idx_{idx}) = 0d0\n"
				)
			
				
		emis_spec = ""		
		for idx, val in options.emission_species.items():
			emis_spec +="dn(1,patmo_idx_"+idx+") = "+"dn(1,patmo_idx_"+idx+") + "+val+"\n"
		
        #replace pragma
		pragmaList = ["#PATMO_ODE","#PATMO_constantspecies","#PATMO_drydeppecies","#PATMO_emissionspecies"]
		replaceList = [fullODE,const_spec,drydep,emis_spec]

		#condition pragmas
		ifPragmas = ["#IFPATMO_useHescape","#IFPATMO_useGravitySettling","#IFPATMO_useAerosolformation","#IFPATMO_useWaterRemoval"]
		ifConditions = [options.useHescape,options.useGravitySettling,options.useAerosolformation,options.useWaterRemoval]
		
		patmo_string.fileReplaceBuild("src_f90/patmo_ode.f90", "build/patmo_ode.f90", \
			pragmaList, replaceList, ifPragmas, ifConditions)

		#create verbatim file
		self.createVerbatimFile()

	#********************
	#creates a file in build directory with the list of the reactions
	def createVerbatimFile(self):
		fout = open("build/reactionsVerbatim.dat","w")
		for reaction in (self.reactions+self.photoReactions+self.reverseReactions):
			fout.write(reaction.getVerbatim()+"\n")
		fout.close()

	#**********************
	#create a topology file with note position (using graphviz), and edges connections map
	def createTopology(self):
		connections = []
		#loop on reactions
		for reaction in self.reactions:
			for reactant in reaction.reactants:
				for product in reaction.products:
					#create edge
					edge = [reactant.name, product.name]
					#if edge known, skip
					#if(not(edge in connections)):
					connections.append(edge)

		#prepare digraph file for graphviz
		dotFile = "digraph g{"
		for edge in connections:
			dotFile += "\""+edge[0]+"\" -> \""+edge[1]+"\";\n"
		dotFile += "}\n"
		#write dot file
		patmo_string.writeFile("tmp.dot", dotFile)

		#call graphviz as external program
		aCall = ["circo","tmp.dot","-Tplain"]
		process = Popen(aCall, stdout=PIPE)
		(output, err) = process.communicate()
		# Decode the output from bytes to a string
		output = output.decode('utf-8')  # or another encoding if necessary
		xposList = []
		yposList = []
		labelList = []
		idxMap = dict()
		#process graphviz output
		for row in output.split("\n"):
			#only use node positions
			if("node" in row):
				arow = [x.strip() for x in row.split(" ") if x.strip()!=""]
				(xpos, ypos, label) = [float(x) for x in arow[2:4]]+[arow[1]]
				#store positions and labels
				xposList.append(xpos)
				yposList.append(ypos)
				label = "\""+label.replace("\"","")+"\""
				labelList.append(label)
				idxMap[label] = len(idxMap)

		#normalize positions
		xposList = [(x-min(xposList))/(max(xposList)-min(xposList)) for x in xposList]
		yposList = [(x-min(yposList))/(max(yposList)-min(yposList)) for x in yposList]

		fout = open("map.dat","w")
		#write positions in map file
		for i in range(len(xposList)):
			fout.write((", ".join([str(x) for x in [xposList[i],yposList[i],labelList[i]]]))+"\n")
		fout.write("\nmap\n")
		icount = 0
		#write connections map in map file
		for edge in connections:
			fout.write(str(icount+1) + ", " + (", ".join([str(idxMap["\""+x+"\""]) for x in edge]))+"\n")
			icount += 1
		fout.close()

		#remove temporary dot file
		os.remove("tmp.dot")



