import os
import patmo_string
import patmo_reaction
import patmo_species
from math import log10,exp
from scipy.interpolate import interp1d

class photochemistry:


	#*************************
	#constructor
	def __init__(self,network,options):
		self.network = network
		self.options = options
		self.energyMetric = dict()

		#skip this if no photochemistry is required
		if(not(self.options.usePhotochemistry)): return

		folder = patmo_string.pathFormat("data/xsecs_swri/")
		allSpecies = network.getSpecies()
		allSpeciesNames = network.getSpeciesNames()

		#loop on species
		for species in allSpecies:
			myReactions = self.findReactions(species,folder)
			#loop on reactions found in xsecs databases
			for reaction in myReactions:
				productsOK = True
				#check if products are present in species list
				for product in reaction.products:
					if(not(product.name in allSpeciesNames)):
						productsOK = False
						break
				#append reaction to network if all products are present
				if(productsOK):
					print ("Xsec found for: "+reaction.getVerbatim())
					network.photoReactions.append(reaction)

		notFoundXsec = []
		#check missing xsec filenames to print warning
		for species in allSpecies:
			foundSpecies = False
			#loop on file names
			for fname in os.listdir(folder):
				if(not(fname.endswith(".dat"))): continue
				if(species.name==fname.split(".")[0]): foundSpecies = True
			#append missing species
			if(not(foundSpecies)): notFoundXsec.append(species.name)
		#print warning
		if(len(notFoundXsec)>0):
			print ("Missing xsec files for "+(", ".join(sorted(notFoundXsec))))

	#*************************
	#parse swri header elements
	def parseSwriHeader(self,headerElement):
		headerElement = headerElement.replace("O1D","O(1D)")
		headerElement = headerElement.replace("C1D","C(1D)")
		headerElement = headerElement.replace("O2d","O2")
		headerElement = headerElement.replace("O1S","O(1S)")
		headerElement = headerElement.replace("sCO","CO")
		headerElement = headerElement.replace("sCH2","CH2")
		headerElement = headerElement.replace("+","+/")
		if(headerElement.endswith("+/")): headerElement += "e-"
		products = []
		for species in headerElement.split("/"):
			mySpecies = patmo_species.species()
			mySpecies.parse(species)
			products.append(mySpecies)
		return products

	#*************************
	#search xsecs database for reactions
	def findReactions(self,species,folder):
		hplanck_eV = 4.135667662e-15 #eV*s
		clight = 2.99792458e10 #cm/s
		skipHead = ["Lambda","Total"]
		myReactions = []
		#loop on swri xsecs
		for fname in os.listdir(folder):
			if(not(fname.endswith(".dat"))): continue
			#check if species name matches the file name
			if(species.name==fname.split(".")[0]):
				readData = isHeader = False
				#open file
				fh = open(folder+fname,"r")
				#read file
				for row in fh:
					srow = row.strip()
					if(srow==""): continue
					#after this line next is header
					if(srow.startswith("0 Branching ratio for")):
						isHeader = True
						continue
					#searh for header
					if(isHeader):
						readData = True
						isHeader = False
						header = [x.strip() for x in srow.split(" ") if x!=""]
						myReactions = dict()
						#loop on header columns to skip
						for head in header:
							if(head in skipHead): continue
							if("tCO" in head): continue
							myReactions[head] = patmo_reaction.reaction()
							myReactions[head].products = self.parseSwriHeader(head)
							myReactions[head].reactants = [species]
							myReactions[head].xsec = []
							myReactions[head].xsecEnergy = []
							#print myReactions[head].getVerbatim()
						continue
					if(not(readData)): continue
					#read data
					arow = [float(x) for x in srow.split(" ") if x!=""]
					for ihead in range(len(header)):
						head = header[ihead]
						if(head in myReactions):
							wl = arow[0]*1e-8 #AA->cm
							myReactions[head].xsec.append(arow[ihead])
							myReactions[head].xsecEnergy.append(hplanck_eV*clight/wl) #eV

				#reverse data to have increasing energy
				for head in myReactions.keys():
					myReactions[head].xsec = myReactions[head].xsec[::-1]
					myReactions[head].xsecEnergy = myReactions[head].xsecEnergy[::-1]
				#copy reactions (dict->list values)
				myReactions = myReactions.values()
				#reaction found in the file can break
				break
		return myReactions

	#*****************
	def createLogMetric(self):

		#skip this if no photochemistry is required
		if(not(self.options.usePhotochemistry)): return

		energyMin = self.options.energyMin
		energyMax = self.options.energyMax
		binsNumber = self.options.photoBinsNumber

		#if optional argument missing find limits automatically
		if(energyMin==None):
			emin = 9e99
			#loop on rates to find limits (min/max)
			for reaction in self.network.photoReactions:
				emin = min(min(reaction.xsecEnergy),emin)
		else:
			emin = energyMin
		#find max limit automatically if required
		if(energyMax==None):
			emax = 0e0
			#loop on rates to find limits (min/max)
			for reaction in self.network.photoReactions:
				emax = max(max(reaction.xsecEnergy),emax)
		else:
			emax = energyMax

		#print message if auto min/max
		if(energyMin==None):
			print ("Automatic min limit (eV):",emin)

		if(energyMax==None):
			print ("Automatic max limit (eV):",emax)


		#log energy limits
		lemin = log10(emin)
		lemax = log10(emax)

		#init energy metric lists
		self.energyMetric["left"] = []
		self.energyMetric["right"] = []
		self.energyMetric["mid"] = []
		self.energyMetric["span"] = []

		reportFolder = "reports/xsecs/"
		xsecBuildFolder = "build/xsecs/"

		#create xsecs folder in build if not there
		if(not(os.path.exists(xsecBuildFolder))):
			os.makedirs(xsecBuildFolder)

		#loop on bins to create log metric
		for i in range(binsNumber):
			energyLeft = 1e1**(i*(lemax-lemin)/binsNumber+lemin)
			energyRight = 1e1**((i+1)*(lemax-lemin)/binsNumber+lemin)
			self.energyMetric["left"].append(energyLeft)
			self.energyMetric["right"].append(energyRight)
			self.energyMetric["mid"].append((energyLeft+energyRight)/2.)
			self.energyMetric["span"].append(energyRight-energyLeft)

		#write metric to file as mid, span (loaded by f90)
		fout = open(patmo_string.pathFormat(xsecBuildFolder)+"photoMetric.dat","w")
		fout.write("#photochemistry metric, mid, span, left, right (in eV)\n")
		for i in range(binsNumber):
			edata = [self.energyMetric["mid"][i], self.energyMetric["span"][i], \
				self.energyMetric["left"][i], self.energyMetric["right"][i]]
			fout.write((" ".join([str(x) for x in edata]))+"\n")
		fout.close()


		#set the starting index for photo rates
		indexBase = len(self.network.reactions)
		index = 1
		#interpolate rates on metric
		for reaction in self.network.photoReactions:
			#dump xsec loaded from file
			reaction.dumpXsec(reportFolder,postpone="_org.dat")

			#interpolate loaded reaction rate
			finterp = interp1d(reaction.xsecEnergy,reaction.xsec)
			minEnergy = min(reaction.xsecEnergy)
			maxEnergy = max(reaction.xsecEnergy)
			#empty xsec, will be replaced with interpolated
			reaction.xsecEnergy = []
			reaction.xsecEnergySpan = []
			reaction.xsec = []
			#interpolate xsecs on the current metric
			for ienergy in range(len(self.energyMetric["mid"])):
				xL = self.energyMetric["left"][ienergy]
				xR = self.energyMetric["right"][ienergy]
				fR = fL = 0e0
				if(xL>=minEnergy and xL<=maxEnergy):
					fL = finterp(xL)
				if(xR>=minEnergy and xR<=maxEnergy):
					fR = finterp(xR)
				#add energy from metric
				reaction.xsecEnergy.append(self.energyMetric["mid"][ienergy])
				reaction.xsecEnergySpan.append(self.energyMetric["right"][ienergy] \
					- self.energyMetric["left"][ienergy])
				#add xsec from interpolation
				reaction.xsec.append((fL+fR)/2.)

			#define reaction index
			reaction.index = indexBase + index
			reaction.photoIndex = index
			#dump interpolated xsec to file
			reaction.dumpXsec(reportFolder,postpone="_interp.dat")
			reaction.dumpXsec(xsecBuildFolder,postpone=".dat")
			reaction.rate = "integrateXsec("+str(index)+", tau(:,:))"
			#increase counter reaction index
			index += 1


	#**********************
	def buildLoadPhotoRates(self):
		allLoad = ""
		for reaction in self.network.photoReactions:
			fname = "xsecs/"+reaction.getFileSafeVerbatim()+".dat"
			allLoad += "!"+reaction.getVerbatim()+"\n"
			allLoad += "call loadPhotoXsec(\""+fname+"\","+str(reaction.photoIndex)+")\n"

		#replace commons pragma
		pragmaList = ["#PATMO_loadAllPhotoXsecs"]
		replaceList = [allLoad]

		patmo_string.fileReplaceBuild("src_f90/patmo_photo.f90", "build/patmo_photo.f90", \
			pragmaList, replaceList)
