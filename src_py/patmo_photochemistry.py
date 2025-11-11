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
						#header = [x.strip() for x in srow.split(" ") if x!=""]
						header = [x for x in srow.replace("\t", " ").split() if x]
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
					arow = [float(x) for x in srow.replace("\t", " ").split() if x]
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
		# Skip if photochemistry is disabled
		if not self.options.usePhotochemistry:
			return
    
		# Physical constants (Planck constant, speed of light)
		hplanck_eV = 4.135667662e-15   # eV*s
		clight     = 2.99792458e10     # cm/s
		nm_to_cm   = 1e-7              # nm → cm
    
		binsNumber = self.options.photoBinsNumber
    
        # ===== Determine energy range from wavelength limits (nm) =====
		wl_min_nm = self.options.wavelengMin
		wl_max_nm = self.options.wavelengMax
    
		if (wl_min_nm is not None) and (wl_max_nm is not None):
            # Ensure correct order
			if wl_min_nm > wl_max_nm:
				wl_min_nm, wl_max_nm = wl_max_nm, wl_min_nm
    
            # Convert wavelength to energy (E = h * c / λ)
            # Note: shorter wavelength -> higher energy
			emin = hplanck_eV * clight / (wl_max_nm * nm_to_cm)  # energy at λ_max
			emax = hplanck_eV * clight / (wl_min_nm * nm_to_cm)  # energy at λ_min
    
			print("Using wavelength limits (nm): [{}, {}]".format(wl_min_nm, wl_max_nm))
			print("Derived energy limits (eV)  :", emin, emax)
		else:
            # If wavelength range not provided, find limits automatically
			emin = 9e99
			emax = 0.0
			for reaction in self.network.photoReactions:
				if reaction.xsecEnergy:
					emin = min(emin, min(reaction.xsecEnergy))
					emax = max(emax, max(reaction.xsecEnergy))
			if not (emin < emax):
				raise RuntimeError("Cannot determine energy range: no photolysis xsecs loaded.")
			print("Automatic energy limits (eV):", emin, emax)
    
        # ===== Create logarithmic bins in energy space =====
		from math import log10
		lemin = log10(emin)
		lemax = log10(emax)

		self.energyMetric = {"left": [], "right": [], "mid": [], "span": []}

		reportFolder    = "reports/xsecs/"
		xsecBuildFolder = "build/xsecs/"
		import os
		if not os.path.exists(xsecBuildFolder):
			os.makedirs(xsecBuildFolder)

		for i in range(binsNumber):
			energyLeft  = 10.0 ** (i      * (lemax - lemin) / binsNumber + lemin)
			energyRight = 10.0 ** ((i + 1) * (lemax - lemin) / binsNumber + lemin)
			self.energyMetric["left"].append(energyLeft)
			self.energyMetric["right"].append(energyRight)
			self.energyMetric["mid"].append((energyLeft + energyRight) / 2.0)
			self.energyMetric["span"].append(energyRight - energyLeft)

        # ===== Write metric file (in eV) for Fortran side =====
		import patmo_string
		fout = open(patmo_string.pathFormat(xsecBuildFolder) + "photoMetric.dat", "w")
		fout.write("# Photochemistry metric: mid, span, left, right (in eV)\n")
		for i in range(binsNumber):
			edata = [
                self.energyMetric["mid"][i],
                self.energyMetric["span"][i],
                self.energyMetric["left"][i],
                self.energyMetric["right"][i],
            ]
			fout.write((" ".join([str(x) for x in edata])) + "\n")
		fout.close()
    
        # ===== Interpolate all reactions on the new metric =====
		indexBase = len(self.network.reactions)
		index = 1
    
		from scipy.interpolate import interp1d
		for reaction in self.network.photoReactions:
            # Save the original cross-section (as loaded from file)
			reaction.dumpXsec(reportFolder, postpone="_org.dat")
    
            # Interpolator (allow extrapolation outside range -> 0)
			finterp = interp1d(reaction.xsecEnergy, reaction.xsec,
                               bounds_error=False, fill_value=0.0)
			minEnergy = min(reaction.xsecEnergy)
			maxEnergy = max(reaction.xsecEnergy)
    
            # Reset arrays before storing interpolated results
			reaction.xsecEnergy = []
			reaction.xsecEnergySpan = []
			reaction.xsec = []
    
			for ienergy in range(len(self.energyMetric["mid"])):
				xL = self.energyMetric["left"][ienergy]
				xR = self.energyMetric["right"][ienergy]
				fL = finterp(xL) if (minEnergy <= xL <= maxEnergy) else 0.0
				fR = finterp(xR) if (minEnergy <= xR <= maxEnergy) else 0.0
    
				reaction.xsecEnergy.append(self.energyMetric["mid"][ienergy])
				reaction.xsecEnergySpan.append(self.energyMetric["right"][ienergy] - self.energyMetric["left"][ienergy])
				reaction.xsec.append((fL + fR) / 2.0)
    
            # Assign indices and dump interpolated results
			reaction.index = indexBase + index
			reaction.photoIndex = index
			reaction.dumpXsec(reportFolder, postpone="_interp.dat")
			reaction.dumpXsec(xsecBuildFolder, postpone=".dat")
			reaction.rate = "integrateXsec(" + str(index) + ", tau(:,:))"
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
