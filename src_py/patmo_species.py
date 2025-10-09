import sys,itertools
import patmo_utils
import patmo_string

class species:

	ref = []
	massDict = dict()
	#load data and mass dictionary from file
	fh = open("data/parseReference/reference.dat","r")
	for row in fh:
		srow = row.strip().replace("\t"," ")
		if(srow==""): continue
		if(srow.startswith("#")): continue
		(reference,mass) = [x.strip() for x in srow.split(" ") if x!=""]
		ref.append(reference)
		massDict[reference] = patmo_string.evalMass(mass)
	fh.close()

	#create XYZ hash for each element
	href = list(itertools.product('XYZ', repeat=3))
	href = [("".join(x)) for x in href][:len(ref)]
	#reference dictionary element=>hash
	dref = dict()
	#backward reference dictionary hash=>element
	drefB = dict()
	#create dictionaries
	for i in range(len(ref)):
		dref[ref[i]] = href[i]
		drefB[href[i]] = ref[i]
	#sort by element length (longest first)
	ref = sorted(ref, key=lambda x:len(x),reverse=True)

	#********************
	#constructor
	def __init__(self):
		self.atomDict = dict()
		self.name = None
		self.charge = 0
		self.label = None
		self.mass = 0e0
		self.index = None
		self.htmlName = None

	#*****************
	#make an html page for this reaction
	def makeHtmlPage(self,network):
		fout = open("htmlDocs/species_"+self.name+".html","w")
		fout.write(patmo_string.readFile("htmlDocs/header.src")+"\n")
		fout.write(self.getHtmlName()+"<br><br>\n")

		fout.write("<p>Formation Reactions:</p>")
		fout.write("<table width=\"40%\">")
		for reaction in self.getFormationReactions(network):
			fout.write(reaction.getHtmlTableRow(highlightSpecies=self.name)+"<br>\n")
		fout.write("</table><br><br>")

		fout.write("<p>Destruction Reactions:</p>")
		fout.write("<table width=\"40%\">")
		for reaction in self.getDestructionReactions(network):
			fout.write(reaction.getHtmlTableRow(highlightSpecies=self.name)+"<br>\n")
		fout.write("</table>")

		fout.write(patmo_string.readFile("htmlDocs/footer.src")+"\n")
		fout.close()

	#******************
	#get the list of reactions that form the self species in the network
	def getFormationReactions(self,network):
		formationReactions = []
		#loop on reactions
		for reaction in network.reactions:
			#if species found in products append reaction
			if(reaction.isContainingSpeciesNameProducts(self.name)):
				formationReactions.append(reaction)

		#return reactions list
		return formationReactions

	#******************
	#get the list of reactions that destroy the self species in the network
	def getDestructionReactions(self,network):
		destructionReactions = []
		#loop on reactions
		for reaction in network.reactions:
			#if species found in products append reaction
			if(reaction.isContainingSpeciesNameReactants(self.name)):
				destructionReactions.append(reaction)

		#return reactions list
		return destructionReactions

	#*******************
	#get the html name of this species
	def getHtmlName(self):
		if(self.htmlName!=None): return self.htmlName
		self.htmlName = ""
		#character after which no subscript number is required, e.g. O(1D)
		noSubChar = ["_","("]
		charOld = None
		#loop on name characters
		for char in self.name[:]:
			#subscript if number and not special number, e.g. O(1D)
			if(patmo_string.isInteger(char) and not(charOld in noSubChar)):
				char = "<sub>"+char+"</sub>"
			self.htmlName += char
			charOld = char
		return self.htmlName

	#********************
	#parse species
	def parse(self,spec):
		self.name = spec
		cspec = spec #copy to replace
		rspec = spec #copy to remove

		self.label = spec.replace("+","j")
		self.label = self.label.replace("-","w")
		self.label = self.label.replace(")","")
		self.label = self.label.replace("(","_")
		self.label = "patmo_idx_"+self.label.strip()

		#check for cations
		if("+" in spec):
			cspec = cspec.replace("+","")
			rspec = rspec.replace("+","")
			self.charge = spec.count("+")
		#check for anions
		if("-" in spec):
			cspec = cspec.replace("-","")
			rspec = rspec.replace("-","")
			self.charge = spec.count("-")

		#loop on elements to replace hash and delete found
		for r in self.ref:
			cspec = cspec.replace(r,"#"+self.dref[r]+"#")
			rspec = rspec.replace(r,"")
		#replace number from remove copy (up to 10)
		for r in range(10):
			rspec = rspec.replace(str(r),"")
		#if something remains in the remove copy rise error
		if(rspec.strip()!=""):
			print ("original:", spec)
			print ("hashed  :", cspec)
			print ("residual:", rspec)
			sys.exit("ERROR: parsing found unknown parts")
		#divide hashes
		parts = cspec.split("#")
		#loop on parts (hashes)
		for i in range(len(parts)):
			part = parts[i]
			#skip empty parts
			if(part.strip()==""): continue
			#skip numbers
			if(patmo_utils.isNumber(parts[i])): continue
			multi = 1 #multiplication factor (default)
			#if it's not the last part check next for multiplier
			if(i<len(parts)-1):
				#check if next part is a number, if so consider as multiplier
				if(patmo_utils.isNumber(parts[i+1])): multi = int(parts[i+1])
			self.mass += self.massDict[self.drefB[part]]*multi
			#store into dictionary (e.g. H2O is {H:2,O:1}
			self.atomDict[self.drefB[part]] = multi

