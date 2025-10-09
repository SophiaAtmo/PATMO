import os,sys


#************************
def evalMass(massString):
	me = 9.10938356e-28
	mp = 1.6726219e-24
	mn = 1.674927471e-24
	mstr = massString
	mDict = {"me":me,"mp":mp,"mn":mn,"mep":me+mp,"mpne":me+mp+mn}
	labels = sorted(mDict.keys(),key=lambda x:len(x),reverse=True)
	for mlabel in labels:
		mstr = mstr.replace(mlabel,str(mDict[mlabel]))
	return eval(mstr)

#*********************
def readFile(fname):
	fileContent = ""
	fh = open(fname,"r")
	for row in fh:
		fileContent += row
	fh.close()
	return fileContent

#*********************
def writeFile(fname,fileContent):
	fh = open(fname,"w")
	fh.write(fileContent)
	fh.close()

#************************
#load source file and prepare the build file replacing prgamas with replacements
def fileReplaceBuild(sourceFile,destinationFile,pragmaList=[],replaceList=[],pragmaIfList=[],ifList=[]):
	#check for errors
	if(len(pragmaList)!=len(replaceList)): sys.exit("ERROR: pragma list lenght is different from replace list!")
	if(len(pragmaIfList)!=len(ifList)): sys.exit("ERROR: pragma IF list lenght is different from replace list!")
	if(not(os.path.isfile(sourceFile))): sys.exit("ERROR: source file "+sourceFile+" not found!")
	if(("src_f90" in destinationFile) or ("src_py" in destinationFile)): sys.exit("ERROR: you are writing in src_* folders!")

	#sort pragma by length in order to have the longer first (to avoid misreplacing)
	pragamaIndexes = [i[0] for i in sorted(enumerate(pragmaList), key=lambda x:len(x[1]),reverse=True)]
	pragamaIfIndexes = [i[0] for i in sorted(enumerate(pragmaIfList), key=lambda x:len(x[1]),reverse=True)]

	#preapre a found map to detect non-found pragmas
	foundMapPragma = [False for i in range(len(pragmaList))]
	foundMapPragmaIf = [False for i in range(len(pragmaIfList))]

	#open source file and destination file
	fh = open(sourceFile,"r")
	fout = open(destinationFile,"w")
	skipDepth = 0 #depth in the IF tree
	skip = dict() #store the last skip conditions found
	skip[0] = False #root never skips
	#loop on the source lines
	for row in fh:
		#if IFPATMO statement then evaluates
		if("#IFPATMO_" in row.strip()):
			foundIf = False
			#loop on if condition list
			for icount in pragamaIfIndexes:
				pragmaIf = pragmaIfList[icount] #IF pragma
				ifCondition = ifList[icount] #IF condition
				#if the line corresponds to the selected IF pragma
				if(row.strip()==pragmaIf):
					skipDepth += 1 #increase skip depth
					foundMapPragmaIf[icount] = True #store found PragmaIf
					skip[skipDepth] = not(ifCondition) #store skip condition for the depth level
					#if parent level is set to skip no matter if the current skips or not
					if(skip[skipDepth-1]): skip[skipDepth] = True
					foundIf = True
					break #once found no need to go on
			#rise error if pragma IF not found in the list
			if(not(foundIf)):
				print ("ERROR: in file " + sourceFile)
				print ("pragma "+row.strip()+" not present in the list!")
				sys.exit()
			continue

		#if ELSE statement found, changes the skip behaviour
		if(row.strip()=="#ELSEPATMO" or row.strip()=="#ELSE_PATMO"):
			skip[skipDepth] = not(skip[skipDepth])
			#if parent depth is skip, ELSE not matters and should skip anyway
			if(skip[skipDepth-1]): skip[skipDepth] = True
			continue

		#go back to the upper depth level
		if(row.strip()=="#ENDIF_PATMO" or row.strip()=="#ENDIFPATMO"):
			skipDepth -= 1
			#check if there are more ENDIFs than IFs (depth cannot be negative)
			if(skipDepth<0):
				print ("ERROR: too many ENDIFs in "+sourceFile)
				sys.exit()
			continue

		#loop over pragmas to replace
		for icount in pragamaIndexes:
			foundMapPragma[icount] = True
			row = row.replace(pragmaList[icount],str(replaceList[icount]))
		#check if pragma is contained in the non-comment region of the line
		rowHasPragma = ("#PATMO" in row.split("!")[0])
		#if row has still a pragma means that substitution didn't worked and rise error
		if(rowHasPragma):
			print ("ERROR: in file " + sourceFile)
			print ("pragma "+row.strip()+" not present in the list!")
			sys.exit()

		#if it has to skip goes to the next line
		if(skip[skipDepth]): continue

		#write the replaced line
		fout.write(row)
	fout.close()
	fh.close()

	#check if IF blocks are open
	if(skipDepth>0):
		print ("ERROR: some IF blocks are still open in "+sourceFile)
		sys.exit()



	#check if all pragmas have been found
	if(False in foundMapPragma):
		print ("ERROR: some pragmas have not been found in file "+sourceFile)
		for icount in pragamaIndexes:
			print (pragmaList[icount], foundMapPragma[icount])
		sys.exit()

	#check if all IF pragmas have been found
	if(False in foundMapPragmaIf):
		print ("ERROR: some pragmas have not been found in file "+sourceFile)
		for icount in pragamaIfIndexes:
			print (pragmaIfList[icount], foundMapPragmaIf[icount])
		sys.exit()

#******************************
#convert an hash in the format e.g. R|R||P|P to R+R->P+P
def hash2verbatim(reactionHash):
	return reactionHash.replace("||"," -> ").replace("|"," + ")

#******************************
def shortcutReplace(argin):

	#shortcuts for temperature
	#note: sharp # is a signpost for (i)
	shortcuts = {"invT":"1d0/Tgas", \
		"T":"Tgas",\
		"invTgas":"1e0/Tgas",\
		"lnTgas#":"log(Tgas)",\
		"Tgas#":"Tgas",\
		"Tgas2#":"Tgas**2",\
		"Tgas3#":"Tgas**3",\
		"Tgas4#":"Tgas**4",\
		"invTgas#":"1e0/Tgas",\
		}

	#operators to replace shortcuts, e.g. *Tgas/)
	maths = ["+","-","*","/","(",")"]


	arg = argin.strip().replace(" ","").lower().replace("(i)","#")
	#loop on shortucts
	for k,v in shortcuts.items():
		#left operator
		for mL in maths:
			#check if rate starts with shortcut+operator
			if(arg.startswith(k.lower()+mL)):
				arg = arg.replace(k.lower()+mL,"("+v.lower()+")"+mL)
			#check if rate ends with operator+shortcut
			if(arg.endswith(mL+k.lower())):
				arg = arg.replace(mL+k.lower(),mL+"("+v.lower()+")")
			#add right operator and replace if necessary
			for mR in maths:
				arg = arg.replace(mL+k.lower()+mR,mL+"("+v.lower()+")"+mR)

	return arg.replace("d","e")


#*******************************
#convert a floating into a F90 double exponential
# e.g. "1.23" -> 1.230000d+00
def f90Double(floating):
	ff = '%e' % floating
	return ff.replace("e","d")

#*******************************
#convert a floating into a F90 double exponential (compact)
# e.g. "0.1" -> 1.000000d-01 -> 1d-1
def f90DoubleCompact(floating):
	ff = '%e' % floating
	fs = ff.replace("e","d")
	fs = fs.replace("d+0","d")
	fs = fs.replace("d-0","d-")
	while("0d" in fs):
		fs = fs.replace("0d","d")
	return fs.replace(".d","d")

#*******************************
#add a slash to the path if not present and remove double slashes
def pathFormat(path):
	spath = path.strip()
	if(spath[-1]=="/"): return spath.replace("//","/")
	return spath.replace("//","/") + "/"

#*******************************
#this function returns if the string line
# starts with one of the items in the array 
# of string aarg
def lbeg(aarg,line):
	for arg in aarg:
		if(line[:len(arg)]==arg): return True
	return False

#*******************************
#this function returns if the string line
# ends with one of the items in the array 
# of string aarg
def lend(aarg,line):
	for arg in aarg:
		if(line[len(line)-len(arg):]==arg): return True
	return False

#*******************************
#check if a string is an integer or not
def isInteger(arg):
	try:
		int(arg)
		return True
	except:
		return False

#*******************************
#check if a string is a number or not
def isNumber(arg):
	try:
		float(arg)
		return True
	except:
		return False

#*******************************
#this function indent f90 file and remove multiple blank lines
def indentF90(filename):
	import os
	#check if the file exists else return
	if(not(os.path.isfile(filename))): return

	#open file for indent
	fh = open(filename,"r")
	arow = [] #array of the lines of the indented file
	is_blank = is_amper = False #flags
	nind = 0 #number indent level
	nspace = 2 #number of space for indent
	tokenclose = ["end do","end if","end function","end subroutine","else if","elseif","else","enddo","end module",\
		"endif","end type", "endtype", "contains","endfunction","endsubroutine","endmodule","end program", "endprogram",\
		"end interface","endinterface","module procedure"]
	tokenopen = ["do ","function","subroutine","contains","else","else if","elseif","module","program", "type,","interface",\
		"module procedure"]
	module_head = "!############### MODULE ##############" #module header comment
	module_head_found = False
	for row in fh:
		srow = row.strip() #trim the row
		#if(module_head in srow): module_head_found = True #do not duplicate module header
		#check module begin
		#if(lbeg(["module"], srow) and not(module_head_found)): 
		#	#arow.append("\n") #blank line
		#	#arow.append(module_head+"\n") #comment
		if(lbeg(tokenclose, srow)): nind -= 1 #check if the line ends with one of tokenclose
		indent = (" "*(nind*nspace)) #compute number of spaces for indent
		if(is_amper): indent = (" "*(2*nspace)) + indent #increas indent in case of previous &
		if(srow.startswith("#")): indent = "" #no indent for pragmas
		if(not(srow=="" and is_blank)): arow.append(indent+srow+"\n") #append indented line to array of rows
		is_amper = False #is a line after ampersend flag
		if(lend(["&"], srow)): is_amper = True #check if the line ends with &
		is_blank = (srow=="") #flag for blank line mode
		if(lbeg(tokenopen, srow)): nind += 1 #check if the line ends with one of tokenclose
		if(lbeg(["if"],srow) and "then" in srow): nind += 1 #check if line stats with if and has then
		if(srow=="do"): nind += 1
	fh.close()

	arowl = arow[:]

	#write the new file
	fh = open(filename,"w")
	for x in arowl:
		fh.write(x.rstrip()+"\n")
	fh.close()

#*******************************
#indent a list of files
def indentFileList(fileList):
	#loop on the list of files
	for fileName in fileList:
		#check if the file exists
		if(not(os.path.isfile(fileName))):
			sys.exit("ERROR: "+fileName+" is not present. Cannot indent!")
		#do intent
		indentF90(fileName)
