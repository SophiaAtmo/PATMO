import shutil,os,sys,inspect
import patmo_string
import patmo_error

def buildTest(testName):

	sourceFolder = patmo_string.pathFormat("tests/"+testName)
	if(not(os.path.isdir("tests/"+testName)) or (testName.strip()=="")):
		print ("ERROR: unknown test name "+testName)
		folders = next(os.walk('tests/'))[1]
		print ("Tests are: "+(", ".join(folders)))
		patmo_error.trigError(__file__,inspect.currentframe())

	destinationFolder = "build/"
	copyList = []
	#read list of files to copy
	fh = open(sourceFolder+"copylist.pcp","r")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		copyList.append(srow)
	fh.close()

	#copy files found in copylist file
	for fname in copyList:
		print ("copy "+(sourceFolder+fname)+" -> "+(patmo_string.pathFormat(destinationFolder)+fname))
		shutil.copyfile(sourceFolder+fname, \
				patmo_string.pathFormat(destinationFolder)+fname)
