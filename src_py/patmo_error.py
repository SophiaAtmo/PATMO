import sys,inspect

#error trigger, requirese fileName as __file__
# and frame as inspect.currentframe()
#Example call is
#   patmo_error.trigError(__file__,inspect.currentframe())
def trigError(fileName,frame):
	srcFileName = fileName.replace(".pyc",".py")
	print ("(Triggered by source file: "+srcFileName+" [" \
		+ str(inspect.getlineno(frame))+"])")
	#stop program execution
	sys.exit()

