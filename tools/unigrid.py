from scipy import interpolate
import sys

#check arguments
if(len(sys.argv)!=3):
	print "Usage is"
	print " python "+sys.argv[0]+" input output"
	sys.exit()

#raead arguments
fileInput = sys.argv[1]
fileOutput = sys.argv[2]

#open file to read
fh = open(fileInput,"rb")
newFile = "#this file has been interpolated to have an equally spaced grid\n"
rowcount = 0 #count read rows
alt = [] #store height
print "reading from "+fileInput
#loop on file rows
for row in fh:
	srow = row.strip()
	if(srow==""): continue
	#store header
	if(srow.startswith("#")):
		newFile += row
		continue
	rowcount += 1
	#first line indicates non-species columns and non-zero columns
	if(rowcount==1):
		newFile += row
		continue
	#second is position header
	if(rowcount==2):
		newFile += row
		#read position header
		head = [x.strip().lower() for x in srow.replace(" ","\t").split("\t") if x!=""]
		iidx = -1
		#if index column is present store it
		if("index" in head): iidx = head.index("index")
		#look for height column, if not present rise error
		try:
			ialt = head.index("alt")
		except:
			print "ERROR: missing \"alt\" in header!"
			sys.exit()
		#init read data matrix
		data = [[] for i in range(len(head))]
		#init interpolated data matrix
		data_int = [[] for i in range(len(head))]
		continue

	#rows larger than second are data
	if(rowcount>2):
		#read row
		arow = [float(x) for x in srow.replace(" ","\t").split("\t") if x!=""]
		#store height
		alt.append(arow[ialt])
		#put data in the columns
		for i in range(len(arow)):
			data[i].append(arow[i])
fh.close()

#prepare equally spaced grid (assuming that spacing starts from zero).
# it uses zero as an additional value, but cut it's not necessary so [1:]
xnew = [round(i*max(alt)/(len(alt)),5) for i in range(len(alt)+1)][1:]

print "interpolating..."

#loop on columns
for i in range(len(data)):
	col = data[i]
	#create interpolator for that column
	f = interpolate.interp1d(alt,col)
	#row counter
	jcount = 0
	#loop on values to interpolate on
	for vx in xnew:
		jcount += 1
		if(i==iidx):
			#index column it's just an integer (no need to interpolate)
			data_int[i].append(jcount)
		else:
			#calc and round interpolated data
			data_int[i].append(round(f(vx),5))

#loop on interpolated data and write tab-spaced columns
for j in range(len(data_int[0])):
	newFile += "\t".join([str(data_int[i][j]) for i in range(len(data))]) + "\n"

#write new file
print "saving to "+fileOutput
fout = open(fileOutput,"w")
for row in newFile:
	fout.write(row)
fout.close()
print "Everything done!"
