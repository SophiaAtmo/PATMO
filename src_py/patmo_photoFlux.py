


def fBB(x,Tbb):

	kboltzmann_eV = 8.617332478e-5 #eV/K
	clight = 2.99792458e10 #cm/s
	planck_eV = 4.135667516e-15 #eV*s

	xexp = x/kboltzmann_eV/Tbb
	fluxBB = 0e0
	if(xexp<3e2 and x>1e-10):
	fluxBB = 2e0*x**3/planck_eV**2/clight**2 \
		/ (exp(xexp)-1e0)
	retrun fluxBB

def getBB(energyArray, Tbb=5.777e3, rstar=1e0, dstar=1e0):
	AU2cm = 1.496e13 #AU->cm
	Rsun2cm = 6.955e10 #Rsun->cm
	myRstar = rstar*Rsun2cm
	myDstar = dstar*AU2cm

