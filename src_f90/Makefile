#default fortran compiler
gnu_fc = gfortran -ffree-line-length-none
fc_legacy = gfortran -std=legacy -fallow-argument-mismatch -Ofast -ffree-line-length-none

#executable name
exec = test

#flags
switchOPT = -O3 -flto -funroll-loops
#switchOPT = -O3 -ipo -ip -unroll
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused -g
switchPRO = $(switchOPT) -pg -traceback -g
switchOMP = $(switchOPT) -openmp

switch = $(switchOPT)


#objects that need fortran, mind the sequence
gnu_objs = opkda1.o opkda2.o opkdmain.o patmo_commons.o patmo_constants.o \
       patmo_parameters.o patmo_utils.o patmo_rates.o patmo_reverseRates.o \
       patmo_photo.o patmo_photoRates.o patmo_sparsity.o patmo_jacobian.o \
    	patmo_ode.o patmo.o test.o\
	   #patmo_write.o\
	   

# default target
all: $(gnu_objs)
	$(gnu_fc) $(switch) $(gnu_objs) -o $(exec)

# Rules for Fortran 77 files with legacy support
opkda1.o: opkda1.f
	$(fc_legacy) -c opkda1.f -o opkda1.o

opkda2.o: opkda2.f
	$(fc_legacy) -c opkda2.f -o opkda2.o

opkdmain.o: opkdmain.f
	$(fc_legacy) -c opkdmain.f -o opkdmain.o

debug: switch = $(switchDBG)
debug: nowarn = -nowarn
debug: all

profile: switch = $(switchPRO)
profile: nowarn = -nowarn
profile: all

omp: openmp
openmp: switch = $(switchOMP)
openmp: all

gfortran: fc = gfortran
gfortran: switch = -ffree-line-length-none
gfortran: all

gfortran_dbg: fc = gfortran
gfortran_dbg: switch = -fbacktrace -g
gfortran_dbg: switch += -ffpe-trap=zero,overflow,invalid
gfortran_dbg: switch += -fbounds-check -ffree-line-length-none
gfortran_dbg: all

coverage: fc = gfortran
coverage: switch = -ffree-line-length-none
coverage: switch += -fprofile-arcs -ftest-coverage
coverage: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f
%.o: %.f90
	$(gnu_fc) $(switch) -c $< -o $@

# Compile Fortran 77 files (legacy support)
%.o: %.f
	$(fc_legacy) -c $< -o $@