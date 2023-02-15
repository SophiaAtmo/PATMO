#default fortran compiler
fc = ifort

#executable name
exec = test

#flags
switchOPT = -O3 -ipo -ip -unroll
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused -g
switchPRO = $(switchOPT) -pg -traceback -g
switchOMP = $(switchOPT) -openmp

switch = $(switchOPT)

#objects
objs = opkda1.o
objs += opkda2.o
objs += opkdmain.o
objs += patmo_commons.o
objs += patmo_constants.o
objs += patmo_parameters.o
objs += patmo_utils.o
objs += patmo_rates.o
objs += patmo_reverseRates.o
objs += patmo_photo.o
objs += patmo_photoRates.o
objs += patmo_sparsity.o
objs += patmo_jacobian.o
objs += patmo_ode.o
objs += patmo.o
objs += test.o

#default target
all: 	$(objs)
	$(fc) $(objs) -o $(exec) $(switch)

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

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) -c $^ -o $@ $(nowarn)