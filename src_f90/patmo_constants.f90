module patmo_constants
  implicit none
  !boltzmann constant, erg/K
  real*8,parameter::kboltzmann = 1.38064852d-16
  !boltzmann constant, eV/K
  real*8,parameter::kboltzmann_eV = 8.617332478d-5
  !seconds per year, s
  real*8,parameter::secondsPerYear = 365d0*24d0*3600d0
  !pi
  real*8,parameter::pi = 4d0*atan(1d0)
  !speed of light, cm/s
  real*8,parameter::clight = 2.99792458d10
  !plank, eV*s
  real*8,parameter::planck_eV = 4.135667516d-15


end module patmo_constants
