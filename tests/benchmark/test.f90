program test
  use patmo
  use patmo_commons
  use patmo_constants
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass
  integer::icell

  !init photochemistry
  call patmo_init()

  call patmo_setEddyKzzAll(1d5)
  call patmo_setGravity(9.8d2)
  call patmo_setDiffusionDzzAll(0d0)
  call patmo_setTgasAll(5d2)
  call patmo_setGridSpacing(1d6)

  x(:) = 0d0
  x(patmo_idx_H2O) = 6.0618d-4
  x(patmo_idx_CH4) = 2.7761d-4
  x(patmo_idx_He) = 0.09691d0
  x(patmo_idx_H2) = 1d0-2d0*x(patmo_idx_H2O)-4d0*x(patmo_idx_CH4)

  !normalize
  x(:) = x(:)/sum(x)

  call patmo_setChemistryAll(x(:))

  !compute hydrostatic equilibrium
  call patmo_hydrostaticEquilibrium(1d0,unitP="bar")
  call patmo_dumpHydrostaticProfile("hydrostat.out")


  !get initial mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass

  dt = secondsPerYear * 1d-3
  tend = secondsPerYear*1d6
  t = 0d0

  !loop on time
  do
     call patmo_printBestFluxes(1,5)
     dt = dt * 1.1
     call patmo_run(dt)

     t = t + dt
     print *,t/tend,t/secondsPerYear,abs(patmo_getTotalMass()-imass)/imass
     call patmo_dumpMixingRatioToFile(33,t,patmo_idx_H2)
     call patmo_dumpMixingRatioToFile(34,t,patmo_idx_CH4)
     if(t>=tend) exit
  end do

  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")

  print *,"mass:",patmo_getTotalMass()


end program test
