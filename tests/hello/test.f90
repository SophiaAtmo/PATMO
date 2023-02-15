program test
  use patmo
  use patmo_commons
  use patmo_constants
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass
  integer::icell

  !init photochemistry
  call patmo_init()

  !load temperature and density profile
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="ppbv")

  !set BB flux (default Sun@1AU)
  call patmo_setFluxBB()

  call patmo_setGravity(9.8d2)

  !compute hydrostatic equilibrium
  call patmo_hydrostaticEquilibrium(9.57d0,unitP="mbar")
  call patmo_dumpHydrostaticProfile("hydrostat.out")
  !set diffusion, same for every layer
  call patmo_setDiffusionDzzAll(1d5)

  !get initial mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass

  dt = secondsPerYear * 1d-3
  tend = secondsPerYear
  t = 0d0

  !loop on time
  do
     dt = dt * 1.1
     call patmo_run(dt)
     t = t + dt
     print *,t/tend,abs(patmo_getTotalMass()-imass)/imass
     call patmo_dumpOpacity("opacity.dat",unitEnergy="micron")
     call patmo_dumpMixingRatioToFile(33,t,patmo_idx_H2)
     call patmo_dumpMixingRatioToFile(34,t,patmo_idx_CH4)
     if(t>=tend) exit
  end do

  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")

  print *,"mass:",patmo_getTotalMass()


end program test
