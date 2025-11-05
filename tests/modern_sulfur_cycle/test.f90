program test
  use patmo
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  use patmo_ode
  use patmo_utils
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass,one_year
  real*8::heff(chemSpeciesNumber)
  real*8::dep(chemSpeciesNumber) !cm/s
  integer::icell,i,j
  real*8::convergence = 100.0

  !init photochemistry
  call patmo_init()

  !load temperature and density profile
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")

  !set BB flux (default Sun@1AU)
  call patmo_setFluxBB()
  call patmo_setGravity(9.8d2)
  !call patmo_setEddyKzzAll(1.0d5)
  
  !read initial value 
  call patmo_dumpHydrostaticProfile("hydrostat.out")
  wetdep(:,:) = 0d0 

  !calculate wet deposition
  call computewetdep(1,2.0d-2)   !OCS
  call computewetdep(11,5.0d-2)  !CS2
	call computewetdep(28,1.0d-1)  !H2S
	call computewetdep(18,4.0d3)  !SO2
  ! call computewetdep(26,5d14)  !H2SO4
  call computewetdep(19,5d14)  !SO4

  !Turco H2SO4
  wetdep(12,26) = 1.77E-06
  wetdep(11,26) = 3.54E-06
  wetdep(10,26) = 5.31E-06
  wetdep(9,26) = 7.08E-06
  wetdep(8,26) = 8.85E-06
  wetdep(7,26) =  1.06E-05
  wetdep(6,26) = 1.24E-05
  wetdep(5,26) = 1.42E-05
  wetdep(4,26) = 1.59E-05
  wetdep(3,26) = 1.77E-05
  wetdep(2,26) = 1.95E-05
  wetdep(1,26) = 2.12E-05

  va(:) = 0d0  
  pa(:) = 0d0
  
  open(60,file="vapor_H2SO4.txt",status="old")  
    do i=13,34
      read(60,*) va(i)
	  end do
  close(60)

  open(61,file="partial_H2SO4.txt",status="old") 
    do i=13,34
      read(61,*) pa(i)
	  end do
  close(61)

  gd(:) = 0d0

  open(62,file="SO4_deposition_rate.txt",status="old")  
    do i=1,60
      read(62,*) gd(i)
	  end do
  close(62)

  !get initial mass, g/cm3
    imass = patmo_getTotalMass()
    print*,"mass:",imass
 
  !loop on time 
  dt = secondsPerDay 
  tend = secondsPerDay*365*10d0 
  one_year = secondsPerDay*365
  t = 0d0

  !loop on time
  do
    dt = dt
    call patmo_run(dt,convergence)
    t = t + dt
    if (t==tend .or. abs(convergence) < 1e-10) then
        
    call patmo_dumpDensityToFile(35,t,patmo_idx_COS)
 	call patmo_dumpDensityToFile(36,t,patmo_idx_S2)
	call patmo_dumpDensityToFile(37,t,patmo_idx_HSO)
	call patmo_dumpDensityToFile(38,t,patmo_idx_HSO2)
 	call patmo_dumpDensityToFile(39,t,patmo_idx_HSO3)
  	call patmo_dumpDensityToFile(40,t,patmo_idx_CS2)
    call patmo_dumpDensityToFile(41,t,patmo_idx_SO3)
    call patmo_dumpDensityToFile(42,t,patmo_idx_CH4O3S)
    call patmo_dumpDensityToFile(43,t,patmo_idx_SO4)
	call patmo_dumpDensityToFile(44,t,patmo_idx_S)
	call patmo_dumpDensityToFile(45,t,patmo_idx_SO2)
	call patmo_dumpDensityToFile(46,t,patmo_idx_SO4)
	call patmo_dumpDensityToFile(47,t,patmo_idx_CS)
	call patmo_dumpDensityToFile(48,t,patmo_idx_SCSOH)
	call patmo_dumpDensityToFile(49,t,patmo_idx_H2SO4)
	call patmo_dumpDensityToFile(50,t,patmo_idx_SO3)
	call patmo_dumpDensityToFile(51,t,patmo_idx_H2S)
	call patmo_dumpDensityToFile(52,t,patmo_idx_SH)
	call patmo_dumpDensityToFile(53,t,patmo_idx_SO)
	call patmo_dumpDensityToFile(55,t,patmo_idx_CH3SCH3)
    
    endif
 
  print '(F11.2,a2)',t/tend*1d2," %"
    if(t>=tend) exit
  end do
  
 !get mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass
  
  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")
  call patmo_dumpJValue("jvalue.dat")
  call patmo_dumpOpacity("opacity.dat")
  call patmo_dumpAllRates("rates.dat")
  call patmo_dumpAllMixingRatioToFile("allNDs.dat")

end program test
!**************

subroutine computewetdep(i,heff)
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::gamma(cellsNumber)  ! Precipation and Nonprecipitation time 
  real*8::wh2o(cellsNumber)   ! Rate of wet removal
  real*8::rkj(cellsNumber,chemSpeciesNumber)    ! Average Remeval Frequency	  
  real*8::y(cellsNumber),fz(cellsNumber),wl,qj(cellsNumber,chemSpeciesNumber)
	real*8::heff !Henry's Constant
  real*8::zkm(cellsNumber),temp(cellsNumber),gam15,gam8
	integer::i,j
	!Gas constant
  real*8,parameter::Rgas_atm = 1.36d-22 !Latm/K/mol

  gam15 = 8.64d5/2d0
  gam8 = 7.0d6/2d0
  wl = 1d0
  zkm(:) = height(:)/1d5
	temp(:) = TgasAll(:)
	 
  do j=1,12

	!FIND APPROPRIATE GAMMA
  if (zkm(j).LE.1.51d0) then
    gamma(j) = gam15
  else if (zkm(j).LT.8d0) then
    gamma(j) = gam15 + (gam8-gam15)*((zkm(j)-1.5d0)/6.5d0)
  else
    gamma(j) = gam8
  end if 

  !FIND WH2O
  if (zkm(j).LE.1d0) then
    y(j) = 11.35d0 + 1d-1*zkm(j)
  else
    y(j) = 11.5444d0 - 0.085333d0*zkm(j) - 9.1111d-03*zkm(j)*zkm(j)
  end if
  wh2o(j) = 10d0**y(j)

	!FIND F(Z)
  if (zkm(j).LE.1.51d0) then
     fz(j) = 1d-1
  else
     fz(j) = 0.16615d0 - 0.04916d0*zkm(j) + 3.37451d-3*zkm(j)*zkm(j)
  end if
		
	   !Raintout rates
        rkj(j,i) = wh2o(j)/55d0 /(av*wl*1.0d-9 + 1d0/(heff*Rgas_atm*temp(j)))
	    qj(j,i) = 1d0 - fz(j) + fz(j)/(gamma(j)*rkj(j,i)) * (1d0 - EXP(-rkj(j,i)*gamma(j)))
        wetdep(j,i) = (1d0 - EXP(-rkj(j,i)*gamma(j)))/(gamma(j)*qj(j,i))
	end do

	!Output 
   if (i==1) then
  	do j=1,12
      open  (20,file="Rainout-OCS.txt")
      write (20,*) 'GAMMA', gamma(j)
      write (20,*) 'ZKM', zkm(j)
      write (20,*) 'WH2O', wh2o(j)
	  write (20,*) 'RKJ', rkj(j,i)
 	  write (20,*) 'QJ', qj(j,i)
	  write (20,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if
	 	 
   if (i==11) then
	do j=1,12
      open  (25,file="Rainout-CS2.txt")
      write (25,*) 'ZKM', zkm(j)
	  write (25,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 
	 	 
   if (i==28) then
 	do j=1,12
      open  (26,file="Rainout-H2S.txt")
      write (26,*) 'ZKM', zkm(j)
	  write (26,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 	
   
   if (i==18) then
   	do j=1,12
      open  (24,file="Rainout-SO2.txt")
      write (24,*) 'ZKM', zkm(j)
	  write (24,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 
   
   if (i==26) then
   	do j=1,12
      open  (22,file="Rainout-H2SO4.txt")
      write (22,*) 'ZKM', zkm(j)
	  write (22,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 
   
   if (i==19) then
	do j=1,12
      open  (21,file="Rainout-SO4.txt")
      write (21,*) 'ZKM', zkm(j)
	  write (21,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if

   end subroutine computewetdep
	!**************


