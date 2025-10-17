module patmo
contains

  !*************
  !initialize all
  subroutine patmo_init()
    use patmo_photo
    use patmo_parameters
    use patmo_utils

#IFPATMO_usePhotochemistry
    !load photo metrics (i.e. binning)
    call loadPhotoMetric("xsecs/photoMetric.dat")
    !load photo cross-sections
    call loadAllPhotoXsecs()

    !init default photoflux
    photoFlux(:) = 0d0
#ENDIFPATMO

    !init cumulative flux for entropy production
    ! (i.e. integrated with time)
    cumulativeFlux(:,:) = 0d0

    !set reactions rate to zero by default
    krate(:,:) = 0d0

    !load verbatim reactions
    call loadReactionsVerbatim()

  end subroutine patmo_init

  !**************
  !run model for a time-step
  subroutine patmo_run(dt)
    !use dvode_f90_m
    use patmo_parameters
    use patmo_commons
    use patmo_ode
    use patmo_jacobian
    use patmo_sparsity
    use patmo_rates
    use patmo_photoRates
    use patmo_reverseRates
    use patmo_utils
    implicit none
    real*8,intent(in)::dt
    real*8::atol(neqAll),rtol(neqAll)
    real*8::tstart,tend,n(neqAll)
    real*8::sumxn(photoBinsNumber),m(speciesNumber)
    !type(VODE_OPTS)::OPTIONS
    integer::istate,itask,i,j

    integer,parameter::meth=2
    integer,parameter::lwm=2*neqAll**2 + 2*neqAll &
         + (neqAll**2+10*neqAll)/2
    integer::itol,iopt,mf,lrw,liw
    integer::iwork(20+9*neqAll+LWM),neq(1)
    real*8::rwork(20+neqAll*6+3*neqAll+lwm)

    lrw = size(rwork)
    liw = size(iwork)

    iwork(:) = 0
    rwork(:) = 0d0

    atol(:) = 1d-10 !absolute tolerances (array)
    rtol(:) = 1d-4 !relative tolerances (array)

    !computes sparsity if not already done
    if(nonZeroElements==0) then
       call computeSparsity()
    end if

    itol = 4
    istate = 1
    itask = 1
    iopt = 0
    MF = 222
    call xsetf(0)

    !set solver options (DVODE_f90)
    !OPTIONS = SET_OPTS(SPARSE_J=.true., ABSERR_VECTOR=ATOL(:), &
    !     RELERR_VECTOR=RTOL(:), MXSTEP=100000, &
    !     USER_SUPPLIED_SPARSITY = .true., &
    !     MA28_RPS = .true., &
    !     USER_SUPPLIED_JACOBIAN = .false.)

    !set the sparsity structure (DVODE_f90)
    !CALL USERSETS_IAJA(iaSparsity, size(iaSparsity), &
    !     jaSparsity, size(jaSparsity))

    tstart = 0d0
    tend = dt

#IFPATMO_use_opacity
    !upper layer opacity is zero
    tauAll(:,cellsNumber) = 0d0
    !loop on cells
    do j=cellsNumber-1,1,-1
       sumxn(:) = 0d0
       !loop on reactions
       do i=1,photoReactionsNumber
          sumxn(:) = sumxn(:) + xsecAll(:,i) * nall(j,photoPartnerIndex(i))
       end do
       tauAll(:,j) = tauAll(:,j+1) + gridSpace(j) * sumxn(:)
    end do
#ENDIFPATMO

    !unroll chemistry
    do i=1,speciesNumber
       n((i-1)*cellsNumber+1:(i*cellsNumber)) = nall(:,i)
    end do
    !unroll Tgas
    n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber)) &
         = TgasAll(:)

    !compute rates
    call computeRates(TgasAll(:))
    call computePhotoRates(tauAll(:,:))
    call computeReverseRates(TgasAll(:))
#IFPATMO_useHescape 
    call computeHescape() 
#ENDIFPATMO
    !compute tot density
    !ntotAll(:) = sum(nall(:,1:chemSpeciesNumber),2) !Trieu commented
    ntotAll(:) = 0.0d0
    do j = 1, chemSpeciesNumber
      if (j /= patmo_idx_M) then
       ntotAll(:) = ntotAll(:) + nAll(:, j)
      end if
    end do                                           !Trieu commented

    !compute mean molecular mass of the whole atmosphere
    ! (averaged between layers)
    m(:) = getSpeciesMass()
    meanMolecularMass = 0d0
    do i=1,cellsNumber
    !   meanMolecularMass = meanMolecularMass &       !Trieu commented
    !        + sum(m(1:chemSpeciesNumber) &
    !        * nAll(i,1:chemSpeciesNumber)) &
    !        / ntotAll(i) / cellsNumber
       do j = 1, chemSpeciesNumber
         if (j /= patmo_idx_M) then
           meanMolecularMass = meanMolecularMass + m(j) * nAll(i,j) / ntotAll(i)
         end if
       end do                                         
    end do

    ! Average across all cells
    meanMolecularMass = meanMolecularMass / cellsNumber  !Trieu commented
    
    !call the solver (DVODE_f90)
    !CALL DVODE_F90(fex, &
    !     neqAll, n(:), &
    !     tstart, tend, ITASK, ISTATE, OPTIONS, &
    !     jex)

    neq(:) = neqAll

    !loop until istate=2 or istate=error
    do
       CALL DLSODES(fex, NEQ(:), n(:), tstart, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
       print *, istate
       !recompute sparsity if required
       if(istate==-5.or.istate==-4) then
          istate = 3
          cycle
       end if
       !loop when max iteration reached
       if(istate/=-1) exit
       istate = 1
    end do

    !check output state
    if(istate/=2) then
       print *,"ERROR: istate=",istate
       stop
    end if

    !avoid negative species
    do i=1,neqAll
       n(i) = max(n(i),0d0)
    end do

    !roll chemistry
    do i=1,speciesNumber
       nall(:,i) = n((i-1)*cellsNumber+1:(i*cellsNumber))
    end do
    !roll Tgas
    TgasAll(:) = n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))

  end subroutine patmo_run

  !****************
  !dump the histogram of the connection degree to ifile, for the
  ! cell icell, at a given time (or any other independent varibles)
  subroutine patmo_dumpWeightedDegreeHistogram(ifile,icell,time)
    use patmo_utils
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::time
    integer,intent(in)::ifile,icell
    integer::hist(speciesNumber),i

    !get histogram
    hist(:) = getDegreeHistogram(nAll(:,:),icell)
    !write histogram to file
    do i=1,speciesNumber
       write(ifile,*) time, i, hist(i)
    end do
    write(ifile,*)

  end subroutine patmo_dumpWeightedDegreeHistogram

  !**************
  subroutine patmo_printBestFluxes(icell,bestFluxesNumber)
    use patmo_commons
    use patmo_utils
    use patmo_parameters
    implicit none
    integer,intent(in)::bestFluxesNumber,icell
    integer::idx(bestFluxesNumber),i
    real*8::flux(reactionsNumber)

    !get fluxes
    flux(:) = getFlux(nAll(:,:),icell)

    idx(:) = getBestFluxIdx(icell,bestFluxesNumber)
    print *,"*************"
    do i=1,bestFluxesNumber
       print *,idx(i),trim(reactionsVerbatim(idx(i))),flux(idx(i))
    end do

  end subroutine patmo_printBestFluxes

  !**************
  !compute cumulative flux for entropy production
  subroutine patmo_computeEntropyProductionFlux(dt)
    use patmo_utils
    implicit none
    real*8,intent(in)::dt

    call computeEntropyProductionFlux(dt)

  end subroutine patmo_computeEntropyProductionFlux

  !***************
  function patmo_getEntropyProduction(timeInterval)
    use patmo_utils
    implicit none
    real*8,intent(in)::timeInterval
    real*8::patmo_getEntropyProduction

    patmo_getEntropyProduction = getEntropyProduction(timeInterval)

  end function patmo_getEntropyProduction

  !***************
  !Assume a black-body flux, with starTbb (K), starRadius (Rsun)
  ! starDistance (AU). Default is Sun at 1AU
  subroutine patmo_setFluxBB(starTbb,starRadius,starDistance)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    use patmo_photo
    implicit none
    real*8,optional,intent(in)::starTbb,starRadius,starDistance
    real*8,parameter::AU2cm=1.496d13 !AU->cm
    real*8,parameter::Rsun2cm=6.963d10 !Rsun->cm
    real*8::Tbb,rstar,dstar
    integer::i

    !default is Sun
    Tbb = 5.777d3 !K
    rstar = Rsun2cm !cm
    dstar = AU2cm !cm

    !check optional parameters
    if(present(starTbb)) Tbb = starTbb
    if(present(starRadius)) rstar = starRadius*Rsun2cm
    if(present(starTbb)) dstar = starDistance*AU2cm

    !integrate flux
    do i=1,photoBinsNumber
       photoFlux(i) = (fluxBB(energyLeft(i),Tbb)+fluxBB(energyRight(i),Tbb)) &
            * energySpan(i)/2d0
    end do

    !scale geometric flux
    photoFlux(:) =  pi*rstar**2/dstar**2 * photoFlux(:)
    open(58,file="solar_flux.txt",status="old")
    do i=1,photoBinsNumber
        read(58,*) photoFlux(i)
    end do
    close(58)

  end subroutine patmo_setFluxBB

  !***************
  !dump opacity to file fname using unitEenergy as unit for
  ! energy (eV or mircron, eV default).
  ! File format is energy,layer,opacity
  subroutine patmo_dumpOpacity(fname,unitEnergy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    character(len=*),optional,intent(in)::unitEnergy
    character(len=*),intent(in)::fname
    character(len=100)::unitE
    integer::i,j

    unitE = "eV"
    if(present(unitEnergy)) unitE = trim(unitEnergy)

    open(22,file=trim(fname),status="replace")
    if(trim(unitE)=="eV") then
       !loop on energy
       do i=1,photoBinsNumber
          !loop on cells
          do j=1,cellsNumber
             write(22,*) energyMid(i),height(j)/1d5,tauAll(i,j)
          end do
          write(22,*)
       end do
    else if(trim(unitE)=="micron") then
       !loop on energy
       do i=photoBinsNumber,1,-1
          !loop on cells
          do j=1,cellsNumber
             !h*c/E -> cm -> micron
             write(22,*) 1d4*planck_eV*clight/energyMid(i),height(j)/1d5,&
                  tauAll(i,j)
          end do
          write(22,*)
       end do
    else
       print *,"ERROR: unknown unit "//trim(unitEnergy)
       stop
    end if
    close(22)

    print *, "Opacity dumped in ", trim(fname)

  end subroutine patmo_dumpOpacity

  !***************
  !find hydrostatic equilbrium knowing the pressure at ground
  ! (pground), using dp/dz = -mu*p*g/k/T
  ! Pressure unit is defined in unitP, default dyne
  subroutine patmo_hydrostaticEquilibrium(pground,unitP)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    real*8,intent(in)::pground
    real*8::p,zold,dz,ntot,n(speciesNumber)
    integer::i
    character(len=*),optional,intent(in)::unitP
    character(len=50)::units

    !optional argument units
    if(present(unitP)) then
       units = trim(unitP)
    else
       !default
       units = "dyne"
    end if

    !convert initial pressure to dyne/cm2 if necessary
    if(trim(units)=="dyne") then
       p = pground
    elseif(trim(units)=="dyne/cm2") then
       p = pground
    elseif(trim(units)=="atm") then
       p = pground*1.013250d6
    elseif(trim(units)=="mbar") then
       p = pground*1d3
    elseif(trim(units)=="bar") then
       p = pground*1d6
    else
       !error if units unknonw
       print *,"ERROR: unknown pressure unit for hydrostatic eq ",trim(units)
       stop
    end if

    !initial conditions
    zold = 0d0
    !loop on cells
    do i=1,cellsNumber
       !difference in height
       dz = height(i)-zold
       !temp array
       n(:) = nall(i,:)
       !compute difference in pressure
       p = p - getMeanMass(n(:)) * p / kboltzmann &
            / TgasAll(i) * gravity * dz
       !total number density p=n*k*T
       ntot = p/kboltzmann/TgasAll(i)
       !resacale abundances depending on total pressure
       nall(i,1:chemSpeciesNumber) = &
            n(1:chemSpeciesNumber) &
            / sum(n(1:chemSpeciesNumber))*ntot
       !store old height
       zold = height(i)
    end do

  end subroutine patmo_hydrostaticEquilibrium

  !***************
  !dump hydrostatic pressure profile to fname.
  ! Format: h/km, p/mbar, Tgas/K
  subroutine patmo_dumpHydrostaticProfile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer::i
    character(len=*),intent(in)::fname
    real*8::ntot

    open(22,file=trim(fname),status="replace")
    write(22,*) "#hydrostatic equilibrium dump"
    write(22,*) "#alt/km p/mbar Tgas/K"
    !loop on cells
    do i=1,cellsNumber
       ntot = sum(nall(i,1:chemSpeciesNumber))
       write(22,*) height(i)/1d5,ntot*kboltzmann*TgasAll(i)/1d3,TgasAll(i)
    end do
    close(33)
    print *,"Hydrostatic equilibrium dumped in ",trim(fname)

  end subroutine patmo_dumpHydrostaticProfile

  !****************
  !load initial profile (density, etc...) from fname.
  ! Height in unitH, species in unitX
  subroutine patmo_loadInitialProfile(fname,unitH,unitX,defaultDensity)
    use patmo_parameters
    use patmo_commons
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=*),optional,intent(in)::unitH,unitX
    character(len=50)::units,unitsX
    character(len=50),allocatable::arow(:)
    real*8,optional,intent(in)::defaultDensity
    real*8,allocatable::x(:),rout(:)
    real*8::zold,defaultN
    integer::ios,i,idx,j,nonZero,offset
    logical::firstRow

    units = "cm"
    !optional argument units height
    if(present(unitH)) units = trim(unitH)

    unitsX = "1/cm3"
    !optional argument units chemical species
    if(present(unitX)) unitsX = trim(unitX)

    defaultN = 0d0
    !optional argument for default density
    if(present(defaultDensity)) defaultN = defaultDensity

    !read file
    print *,"reading ",trim(fname)
    open(22,file=trim(fname),status="old",iostat=ios)
    !check for file opening
    if(ios/=0) then
       print *,"ERROR: problem while opening ",trim(fname)
       stop
    end if

    !read until comment is found
    do
       read(22,*,iostat=ios) offset,nonZero
       if(ios==0) exit
    end do
    allocate(arow(offset+nonZero))
    allocate(x(nonZero))
    allocate(rout(offset))
    read(22,*) arow(:)

    !set default abundance
    nall(:,1:chemSpeciesNumber) = defaultN
    !set not chemial species to zero
    nall(:,chemSpeciesNumber+1:speciesNumber) = 0d0
    !loop on cells (file lines have to be the same number)
    do j=1,cellsNumber
       !read data+chemistry
       read(22,*,iostat=ios) rout(:),x(:)
       if(ios/=0) then
          print *,"ERROR: problem while reading ",trim(fname)
          if(j>1) print *,&
               "(could be less file lines than declared lines number)"
          stop
       end if

       !loop on data accoding to header
       do i=1,offset
          if(trim(arow(i))=="alt") then
             height(j) = rout(i)
          elseif(trim(arow(i))=="Tgas") then
             TgasAll(j) = rout(i)
          elseif(trim(arow(i))=="Dzz") then
             diffusionDzz(j) = rout(i)
          elseif(trim(arow(i))=="Kzz") then
             eddyKzz(j) = rout(i)
          elseif(trim(arow(i))=="index") then
             continue
          elseif(trim(arow(i))=="idx") then
             continue
          elseif(trim(arow(i))=="dummy") then
             continue
          else
             print *,"ERROR: unknown header element: ",trim(arow(i))
             stop
          end if
       end do

       !load species into common array
       do i=1,nonZero
          idx = getSpeciesIndex(arow(offset+i),error=.false.)
          if(idx/=-1) nall(j,idx) = x(i)
       end do

       !convert units if necessary
       if(trim(unitsX)=="ppbv") then
          nall(j,:) = nall(j,:)/sum(nall(j,1:chemSpeciesNumber))
       elseif(trim(unitsX)=="1/cm3") then
          continue
       else
          print *,&
               "ERROR: unknown chemical abundance units while reading profile",&
               trim(fname),trim(units)
          stop
       end if

    end do
    close(22)
    deallocate(arow)
    deallocate(x)
    deallocate(rout)

    !convert units if necessary
    if(trim(units)=="km") then
       height(:) = height(:)*1d5
    elseif(trim(units)=="cm") then
       continue
    else
       print *,"ERROR: unknown units while reading profile", &
            trim(fname),trim(units)
       stop
    end if

    !store inverse grid space squared, 1/dz**2, and dz
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/(height(j)-zold)**2
       gridSpace(j) = (height(j)-zold)
       zold = height(j)
    end do

  end subroutine patmo_loadInitialProfile

  !****************
  !return total mass in g/cm3
  function patmo_getTotalMass()
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    integer::icell
    real*8::patmo_getTotalMass
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    patmo_getTotalMass = 0d0
    do icell=1,cellsNumber
       patmo_getTotalMass = patmo_getTotalMass &
            + sum(m(1:chemSpeciesNumber) &
            * nall(icell,1:chemSpeciesNumber))
    end do

  end function patmo_getTotalMass

#PATMO_massNucleiFunctions

  !***************
  !set uniform grid spacing, cm
 subroutine patmo_setGridSpacing(dz)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::dz
    real*8::zold
    integer::j

    gridSpace(:) = dz
    !store inverse grid space squared, 1/dz**2, and height
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/gridSpace(j)**2
       height(j) = zold
       zold = zold + gridSpace(j)
    end do

 end subroutine patmo_setGridSpacing

  !***************
  !set thermal diffusion
 subroutine patmo_setThermalDiffusion(alpha)
    use patmo_parameters
    implicit none
    real*8,intent(in)::alpha

    thermalDiffusionFactor = alpha

 end subroutine patmo_setThermalDiffusion

  !***************
  !set eddy Kzz coefficient of icell layer
 subroutine patmo_setEddyKzz(icell,kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz
    integer,intent(in)::icell

    eddyKzz(icell) = kzz

 end subroutine patmo_setEddyKzz

  !***************
  !set eddy Kzz, same for all layers
 subroutine patmo_setEddyKzzAll(kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz

    eddyKzz(:) = kzz

 end subroutine patmo_setEddyKzzAll

  !***************
  !set diffusion Dzz for layer icell
 subroutine patmo_setDiffusionDzz(icell,dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz
    integer,intent(in)::icell

    diffusionDzz(icell) = dzz

 end subroutine patmo_setDiffusionDzz

  !***************
  !set diffusion Dzz, same for all layers
 subroutine patmo_setDiffusionDzzAll(dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz

    diffusionDzz(:) = dzz

 end subroutine patmo_setDiffusionDzzAll

  !***************
  !append density of species idx to file number ifile
 subroutine patmo_dumpDensityToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx)
    end do
    write(ifile,*)

 end subroutine patmo_dumpDensityToFile

  !****************
  !dump all mixing rations to file (one column one species)
  ! first column is layer number
 subroutine patmo_dumpAllMixingRatioToFile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=500)::header
    character(len=maxNameLength)::names(speciesNumber)
    integer::i

    names(:) = getSpeciesNames()
    !prepare header (species names)
    header = "#layer"
    do i=1,chemSpeciesNumber
       header = trim(header)//" "//names(i)
    end do

    !open file to dump mixing ratios
    open(67,file=trim(fname),status="replace")
    !write file header (species names)
    write(67,*) trim(header)
    !write mixing ratios
    do i=1,cellsNumber
       write(67,'(I5,999E17.8e3)') i,nall(i,1:chemSpeciesNumber)
    end do
    close(67)

 end subroutine patmo_dumpAllMixingRatioToFile

  !***************
  !append mixing ration of species idx to file number ifile
 subroutine patmo_dumpMixingRatioToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx) &
            / sum(nall(i,1:chemSpeciesNumber))
    end do
    write(ifile,*)

 end subroutine patmo_dumpMixingRatioToFile

  !****************
  !set gravity in cm/s2
 subroutine patmo_setGravity(g)
    use patmo_parameters
    implicit none
    real*8,intent(in)::g

    gravity = g

 end subroutine patmo_setGravity

  !****************
  !set chemistry of layer icell
 subroutine patmo_setChemistry(icell,n)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::n(speciesNumber)

    nall(icell,:) = n(:)

 end subroutine patmo_setChemistry

  !****************
  !set the same chemistry for all the layers
 subroutine patmo_setChemistryAll(n)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::n(speciesNumber)
    integer::icell

    do icell=1,cellsNumber
       nall(icell,:) = n(:)
    end do

 end subroutine patmo_setChemistryAll

  !**************
  !set Tgas for layer icell
 subroutine patmo_setTgas(icell,Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::Tgas

    TgasAll(icell) = Tgas

 end subroutine patmo_setTgas

  !**************
  !set the same Tgas for all layers
 subroutine patmo_setTgasAll(Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::Tgas

    TgasAll(:) = Tgas

 end subroutine patmo_setTgasAll

  !**************
  !get density of species idx_species at layer icell
 function patmo_getDensity(icell,idx_species)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell,idx_species
    real*8::patmo_getDensity

    patmo_getDensity = nall(icell,idx_species)

 end function patmo_getDensity

  !**************
  !get Tgas of layer icell
 function patmo_getTgas(icell)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8::patmo_getTgas

    patmo_getTgas = TgasAll(icell)

 end function patmo_getTgas

  !**************
  !dump J-Values
 subroutine patmo_dumpJValue(fname)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    integer::i

    open(22,file=trim(fname),status="replace")
    write(22,*) "altitude/km#PATMO_JValueReactions"
    !loop on cells
    do i=1,cellsNumber
        #PATMO_JValues
    end do
    write(22,*)
    close(22)

 end subroutine patmo_dumpJValue

  !**************
  !dump all reaction rates
 subroutine patmo_dumpAllRates(fname)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_rates  !Trieu added
    implicit none
    character(len=*),intent(in)::fname
    integer::i
    real*8::n(cellsNumber,speciesNumber)  !Trieu added

    call computeTotalDensityExcludingM(n) !Trieu added
    nall(:,patmo_idx_M) = n(:, patmo_idx_M)  !Trieu added
    
    open(22,file=trim(fname),status="replace")
    write(22,*) "altitude/km#PATMO_DumpReactions"
    !loop on cells
    do i=1,cellsNumber
        #PATMO_DumpAllReactionRates
    end do
    write(22,*)
    close(22)

 end subroutine patmo_dumpAllRates

 subroutine patmo_dumpAllNumberDensityDifference(ifile,nb,na)
  use patmo_commons
  use patmo_parameters
  implicit none
  integer,intent(in)::ifile
  real*8,intent(in)::nb(neqAll), na(cellsNumber, speciesNumber)
  real*8::deltaNAll(cellsNumber, speciesNumber)
  integer::i,j

  !compute deference
  do i = 1, speciesNumber
      deltaNAll(:, i) = nb((i - 1) * cellsNumber + 1 : (i * cellsNumber)) - na(:, i)
  end do

  ! do j = 1, cellsNumber
  !     do i = 1, speciesNumber
  !         write(ifile, *) j, i, deltaNAll(j, i)
  !     end do
  ! end do
  do i = 1, speciesNumber
      write (ifile, *) i, deltaNAll(:, i)
  end do
  write(ifile,*)

 end subroutine patmo_dumpAllNumberDensityDifference

#IFPATMO_useHescape_dump
 subroutine computeHescape()
   use patmo_constants
   use patmo_parameters
   implicit none
    
   Hesc = 2.5d8 * &
            (nAll(cellsNumber, patmo_idx_H) &
            + 2d0 * nAll(cellsNumber, patmo_idx_H2) &
            + 2d0 * nAll(cellsNumber, patmo_idx_H2O) &
            + 4d0 * nAll(cellsNumber, patmo_idx_CH4)) &
            / sum(nAll(cellsNumber,1:chemSpeciesNumber))
   H2esc = 2.5d8 * &
            (nAll(cellsNumber, patmo_idx_H2) &
            + nAll(cellsNumber, patmo_idx_H2O)&
            + 2d0 * nAll(cellsNumber, patmo_idx_CH4))&
            / sum(nAll(cellsNumber,1:chemSpeciesNumber))
 end subroutine computeHescape

 subroutine patmo_dumpHescape(ifile, time)
   use patmo_parameters
   implicit none
   integer,intent(in)::ifile
   real*8, intent(in)::time
   write (ifile, *) time, Hesc, H2esc
 end subroutine patmo_dumpHescape
#ENDIFPATMO

end module patmo
