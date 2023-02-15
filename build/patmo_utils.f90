module patmo_utils
contains

  !NOTE
  !references for some of these functions are in
  !https://bitbucket.org/tgrassi/planatmo/wiki/user_functions

  !**************
  subroutine computeEntropyProductionFlux(dt)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::dt
    real*8::flux(reactionsNumber)
    integer::icell

    !loop on cells
    do icell=1,cellsNumber
      !get reaction fluxes
      flux(:) = getFlux(nAll(:,:),icell)
      cumulativeFlux(icell,:) = &
          cumulativeFlux(icell,:) + flux(:) * dt
    end do

  end subroutine computeEntropyProductionFlux

  !*************
  function getEntropyProduction(timeInterval)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::timeInterval
    real*8::getEntropyProduction,ds,dsp,ratio
    integer::i,icell,lowReverse(chemReactionsNumber)

    !small Reverse/Forward ratio flags (to rise warning message)
    lowReverse(:) = 0

    open(22,file="entropyReport.out",status="replace")
    write(22,*) "#cell_number  reaction_index   contribution_to_entropy   Rb/Rf"
    ds = 0d0
    !loop on cells
    do icell=1,cellsNumber
      !get reaction fluxes
      do i=1,chemReactionsNumber
        !check if Reverse/Forward ratio is small
        ratio = cumulativeFlux(icell,reverseReactionsOffset+i) &
            / cumulativeFlux(icell,i)
        if(ratio<1d-10) then
          lowReverse(i) = 1
        end if
        !sum entropy production
        dsp = (cumulativeFlux(icell,i) &
            - cumulativeFlux(icell,reverseReactionsOffset+i)) &
            * log(cumulativeFlux(icell,i) &
            / cumulativeFlux(icell,reverseReactionsOffset+i))
        ds = ds + dsp
        write(22,'(2I5,2E17.8e3)') icell,i,dsp/timeInterval,ratio
      end do
      write(22,*)
    end do
    close(22)

    !rise warning if Reverse/Forward ratio found
    if(sum(lowReverse)>0) then
      print *,"WARNING: suspiciously low Reverse/Forward flux ratio for rates:"
      do i=1,chemReactionsNumber
        if(lowReverse(i)==1) print *,i
      end do
      print *,"check entropyReport.out file"
    end if

    !divide by total time
    getEntropyProduction = ds / timeInterval

  end function getEntropyProduction

  !**************
  !get the best nbest reaction fluxes indexes
  function getBestFluxIdx(icell,nbest)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8::flux(reactionsNumber),maxFlux
    integer,intent(in)::icell,nbest
    integer::i,j,idxList(nbest),getBestFluxIdx(nbest)

    !get fluxes
    flux(:) = getFlux(nAll(:,:),icell)

    !loop on number of best flux required
    do j=1,nbest
      maxFlux = 0d0
      !loop on reactions to find max
      do i=1,reactionsNumber
        if(flux(i).ge.maxFlux) then
          idxList(j) = i
          maxFlux = flux(i)
        end if
      end do
      !set max flux to zero to find next best
      flux(idxList(j)) = 0d0
    end do

    getBestFluxIdx(:) = idxList(:)

  end function getBestFluxIdx

  !**************
  !get reaction fluxes for cell icell
  function getFlux(nin,icell)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::nin(cellsNumber,speciesNumber)
    integer,intent(in)::icell
    integer::i
    real*8::n(speciesNumber)
    real*8::getFlux(reactionsNumber)

    !local copy for the given cell
    n(:) = nin(icell,:)
    n(positionDummy) = 1d0

    !loop on reactions
    do i=1,reactionsNumber
      getFlux(i) = krate(icell,i) &
          * n(indexReactants1(i)) &
          * n(indexReactants2(i))
    end do

  end function getFlux

  !***************
  !get the degree weighted by flux cut at a fraction
  ! of the maximum
  function getDegreeFluxWeighted(n,icell,fraction)
    use patmo_commons
    implicit none
    real*8,intent(in)::n(cellsNumber,speciesNumber),fraction
    integer,intent(in)::icell
    integer::degree(speciesNumber),i
    integer::getDegreeFluxWeighted(speciesNumber)
    real*8::flux(reactionsNumber),fluxFrac

    !get flux
    flux(:) = getFlux(n(:,:),icell)
    !fraction of the maximum flux
    fluxFrac = maxval(flux(:))*fraction

    degree(:) = 0
    do i=1,reactionsNumber
      if(flux(i)>=fluxFrac) then
        degree(indexReactants1(i)) = degree(indexReactants1(i)) + 1
        degree(indexReactants2(i)) = degree(indexReactants2(i)) + 1
        degree(indexProducts1(i)) = degree(indexProducts1(i)) + 1
        degree(indexProducts2(i)) = degree(indexProducts2(i)) + 1
        degree(indexProducts3(i)) = degree(indexProducts3(i)) + 1

      end if
    end do

    getDegreeFluxWeighted(:) = degree(:)

  end function getDegreeFluxWeighted

  !********************
  function getDegreeHistogram(n,icell)
    use patmo_commons
    implicit none
    real*8,intent(in)::n(cellsNumber,speciesNumber)
    integer,intent(in)::icell
    integer::i,degree(speciesNumber)
    integer::hist(speciesNumber)
    integer::getDegreeHistogram(speciesNumber)
    real*8::fraction

    fraction = 1d-1

    degree(:) = getDegreeFluxWeighted(n(:,:),icell,fraction)
    hist(:) = 0
    do i=1,speciesNumber
      hist(degree(i)) = hist(degree(i)) + 1
    end do

    getDegreeHistogram(:) = hist(:)

  end function getDegreeHistogram

  !**************
  !get an array with all the masses
  function getSpeciesMass()
    use patmo_commons
    implicit none
    real*8::getspeciesMass(speciesNumber)

    getSpeciesMass(patmo_idx_COS) = 9.040424d-23
    getSpeciesMass(patmo_idx_S2) = 1.071507d-22
    getSpeciesMass(patmo_idx_O3) = 8.036305d-23
    getSpeciesMass(patmo_idx_O2) = 5.357536d-23
    getSpeciesMass(patmo_idx_NO2) = 7.701459d-23
    getSpeciesMass(patmo_idx_NO) = 5.02269d-23
    getSpeciesMass(patmo_idx_HSO) = 8.203658d-23
    getSpeciesMass(patmo_idx_HO2) = 5.52489d-23
    getSpeciesMass(patmo_idx_HSO2) = 1.088243d-22
    getSpeciesMass(patmo_idx_HSO3) = 1.356119d-22
    getSpeciesMass(patmo_idx_CS2) = 1.171919d-22
    getSpeciesMass(patmo_idx_CH4O3S) = 1.506737d-22
    getSpeciesMass(patmo_idx_CO) = 3.682888d-23
    getSpeciesMass(patmo_idx_H) = 1.673533d-24
    getSpeciesMass(patmo_idx_O) = 2.678768d-23
    getSpeciesMass(patmo_idx_N) = 2.343922d-23
    getSpeciesMass(patmo_idx_S) = 5.357536d-23
    getSpeciesMass(patmo_idx_SO2) = 1.071507d-22
    getSpeciesMass(patmo_idx_SO4) = 1.607261d-22
    getSpeciesMass(patmo_idx_CS) = 6.361656d-23
    getSpeciesMass(patmo_idx_N2) = 4.687844d-23
    getSpeciesMass(patmo_idx_CO2) = 6.361656d-23
    getSpeciesMass(patmo_idx_SCSOH) = 1.456531d-22
    getSpeciesMass(patmo_idx_OH) = 2.846122d-23
    getSpeciesMass(patmo_idx_H2) = 3.347066d-24
    getSpeciesMass(patmo_idx_H2SO4) = 1.640732d-22
    getSpeciesMass(patmo_idx_SO3) = 1.339384d-22
    getSpeciesMass(patmo_idx_H2S) = 5.692243d-23
    getSpeciesMass(patmo_idx_H2O) = 3.013475d-23
    getSpeciesMass(patmo_idx_SH) = 5.52489d-23
    getSpeciesMass(patmo_idx_SO) = 8.036305d-23
    getSpeciesMass(patmo_idx_CS2E) = 1.272415d-22
    getSpeciesMass(patmo_idx_CH3SCH3) = 8.369896d-23

  end function getSpeciesMass

  !*****************
  !get an array with all the species names
  function getSpeciesNames()
    use patmo_commons
    implicit none
    character(len=maxNameLength)::getSpeciesNames(speciesNumber)

    getSpeciesNames(patmo_idx_COS) = "COS"
    getSpeciesNames(patmo_idx_S2) = "S2"
    getSpeciesNames(patmo_idx_O3) = "O3"
    getSpeciesNames(patmo_idx_O2) = "O2"
    getSpeciesNames(patmo_idx_NO2) = "NO2"
    getSpeciesNames(patmo_idx_NO) = "NO"
    getSpeciesNames(patmo_idx_HSO) = "HSO"
    getSpeciesNames(patmo_idx_HO2) = "HO2"
    getSpeciesNames(patmo_idx_HSO2) = "HSO2"
    getSpeciesNames(patmo_idx_HSO3) = "HSO3"
    getSpeciesNames(patmo_idx_CS2) = "CS2"
    getSpeciesNames(patmo_idx_CH4O3S) = "CH4O3S"
    getSpeciesNames(patmo_idx_CO) = "CO"
    getSpeciesNames(patmo_idx_H) = "H"
    getSpeciesNames(patmo_idx_O) = "O"
    getSpeciesNames(patmo_idx_N) = "N"
    getSpeciesNames(patmo_idx_S) = "S"
    getSpeciesNames(patmo_idx_SO2) = "SO2"
    getSpeciesNames(patmo_idx_SO4) = "SO4"
    getSpeciesNames(patmo_idx_CS) = "CS"
    getSpeciesNames(patmo_idx_N2) = "N2"
    getSpeciesNames(patmo_idx_CO2) = "CO2"
    getSpeciesNames(patmo_idx_SCSOH) = "SCSOH"
    getSpeciesNames(patmo_idx_OH) = "OH"
    getSpeciesNames(patmo_idx_H2) = "H2"
    getSpeciesNames(patmo_idx_H2SO4) = "H2SO4"
    getSpeciesNames(patmo_idx_SO3) = "SO3"
    getSpeciesNames(patmo_idx_H2S) = "H2S"
    getSpeciesNames(patmo_idx_H2O) = "H2O"
    getSpeciesNames(patmo_idx_SH) = "SH"
    getSpeciesNames(patmo_idx_SO) = "SO"
    getSpeciesNames(patmo_idx_CS2E) = "CS2E"
    getSpeciesNames(patmo_idx_CH3SCH3) = "CH3SCH3"

  end function getSpeciesNames

  !***************************
  function getTotalMassNuclei_C()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::icell
    real*8::getTotalMassNuclei_C
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    getTotalMassNuclei_C = 0d0

    do icell=1,cellsNumber
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_COS) * nall(icell,patmo_idx_COS)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_CS2) * nall(icell,patmo_idx_CS2)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_CH4O3S) * nall(icell,patmo_idx_CH4O3S)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_CO) * nall(icell,patmo_idx_CO)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_CS) * nall(icell,patmo_idx_CS)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_CO2) * nall(icell,patmo_idx_CO2)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_SCSOH) * nall(icell,patmo_idx_SCSOH)
      getTotalMassNuclei_C = getTotalMassNuclei_C + m(patmo_idx_CH3SCH3) * nall(icell,patmo_idx_CH3SCH3)
    end do

  end function

  !***************************
  function getTotalMassNuclei_H()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::icell
    real*8::getTotalMassNuclei_H
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    getTotalMassNuclei_H = 0d0

    do icell=1,cellsNumber
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_HSO) * nall(icell,patmo_idx_HSO)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_HO2) * nall(icell,patmo_idx_HO2)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_HSO2) * nall(icell,patmo_idx_HSO2)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_HSO3) * nall(icell,patmo_idx_HSO3)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_CH4O3S) * nall(icell,patmo_idx_CH4O3S) * 4d0
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_H) * nall(icell,patmo_idx_H)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_SCSOH) * nall(icell,patmo_idx_SCSOH)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_OH) * nall(icell,patmo_idx_OH)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_H2) * nall(icell,patmo_idx_H2) * 2d0
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_H2SO4) * nall(icell,patmo_idx_H2SO4) * 2d0
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_H2S) * nall(icell,patmo_idx_H2S) * 2d0
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_H2O) * nall(icell,patmo_idx_H2O) * 2d0
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_SH) * nall(icell,patmo_idx_SH)
      getTotalMassNuclei_H = getTotalMassNuclei_H + m(patmo_idx_CH3SCH3) * nall(icell,patmo_idx_CH3SCH3) * 3d0
    end do

  end function

  !***************************
  function getTotalMassNuclei_O()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::icell
    real*8::getTotalMassNuclei_O
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    getTotalMassNuclei_O = 0d0

    do icell=1,cellsNumber
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_COS) * nall(icell,patmo_idx_COS)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_O3) * nall(icell,patmo_idx_O3) * 3d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_O2) * nall(icell,patmo_idx_O2) * 2d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_NO2) * nall(icell,patmo_idx_NO2) * 2d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_NO) * nall(icell,patmo_idx_NO)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_HSO) * nall(icell,patmo_idx_HSO)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_HO2) * nall(icell,patmo_idx_HO2) * 2d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_HSO2) * nall(icell,patmo_idx_HSO2) * 2d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_HSO3) * nall(icell,patmo_idx_HSO3) * 3d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_CH4O3S) * nall(icell,patmo_idx_CH4O3S) * 3d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_CO) * nall(icell,patmo_idx_CO)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_O) * nall(icell,patmo_idx_O)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_SO2) * nall(icell,patmo_idx_SO2) * 2d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_SO4) * nall(icell,patmo_idx_SO4) * 4d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_CO2) * nall(icell,patmo_idx_CO2) * 2d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_SCSOH) * nall(icell,patmo_idx_SCSOH)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_OH) * nall(icell,patmo_idx_OH)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_H2SO4) * nall(icell,patmo_idx_H2SO4) * 4d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_SO3) * nall(icell,patmo_idx_SO3) * 3d0
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_H2O) * nall(icell,patmo_idx_H2O)
      getTotalMassNuclei_O = getTotalMassNuclei_O + m(patmo_idx_SO) * nall(icell,patmo_idx_SO)
    end do

  end function

  !***************************
  function getTotalMassNuclei_N()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::icell
    real*8::getTotalMassNuclei_N
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    getTotalMassNuclei_N = 0d0

    do icell=1,cellsNumber
      getTotalMassNuclei_N = getTotalMassNuclei_N + m(patmo_idx_NO2) * nall(icell,patmo_idx_NO2)
      getTotalMassNuclei_N = getTotalMassNuclei_N + m(patmo_idx_NO) * nall(icell,patmo_idx_NO)
      getTotalMassNuclei_N = getTotalMassNuclei_N + m(patmo_idx_N) * nall(icell,patmo_idx_N)
      getTotalMassNuclei_N = getTotalMassNuclei_N + m(patmo_idx_N2) * nall(icell,patmo_idx_N2) * 2d0
    end do

  end function

  !***************************
  function getTotalMassNuclei_S()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::icell
    real*8::getTotalMassNuclei_S
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    getTotalMassNuclei_S = 0d0

    do icell=1,cellsNumber
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_COS) * nall(icell,patmo_idx_COS)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_S2) * nall(icell,patmo_idx_S2) * 2d0
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_HSO) * nall(icell,patmo_idx_HSO)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_HSO2) * nall(icell,patmo_idx_HSO2)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_HSO3) * nall(icell,patmo_idx_HSO3)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_CS2) * nall(icell,patmo_idx_CS2) * 2d0
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_CH4O3S) * nall(icell,patmo_idx_CH4O3S)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_S) * nall(icell,patmo_idx_S)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_SO2) * nall(icell,patmo_idx_SO2)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_SO4) * nall(icell,patmo_idx_SO4)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_CS) * nall(icell,patmo_idx_CS)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_SCSOH) * nall(icell,patmo_idx_SCSOH)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_H2SO4) * nall(icell,patmo_idx_H2SO4)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_SO3) * nall(icell,patmo_idx_SO3)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_H2S) * nall(icell,patmo_idx_H2S)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_SH) * nall(icell,patmo_idx_SH)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_SO) * nall(icell,patmo_idx_SO)
      getTotalMassNuclei_S = getTotalMassNuclei_S + m(patmo_idx_CH3SCH3) * nall(icell,patmo_idx_CH3SCH3)
    end do

  end function

  !***************************
  function getTotalMassNuclei_CS2E()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::icell
    real*8::getTotalMassNuclei_CS2E
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    getTotalMassNuclei_CS2E = 0d0

    do icell=1,cellsNumber
      getTotalMassNuclei_CS2E = getTotalMassNuclei_CS2E + m(patmo_idx_CS2E) * nall(icell,patmo_idx_CS2E)
    end do

  end function

  !*******************
  !returns an array with the list of the reactions verbatims
  function getReactionsVerbatim()
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=maxNameLength)::getReactionsVerbatim(reactionsNumber)

    getReactionsVerbatim(:) = reactionsVerbatim(:)

  end function getReactionsVerbatim

  !*****************
  !load reaction verbatim from file
  subroutine loadReactionsVerbatim()
    use patmo_commons
    use patmo_parameters
    implicit none
    integer::ios,i
    character(len=maxNameLength)::cout

    open(33,file="reactionsVerbatim.dat",iostat=ios,status="old")
    if(ios/=0) then
      print *,"ERROR: problem load in verbatim reaction names"
      stop
    end if

    !loop on file to read verbatims
    do i=1,reactionsNumber
      read(33,'(a)',iostat=ios) cout
      reactionsVerbatim(i) = trim(cout)
    end do
    close(33)

  end subroutine loadReactionsVerbatim

  !*****************
  !given a species name returns the index
  function getSpeciesIndex(name,error)
    use patmo_commons
    implicit none
    integer::getSpeciesIndex,i
    logical,optional,intent(in)::error
    logical::riseError
    character(len=*),intent(in)::name
    character(len=maxNameLength)::names(speciesNumber)

    if(present(error)) then
      riseError = error
    else
      riseError = .true.
    end if

    names(:) = getSpeciesNames()

    !loop on species names to find index
    do i=1,speciesNumber
      if(trim(name)==trim(names(i))) then
        getSpeciesIndex = i
        return
      end if
    end do

    if(riseError) then
      print *,"ERROR: index not found for species: ",trim(name)
      print *," Available species are:"
      do i=1,speciesNumber
        print *," "//trim(names(i))
      end do
      stop
    else
      getSpeciesIndex = -1
    end if

  end function getSpeciesIndex

  !****************
  !mean molecular weight in g
  function getMeanMass(n)
    use patmo_commons
    implicit none
    real*8::getMeanMass,m(speciesNumber)
    real*8,intent(in)::n(speciesNumber)

    m(:) = getSpeciesMass()
    getMeanMass = sum(n(1:chemSpeciesNumber) &
        * m(1:chemSpeciesNumber)) &
        / sum(n(1:chemSpeciesNumber))

  end function getMeanMass

end module patmo_utils
