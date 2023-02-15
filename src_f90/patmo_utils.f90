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
            * #PATMO_flux
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
#PATMO_weightedDegree
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

#PATMO_mass

  end function getSpeciesMass

  !*****************
  !get an array with all the species names
  function getSpeciesNames()
    use patmo_commons
    implicit none
    character(len=maxNameLength)::getSpeciesNames(speciesNumber)

#PATMO_speciesNames

  end function getSpeciesNames

#PATMO_massNucleiFunctions

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
