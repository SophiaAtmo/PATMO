module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

#PATMO_photoRates

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber)
    integer::j

    !loop on cells (stride photobins)
    do j=1,cellsNumber
       integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:) &
            /energyMid(:)*energySpan(:)*exp(-tau(:,j))) / planck_eV
    end do

  end function integrateXsec

end module patmo_photoRates
