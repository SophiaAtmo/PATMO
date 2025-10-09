module patmo_reverseRates
contains

  !compute reverse rates using thermochemistry polynomials
  subroutine computeReverseRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(:)
    real*8::Tgas(cellsNumber)
    real*8::lnTgas(cellsNumber)
    real*8::Tgas2(cellsNumber)
    real*8::Tgas3(cellsNumber)
    real*8::Tgas4(cellsNumber)
    real*8::invTgas(cellsNumber)
    real*8::ntot(cellsNumber)
    integer::i

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

    !extrapolate lower and upper limits
    do i=1,cellsNumber
       Tgas(i) = max(inTgas(i),#PATMO_Tmin)
       Tgas(i) = min(Tgas(i),#PATMO_Tmax)
    end do

    lnTgas(:) = log(Tgas(:))
    Tgas2(:) = Tgas(:)**2
    Tgas3(:) = Tgas(:)**3
    Tgas4(:) = Tgas(:)**4
    invTgas(:) = 1d0/Tgas(:)

#PATMO_reverseRates
    do i=1,cellsNumber
      krate(i,473) = 1.57d13*(Tgas(i)/298)*exp(-4883/Tgas(i))
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates
