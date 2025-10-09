module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,M(cellsNumber)
    integer::icell
    real*8::n(cellsNumber,speciesNumber)
#IFPATMO_has3body
    real*8::klow(cellsNumber),kinf(cellsNumber)
    real*8::Pr(cellsNumber),Fc(cellsNumber)
    real*8::y(cellsNumber),F(cellsNumber)
#ENDIFPATMO

    call computeTotalDensityExcludingM(n)
    !total density per layer
    M(:) = n(:, patmo_idx_M)

    !loop on cells
    do icell=1,cellsNumber
       Tgas = inTgas(icell)
       T = Tgas
       invT = 1d0/Tgas
#PATMO_rates
    end do
  
  end subroutine computeRates
   
  subroutine computeTotalDensityExcludingM(n)
    use patmo_commons
    use patmo_parameters
    implicit none

    real*8, intent(inout) :: n(:,:)  ! n(cellsNumber, speciesNumber)
    integer :: j

    n(:, patmo_idx_M) = 0.0d0
    do j = 1, chemSpeciesNumber
      if (j /= patmo_idx_M) then
        n(:, patmo_idx_M) = n(:, patmo_idx_M) + nAll(:, j)
      end if
    end do
  end subroutine computeTotalDensityExcludingM

end module patmo_rates
