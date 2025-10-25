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

 !total density per layer
    M(:) = 0.5*sum(nAll(:,1:chemSpeciesNumber),2)
    !loop on cells
    do icell=1,cellsNumber
       Tgas = inTgas(icell)
       T = Tgas
       invT = 1d0/Tgas
#PATMO_rates
    end do
  
  end subroutine computeRates
   
end module patmo_rates
