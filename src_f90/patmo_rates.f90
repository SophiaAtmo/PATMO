module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,ntot(cellsNumber),kmax
    integer::icell, i
#IFPATMO_has3body
    real*8::klow(cellsNumber),kinf(cellsNumber)
    real*8::Pr(cellsNumber),Fc(cellsNumber)
    real*8::y(cellsNumber),F(cellsNumber)
#ENDIFPATMO

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

    !loop on cells
    do icell=1,cellsNumber
       Tgas = inTgas(icell)
       T = Tgas
       invT = 1d0/Tgas
#PATMO_rates

!!$       kmax = 1d-2
!!$       if(icell==1) print '(2a7,2a17)', "icell", "irate", "Tgas", "krate"
!!$       do i=1,size(krate,2)
!!$          if(krate(icell,i)>kmax) then
!!$             print '(2I7,2E17.8e3)', i, icell, T, krate(icell,i)
!!$          end if
!!$       end do
    end do

  end subroutine computeRates

end module patmo_rates
