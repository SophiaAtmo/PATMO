module patmo_thermo
contains

  ! ********************
  subroutine computeTgas(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber, cellsNumber)

    TgasAll(:) = TgasAll(:)

  end subroutine computeTgas

end module patmo_thermo
