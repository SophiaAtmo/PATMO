program test
  use patmo
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass
  real*8::convergence = 100.0
  integer::icell

  !init photochemistry
  call patmo_init()

  !load temperature and density profile
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")

  !set BB flux (default Sun@1AU)
  call patmo_setFluxBB()

  call patmo_setGravity(9.8d2)

  !compute hydrostatic equilibrium
  ! call patmo_hydrostaticEquilibrium(9.57d0,unitP="mbar")
  call patmo_dumpHydrostaticProfile("hydrostat.out")
  !set diffusion, same for every layer
  ! call patmo_setDiffusionDzzAll(1d5)
  ! call patmo_setEddyKzzAll(1.0d5)

  ! Wet Deposition
  wetdep(:,:) = 0d0

  ! Turco H2SO4 Wet Deposition
 ! wetdep(12,patmo_idx_H2SO4) = 1.77d-6
 ! wetdep(11,patmo_idx_H2SO4) = 3.54d-6
 ! wetdep(10,patmo_idx_H2SO4) = 5.31d-6
 ! wetdep(9,patmo_idx_H2SO4) = 7.08d-6
 ! wetdep(8,patmo_idx_H2SO4) = 8.85d-6
 ! wetdep(7,patmo_idx_H2SO4) = 1.06d-5
 ! wetdep(6,patmo_idx_H2SO4) = 1.24d-5
 ! wetdep(5,patmo_idx_H2SO4) = 1.42d-5
 ! wetdep(4,patmo_idx_H2SO4) = 1.59d-5
 ! wetdep(3,patmo_idx_H2SO4) = 1.77d-5
 ! wetdep(2,patmo_idx_H2SO4) = 1.95d-5
 ! wetdep(1,patmo_idx_H2SO4) = 2.12d-5

  call computeSulfuricAcidAerosol()

  !get initial mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass

  dt = secondsPerDay
  tend = secondsPerYear * 10d1
  t = 0d0

  !loop on time
  do
     dt = dt
     !call computeWetDep(patmo_idx_COS, 2.0d-2, .FALSE.)
     !call computeWetDep(patmo_idx_CS2, 5.0d-2, .FALSE.)
     !call computeWetDep(patmo_idx_H2S, 1.0d-1, .FALSE.)
     !call computeWetDep(patmo_idx_SO2, 4.0d3, .FALSE.)
     !call computeWetDep(patmo_idx_SO4, 5d14, .FALSE.)
     !call computeWetDep(patmo_idx_H2O2, 2d5, .FALSE.)
     !call computeWetDep(patmo_idx_HO2, 3.3d4, .FALSE.)
     !call computeWetDep(patmo_idx_HNO3, 7d11, .FALSE.)
     !call computeWetDep(patmo_idx_CH2O, 1.3d4, .FALSE.)
     call computeH2OCondensation((1.8d1 / Av), 1d0, 1d-4, .FALSE.)
     call patmo_run(dt,convergence)
     t = t + dt
     if(t==tend .or. abs(convergence) < 1d-10) then
       !  print *,t/tend,abs(patmo_getTotalMass()-imass)/imass
       !call patmo_dumpDensityToFile(30, t, patmo_idx_CH4)
       !call patmo_dumpDensityToFile(31, t, patmo_idx_H2)
       !call patmo_dumpDensityToFile(32, t, patmo_idx_H2O)
       !call patmo_dumpDensityToFile(33, t, patmo_idx_O2)
       !call patmo_dumpDensityToFile(34, t, patmo_idx_O3)
       !call patmo_dumpDensityToFile(35, t, patmo_idx_O)
       !call patmo_dumpDensityToFile(36, t, patmo_idx_O_1D)
       !call patmo_dumpDensityToFile(37, t, patmo_idx_NO)
       !call patmo_dumpDensityToFile(38, t, patmo_idx_NO2)
       !call patmo_dumpDensityToFile(39, t, patmo_idx_N2O)
       !call patmo_dumpDensityToFile(40, t, patmo_idx_NO3)
       !call patmo_dumpDensityToFile(41, t, patmo_idx_N2O5)
       !call patmo_dumpDensityToFile(42, t, patmo_idx_H)
       !call patmo_dumpDensityToFile(43, t, patmo_idx_OH)
       !call patmo_dumpDensityToFile(44, t, patmo_idx_HO2)
       !call patmo_dumpDensityToFile(45, t, patmo_idx_COS)
       !call patmo_dumpDensityToFile(46, t, patmo_idx_SO2)
       !call patmo_dumpDensityToFile(47, t, patmo_idx_H2SO4)
       !call patmo_dumpDensityToFile(48, t, patmo_idx_CS2)
       !call patmo_dumpDensityToFile(49, t, patmo_idx_H2S)
       !call patmo_dumpDensityToFile(50, t, patmo_idx_CH3SCH3)
       !call patmo_dumpDensityToFile(51, t, patmo_idx_SO3)
       !call patmo_dumpDensityToFile(52, t, patmo_idx_SO)
       !call patmo_dumpDensityToFile(53, t, patmo_idx_SO4)
       !call patmo_dumpDensityToFile(54, t, patmo_idx_H2O2)
       !call patmo_dumpDensityToFile(55, t, patmo_idx_CO)
       !call patmo_dumpDensityToFile(56, t, patmo_idx_CH3)
       !call patmo_dumpDensityToFile(57, t, patmo_idx_CH3O2)
       !call patmo_dumpDensityToFile(58, t, patmo_idx_CH3OOH)
       !call patmo_dumpDensityToFile(59, t, patmo_idx_CH2O)
       !call patmo_dumpDensityToFile(60, t, patmo_idx_CHO)
       !call patmo_dumpDensityToFile(61, t, patmo_idx_S)
       !call patmo_dumpDensityToFile(62, t, patmo_idx_S2)
       !call patmo_dumpDensityToFile(63, t, patmo_idx_S3)
       !call patmo_dumpDensityToFile(64, t, patmo_idx_S4)
       !call patmo_dumpDensityToFile(65, t, patmo_idx_S8)
       !call patmo_dumpDensityToFile(66, t, patmo_idx_CH2)
       !call patmo_dumpDensityToFile(67, t, patmo_idx_CH)
       !call patmo_dumpDensityToFile(68, t, patmo_idx_NH3)
       !call patmo_dumpDensityToFile(69, t, patmo_idx_NH2)
       !call patmo_dumpDensityToFile(70, t, patmo_idx_NH)
       !call patmo_dumpDensityToFile(71, t, patmo_idx_N)
       !call patmo_dumpDensityToFile(72, t, patmo_idx_N2H4)
       !call patmo_dumpDensityToFile(73, t, patmo_idx_N2H3)
       !call patmo_dumpDensityToFile(74, t, patmo_idx_CO2)
       !call patmo_dumpDensityToFile(75, t, patmo_idx_HOCO)
     endif
     
     print '(F20.1, F20.2, A1, F20.10)', t, t / tend * 1d2, "%", (patmo_getTotalMass()-imass)/imass
     if(t>=tend) exit
  end do

  call patmo_dumpOpacity("opacity.dat")
  call patmo_dumpJValue("jvalue.dat")
  call patmo_dumpAllRates("rates.dat")
  call patmo_dumpAllMixingRatioToFile("allNDs.dat")
  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")

  !call computeWetDep(patmo_idx_COS, 2.0d-2, .TRUE.)
  !call computeWetDep(patmo_idx_CS2, 5.0d-2, .TRUE.)
  !call computeWetDep(patmo_idx_H2S, 1.0d-1, .TRUE.)
  !call computeWetDep(patmo_idx_SO2, 4.0d3, .TRUE.)
  !call computeWetDep(patmo_idx_SO4, 5d14, .TRUE.)
  !call computeWetDep(patmo_idx_H2O2, 2d5, .TRUE.)
  !call computeWetDep(patmo_idx_HO2, 3.3d4, .TRUE.)
  !call computeWetDep(patmo_idx_HNO3, 7d11, .TRUE.)
  !call computeWetDep(patmo_idx_CH2O, 1.3d4, .TRUE.)
  call computeH2OCondensation((1.8d1 / Av), 1d0, 1d-4, .True.)

  print *,"mass:",patmo_getTotalMass()


end program test

subroutine computeWetDep(idx, Heff, isWrite)

    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    integer, parameter :: tropos = 12
    real(8), intent(in) :: Heff
    integer, intent(in) :: idx
    logical, intent(in) :: isWrite
    real(8), parameter :: R = 1.36d-22
    integer :: icell
    real(8) :: z(cellsNumber)
    real(8) :: WH2O(cellsNumber), L(cellsNumber), rkj(cellsNumber), fz(cellsNumber), Qj(cellsNumber), gamma(cellsNumber)
    real(8), parameter :: boundBottom = 1.51d0, boundTop = 8d0
    real(8), parameter :: gammaBottom = 8.65d5/2d0, gammaTop = 7.0d6/2d0
    real(8), parameter :: fzBottom = 1d-1, fzTop = 1.24d-2

    ! INITIALIZE
    z(:)        = height(:) / 1d5
    WH2O(:)     = 0d0
    L(:)        = 0d0
    rkj(:)      = 0d0
    fz(:)       = 0d0
    Qj(:)       = 0d0
    gamma(:)    = 0d0

    do icell = 1, tropos

        ! GAMMA
        if (z(icell) <= boundBottom) then
            gamma(icell) = gammaBottom
        else if (z(icell) >= boundTop) then
            gamma(icell) = gammaTop
        else
            gamma(icell) = &
                  gammaBottom &
                + ((gammaTop - gammaBottom) / (boundTop - boundBottom)) &
                * (z(icell) - boundBottom)
        end if

        ! f(z)
        if (z(icell) <= boundBottom) then
            fz(icell) = fzBottom
        else if (z(icell) >= boundTop) then
            fz(icell) = fzTop
        else
            ! ###PENDING###
            ! fz(icell) = 0.16615d0 - 0.04916d0 * z(icell) + 3.37451d-3 * z(icell) * z(icell)
            fz(icell) = 0.1284d0 * z(icell) ** (-1.128d0)
        end if

        ! WH2O ###PENDING###
        WH2O(icell) = (nAll(icell, patmo_idx_H2O) - nAll(icell + 1, patmo_idx_H2O)) / gridSpace(icell)

        ! Rainout Rates
        if (ABS(WH2O(icell)) >= 1d-15) then
            rkj(icell) = WH2O(icell) / 5.5d1 / (Av * 1d-9 + 1d0 / (Heff * R * TgasAll(icell)))
            Qj(icell)  = 1d0 - fz(icell) + fz(icell) / (gamma(icell) * rkj(icell)) * (1d0 - exp(-rkj(icell) * gamma(icell)))
            wetdep(icell, idx) = (1d0 - exp(-rkj(icell) * gamma(icell))) / (gamma(icell) * Qj(icell))
        else
            wetdep(icell, idx) = 0d0
        end if
    end do

    if (isWrite) then
        write(1000 + idx, *) "icell, WH2O(icell), fz(icell), gamma(icell), wetdep(icell,idx)"
        do icell=1, tropos
            write(1000 + idx, *) icell, WH2O(icell), fz(icell), gamma(icell), wetdep(icell, idx)
        end do
    end if

end subroutine computeWetDep

subroutine computeSulfuricAcidAerosol()
    use patmo_parameters
    implicit none
    integer :: icell
    integer, parameter :: initialStrat = 13, finalStrat = 34

    va(:) = 0d0
    pa(:) = 0d0
    gd(:) = 0d0

    open(60, file="vapor_H2SO4.txt", status="old")
    open(61, file="partial_H2SO4.txt", status="old")
    open(62, file="SO4_deposition_rate.txt", status="old")
    do icell = initialStrat, finalStrat
        read(60, *) va(icell)
        read(61, *) pa(icell)
    end do
    do icell = 1, cellsNumber
        read(62, *) gd(icell)
    end do
    close(60)
    close(61)
    close(62)

end subroutine computeSulfuricAcidAerosol

subroutine computeH2OCondensation(m, rho, r_a, isWrite)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    intrinsic sqrt
    integer :: icell
    real(8), intent(IN) :: m, rho, r_a
    logical, intent(IN) :: isWrite
    real(8), dimension(cellsNumber) :: n_v, p_v, t_c
    real(8), parameter :: A = 4.6543d0, B = 1.435264d3, C = -6.4848d1
            ! Stull (1947), taken from NIST chemistry webbook.

    p_v(:) = (1d1 ** (A - (B / (C + TgasAll(:))))) * 1d-6 ! dyne cm-2
    n_v(:) = p_v(:) / (kboltzmann * TgasAll(:))
    t_c(:) = m / (4d0 * rho) &
            * (8d0 * kboltzmann * TgasAll(:) / (pi * m)) ** (5d-1) &
            * (nAll(:, patmo_idx_H2O) - n_v(:) / r_a)

    do icell = 1, cellsNumber
        if (t_c(icell) > 0) then
            condenseH2O(icell) = t_c(icell)
        else
            condenseH2O(icell) = 0d0
        end if
    end do

    if (isWrite) then
        open(100, file="condense.out")
        write(100, *) "icell, T, n(H2O), p_v, n_v, t_c"
        do icell=1, cellsNumber
            write(100, *) icell, TgasAll(icell), nAll(icell, patmo_idx_H2O), p_v(icell) * 1d6, n_v(icell), t_c(icell)
        end do
        close(100)
    end if

end subroutine computeH2OCondensation
