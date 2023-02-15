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
      Tgas(i) = max(inTgas(i),2d2)
      Tgas(i) = min(Tgas(i),5d3)
    end do

    lnTgas(:) = log(Tgas(:))
    Tgas2(:) = Tgas(:)**2
    Tgas3(:) = Tgas(:)**3
    Tgas4(:) = Tgas(:)**4
    invTgas(:) = 1d0/Tgas(:)

    !CO2 + SH -> COS + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,57) = krate(i,1)*exp(-2.775076d-1*(lnTgas(i)-1d0) &
            + 1.258706d-3*Tgas(i) &
            - 4.510004d-7*Tgas2(i) &
            - 6.102868d-11*Tgas3(i) &
            + 6.191498d-14*Tgas4(i) &
            - 1.770436d4*invTgas(i) &
            + 1.658291d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,57) = krate(i,1)*exp(5.450499d-1*(lnTgas(i)-1d0) &
            - 3.939921d-4*Tgas(i) &
            + 5.516733d-8*Tgas2(i) &
            - 4.211364d-12*Tgas3(i) &
            + 1.161815d-16*Tgas4(i) &
            - 1.759891d4*invTgas(i) &
            - 2.155117d0)
      else
        krate(i,57) = 0d0
      end if
    end do

    !CO + SO -> COS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,58) = krate(i,2)*exp(-2.257872d0*(lnTgas(i)-1d0) &
            + 8.400735d-3*Tgas(i) &
            - 5.554705d-6*Tgas2(i) &
            + 2.47746d-9*Tgas3(i) &
            - 4.967152d-13*Tgas4(i) &
            - 2.581411d4*invTgas(i) &
            + 5.859493d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,58) = krate(i,2)*exp(9.007697d-1*(lnTgas(i)-1d0) &
            + 1.738856d-4*Tgas(i) &
            - 5.041413d-8*Tgas2(i) &
            + 5.80007d-12*Tgas3(i) &
            - 2.538475d-16*Tgas4(i) &
            - 2.530287d4*invTgas(i) &
            - 8.614472d0)
      else
        krate(i,58) = 0d0
      end if
    end do

    !SH + COS -> CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,59) = krate(i,3)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,59) = krate(i,3)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,59) = 0d0
      end if
    end do

    !CS + SO -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,60) = krate(i,4)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,60) = krate(i,4)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,60) = 0d0
      end if
    end do

    !COS + S -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,61) = krate(i,5)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,61) = krate(i,5)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,61) = 0d0
      end if
    end do

    !CO + S2 -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,62) = krate(i,6)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,62) = krate(i,6)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,62) = 0d0
      end if
    end do

    !SCSOH -> CS2 + OH
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,63) = krate(i,7)*exp(8.837417d-1*(lnTgas(i)-1d0) &
!            + 6.639121d-3*Tgas(i) &
!            - 1.182399d-5*Tgas2(i) &
!            + 6.691894d-9*Tgas3(i) &
!            - 1.506943d-12*Tgas4(i) &
!            - 2.240521d4*invTgas(i) &
!            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,63) = krate(i,7)*exp(2.321241d0*(lnTgas(i)-1d0) &
!            - 6.394723d-3*Tgas(i) &
!            + 7.538362d-7*Tgas2(i) &
!            - 5.96604d-11*Tgas3(i) &
!            + 2.152074d-15*Tgas4(i) &
!            - 2.306419d4*invTgas(i) &
!            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
!      else
        krate(i,63) = 0d0
!      end if
    end do

    !COS + HSO2 -> SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,64) = krate(i,8)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,64) = krate(i,8)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,64) = 0d0
      end if
    end do

    !COS + O -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,65) = krate(i,9)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,65) = krate(i,9)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,65) = 0d0
      end if
    end do

    !COS + O2 -> CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,66) = krate(i,10)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,66) = krate(i,10)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,66) = 0d0
      end if
    end do

    !CO + S -> CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,67) = krate(i,11)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,67) = krate(i,11)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,67) = 0d0
      end if
    end do

    !H2O + SH -> H2S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,68) = krate(i,12)*exp(2.289248d-1*(lnTgas(i)-1d0) &
            - 2.744914d-3*Tgas(i) &
            + 3.195679d-6*Tgas2(i) &
            - 1.867226d-9*Tgas3(i) &
            + 4.404765d-13*Tgas4(i) &
            - 1.407682d4*invTgas(i) &
            + 2.589349d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,68) = krate(i,12)*exp(1.087555d-1*(lnTgas(i)-1d0) &
            + 2.368912d-4*Tgas(i) &
            - 5.712312d-8*Tgas2(i) &
            + 6.875396d-12*Tgas3(i) &
            - 3.158447d-16*Tgas4(i) &
            - 1.386166d4*invTgas(i) &
            - 4.086126d-1)
      else
        krate(i,68) = 0d0
      end if
    end do

    !OH + SH -> H2S + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,69) = krate(i,13)*exp(-3.881414d-1*(lnTgas(i)-1d0) &
            - 3.001707d-3*Tgas(i) &
            + 3.8507d-6*Tgas2(i) &
            - 2.188698d-9*Tgas3(i) &
            + 4.983884d-13*Tgas4(i) &
            - 6.167561d3*invTgas(i) &
            + 1.669856d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,69) = krate(i,13)*exp(-3.476293d-1*(lnTgas(i)-1d0) &
            + 6.02411d-4*Tgas(i) &
            - 8.878291d-8*Tgas2(i) &
            + 8.146112d-12*Tgas3(i) &
            - 3.109825d-16*Tgas4(i) &
            - 5.80616d3*invTgas(i) &
            - 2.936611d-1)
      else
        krate(i,69) = 0d0
      end if
    end do

    !H2 + SH -> H2S + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,70) = krate(i,14)*exp(5.912447d-1*(lnTgas(i)-1d0) &
            - 6.552842d-3*Tgas(i) &
            + 6.759321d-6*Tgas2(i) &
            - 3.681057d-9*Tgas3(i) &
            + 8.29721d-13*Tgas4(i) &
            - 6.805796d3*invTgas(i) &
            - 1.615769d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,70) = krate(i,14)*exp(-4.856017d-1*(lnTgas(i)-1d0) &
            + 7.564716d-4*Tgas(i) &
            - 1.126842d-7*Tgas2(i) &
            + 9.954863d-12*Tgas3(i) &
            - 3.737095d-16*Tgas4(i) &
            - 6.564682d3*invTgas(i) &
            + 1.206637d0)
      else
        krate(i,70) = 0d0
      end if
    end do

    !H2O + HSO -> H2S + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,71) = krate(i,15)*exp(8.775649d-2*(lnTgas(i)-1d0) &
            - 4.496811d-4*Tgas(i) &
            + 3.892062d-7*Tgas2(i) &
            - 1.499713d-10*Tgas3(i) &
            + 2.458951d-14*Tgas4(i) &
            - 3.069931d4*invTgas(i) &
            + 2.097107d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,71) = krate(i,15)*exp(1.268016d-1*(lnTgas(i)-1d0) &
            - 1.406198d-5*Tgas(i) &
            + 2.514851d-8*Tgas2(i) &
            - 3.019166d-12*Tgas3(i) &
            + 1.213293d-16*Tgas4(i) &
            - 3.060978d4*invTgas(i) &
            - 3.04536d-1)
      else
        krate(i,71) = 0d0
      end if
    end do

    !H + SO -> SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,72) = krate(i,16)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,72) = krate(i,16)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,72) = 0d0
      end if
    end do

    !OH + SO -> SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,73) = krate(i,17)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,73) = krate(i,17)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,73) = 0d0
      end if
    end do

    !HSO + O2 -> SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,74) = krate(i,18)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,74) = krate(i,18)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,74) = 0d0
      end if
    end do

    !SO2 + O2 -> SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,75) = krate(i,19)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,75) = krate(i,19)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,75) = 0d0
      end if
    end do

    !SO2 + O -> SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,76) = krate(i,20)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,76) = krate(i,20)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,76) = 0d0
      end if
    end do

    !SO2 + H -> SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,77) = krate(i,21)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,77) = krate(i,21)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,77) = 0d0
      end if
    end do

    !SO2 + NO -> SO + NO2
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,78) = krate(i,22)*exp(-3.307801d-1*(lnTgas(i)-1d0) &
!            - 7.751533d-4*Tgas(i) &
!            + 1.465141d-6*Tgas2(i) &
!            - 9.739932d-10*Tgas3(i) &
!            + 2.437221d-13*Tgas4(i) &
!            - 2.95164d4*invTgas(i) &
!            + 2.427759d0)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,78) = krate(i,22)*exp(2.087491d-1*(lnTgas(i)-1d0) &
!            - 1.603123d-4*Tgas(i) &
!            + 4.013125d-8*Tgas2(i) &
!            - 2.849403d-12*Tgas3(i) &
!            + 7.811797d-17*Tgas4(i) &
!            - 2.92732d4*invTgas(i) &
!            - 9.201861d-1)
!      else
        krate(i,78) = 0d0
!      end if
    end do

    !SO + O -> S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,79) = krate(i,23)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,79) = krate(i,23)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,79) = 0d0
      end if
    end do

    !O2 + SO -> S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,80) = krate(i,24)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,80) = krate(i,24)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,80) = 0d0
      end if
    end do

    !H + SO -> S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,81) = krate(i,25)*exp(1.906453d-1*(lnTgas(i)-1d0) &
            + 2.350427d-3*Tgas(i) &
            - 3.539649d-6*Tgas2(i) &
            + 2.165587d-9*Tgas3(i) &
            - 5.105171d-13*Tgas4(i) &
            - 1.088276d4*invTgas(i) &
            + 4.006756d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,81) = krate(i,25)*exp(-7.510469d-1*(lnTgas(i)-1d0) &
            + 1.095328d-4*Tgas(i) &
            - 7.977419d-9*Tgas2(i) &
            + 9.398975d-13*Tgas3(i) &
            - 5.635597d-17*Tgas4(i) &
            - 1.14541d4*invTgas(i) &
            + 6.538105d0)
      else
        krate(i,81) = 0d0
      end if
    end do

    !OH + SO3 -> SO2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,82) = krate(i,26)*exp(1.61001d0*(lnTgas(i)-1d0) &
            - 8.009682d-3*Tgas(i) &
            + 6.273806d-6*Tgas2(i) &
            - 3.072642d-9*Tgas3(i) &
            + 6.534202d-13*Tgas4(i) &
            - 8.876536d3*invTgas(i) &
            - 1.31498d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,82) = krate(i,26)*exp(-5.787838d-1*(lnTgas(i)-1d0) &
            - 1.413477d-4*Tgas(i) &
            + 6.323959d-8*Tgas2(i) &
            - 7.909567d-12*Tgas3(i) &
            + 3.569003d-16*Tgas4(i) &
            - 9.036184d3*invTgas(i) &
            + 7.706091d0)
      else
        krate(i,82) = 0d0
      end if
    end do

    !SO3 + O2 -> SO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,83) = krate(i,27)*exp(9.251222d-1*(lnTgas(i)-1d0) &
            - 4.310392d-3*Tgas(i) &
            + 4.183748d-6*Tgas2(i) &
            - 2.427028d-9*Tgas3(i) &
            + 5.828176d-13*Tgas4(i) &
            - 2.890986d4*invTgas(i) &
            - 5.108402d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,83) = krate(i,27)*exp(6.75679d0*(lnTgas(i)-1d0) &
            - 6.822651d-3*Tgas(i) &
            + 1.426713d-6*Tgas2(i) &
            - 1.55403d-10*Tgas3(i) &
            + 6.59569d-15*Tgas4(i) &
            - 2.659453d4*invTgas(i) &
            - 3.370434d1)
      else
        krate(i,83) = 0d0
      end if
    end do

    !SO2 + OH -> HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,84) = krate(i,28)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,84) = krate(i,28)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,84) = 0d0
      end if
    end do

    !O2 + O2 + SH -> HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,85) = krate(i,29)*exp(-3.706548d0*(lnTgas(i)-1d0) &
            + 5.543697d-4*Tgas(i) &
            + 4.589081d-6*Tgas2(i) &
            - 3.664404d-9*Tgas3(i) &
            + 9.760036d-13*Tgas4(i) &
            + 1.734979d3*invTgas(i) &
            + 4.837009d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,85) = krate(i,29)*exp(6.324081d0*(lnTgas(i)-1d0) &
            - 5.98477d-3*Tgas(i) &
            + 1.287301d-6*Tgas2(i) &
            - 1.430812d-10*Tgas3(i) &
            + 6.126488d-15*Tgas4(i) &
            + 5.307368d3*invTgas(i) &
            - 5.070441d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,85) = 0d0
      end if
    end do

    !HO2 + SO2 -> HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,86) = krate(i,30)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,86) = krate(i,30)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,86) = 0d0
      end if
    end do

    !HO2 + SO3 -> HSO3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,87) = krate(i,31)*exp(2.907618d-1*(lnTgas(i)-1d0) &
            + 2.512911d-3*Tgas(i) &
            - 2.879387d-6*Tgas2(i) &
            + 1.644321d-9*Tgas3(i) &
            - 3.807075d-13*Tgas4(i) &
            - 3.383301d3*invTgas(i) &
            - 1.13614d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,87) = krate(i,31)*exp(3.861018d-1*(lnTgas(i)-1d0) &
            - 9.114733d-5*Tgas(i) &
            - 1.33545d-8*Tgas2(i) &
            + 3.127919d-12*Tgas3(i) &
            - 1.627562d-16*Tgas4(i) &
            - 3.561354d3*invTgas(i) &
            + 4.998669d-1)
      else
        krate(i,87) = 0d0
      end if
    end do

    !SO3 -> SO2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,88) = krate(i,32)*exp(4.468463d0*(lnTgas(i)-1d0) &
            - 8.475314d-3*Tgas(i) &
            + 4.624041d-6*Tgas2(i) &
            - 1.883545d-9*Tgas3(i) &
            + 3.626003d-13*Tgas4(i) &
            - 4.110368d4*invTgas(i) &
            - 3.083707d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,88) = krate(i,32)*exp(6.310961d-1*(lnTgas(i)-1d0) &
            - 5.418875d-4*Tgas(i) &
            + 7.125405d-8*Tgas2(i) &
            - 5.612989d-12*Tgas3(i) &
            + 2.029648d-16*Tgas4(i) &
            - 4.192898d4*invTgas(i) &
            + 1.551566d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,88) = 0d0
      end if
    end do

    !HSO3 -> SO2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,89) = krate(i,33)*exp(4.482077d0*(lnTgas(i)-1d0) &
            - 9.672905d-3*Tgas(i) &
            + 5.280526d-6*Tgas2(i) &
            - 2.124199d-9*Tgas3(i) &
            + 4.034084d-13*Tgas4(i) &
            - 1.063905d4*invTgas(i) &
            - 5.185011d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,89) = krate(i,33)*exp(2.856084d-2*(lnTgas(i)-1d0) &
            - 4.957808d-4*Tgas(i) &
            + 7.049484d-8*Tgas2(i) &
            - 5.55515d-12*Tgas3(i) &
            + 1.947841d-16*Tgas4(i) &
            - 1.159243d4*invTgas(i) &
            + 1.639613d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,89) = 0d0
      end if
    end do

    !H2SO4 -> SO3 + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,90) = krate(i,34)*exp(2.039365d0*(lnTgas(i)-1d0) &
            - 8.55842d-3*Tgas(i) &
            + 5.821711d-6*Tgas2(i) &
            - 2.687074d-9*Tgas3(i) &
            + 5.604965d-13*Tgas4(i) &
            - 1.132526d4*invTgas(i) &
            + 8.315681d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,90) = krate(i,34)*exp(-1.361725d0*(lnTgas(i)-1d0) &
            + 5.032744d-5*Tgas(i) &
            + 1.803258d-8*Tgas2(i) &
            - 2.660542d-12*Tgas3(i) &
            + 1.317101d-16*Tgas4(i) &
            - 1.191318d4*invTgas(i) &
            + 2.406728d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,90) = 0d0
      end if
    end do

    !SO2 + OH + OH -> H2SO4
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,91) = krate(i,35)*exp(-7.124894d0*(lnTgas(i)-1d0) &
    !        + 1.677694d-2*Tgas(i) &
    !        - 9.790732d-6*Tgas2(i) &
    !        + 4.249147d-9*Tgas3(i) &
    !        - 8.65185d-13*Tgas4(i) &
    !        + 6.03382d4*invTgas(i) &
    !        - 3.821053d0)*(1.3806488d-22*Tgas(i))**(2)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,91) = krate(i,35)*exp(2.742437d-1*(lnTgas(i)-1d0) &
    !        + 8.570799d-4*Tgas(i) &
    !        - 1.209464d-7*Tgas2(i) &
    !        + 9.544247d-12*Tgas3(i) &
    !        - 3.298127d-16*Tgas4(i) &
    !        + 6.189766d4*invTgas(i) &
    !        - 3.946799d1)*(1.3806488d-22*Tgas(i))**(2)
    !  else
        krate(i,91) = 0d0
    !  end if
    end do

    !O3 -> O + O2
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,92) = krate(i,36)*exp(3.543341d0*(lnTgas(i)-1d0) &
    !        - 4.164922d-3*Tgas(i) &
    !        + 4.402935d-7*Tgas2(i) &
    !        + 5.434827d-10*Tgas3(i) &
    !        - 2.202172d-13*Tgas4(i) &
    !        - 1.219382d4*invTgas(i) &
    !        - 2.572867d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,92) = krate(i,36)*exp(-6.125694d0*(lnTgas(i)-1d0) &
    !        + 6.280764d-3*Tgas(i) &
    !        - 1.355459d-6*Tgas2(i) &
    !        + 1.4979d-10*Tgas3(i) &
    !        - 6.392726d-15*Tgas4(i) &
    !        - 1.533445d4*invTgas(i) &
    !        + 4.921999d1)*(1.3806488d-22*Tgas(i))**(-1)
    !  else
        krate(i,92) = 0d0
   !   end if
    end do

    !N + N -> N2
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,93) = krate(i,37)*exp(-1.468995d0*(lnTgas(i)-1d0) &
    !        - 6.183049d-5*Tgas(i) &
    !        - 8.383324d-8*Tgas2(i) &
    !        + 2.029422d-10*Tgas3(i) &
    !        - 7.044062d-14*Tgas4(i) &
    !        + 1.132563d5*invTgas(i) &
    !        - 5.420347d0)*(1.3806488d-22*Tgas(i))**(1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,93) = krate(i,37)*exp(-1.879309d0*(lnTgas(i)-1d0) &
    !        + 5.235596d-4*Tgas(i) &
    !        - 4.24307d-8*Tgas2(i) &
    !        + 1.512378d-12*Tgas3(i) &
    !        - 2.676777d-17*Tgas4(i) &
    !        + 1.131915d5*invTgas(i) &
    !        - 3.427331d0)*(1.3806488d-22*Tgas(i))**(1)
    !  else
        krate(i,93) = 0d0
    !  end if
    end do

    !SO4 -> SO2
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,94) = krate(i,38)*exp(2.67174d0*(lnTgas(i)-1d0) &
    !        - 1.667105d-2*Tgas(i) &
    !        + 9.453086d-6*Tgas2(i) &
    !        - 3.840731d-9*Tgas3(i) &
    !        + 7.244981d-13*Tgas4(i) &
    !        + 6.428951d3*invTgas(i) &
    !        - 1.111545d1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,94) = krate(i,38)*exp(-4.700143d0*(lnTgas(i)-1d0) &
    !        - 6.490479d-4*Tgas(i) &
    !        + 8.865994d-8*Tgas2(i) &
    !        - 7.503263d-12*Tgas3(i) &
    !        + 2.833749d-16*Tgas4(i) &
    !        + 4.953039d3*invTgas(i) &
    !        + 2.414516d1)
    !  else
        krate(i,94) = 0d0
    !  end if
    end do

    !SO2 -> CH3SCH3 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,95) = krate(i,39)*exp(4.774011d0*(lnTgas(i)-1d0) &
            - 1.557652d-3*Tgas(i) &
            + 7.154449d-6*Tgas2(i) &
            - 4.178225d-9*Tgas3(i) &
            + 9.461989d-13*Tgas4(i) &
            - 5.983783d4*invTgas(i) &
            - 3.866955d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,95) = krate(i,39)*exp(3.625742d0*(lnTgas(i)-1d0) &
            + 6.941559d-3*Tgas(i) &
            - 8.105797d-7*Tgas2(i) &
            + 6.424539d-11*Tgas3(i) &
            - 2.282433d-15*Tgas4(i) &
            - 5.948346d4*invTgas(i) &
            - 1.265807d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,95) = 0d0
      end if
    end do

    !SO2 -> CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,96) = krate(i,40)*exp(5.597728d0*(lnTgas(i)-1d0) &
            - 1.118526d-3*Tgas(i) &
            + 6.816712d-6*Tgas2(i) &
            - 3.990816d-9*Tgas3(i) &
            + 9.087257d-13*Tgas4(i) &
            - 3.408447d4*invTgas(i) &
            - 6.022887d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,96) = krate(i,40)*exp(3.920635d0*(lnTgas(i)-1d0) &
            + 7.508924d-3*Tgas(i) &
            - 8.588814d-7*Tgas2(i) &
            + 6.733831d-11*Tgas3(i) &
            - 2.3796d-15*Tgas4(i) &
            - 3.395525d4*invTgas(i) &
            - 3.431551d-1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,96) = 0d0
      end if
    end do

    !SO2 + CH4O3S -> CH3SCH3 + OH
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,97) = krate(i,41)*exp(3.171767d-1*(lnTgas(i)-1d0) &
    !        - 2.342044d-3*Tgas(i) &
    !        - 6.420484d-7*Tgas2(i) &
    !        + 8.147534d-10*Tgas3(i) &
    !        - 2.199762d-13*Tgas4(i) &
    !        - 4.031441d4*invTgas(i) &
    !        - 8.072663d0)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,97) = krate(i,41)*exp(-2.545704d0*(lnTgas(i)-1d0) &
    !        - 2.859464d-4*Tgas(i) &
    !        + 5.634379d-8*Tgas2(i) &
    !        - 5.532956d-12*Tgas3(i) &
    !        + 2.123002d-16*Tgas4(i) &
    !        - 4.130451d4*invTgas(i) &
    !        + 7.676252d0)
    !  else
        krate(i,97) = 0d0
   !   end if
    end do

    !CS2 -> CS2E + O2
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,98) = krate(i,42)*exp(3.782456d0*(lnTgas(i)-1d0) &
    !        - 1.498367d-3*Tgas(i) &
    !       + 1.641217d-6*Tgas2(i) &
    !        - 8.067746d-10*Tgas3(i) &
    !        + 1.621864d-13*Tgas4(i) &
    !        + 1.063944d3*invTgas(i) &
    !        + 3.657676d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,98) = krate(i,42)*exp(3.660961d0*(lnTgas(i)-1d0) &
    !        + 3.281829d-4*Tgas(i) &
    !        - 2.352494d-8*Tgas2(i) &
    !        + 1.714983d-12*Tgas3(i) &
    !        - 6.495672d-17*Tgas4(i) &
    !        + 1.215977d3*invTgas(i) &
    !        + 3.415363d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  else
        krate(i,98) = 0d0
    !  end if
    end do

    !CS2 -> CS2E + N2
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,99) = krate(i,43)*exp(3.531005d0*(lnTgas(i)-1d0) &
    !        - 6.183049d-5*Tgas(i) &
    !        - 8.383324d-8*Tgas2(i) &
    !        + 2.029422d-10*Tgas3(i) &
    !        - 7.044062d-14*Tgas4(i) &
    !        + 1.046976d3*invTgas(i) &
    !        + 2.96747d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,99) = krate(i,43)*exp(2.952576d0*(lnTgas(i)-1d0) &
    !        + 6.984502d-4*Tgas(i) &
    !        - 8.210527d-8*Tgas2(i) &
    !        + 6.550085d-12*Tgas3(i) &
    !        - 2.303776d-16*Tgas4(i) &
    !        + 9.239487d2*invTgas(i) &
    !        + 5.871888d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  else
        krate(i,99) = 0d0
    !  end if
    end do

    !CS + SO2 -> CS2E + O2
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,100) = krate(i,44)*exp(-1.451291d0*(lnTgas(i)-1d0) &
    !        + 7.972313d-3*Tgas(i) &
    !        - 6.985411d-6*Tgas2(i) &
    !        + 3.721109d-9*Tgas3(i) &
    !        - 8.390358d-13*Tgas4(i) &
    !        - 1.624584d4*invTgas(i) &
    !        + 3.123155d0)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,100) = krate(i,44)*exp(4.561796d-1*(lnTgas(i)-1d0) &
    !        - 3.051946d-5*Tgas(i) &
    !        + 9.916702d-9*Tgas2(i) &
    !        - 1.934438d-14*Tgas3(i) &
    !        - 4.758353d-17*Tgas4(i) &
    !        - 1.615798d4*invTgas(i) &
    !        - 4.34393d0)
    !  else
        krate(i,100) = 0d0
    !  end if
    end do

    !S + SO -> S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,101) = krate(i,45)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,101) = krate(i,45)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,101) = 0d0
      end if
    end do

    !CO + SO -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,102) = krate(i,46)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,102) = krate(i,46)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,102) = 0d0
      end if
    end do

    !CS2 + OH -> SCSOH
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,103) = krate(i,47)*exp(-8.837417d-1*(lnTgas(i)-1d0) &
!            - 6.639121d-3*Tgas(i) &
!            + 1.182399d-5*Tgas2(i) &
!            - 6.691894d-9*Tgas3(i) &
!            + 1.506943d-12*Tgas4(i) &
!            + 2.240521d4*invTgas(i) &
!            - 9.828921d0)*(1.3806488d-22*Tgas(i))**(1)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,103) = krate(i,47)*exp(-2.321241d0*(lnTgas(i)-1d0) &
!            + 6.394723d-3*Tgas(i) &
!            - 7.538362d-7*Tgas2(i) &
!            + 5.96604d-11*Tgas3(i) &
!            - 2.152074d-15*Tgas4(i) &
!            + 2.306419d4*invTgas(i) &
!            - 7.693985d0)*(1.3806488d-22*Tgas(i))**(1)
!      else
        krate(i,103) = 0d0
!      end if
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates
