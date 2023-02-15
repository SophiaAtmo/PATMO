module patmo_photo
contains

  !**************
  subroutine loadPhotoMetric(fname)
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    character::aout
    real*8::rout(4)
    integer::ios,i

    !open metric file (energy: mid/eV, span/eV, left/eV, right/eV)
    open(22,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: problem loading "//trim(fname)
      stop
    end if

    read(22,'(a)') aout

    !loop on bins to read (already interpolated by python)
    do i=1,photoBinsNumber
      read(22,*,iostat=ios) rout(:)
      if(ios/=0) then
        print *,"ERROR: problem while reading "//trim(fname)
        stop
      end if
      !load energy mid, span, left, right (eV)
      energyMid(i) = rout(1)
      energySpan(i) = rout(2)
      energyLeft(i) = rout(3)
      energyRight(i) = rout(4)
    end do

    close(22)

  end subroutine loadPhotoMetric

  !***************
  subroutine loadAllPhotoXsecs()
    use patmo_parameters
    implicit none

    !COS -> CO + S
    call loadPhotoXsec("xsecs/COS__CO_S.dat",1)
    !O3 -> O2 + O
    call loadPhotoXsec("xsecs/O3__O2_O.dat",2)
    !O2 -> O + O
    call loadPhotoXsec("xsecs/O2__O_O.dat",3)
    !CS2 -> CS + S
    call loadPhotoXsec("xsecs/CS2__CS_S.dat",4)
    !CS2 -> CS2E
    call loadPhotoXsec("xsecs/CS2__CS2E.dat",5)
    !SO2 -> SO + O
    call loadPhotoXsec("xsecs/SO2__SO_O.dat",6)
    !SO3 -> SO2 + O
    call loadPhotoXsec("xsecs/SO3__SO2_O.dat",7)
    !H2S -> SH + H
    call loadPhotoXsec("xsecs/H2S__SH_H.dat",8)
    !SO -> S + O
    call loadPhotoXsec("xsecs/SO__S_O.dat",9)

  end subroutine loadAllPhotoXsecs

  !***************
  subroutine loadPhotoXsec(fname,index)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    character(len=100)::aout
    integer,intent(in)::index
    integer::ios,i
    real*8::rout(2)

    !open xsec file
    open(22,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: problem loading "//trim(fname)
      stop
    end if

    !skip header (one line)
    read(22,'(a)') aout

    !loop on bins to read (already interpolated by python)
    do i=1,photoBinsNumber
      read(22,*,iostat=ios) rout(:)
      if(ios/=0) then
        print *,"ERROR: problem while reading "//trim(fname)
        stop
      end if
      !load xsecs for all cells (cm2)
      xsecAll(i,index) = rout(2)
    end do

    close(22)
  end subroutine loadPhotoXsec

  !**************************
  function fluxBB(x,Tbb)
    use patmo_constants
    implicit none
    real*8,intent(in)::x,Tbb
    real*8::fluxBB,xexp

    !exponent
    xexp = x/kboltzmann_eV/Tbb

    !default value
    fluxBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      fluxBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function fluxBB

end module patmo_photo
