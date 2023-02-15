module patmo_jacobian
contains

  !*********************
  !jacobian for chemistry + diffusion (DVODE_f90)
  subroutine jex(neq, t, nin, ian, jan, nz, pd)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    integer::neq,ian(size(iaSparsity(:)))
    integer::jan(size(jaSparsity(:))),nz,i,j,ii,k
    real*8::pd(nz), t, nin(neqAll)
    real*8::pd_vec(neqAll,neqAll)
    real*8::flux1(cellsNumber)
    real*8::flux2(cellsNumber)
    real*8::flux3(cellsNumber)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::mu(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)

    !when the function is called with nz==0 returns the non-zero elements
    ! known from sparsity calculation in patmo.f90
    if(nz==0) then
       nz = nonZeroElements
       return
    end if

    !get the IAN/JAN sparsity structure arrays from the bundle
    ian(:) = iaSparsity(:)
    jan(:) = jaSparsity(:)

    !init jacobian matric
    pd_vec(:,:) = 0d0

    !roll chemistry
    do i=1,speciesNumber
       n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do
    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))

    !chemistry jacobian
    do i=1,reactionsNumber
       do j=1,cellsNumber
#PATMO_jacobian
       end do
    end do

    m(:) = getSpeciesMass()

    !prepare diffusion elements
    do j=1,cellsNumber
       ngas(j) = sum(n(j,:))
       mu(j) = sum(m(:)*n(j,:))/ngas(j)
    end do

    !forward grid points
    do j=1,cellsNumber-1
       dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
       kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
       Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
       Tgas_p(j) = Tgas(j+1)
       ngas_p(j) = ngas(j+1)
       ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
       n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
       dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
       kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
       Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
       Tgas_m(j) = Tgas(j-1)
       ngas_m(j) = ngas(j-1)
       ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
       n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)


    !compute some prefactor to save time in the next loop
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    do i=1,speciesNumber
       prem(:) = (mu(:)-m(i))*gravity/kboltzmann*gridSpace(:)
       d_hp(:,i) =  dzzh_hp(:) &
            * (prem(:)*iTgas_hp(:) &
            - therm_hp(:))
       d_hm(:,i) = dzzh_hm(:) &
            * (prem(:)*iTgas_hm(:) &
            - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    !total density ratios
    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

    !diffusion jacobian
    do i=1,speciesNumber
       !dndot/dnij
       do j=1,cellsNumber
          pd_vec(cellsNumber*(i-1)+j,cellsNumber*(i-1)+j) = &
               pd_vec(cellsNumber*(i-1)+j,cellsNumber*(i-1)+j) &
               - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
               + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j))
       end do

       !dndot/dni(j+1)
       do j=1,cellsNumber-1
          pd_vec(cellsNumber*(i-1)+j,cellsNumber*(i-1)+j+1) = &
               pd_vec(cellsNumber*(i-1)+j,cellsNumber*(i-1)+j+1) &
               + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j)
       end do

       !dndot/dni(j-1)
       do j=2,cellsNumber
          pd_vec(cellsNumber*(i-1)+j,cellsNumber*(i-1)+j-1) = &
               pd_vec(cellsNumber*(i-1)+j,cellsNumber*(i-1)+j-1) &
               + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j)
       end do

    end do

    !unpack using the Yale sparsity format (IAN/JAN)
    ii = 0 !unrolled Jacobian index
    !loop on columns
    do j=1,neq
       !loop on number of non-zero elements of the jth column
       do k=iaSparsity(j),iaSparsity(j+1)-1
          !increase the index of the unrolled Jacobian
          ii = ii + 1
          !unroll the dense Jacobian (ja gives the row index)
          pd(ii) = pd_vec(jaSparsity(k),j)
       end do
    end do

  end subroutine jex


  !********************
  !Dummy jacobian for DLSODES
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    implicit none
    integer::neq, j, ian, jan
    real*8::tt, n(neq), pdj(neq)

  end subroutine jes
end module patmo_jacobian

