module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
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
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
       n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

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

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
       prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
       d_hp(:,i) =  dzzh_hp(:) &
            * (prem(:)*iTgas_hp(:) &
            - therm_hp(:))
       d_hm(:,i) = dzzh_hm(:) &
            * (prem(:)*iTgas_hm(:) &
            - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0
#PATMO_ODE

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

    do i=1,chemSpeciesNumber
       dn(:,i) = dn(:,i) &
            + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
            - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
            + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
            + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
    end do
    
    !Chemical Species with constant concentration
#PATMO_constantspecies
    
    ! Gravity Settling
#PATMO_gravitysettling
   
    ! Dry Deposition: assumed a deposition rate of 0.1 cm/s 
    !dn(1,patmo_idx_A)=dn(1,patmo_idx_A) - 0.1/(layer_thickness(in cm))*n(1,patmo_idx_A)
    ! Fix the mixing ratio of CH4 and O2 at the bottom layer as a constant (Claire et al., 2014; Zahnle et al., 2006)
#PATMO_drydeppecies
      
    !Volcanic emission
    !The release of 1 Tmol/year from Claire et al., 2014, with an H2S:SO2 ratio of 1:10
    !The release of molecular hydrogen 3 Tmol/year from Claire et al., 2014
#PATMO_emissionspecies

#IFPATMO_useWaterRemoval    
    ! Water Removal
    dn(:,patmo_idx_H2O) = dn(:,patmo_idx_H2O) - n(:,patmo_idx_H2O) * condenseH2O(:)
#ENDIFPATMO
    ! Wet Deposition
    do j=12, 2, -1
        do i = 1, chemSpeciesNumber
            dn(j,     i) = dn(j,     i) - wetdep(j, i) * n(j, i)
            dn(j - 1, i) = dn(j - 1, i) + wetdep(j, i) * n(j, i)
        end do
    end do
    do i = 1, chemSpeciesNumber
        dn(1, i) = dn(1, i) - wetdep(1, i) * n(1, i)
    end do
#IFPATMO_useAerosolformation   
    !aerosol formation
    do i=13,34
         if (va(i) <= n(i, patmo_idx_H2SO4) .AND. pa(i) >= n(i, patmo_idx_H2SO4)) then
            dn(i, patmo_idx_H2SO4) = dn(i, patmo_idx_H2SO4) - (n(i, patmo_idx_H2SO4) - va(i))
            dn(i, patmo_idx_SO4)   = dn(i, patmo_idx_SO4)   + (n(i, patmo_idx_H2SO4) - va(i))
         end if
    end do	
#ENDIFPATMO

    
#IFPATMO_useHescape
    ! Hydrogen Escape
    if (n(cellsNumber, patmo_idx_H) > Hesc) then
        dn(cellsNumber, patmo_idx_H) = dn(cellsNumber, patmo_idx_H) - Hesc
        !print *, "triggered H escape"
    else
        n(cellsNumber, patmo_idx_H) = 0d0
        dn(cellsNumber, patmo_idx_H) = 0d0
    end if
    
    if (n(cellsNumber, patmo_idx_H2) > H2esc) then
        dn(cellsNumber, patmo_idx_H2) = dn(cellsNumber, patmo_idx_H2) - H2esc
        !print *, "triggered H2 escape"
    else
        n(cellsNumber, patmo_idx_H2) = 0d0
        dn(cellsNumber, patmo_idx_H2) = 0d0
    end if
#ENDIFPATMO
    
    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
       dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
