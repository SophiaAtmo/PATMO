module patmo_parameters
  use patmo_commons
  implicit none
  real*8::thermalDiffusionFactor
  real*8::diffusionDzz(cellsNumber)
  real*8::eddyKzz(cellsNumber)
  real*8::gravity
  real*8::krate(cellsNumber,reactionsNumber)
  real*8::nall(cellsNumber,speciesNumber)
  real*8::TgasAll(cellsNumber)
  real*8::ntotAll(cellsNumber)
  real*8::meanMolecularMass
  real*8::height(cellsNumber)
  real*8::idh2(cellsNumber)
  real*8::gridSpace(cellsNumber)
  real*8::xsecAll(photoBinsNumber,photoReactionsNumber)
  real*8::photoFlux(photoBinsNumber)
  real*8::energyMid(photoBinsNumber)
  real*8::energySpan(photoBinsNumber)
  real*8::energyLeft(photoBinsNumber)
  real*8::energyRight(photoBinsNumber)
  real*8::tauAll(photoBinsNumber,cellsNumber)
  real*8::cumulativeFlux(cellsNumber,reactionsNumber)
  real*8::wetdep(cellsNumber,chemSpeciesNumber)
  real*8::drydep(cellsNumber,chemSpeciesNumber)
  real*8::va(cellsNumber), pa(cellsNumber), gd(cellsNumber)
  real*8::condenseH2O(cellsNumber)
  real*8::Hesc, H2esc
  real*8::O2Flux, CH4Flux, S8SurFall, SO4SurFall
  integer::iaSparsity(neqAll+1)
  integer,allocatable::jaSparsity(:)
  integer::nonZeroElements
  logical::ratesUpdate
  character(len=maxNameLength)::reactionsVerbatim(reactionsNumber)
end module patmo_parameters
