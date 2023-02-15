module patmo_commons
  implicit none

#PATMO_commons

  integer,parameter::chemReactionsOffset = 0
  integer,parameter::photoReactionsOffset = chemReactionsNumber
  integer,parameter::reverseReactionsOffset = &
       photoReactionsOffset + photoReactionsNumber

  integer,parameter::neqAll = speciesNumber*cellsNumber
  integer,parameter::maxNameLength = 50

#PATMO_photochemPartners

#PATMO_reaction_arrays

end module patmo_commons
