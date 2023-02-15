module patmo_commons
  implicit none

  integer,parameter::reactionsNumber = 103
  integer,parameter::chemReactionsNumber = 47
  integer,parameter::photoReactionsNumber = 9
  integer,parameter::reverseReactionsNumber = 47
  integer,parameter::chemSpeciesNumber = 33
  integer,parameter::speciesNumber = 35
  integer,parameter::positionTgas = 34
  integer,parameter::positionDummy = 35
  integer,parameter::cellsNumber = 60
  integer,parameter::photoBinsNumber = 4440
  integer,parameter::patmo_idx_COS = 1
  integer,parameter::patmo_idx_S2 = 2
  integer,parameter::patmo_idx_O3 = 3
  integer,parameter::patmo_idx_O2 = 4
  integer,parameter::patmo_idx_NO2 = 5
  integer,parameter::patmo_idx_NO = 6
  integer,parameter::patmo_idx_HSO = 7
  integer,parameter::patmo_idx_HO2 = 8
  integer,parameter::patmo_idx_HSO2 = 9
  integer,parameter::patmo_idx_HSO3 = 10
  integer,parameter::patmo_idx_CS2 = 11
  integer,parameter::patmo_idx_CH4O3S = 12
  integer,parameter::patmo_idx_CO = 13
  integer,parameter::patmo_idx_H = 14
  integer,parameter::patmo_idx_O = 15
  integer,parameter::patmo_idx_N = 16
  integer,parameter::patmo_idx_S = 17
  integer,parameter::patmo_idx_SO2 = 18
  integer,parameter::patmo_idx_SO4 = 19
  integer,parameter::patmo_idx_CS = 20
  integer,parameter::patmo_idx_N2 = 21
  integer,parameter::patmo_idx_CO2 = 22
  integer,parameter::patmo_idx_SCSOH = 23
  integer,parameter::patmo_idx_OH = 24
  integer,parameter::patmo_idx_H2 = 25
  integer,parameter::patmo_idx_H2SO4 = 26
  integer,parameter::patmo_idx_SO3 = 27
  integer,parameter::patmo_idx_H2S = 28
  integer,parameter::patmo_idx_H2O = 29
  integer,parameter::patmo_idx_SH = 30
  integer,parameter::patmo_idx_SO = 31
  integer,parameter::patmo_idx_CS2E = 32
  integer,parameter::patmo_idx_CH3SCH3 = 33

  integer,parameter::chemReactionsOffset = 0
  integer,parameter::photoReactionsOffset = chemReactionsNumber
  integer,parameter::reverseReactionsOffset = &
      photoReactionsOffset + photoReactionsNumber

  integer,parameter::neqAll = speciesNumber*cellsNumber
  integer,parameter::maxNameLength = 50

  integer,dimension(photoReactionsNumber)::photoPartnerIndex = (/patmo_idx_COS,patmo_idx_O3,patmo_idx_O2,patmo_idx_CS2,patmo_idx_CS2,patmo_idx_SO2,patmo_idx_SO3,patmo_idx_H2S,patmo_idx_SO/)

  integer,parameter,dimension(reactionsNumber)::indexReactants2 = (/patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_HO2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_NO2,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_OH,&
      patmo_idx_HO2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_H2O,&
      positionDummy,&
      patmo_idx_O2,&
      positionDummy,&
      positionDummy,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_N2,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_SH,&
      patmo_idx_SO,&
      patmo_idx_COS,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_S2,&
      positionDummy,&
      patmo_idx_HSO2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_S,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_HSO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_NO,&
      patmo_idx_O,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO3,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_SO2,&
      patmo_idx_SO3,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_OH,&
      positionDummy,&
      patmo_idx_N,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_CH4O3S,&
      positionDummy,&
      positionDummy,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_OH/)
  integer,parameter,dimension(reactionsNumber)::indexReactants1 = (/patmo_idx_COS,&
      patmo_idx_COS,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_SCSOH,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_H2S,&
      patmo_idx_H2S,&
      patmo_idx_H2S,&
      patmo_idx_H2S,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_HSO,&
      patmo_idx_HSO,&
      patmo_idx_HSO2,&
      patmo_idx_HSO3,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO3,&
      patmo_idx_H2SO4,&
      patmo_idx_O,&
      patmo_idx_N2,&
      patmo_idx_SO2,&
      patmo_idx_CH3SCH3,&
      patmo_idx_CH3SCH3,&
      patmo_idx_CH3SCH3,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_S2,&
      patmo_idx_CS,&
      patmo_idx_SCSOH,&
      patmo_idx_COS,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_SO2,&
      patmo_idx_SO3,&
      patmo_idx_H2S,&
      patmo_idx_SO,&
      patmo_idx_CO2,&
      patmo_idx_CO,&
      patmo_idx_SH,&
      patmo_idx_CS,&
      patmo_idx_COS,&
      patmo_idx_CO,&
      patmo_idx_SCSOH,&
      patmo_idx_COS,&
      patmo_idx_COS,&
      patmo_idx_COS,&
      patmo_idx_CO,&
      patmo_idx_H2O,&
      patmo_idx_OH,&
      patmo_idx_H2,&
      patmo_idx_H2O,&
      patmo_idx_H,&
      patmo_idx_OH,&
      patmo_idx_HSO,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_O2,&
      patmo_idx_H,&
      patmo_idx_OH,&
      patmo_idx_SO3,&
      patmo_idx_SO2,&
      patmo_idx_O2,&
      patmo_idx_HO2,&
      patmo_idx_HO2,&
      patmo_idx_SO3,&
      patmo_idx_HSO3,&
      patmo_idx_H2SO4,&
      patmo_idx_SO2,&
      patmo_idx_O3,&
      patmo_idx_N,&
      patmo_idx_SO4,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS,&
      patmo_idx_S,&
      patmo_idx_CO,&
      patmo_idx_CS2/)
  integer,parameter,dimension(reactionsNumber)::indexProducts2 = (/patmo_idx_SH,&
      patmo_idx_SO,&
      patmo_idx_COS,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_S2,&
      positionDummy,&
      patmo_idx_HSO2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_S,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_HSO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_NO,&
      patmo_idx_O,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO3,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_SO2,&
      patmo_idx_SO3,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_OH,&
      positionDummy,&
      patmo_idx_N,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_CH4O3S,&
      positionDummy,&
      positionDummy,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_OH,&
      patmo_idx_S,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_S,&
      positionDummy,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_HO2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_OH,&
      patmo_idx_NO2,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_OH,&
      patmo_idx_HO2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_H2O,&
      positionDummy,&
      patmo_idx_O2,&
      positionDummy,&
      positionDummy,&
      patmo_idx_O,&
      patmo_idx_OH,&
      patmo_idx_OH,&
      patmo_idx_O2,&
      patmo_idx_N2,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      positionDummy/)
  integer,parameter,dimension(reactionsNumber)::indexProducts3 = (/positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_SH,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_OH,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy/)
  integer,parameter,dimension(reactionsNumber)::indexProducts1 = (/patmo_idx_CO2,&
      patmo_idx_CO,&
      patmo_idx_SH,&
      patmo_idx_CS,&
      patmo_idx_COS,&
      patmo_idx_CO,&
      patmo_idx_SCSOH,&
      patmo_idx_COS,&
      patmo_idx_COS,&
      patmo_idx_COS,&
      patmo_idx_CO,&
      patmo_idx_H2O,&
      patmo_idx_OH,&
      patmo_idx_H2,&
      patmo_idx_H2O,&
      patmo_idx_H,&
      patmo_idx_OH,&
      patmo_idx_HSO,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO,&
      patmo_idx_O2,&
      patmo_idx_H,&
      patmo_idx_OH,&
      patmo_idx_SO3,&
      patmo_idx_SO2,&
      patmo_idx_O2,&
      patmo_idx_HO2,&
      patmo_idx_HO2,&
      patmo_idx_SO3,&
      patmo_idx_HSO3,&
      patmo_idx_H2SO4,&
      patmo_idx_SO2,&
      patmo_idx_O3,&
      patmo_idx_N,&
      patmo_idx_SO4,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS,&
      patmo_idx_S,&
      patmo_idx_CO,&
      patmo_idx_CS2,&
      patmo_idx_CO,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_CS,&
      patmo_idx_CS2E,&
      patmo_idx_SO,&
      patmo_idx_SO2,&
      patmo_idx_SH,&
      patmo_idx_S,&
      patmo_idx_COS,&
      patmo_idx_COS,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_CS2,&
      patmo_idx_SCSOH,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_CS,&
      patmo_idx_H2S,&
      patmo_idx_H2S,&
      patmo_idx_H2S,&
      patmo_idx_H2S,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SH,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_SO,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_S,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_HSO,&
      patmo_idx_HSO,&
      patmo_idx_HSO2,&
      patmo_idx_HSO3,&
      patmo_idx_SO2,&
      patmo_idx_SO2,&
      patmo_idx_SO3,&
      patmo_idx_H2SO4,&
      patmo_idx_O,&
      patmo_idx_N2,&
      patmo_idx_SO2,&
      patmo_idx_CH3SCH3,&
      patmo_idx_CH3SCH3,&
      patmo_idx_CH3SCH3,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_CS2E,&
      patmo_idx_S2,&
      patmo_idx_CS,&
      patmo_idx_SCSOH/)

end module patmo_commons
