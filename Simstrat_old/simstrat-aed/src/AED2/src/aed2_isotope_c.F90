!###############################################################################
!#                                                                             #
!# aed2_isotope_c.F90                                                          #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created 9 May 2013                                                          #
!###############################################################################

#include "aed2.h"
#define PURE
!
MODULE aed2_isotope_c
!-------------------------------------------------------------------------------
! aed2_isotope_c --- isotope_c biogeochemical model
!
! Created for Caboolture estuary to test O, N and C isotope_c submodel
! Matt Hipsey
! Modifications:
! - 07Nov2013 Inclusion of diagnostic variables for fluxes between pools (SA)
! - 05Feb2014 Inclusion of source and sink of DIC, DOC and POC pool (SA)
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,  ONLY: aed2_gas_piston_velocity, aed2_oxygen_sat

   IMPLICIT NONE

   PRIVATE    ! By default make everything private
!
   PUBLIC aed2_isotope_c_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_isotope_c_data_t
      !# Variable identifiers
      INTEGER  :: id_oxy, id_oxy18
      INTEGER  :: id_nit, id_nit15
      INTEGER  :: id_amm, id_amm15
      INTEGER  :: id_don, id_don15
      INTEGER  :: id_pon, id_pon15
      INTEGER  :: id_doc, id_doc13
      INTEGER  :: id_dic, id_dic13
      INTEGER  :: id_pH
      INTEGER  :: id_Fsed_amm,id_Fsed_nit,id_Fsed_oxy
      INTEGER  :: id_ta
      !INTEGER  :: id_ret  ! retention time
      INTEGER  :: id_chl  ! chlorophyll-a
      INTEGER  :: id_pCO2m ! pCO2 measured

      INTEGER  :: id_temp, id_salt
      INTEGER  :: id_wind
      INTEGER  :: id_nitrif,id_denit,id_sed_amm,id_sed_nit
      INTEGER  :: id_del13CDOC,id_del13CDIC
      INTEGER  :: id_oxy_sat
! In/out Flux
      INTEGER  :: id_DOCsed2wat, id_DICsed2wat
      INTEGER  :: id_DOC13sed2wat, id_DIC13sed2wat
      INTEGER  :: id_DICwat2atm, id_DICatm2wat
      INTEGER  :: id_DIC13wat2atm, id_DIC13atm2wat

      INTEGER  :: id_DIC2POC, id_DIC132POC13     ! Primary Production
      INTEGER  :: id_POC2DIC, id_POC132DIC13     ! Respiration
      INTEGER  :: id_DOC2DIC, id_DOC132DIC13     ! Mineralisation
      INTEGER  :: id_PhyMor2DOC, id_PhyMor2DOC13 ! PhytoMortality
      INTEGER  :: id_PhyEx2DOC, id_PhyEx2DOC13   ! Exudation
      INTEGER  :: id_POC2DOCsed, id_POC132DOC13sed! Sedimentation
      INTEGER  :: id_Resp2DIC, id_Resp2DIC13     ! Respiration to DIC
      INTEGER  :: id_Det2DOC, id_Det2DOC13       ! Decomposition Det to DOC
      INTEGER  :: id_PPDIC, id_PPDIC13           ! PP DIC uptake

      INTEGER  :: id_oxy_exch
      INTEGER  :: id_totC,id_totO


      !# Model parameters
      AED_REAL :: Rnitrif,Rdenit,Rommin,         &
                  Rmort,Rdecomp,Rresus,          &
                  fdoc,fdetr,Cupt,CphyPOC,       &
                  maxU,fexud,                    &
                  Fsed_oxy,Fsed_amm,Fsed_nit,    &
                  Fsed_don,Fsed_doc,Fsed_dic,    &
                  Knitrif,Kdenit,Kommin,         &
                  Ksed_amm,Ksed_nit,Ksed_oxy,    &
                  theta_nitrif,theta_denit,theta_ommin, &
                  theta_sed,atmco2,ionic,        &
                  alpha_doc_sed,alpha_dic_sed,   &
                  alpha_dic_phy,                 &
                  alpha_doc_min,alpha_dic_atm,   &
                  alpha_doc_decomp,alpha_dic_resp, &
                  alpha_doc_mort,alpha_doc_exud

      LOGICAL  :: use_sed_model

     CONTAINS
         PROCEDURE :: define            => aed2_define_isotope_c
         PROCEDURE :: calculate_surface => aed2_calculate_surface_isotope_c
         PROCEDURE :: calculate         => aed2_calculate_isotope_c
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_isotope_c
!        PROCEDURE :: mobility          => aed2_mobility_isotope_c
!        PROCEDURE :: light_extinction  => aed2_light_extinction_isotope_c
!        PROCEDURE :: delete            => aed2_delete_isotope_c


   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_isotope_c(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in)                       :: namlst
!
!LOCALS

   AED_REAL          :: oxy_initial=4.5, oxy18_initial=4.5
   AED_REAL          :: nit_initial=4.5, nit15_initial=4.5
   AED_REAL          :: amm_initial=4.5, amm15_initial=4.5
   AED_REAL          :: don_initial=4.5, don15_initial=4.5
   AED_REAL          :: doc_initial=4.5, doc13_initial=4.5
   AED_REAL          :: dic_initial=4.5, dic13_initial=4.5
   AED_REAL          :: ph_initial=7.5
   AED_REAL          :: alkalinity_initial=1000
   AED_REAL          :: chl_initial=3
   AED_REAL          :: pCO2m_initial=500
   !AED_REAL          :: retention_initial=0.0
   AED_REAL          :: atmco2   = 395.
   AED_REAL          :: ionic    = 1.5
   AED_REAL          :: Fsed_don = 3.5
   AED_REAL          :: Fsed_doc = 3.5
   AED_REAL          :: Fsed_dic = 3.5
   AED_REAL          :: Fsed_oxy = 3.5
   AED_REAL          :: Ksed_oxy = 30.0
   AED_REAL          :: Rnitrif  = 0.01
   AED_REAL          :: Rdenit = 0.01
   AED_REAL          :: Rommin = 0.01
   AED_REAL          :: Rmort = 0.01
   AED_REAL          :: Rdecomp = 0.01
   AED_REAL          :: Rresus = 0.01
   AED_REAL          :: fdoc = 0.2
   AED_REAL          :: fdetr = 0.2
   AED_REAL          :: Cupt = 25
   AED_REAL          :: CphyPOC = 0.13
   AED_REAL          :: maxU = 0.2
   AED_REAL          :: fexud = 0.05
   AED_REAL          :: Fsed_amm = 3.5
   AED_REAL          :: Fsed_nit = 3.5
   AED_REAL          :: Knitrif = 150.0
   AED_REAL          :: Kdenit = 150.0
   AED_REAL          :: Kommin = 150.0
   AED_REAL          :: Ksed_amm = 30.0
   AED_REAL          :: Ksed_nit = 30.0
   AED_REAL          :: theta_nitrif = 1.0
   AED_REAL          :: theta_denit = 1.0
   AED_REAL          :: theta_ommin = 1.0
   AED_REAL          :: theta_sed = 1.04
   CHARACTER(len=64) :: nitrif_reactant_variable=''
   CHARACTER(len=64) :: denit_product_variable=''
   CHARACTER(len=64) :: Fsed_amm_variable=''
   CHARACTER(len=64) :: Fsed_nit_variable=''
   CHARACTER(len=64) :: Fsed_oxy_variable=''
   AED_REAL          :: alpha_doc_sed = 1.03
   AED_REAL          :: alpha_dic_sed = 1.03
   AED_REAL          :: alpha_doc_min = 1.04
   AED_REAL          :: alpha_dic_phy = 1.04
   AED_REAL          :: alpha_dic_atm = 1.04
   AED_REAL          :: alpha_doc_decomp = 1.01
   AED_REAL          :: alpha_dic_resp  = 0.97
   AED_REAL          :: alpha_doc_mort = 1.03
   AED_REAL          :: alpha_doc_exud = 1.03
   INTEGER           :: status

   NAMELIST /aed2_isotope_c/ oxy_initial, oxy18_initial, &
                          nit_initial, nit15_initial,   &
                          amm_initial, amm15_initial,   &
                          don_initial, don15_initial,   &
                          doc_initial, doc13_initial,   &
                          dic_initial, dic13_initial,   &
                          ph_initial,alkalinity_initial,&
                          chl_initial, pCO2m_initial,   &
                          atmco2,ionic,                 &
                          Rnitrif,Rdenit,Rommin,        &
                          Rmort,Rdecomp,Rresus,         &
                          fdoc,fdetr,Cupt,CphyPOC,      &
                          maxU,fexud,                   &
                          Knitrif,Kdenit,Kommin,        &
                          Fsed_oxy,Fsed_amm,Fsed_nit,   &
                          Fsed_don,Fsed_doc,Fsed_dic,   &
                          Ksed_oxy,Ksed_amm,Ksed_nit,   &
                          theta_nitrif,theta_denit,theta_ommin,    &
                          theta_sed,                    &
                          alpha_doc_sed,alpha_dic_sed,  &
                          alpha_doc_min,  &
                          alpha_dic_phy,alpha_dic_atm,  &
                          alpha_doc_decomp,alpha_dic_resp, &
                          alpha_doc_mort,alpha_doc_exud

!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_isotope_c,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist aed2_isotope_c'
      STOP
   ENDIF

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   ! OXYGEN
   data%Fsed_oxy      = Fsed_oxy/secs_per_day
   data%Ksed_oxy      = Ksed_oxy

   ! NITROGEN
   data%Rnitrif  = Rnitrif/secs_per_day
   data%Rdenit   = Rdenit/secs_per_day
   data%Fsed_amm = Fsed_amm/secs_per_day
   data%Fsed_nit = Fsed_nit/secs_per_day
   data%Fsed_don = Fsed_don/secs_per_day
   data%Knitrif  = Knitrif
   data%Kdenit   = Kdenit
   data%Kommin   = Kommin
   data%Ksed_amm  = Ksed_amm
   data%Ksed_nit  = Ksed_nit
   data%theta_nitrif = theta_nitrif
   data%theta_denit  = theta_denit
   data%theta_ommin  = theta_ommin
   data%theta_sed = theta_sed

   ! CARBON
   data%Fsed_doc = Fsed_doc/secs_per_day
   data%Fsed_dic = Fsed_dic/secs_per_day
   data%Rommin   = Rommin/secs_per_day
   data%Rmort    = Rmort/secs_per_day
   data%Rdecomp  = Rdecomp/secs_per_day
   data%Rresus   = Rresus/secs_per_day
   data%Cupt = Cupt/secs_per_day
   data%maxU = maxU/secs_per_day
   data%atmco2 = atmco2
   data%ionic = ionic
   data%fdoc = fdoc
   data%fdetr = fdetr
   data%CphyPOC = CphyPOC
   data%fexud = fexud
   data%alpha_doc_sed = alpha_doc_sed
   data%alpha_dic_sed = alpha_dic_sed
   data%alpha_doc_min = alpha_doc_min
   data%alpha_dic_phy = alpha_dic_phy
   data%alpha_dic_atm = alpha_dic_atm
   data%alpha_doc_decomp = alpha_doc_decomp
   data%alpha_dic_resp = alpha_dic_resp
   data%alpha_doc_mort = alpha_doc_mort
   data%alpha_doc_exud = alpha_doc_exud

   ! Register state variables
   data%id_oxy = aed2_define_variable( 'oxy','mmol/m**3','oxygen',              &
                                    oxy_initial,minimum=zero_)
   data%id_oxy18 = aed2_define_variable( 'oxy18','mmol/m**3','oxygen',          &
                                    oxy_initial,minimum=zero_)
   data%id_nit = aed2_define_variable( 'nit','mmol/m**3','nitrate',             &
                                    nit_initial,minimum=zero_)
   data%id_nit15 = aed2_define_variable( 'nit15','mmol/m**3','nitrate',         &
                                    nit15_initial,minimum=zero_)
   data%id_amm = aed2_define_variable( 'amm','mmol/m**3','ammonium',            &
                                    amm_initial,minimum=zero_)
   data%id_amm15 = aed2_define_variable( 'amm15','mmol/m**3','ammonium',        &
                                    amm15_initial,minimum=zero_)
   data%id_don = aed2_define_variable( 'don','mmol/m**3','don',                 &
                                    don_initial,minimum=zero_)
   data%id_don15 = aed2_define_variable( 'don15','mmol/m**3','don',             &
                                    don15_initial,minimum=zero_)
   data%id_doc = aed2_define_variable( 'doc','mmol/m**3','doc',                 &
                                    doc_initial,minimum=zero_)
   data%id_doc13 = aed2_define_variable( 'doc13','mmol/m**3','doc',             &
                                    doc13_initial,minimum=zero_)
   data%id_dic = aed2_define_variable( 'dic','mmol/m**3','dic',                 &
                                    dic_initial,minimum=zero_)
   data%id_dic13 = aed2_define_variable( 'dic13','mmol/m**3','dic',             &
                                    dic13_initial,minimum=zero_)
   data%id_ph = aed2_define_variable( 'pH','pH unit','pH',                      &
                                    ph_initial,minimum=zero_)
   data%id_ta = aed2_define_variable( 'alkalinity',' mmol/m**3','alkalinity',   &
                                    alkalinity_initial,minimum=zero_)
   !data%id_ret = aed2_define_variable( 'retention','time','water_retention',   &
   !                                 zero_,minimum=zero_)
   data%id_chl = aed2_define_variable( 'chl','ug/L','chlorophyll-a',   &
                                    chl_initial,minimum=zero_)
   data%id_pCO2m = aed2_define_variable( 'pCO2m','uAtm','pCO2 measured',   &
                                    pCO2m_initial,minimum=zero_)

! Register external state variable dependencies
!   data%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
!   IF (data%use_oxy) THEN
!     data%id_oxy = data%register_state_dependency(nitrif_reactant_variable)
!!    data%id_oxy18 = data%register_state_dependency(nitrif_reactant_variable)
!   ENDIF
!   data%use_no2 = denit_product_variable .NE. '' !This means n2 module switched on
!!  IF (data%use_no2) data%id_denit_product = data%register_state_dependency(denit_product_variable)
!
!   data%use_sed_model = Fsed_amm_variable .NE. ''
!   IF (data%use_sed_model) THEN

!     data%id_Fsed_amm = data%register_state_dependency(Fsed_amm_variable,benthic=.true.)
!     data%id_Fsed_nit = data%register_state_dependency(Fsed_nit_variable,benthic=.true.)
!   ENDIF

!
   data%id_del13CDOC = aed2_define_diag_variable( 'del13CDOC','o/oo',        &
                     'del ratio of 13C-DOC')

   data%id_del13CDIC = aed2_define_diag_variable( 'del13CDIC','o/oo',        &
                     'del ratio of 13C-DIC')

! Fluxes

   data%id_DOCsed2wat = aed2_define_diag_variable( 'DOCsed2wat','mmol/m**3',  &
                     'DOC release')
   data%id_DOC13sed2wat = aed2_define_diag_variable( 'DOC13sed2wat','mmol/m**3',  &
                     'DOC13 release')
   data%id_DICsed2wat = aed2_define_diag_variable( 'DICsed2wat','mmol/m**3',  &
                     'DIC release')
   data%id_DIC13sed2wat = aed2_define_diag_variable( 'DIC13sed2wat','mmol/m**3',  &
                     'DIC13 release')
   data%id_DICwat2atm = aed2_define_diag_variable( 'DICwat2atm','mmol/m**3',  &
                     'DIC outflux')
   data%id_DICatm2wat = aed2_define_diag_variable( 'DICatm2wat','mmol/m**3',  &
                     'DIC influx')
   data%id_DIC13wat2atm = aed2_define_diag_variable( 'DIC13wat2atm','mmol/m**3',  &
                     'DIC13 outflux')
   data%id_DIC13atm2wat = aed2_define_diag_variable( 'DIC13atm2wat','mmol/m**3',  &
                     'DIC13 influx')
   data%id_DOC2DIC = aed2_define_diag_variable( 'DOC2DIC','mmol/m**3',  &
                     'DOC convert to DIC')
   data%id_DOC132DIC13 = aed2_define_diag_variable( 'DOC132DIC13','mmol/m**3',  &
                     'DOC13 convert to DIC13')
   data%id_PhyMor2DOC = aed2_define_diag_variable( 'PhyMor2DOC','mmol/m**3',  &
                      'POC Phyto Mortality to DOC')
   data%id_PhyMor2DOC13 = aed2_define_diag_variable( 'PhyMor2DOC13','mmol/m**3',  &
                      'POC Phyto Mortality to DOC13')
   data%id_PhyEx2DOC = aed2_define_diag_variable( 'PhyEx2DOC','mmol/m**3',  &
                      'POC Phyto Excudates to DOC')
   data%id_PhyEx2DOC13 = aed2_define_diag_variable( 'PhyEx2DOC13','mmol/m**3', &
                      'POC Phyto Excudates to DOC13')
   data%id_Det2DOC = aed2_define_diag_variable( 'Det2DOC','mmol/m**3', &
                      'Detritus to DOC')
   data%id_Det2DOC13 = aed2_define_diag_variable( 'Det2DOC13','mmol/m**3', &
                      'Detritus to DOC13')
   data%id_Resp2DIC = aed2_define_diag_variable( 'Resp2DIC','mmol/m**3', &
                      'Resp to DIC')
   data%id_Resp2DIC13 = aed2_define_diag_variable( 'Resp2DIC13','mmol/m**3', &
                      'Resp to DIC13')
   data%id_PPDIC = aed2_define_diag_variable( 'PPDIC','mmol/m**3', &
                      'PP to DIC')
   data%id_PPDIC13 = aed2_define_diag_variable( 'PPDIC13','mmol/m**3', &
                      'PP to DIC13')

   data%id_oxy_sat = aed2_define_diag_variable( 'oxySat','mmol/m2','oxy saturation')

   data%id_oxy_exch = aed2_define_sheet_diag_variable('oxyAtmFlux',     &
                      'mmol/m2','oxy_atm_flux')



   ! Register environmental dependencies

    data%id_temp = aed2_locate_global( 'temperature')
    data%id_salt = aed2_locate_global( 'salinity')
    data%id_wind = aed2_locate_global( 'wind_speed')

END SUBROUTINE aed2_define_isotope_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_isotope_c(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_isotope_c model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: temp,salt,wind   !State variables
   AED_REAL           :: oxy,amm,nit,don,doc,dic !State variables
   AED_REAL           :: doc13,dic13 !State variables
   AED_REAL           :: alkalinity
   !AED_REAL          :: retention
   AED_REAL           :: chl,pCO2m
   AED_REAL           :: Tabs
   AED_REAL           :: pCO2a,pCO2d,Ko,kCO2,pH,HCO3,CO3,windHt
   AED_REAL           :: nitrification,denitrification,mineralisation,resuspension
   AED_REAL           :: DOC2DIC,DOC132DIC13,CPhy,CMort,PhyMor2DOC,PhyMor2DOC13
   AED_REAL           :: PhyEx2DOC,PhyEx2DOC13,Det2DOC,Det2DOC13
   AED_REAL           :: DOC132DIC13_Rat,Det2DOC13_Rat,Resp2DIC13_Rat,PPDIC13_Rat
   AED_REAL           :: PhyEx2DOC13_Rat,PhyMor2DOC13_Rat
   AED_REAL           :: Cdetr,Resp2DIC,Resp2DIC13,PPDIC,PPDIC13,alphaDICpp
   AED_REAL           :: del13CDIC,del13CDOC
   AED_REAL,PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to ammonium utilised during nitrification
   AED_REAL,PARAMETER :: Yoxy_ommin  = 1. !ratio of oxygen to carbion utilised during mineralisation
   AED_REAL,PARAMETER :: CChl = 50. ! C/Chl ratio
   AED_REAL,PARAMETER :: del13CPOC = -23 ! median d13CPOC during CR2
   AED_REAL,PARAMETER :: del13CPhy = -34 ! Hamilton & Lewis (1992)
   AED_REAL,PARAMETER :: del13CDet = -28 ! d13C Detritus (Otero et al, 2000)
   AED_REAL,PARAMETER :: Rf_13C = 0.011180 !(VPDB), Fry(2006)Appendix 1, p. 285
   AED_REAL,PARAMETER :: Rf_15N = 0.0036765 !(VPDB)
   AED_REAL,PARAMETER :: Rf_18O = 0.0020052 !(SMOW)
   AED_REAL,PARAMETER :: Rf_17O = 0.0003799 !(SMOW)
   !AED_REAL,PARAMETER :: ret_dt = 1.0 !per time step
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen
   amm = _STATE_VAR_(data%id_amm)! ammonium
   nit = _STATE_VAR_(data%id_nit)! nitrate
   don = _STATE_VAR_(data%id_don)! dissolved organic nitrogen
   doc = _STATE_VAR_(data%id_doc)! dissolved organic carbon
   dic = _STATE_VAR_(data%id_dic)! dissolved inorganic carbon
   doc13 = _STATE_VAR_(data%id_doc13)! dissolved organic carbon C13
   dic13 = _STATE_VAR_(data%id_dic13)! dissolved inorganic carbon C13
   alkalinity = _STATE_VAR_(data%id_ta)    ! alkalinity
   pH = _STATE_VAR_(data%id_pH)            ! pH
   !retention = _STATE_VAR_(data%id_ret)   ! retention time
   chl = _STATE_VAR_(data%id_chl)          ! chlorophyll-a
   pCO2m = _STATE_VAR_(data%id_pCO2m)      ! measured pCO2


   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp) ! temperature [C]
   salt = _STATE_VAR_(data%id_salt) ! salinity    [psu]
   wind = _STATE_VAR_S_(data%id_wind)    ! Wind speed at 10 m above surface (m/s)

  !-------------------------
   !correct sim temp (overheating about 1 to 2 degC)
    temp=temp-1
   !------------------------
!print *, '========================'
!print *, '=== C.isotope_c_do start=== dic13=',dic13,'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3


  ! Define some intermediate rates units [/day]
   nitrification   = fnitrif(data,oxy,temp)
   denitrification = fdenit(data,oxy,temp)
   mineralisation  = fommin(data,oxy,temp)
   resuspension = fthetased(data,oxy,temp)

   ! Set temporal derivatives
   ! OXY
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + ( -Yoxy_nitrif*amm*nitrification - Yoxy_ommin*doc*mineralisation )
   ! AMM
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + ( -amm*nitrification + don*mineralisation)
   ! NIT
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + ( +amm*nitrification - nit*denitrification)
   ! DON
   _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + ( -don*mineralisation)

!-------------------------------------------------

   del13CDIC=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
   del13CDOC=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3

  ! Mineralisation DOC-->DIC
    DOC2DIC = doc*mineralisation
    del13Cdoc = (((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
    !DOC132DIC13_Rat = (del13Cdoc+1000)/(del13Cdoc+1000+(1000/Rf_13C)/data%alpha_doc_min)
    DOC132DIC13_Rat = (del13Cdoc/1000+1)*Rf_13C/data%alpha_doc_min
    DOC132DIC13 = DOC2DIC-(1/(1+DOC132DIC13_Rat)*DOC2DIC)
    !DOC132DIC13 = doc13*mineralisation*data%alpha_doc_min

  ! Mortality POC phyto-->DOC
    Cphy = CChl*chl/12 !uM
    CMort = data%Rmort*CPhy !uM
    PhyMor2DOC = data%fdoc*CMort !uM
    !PhyMor2DOC13_Rat = (del13CPOC+1000)/(del13CPOC+1000+(1000/Rf_13C)/data%alpha_doc_mort)
    PhyMor2DOC13_Rat = (del13CPOC/1000+1)*Rf_13C/data%alpha_doc_mort
    PhyMor2DOC13 = PhyMor2DOC-(1/(1+PhyMor2DOC13_Rat)*PhyMor2DOC)
   !PhyMor2DOC13 = (del13CPOC+1000)/(del13CPOC+1000+(1000/Rf_13C))*PhyMor2DOC*data%alpha_doc_mort

  ! Exudation POC phyto-->DOC
    PhyEx2DOC = data%fexud*(data%maxU*Cphy)
    !PhyEx2DOC13_Rat = (del13CPOC+1000)/(del13CPOC+1000+(1000/Rf_13C)/data%alpha_doc_exud);
    PhyEx2DOC13_Rat = (del13CPOC/1000+1)*Rf_13C/data%alpha_doc_exud
    PhyEx2DOC13 = PhyEx2DOC-(1/(1+PhyEx2DOC13_Rat)*PhyEx2DOC);
    !PhyEx2DOC13 = (del13CPOC+1000)/(del13CPOC+1000+(1000/Rf_13C))*PhyEx2DOC*data%alpha_doc_exud

  ! Decomposition of POC Detritus -->DOC
    Cdetr = (1-data%CphyPOC)*Cphy/data%CphyPOC ! C in detritus = (Det/POC)*Cphyto/(Cphyto/POC)
    Det2DOC = data%Rdecomp*Cdetr  !uM
    Det2DOC13_Rat = (del13CDet/1000+1)*Rf_13C/data%alpha_doc_decomp;
    Det2DOC13 = Det2DOC-(1/(1+Det2DOC13_Rat)*Det2DOC);
    !Det2DOC13 = (del13CDet+1000)/(del13CDet+1000+(1000/Rf_13C))*Det2DOC*data%alpha_doc_decomp

  ! Resuspension DOC from sediment
  !  resusDOC = docSed*resuspenstion
  !  resusDOC13 = (del13DOCSed+1000)/(del13CDOCSed+1000+(1000/Rf_13C))*resusDOC*data%alpha_doc_sed

  ! DIC due to dead phyto respiration
    Resp2DIC = (1-data%fDetr-data%fdoc)*CMort
    Resp2DIC13_Rat = (del13CPOC/1000+1)*Rf_13C/data%alpha_dic_resp
    !Resp2DIC13_Rat = (del13CPOC+1000)/(del13CPOC+1000+(1000/Rf_13C)/data%alpha_dic_resp)
    Resp2DIC13 = Resp2DIC-(1/(1+Resp2DIC13_Rat)*Resp2DIC);
    !Resp2DIC13 = (del13CPOC+1000)/(del13CPOC+1000+(1000/Rf_13C))*Resp2DIC*data%alpha_dic_resp

  ! -DIC due to DIC uptake------------------
    ! 1st calculate Solubility, Ko (mol/L/atm) (Weiss, 1974)
     Tabs = temp + 273.15
     Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
                  + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
     Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
     Ko = exp(Ko)

    ! 2nd calculate gas piston velocity
    ! windHt = 10.  !reference wind measurement height
     windHt = 48-15.  ! wind measurement height

     kCO2 = aed2_gas_piston_velocity(windHt,wind,temp,salt)
    ! 3rd dissolved pCO2

    pCO2a = aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3) / Ko !pCO2 in atm

    !print *, 'before adjustment, pCO2=', pCO2a, 'pCO2m=', pCO2m

    !temporary to fix the ratio measured and calculated pCO2a
    pCO2a = 1.5*pCO2a

    pCO2d =  kCO2 * (Ko*1e6) * pCO2a ! pCO2 in uM

    ! 4th calculate C uptake, etc
    PPDIC = chl*data%Cupt/12
    alphaDICpp = 1.0007+(0.0283*0.69*(pCO2d-5.6)/pCO2d) ! (Yoshioka,1997)
    PPDIC13_Rat = (del13CDIC/1000+1)*Rf_13C/alphaDICpp
    !PPDIC13_Rat = (del13CDIC+1000)/(del13CDIC+1000+(1000/Rf_13C)/alphaDICpp)
    PPDIC13 = PPDIC-(1/(1+PPDIC13_Rat)*PPDIC)
    !PPDIC13 = (del13CDIC+1000)/(del13CDIC+1000+(1000/Rf_13C))*PPDIC*alphaDICpp

!print *, 'DOC2DIC=',DOC2DIC, 'PhyMor2DOC=',PhyMor2DOC,'PhyEx2DOC=',PhyEx2DOC,'Det2DOC=',Det2DOC
!print *, 'PPDIC=',' DOC2DIC,Resp2DIC,PPDIC=', DOC2DIC,Resp2DIC,PPDIC, DOC2DIC+Resp2DIC-PPDIC
!print *, 'del13CDOC=', (((doc13/(doc-doc13))/Rf_13C)-1.)*1e3, 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3

!--Enabling ------------------------------------
!
   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + ( -DOC2DIC+(PhyMor2DOC+PhyEx2DOC)+Det2DOC)
   _FLUX_VAR_(data%id_doc13) = _FLUX_VAR_(data%id_doc13) + ( -DOC132DIC13+(PhyMor2DOC13+PhyEx2DOC13)+Det2DOC13)

   _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + ( +DOC2DIC+Resp2DIC-PPDIC)
   _FLUX_VAR_(data%id_dic13) = _FLUX_VAR_(data%id_dic13) + ( +DOC132DIC13+Resp2DIC13-PPDIC13)

   _DIAG_VAR_(data%id_del13CDOC) =  del13CDOC
   _DIAG_VAR_(data%id_del13CDIC) =  del13CDIC

!-------------------------------------------------
! temporary to simulate states as tracers
!
!   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + ( 0)
!   _FLUX_VAR_(data%id_doc13) = _FLUX_VAR_(data%id_doc13) + ( 0)

!   _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + ( 0)
!   _FLUX_VAR_(data%id_dic13) = _FLUX_VAR_(data%id_dic13) + ( 0)

!   _DIAG_VAR_(data%id_del13CDOC) =  del13CDOC
!   _DIAG_VAR_(data%id_del13CDIC) =  del13CDIC
!------------------

   !Fluxes
   _DIAG_VAR_(data%id_DOC2DIC) =  DOC2DIC
   _DIAG_VAR_(data%id_DOC132DIC13) =  DOC132DIC13
   _DIAG_VAR_(data%id_PhyMor2DOC) =  PhyMor2DOC
   _DIAG_VAR_(data%id_PhyMor2DOC13) =  PhyMor2DOC13
   _DIAG_VAR_(data%id_PhyEx2DOC) =  PhyEx2DOC
   _DIAG_VAR_(data%id_PhyEx2DOC13) =  PhyEx2DOC13
   _DIAG_VAR_(data%id_Det2DOC) =  Det2DOC
   _DIAG_VAR_(data%id_Det2DOC13) =  Det2DOC13
   _DIAG_VAR_(data%id_Resp2DIC) =  Resp2DIC
   _DIAG_VAR_(data%id_Resp2DIC13) =  Resp2DIC13
   _DIAG_VAR_(data%id_PPDIC) =  PPDIC
   _DIAG_VAR_(data%id_PPDIC13) =  PPDIC13



!print *, '=== C.isotope_c_do ending ==='
!print *, 'oxy=', oxy, 'doc=', doc
!print *, 'pCO2d=', pCO2d, '%diff pCO2 (meas-calc)=', (pCO2m-pCO2a)/pCO2m*100

!print *, 'dic=',dic,'doc=',doc, 'dic13=', dic13, 'doc13=', doc13
!print *, 'del13CDOC=', (((doc13/(doc-doc13))/Rf_13C)-1.)*1e3, 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
!print *, 'dic=',dic, 'dic13=', dic13, 'del13C=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3

!print *, 'DICchange=',+DOC2DIC+Resp2DIC-PPDIC, 'DIC13change=',+DOC132DIC13+Resp2DIC13-PPDIC13
!print *, 'DOCchange=',-DOC2DIC+(PhyMor2DOC+PhyEx2DOC)+Det2DOC, 'DOC13change=',-DOC132DIC13+(PhyMor2DOC13+PhyEx2DOC13)+Det2DOC13

if ((((dic13/(dic-dic13))/Rf_13C)-1.)*1e3<-30) then
 print *,  'aaa...error.. dic=',dic, 'dic13=', dic13, 'del13C_dic=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
 STOP
endif


!if ((((doc13/(doc-doc13))/Rf_13C)-1.)*1e3>-10) then
! print *,  'error.. doc=',doc, 'doc13=', doc13, 'del13C_doc=',(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
! STOP
!endif


if ((((dic13/(dic-dic13))/Rf_13C)-1.)*1e3>25) then
 print *,  'bbb error.. dic=',dic, 'dic13=', dic13, 'del13C_dic=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
 STOP
endif


if (doc>4000) then
 print *,  'error.. doc=',doc, 'doc13=', doc13, 'del13C_doc=',(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
 STOP
endif


END SUBROUTINE aed2_calculate_isotope_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_isotope_c(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED isotope_c.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: amm,nit,oxy,don,doc,dic
   AED_REAL :: doc13,dic13
   AED_REAL :: pH ! for testing only

   ! Temporary variables
   AED_REAL :: oxy_flux, Fsed_oxy
   AED_REAL :: amm_flux,nit_flux,don_flux
   AED_REAL :: Fsed_amm,Fsed_nit,Fsed_don
   AED_REAL :: doc_flux,dic_flux
   AED_REAL :: Fsed_doc,Fsed_dic
   AED_REAL :: del13CDIC,del13CDOC
   AED_REAL :: doc13_flux_Rat,dic13_flux_Rat

   ! Parameters
   AED_REAL,PARAMETER :: d13DOCSed =-28
   AED_REAL,PARAMETER :: Rf_13C = 0.011180 !(VPDB), Fry(2006)Appendix 1, p. 285

!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen
   amm = _STATE_VAR_(data%id_amm)! ammonium
   nit = _STATE_VAR_(data%id_nit)! nitrate
   don = _STATE_VAR_(data%id_don)! dissolved organic nitrogen
   doc = _STATE_VAR_(data%id_doc)! dissolved organic carbon
   dic = _STATE_VAR_(data%id_dic)! dissolved inorganic carbon
   doc13 = _STATE_VAR_(data%id_doc13)! doc13
   dic13 = _STATE_VAR_(data%id_dic13)! dic13
   pH = _STATE_VAR_(data%id_pH)! pH for testing only


   !-----------------------
   !correct sim temp (overheating about 1 to 2 degC)
   ! temp=temp-1
   !------------------------

!print *, 'B1. Do Benthos Start, dic=', dic, 'dic13=', dic13, 'del13C=', (((dic13/(dic-dic13))/0.011180)-1.)*1e3

   Fsed_oxy = data%Fsed_oxy
   Fsed_amm = data%Fsed_amm
   Fsed_nit = data%Fsed_nit
   Fsed_don = data%Fsed_don
   Fsed_doc = data%Fsed_doc
   Fsed_dic = data%Fsed_dic

   ! Sediment flux dependent on oxygen and temperature
   oxy_flux = Fsed_oxy * oxy/(data%Ksed_oxy+oxy) * (data%theta_sed**(temp-20.0))
   amm_flux = Fsed_amm * data%Ksed_amm/(data%Ksed_amm+oxy) * (data%theta_sed**(temp-20.0))
   nit_flux = Fsed_nit * oxy/(data%Ksed_nit+oxy) * (data%theta_sed**(temp-20.0))
   don_flux = Fsed_don * (data%theta_sed**(temp-20.0))
   doc_flux = Fsed_doc * (data%theta_sed**(temp-20.0))
   !dic_flux = -oxy_flux
   dic_flux = -0.8*oxy_flux

   !calculate delC before fluxes
   del13Cdic=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
   del13Cdoc=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (oxy_flux)
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (amm_flux)
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + (nit_flux)
   _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + (don_flux)
   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + (doc_flux)
   _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (dic_flux)

   !Calculate C13 Ratios after fractionation
    doc13_flux_Rat=(del13CDOC/1000+1)*Rf_13C/data%alpha_doc_sed
    dic13_flux_Rat=(del13CDIC/1000+1)*Rf_13C/data%alpha_dic_sed
   _FLUX_VAR_(data%id_doc13) = _FLUX_VAR_(data%id_doc13) + (doc_flux-(1/(1+doc13_flux_Rat)*doc_flux))
   _FLUX_VAR_(data%id_dic13) = _FLUX_VAR_(data%id_dic13) + (dic_flux-(1/(1+dic13_flux_Rat)*dic_flux))

! o.m is degraded, release isotopically light C, lowering delC13 of dic in depth

   _DIAG_VAR_(data%id_DICsed2wat) =  dic_flux
   _DIAG_VAR_(data%id_DOCsed2wat) =  doc_flux
   _DIAG_VAR_(data%id_DIC13sed2wat) =  dic_flux-(1/(1+dic13_flux_Rat)*dic_flux)
   _DIAG_VAR_(data%id_DOC13sed2wat) =  doc_flux-(1/(1+doc13_flux_Rat)*doc_flux)


!print *, '== B2. set bottom exhange'
!print *, 'doc=', doc, 'doc13=', doc13, 'del13CDOC=',del13CDOC
!print *, 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',del13CDIC



!print *, 'Do benthos end dic=', dic, 'dic13=', dic13, 'del13CDIC=',del13Cdic
!print *, 'Do benthos end doc=', doc, 'oxy=', oxy

END SUBROUTINE aed2_calculate_benthic_isotope_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed2_calculate_surface_isotope_c(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed oxygen model
!
! see aed2_oxygen_do_surface_exchange
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, pCO2m

   ! State
   AED_REAL :: oxy,dic,pH,alkalinity,dic13,atm13C

   ! Temporary variables
   AED_REAL :: oxy_atm_flux = zero_
   AED_REAL :: Coxy_air = zero_ !Dissolved oxygen in the air phase
   AED_REAL :: koxy_trans = zero_
   AED_REAL :: f_pres  = 1.0      ! Pressure correction function only applicable at high altitudes
   AED_REAL :: pCO2,FCO2
   AED_REAL :: Ko, KCO2
   AED_REAL :: Tabs,windHt
   AED_REAL :: Faw, Fwa, FC13
   AED_REAL :: HCO3, CO3, CO2
   AED_REAL :: Rdic, RCO2a
   AED_REAL :: eta_dic_g, eta_HCO3_g, eta_CO3_g, eta_aq_g, eta_as
   AED_REAL :: alpha_dic_g, alpha_HCO3, alpha_CO3, alpha_aq_g, alpha_k
   AED_REAL :: FC13wat, FC13air

   AED_REAL, PARAMETER :: Rf_13C = 0.011180 !(VPDB), Fry(2006)Appendix 1, p. 285
   AED_REAL, PARAMETER :: delC13_atm = -8 !o/oo, delC13 of atmospheric CO2

   AED_REAL, PARAMETER :: eta_k = -0.81 !o/oo, kinetic isotope_c fract during CO2 gas transfer  (Zhang, 1995)
   !AED_REAL, PARAMETER :: eta_aq_g = -1.19 !o/oo, equilibrium fractination during CO2 dissolution (Zhang, 1995)

!
!-------------------------------------------------------------------------------
!BEGIN

   !---- First, define environment context
   !Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind = _STATE_VAR_S_(data%id_wind)    ! Wind speed at 10 m above surface (m/s)
   !wind) ! Wind speed at 10 m above surface (m/s = _STATE_VAR_S_(data%id_wind)

   !windHt = 10.  !reference wind measurement height
   windHt = 48-15.  ! wind measurement height

   !-----------------------
   !correct sim temp (overheating about 1 to 2 degC)
    temp=temp-1
   !------------------------

   !---- Second, flux oxygen
    ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! Concentration of oxygen in surface layer

   koxy_trans = aed2_gas_piston_velocity(windHt,wind,temp,salt)

!print *, 'koxy_trans=', koxy_trans, 'wind=', wind, 'temp=', temp, 'salt=', salt

   ! First get the oxygen concentration in the air phase at interface
   ! Taken from Riley and Skirrow (1974)
   f_pres = 1.0
   Coxy_air = f_pres * aed2_oxygen_sat(salt,temp)

   ! Get the oxygen flux
   oxy_atm_flux = koxy_trans * (Coxy_air - oxy)

   ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
   _FLUX_VAR_T_(data%id_oxy) = oxy_atm_flux

!print *, 'Coxy_air=',Coxy_air, 'oxy_atm_flux=',oxy_atm_flux,'oxy=',oxy

   !---- Third, flux CO2 --------
   ! Retrieve current (local) state variable values.
   dic = _STATE_VAR_(data%id_dic)          ! Concentration of inorganic carbon in surface layer
   dic13 = _STATE_VAR_(data%id_dic13)      ! Concentration of inorganic carbon13 in surface layer
   pH = _STATE_VAR_(data%id_pH)            ! pH in surface layer
   alkalinity = _STATE_VAR_(data%id_ta)    ! alkalinity in surface layer

!print *,'======================================'
!print *, '== A. Start surface exchange ', 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',((dic13/(dic-dic13)/Rf_13C)-1.)*1e3


   kCO2 = aed2_gas_piston_velocity(windHt,wind,temp,salt)

   ! Solubility, Ko (mol/L/atm) (Weiss, 1974)
   Tabs = temp + 273.15
   Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
                  + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
   Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
   Ko = exp(Ko)

   ! pCO2 in surface water layer

   !print *, '== A. surface exchange ', 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3, 'pH=', pH

   pCO2 = aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3) / Ko

   !print *, '== aed2_carbon_co2 , pH=',pH,' alkalinity=',alkalinity,' HCO3=',HCO3,' CO3=',CO3,' pCO2=',pCO2,' Ko=',Ko

   ! mmol/m2/s = m/s * mmol/L/atm * atm
   !FCO2 = kCO2 * Ko*1e6 * (pCO2 - data%atmco2)

   !------ Yanti correction (20/5/2013)
   ! pCO2 is actually in uatm (=ppm)
     FCO2 = - kCO2 * Ko*1e6 * (pCO2 - data%atmco2) *1e-6 ! dCO2/dt

   ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
   _FLUX_VAR_T_(data%id_dic) = FCO2

!------------------------------------------------------------------------
!
!   !Lynch-Stieglitz etal (1995)
!   Faw=kCO2 * Ko*1e6 * data%atmco2 *1e-6        !atm to surface
!   Fwa=kCO2 * Ko*1e6 * pCO2 *1e-6               !water to atm
!
!   atm13C = (delC13_atm+1000)/((delC13_atm+1000+(1000/Rf_13C))) !note this is fraction of 13C_atm,13C/C
!
!   FC13=((-0.373/Tabs+1.00019)*data%alpha_dic_atm*Faw*atm13C)-((-9.866/Tabs+1.02412)*data%alpha_dic_atm*Fwa*dic13/dic)
!---------------------------------------------------------------------------
!
!   Based on Zhang (1995)
!   Note, Zhang definition: alpha=k13/k (not k12/k13), R=13C/12C
!
!   First, calculate the individual temperature-dependent fractination of each CO2, HCO3, CO3

    eta_HCO3_g = -(0.1141*temp)+10.78     ! fractination between bicarb and CO2
    eta_CO3_g = -(0.052*temp)+7.22        ! fractination between carbonate and CO2
    eta_aq_g = -(0.0049*temp)-1.31        ! fractination during gass dissolution CO2

!   Second, calculate kinetic isotope_c fractionation during CO2 exhange (alpha_k) and
!   equilibrium fractionation during CO2(g) dissolution (alpha_aq_g)

    alpha_k = (eta_k/1000)+1
    alpha_aq_g = (eta_aq_g/1000)+1

!   Third, calculate total fractination during CO2 expansion
    eta_as = (alpha_k*alpha_aq_g-1)*1000

!   Fourth, calculate the equilibrium fractination factor between DIC and CO2(g)
    CO2=aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3)
    eta_dic_g = (CO2/dic*eta_as)+(HCO3/dic*eta_HCO3_g)+(CO3/dic*eta_CO3_g)
    alpha_dic_g = (eta_dic_g/1000)+1

!   Last, calculate the d13CO2/dt
    Rdic = dic13/(dic-dic13) ! 13C/12C ratio of dic in water
    RCO2a = (delC13_atm/1000+1)*Rf_13C  ! 13C/12C ratio of CO2 in atmosphere

    !print *,'Rdic=',Rdic,'RCO2a=',Rco2a,'kCO2=',kco2

    FC13 = - kCO2*alpha_k*alpha_aq_g*(Ko*1e6)*((pCO2*1e-6)*Rdic/alpha_dic_g-(data%atmco2*1e-6)*RCO2a)

    FC13wat=-kCO2*alpha_k*alpha_aq_g*(Ko*1e6)*((pCO2*1e-6)*Rdic/alpha_dic_g)
    FC13air=-kCO2*alpha_k*alpha_aq_g*(Ko*1e6)*((data%atmco2*1e-6)*RCO2a)

!   Diagnosis

    _DIAG_VAR_(data%id_DICwat2atm) =  - kCO2 * Ko*1e6 * pCO2*1e-6
    _DIAG_VAR_(data%id_DICatm2wat) =  - kCO2 * Ko*1e6 * (-data%atmco2)*1e-6

    _DIAG_VAR_(data%id_DIC13wat2atm) =  FC13wat
    _DIAG_VAR_(data%id_DIC13atm2wat) =  FC13air

!print *, 'pCO2=', pCO2, 'atmco2=', data%atmco2,'(pCO2*1e-6)*Rdic/alpha_dic_g=', (pCO2*1e-6)*Rdic/alpha_dic_g, '(data%atmco2*1e-6)*RCO2a=', (data%atmco2*1e-6)*RCO2a
!print *, 'pCO2=', pCO2, 'atmco2=', data%atmco2

   _FLUX_VAR_T_(data%id_dic13) = FC13

!print *, '== A. End surface exchange ', 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3, 'FC13=',FC13

!----------------------
!   IF(FCO2>0.) THEN
!     ! pCO2 is high and flux is from the water to air
!     _FLUX_VAR_T_(data%id_dic13) = FCO2*(dic13/dic)/data%alpha_dic_atm !loosing
!
!   ELSE !FCO2<0., gaining
!     ! pCO2 is low and flux is from air to water
!     !atm13C = ((delC13_atm/1000.)+1)* Rf_13C !note, this is Rs
!      atm13C = (delC13_atm+1000)/((delC13_atm+1000+(1000/Rf_13C))) !note this is fraction of 13C_atm, 13C/C
!     _FLUX_VAR_T_(data%id_dic13) = FCO2*atm13C
!
!   END IF
!----------------------

!   if (FCO2>0) then
!      print *, 'influx,wtr=', FC13wat,'air=',FC13air,'FC13=',FC13,'Rs=', dic13/(dic-dic13),'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
!   else
!      print *, 'outflux,wtr=', FC13wat,'air=',FC13air,'FC13=',FC13,'Rs=', dic13/(dic-dic13), 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
!   endif

   !---- Fourth, set diagnostics

   ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).
   _DIAG_VAR_S_(data%id_oxy_exch) = oxy_atm_flux
   _DIAG_VAR_(data%id_oxy_sat) =  Coxy_air


END SUBROUTINE aed2_calculate_surface_isotope_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fnitrif(data,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                 :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   fnitrif = data%Rnitrif * oxy/(data%Knitrif+oxy) * (data%theta_nitrif**(temp-20.0))

END FUNCTION fnitrif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdenit(data,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                    :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN

   fdenit = data%Rdenit * data%Kdenit/(data%Kdenit+oxy) * (data%theta_denit**(temp-20.0))

END FUNCTION fdenit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fommin(data,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for organic matter mineralisation
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                   :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN

   fommin = data%Rommin * oxy/(data%Kommin+oxy) * (data%theta_ommin**(temp-20.0))

END FUNCTION fommin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fthetased(data,oxy,temp)
!-------------------------------------------------------------------------------
! Temperature function for resuspension
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                   :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN

   fthetased = data%Rresus * oxy/(data%Kommin+oxy) * (data%theta_sed**(temp-20.0))

END FUNCTION fthetased
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
!PURE AED_REAL FUNCTION aed2_carbon_co2(data,temp,dic,pH,alkalinity)
AED_REAL FUNCTION aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_c_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)                   :: dic, temp, alkalinity
   AED_REAL, INTENT(OUT)                  :: HCO3, CO3
   AED_REAL, INTENT(INOUT)                :: pH

!
!LOCALS
   ! Temporary variables
   AED_REAL :: pKh, pKa1, pKa2, pKw, TA
   AED_REAL :: Kh, Kw, Ka1, Ka2, is_f
   AED_REAL :: alf_CO2, alf_HCO3, alf_CO3
   AED_REAL :: H, CO2
!-------------------------------------------------------------------------------
!BEGIN

   ! Acidity constants temperature dependence (mol/kg atm)at ionic strength=0
   ! Note (Y13/4/13), have checked the values are within the range of values provided
   ! in Butler 1998 Table 10.1 (selected values of equilibrium cosntant extrapoalted to I=0)

   pKh = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
   pKa1 = 0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
   pKa2 = 0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
   pKw = 0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909

   ! Ionic strength dependence

   ! 1st calculate ionic strength function is_f
   is_f = (((SQRT(data%ionic)) / (1+SQRT(data%ionic))) - 0.20*data%ionic) * &
                       (298.0/(temp+273.))**0.666667

   ! pKh = pKh(0) + bI
   ! b = 0.105 (Butler, 1982)
   pKh = pKh + 0.105*data%ionic

   ! pKw = pKw(0) - f
   pKw = pKw - is_f

   ! pKa1 = pKa1(0) - f - bI
   pKa1 = pKa1 - is_f - 0.105*data%ionic

   !pKa2 = pKa2(0) - 2f
   !pKa2 = pKa2 + 2.0*is_f !YA
   pKa2 = pKa2 - 2.0*is_f

   ! Convert from pK etc to Kh, Kw, Ka1, Ka2
   Kh  = 10.**(-pKh)
   Ka1 = 10.**(-pKa1)
   Ka2 = 10.**(-pKa2)
   Kw  = 10.**(-pKw)

   ! Calculate the speciation to know the molar mass of DIC                                                             !
   H    = 10.**(-pH)
   alf_CO2  = (H*H)/(H*H + Ka1*H + Ka1*Ka2)  !alfa0, this is (H2CO3+CO2)
   alf_HCO3 = (Ka1*H)/(H*H + Ka1*H + Ka1*Ka2) !alfa1
   alf_CO3  = (Ka1*Ka2)/(H*H + Ka1*H + Ka1*Ka2) !alfa2

   ! and update speciation (mol C/L),fractination
   CO2  = dic*alf_CO2
   HCO3 = dic*alf_HCO3
   CO3  = dic*alf_CO3

   ! calculate TA
   TA = dic * (Ka1*H + 2.0*Ka1*Ka2) / (H*H + Ka1*H + Ka1*Ka2) + (Kw/H) - H

   aed2_carbon_co2 = CO2

   !print *, 'End aed2_carbon_co2,pH=', pH, '% diff |meas-calc| alk=', (alkalinity-TA)/alkalinity*100

END FUNCTION aed2_carbon_co2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed2_isotope_c
