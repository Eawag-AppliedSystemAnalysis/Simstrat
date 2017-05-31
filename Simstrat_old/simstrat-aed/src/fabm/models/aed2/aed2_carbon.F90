!###############################################################################
!#                                                                             #
!# aed2_carbon.F90                                                             #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_carbon
!-------------------------------------------------------------------------------
! aed2_carbon --- carbon biogeochemical model
!
! The AED module carbon contains equations that describe exchange of
! soluable reactive carbon across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_util,  ONLY: aed2_gas_piston_velocity

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_carbon_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_carbon_data_t
      !# Variable identifiers
      INTEGER  :: id_dic, id_pH, id_ch4, id_oxy
      INTEGER  :: id_Fsed_dic, id_Fsed_ch4
      INTEGER  :: id_temp, id_salt
      INTEGER  :: id_wind
      INTEGER  :: id_ch4ox
      INTEGER  :: id_sed_dic
      INTEGER  :: id_atm_co2_exch, id_atm_ch4_exch

      !# Model parameters
      AED_REAL :: Fsed_dic, Ksed_dic, theta_sed_dic
      AED_REAL :: Fsed_ch4, Ksed_ch4, theta_sed_ch4
      AED_REAL :: Rch4ox, Kch4ox, vTch4ox, atmco2,ionic
      LOGICAL  :: use_oxy, use_sed_model_dic, use_sed_model_ch4
      LOGICAL  :: simDIC, simCH4

     CONTAINS
         PROCEDURE :: define            => aed2_define_carbon
         PROCEDURE :: calculate_surface => aed2_calculate_surface_carbon
         PROCEDURE :: calculate         => aed2_calculate_carbon
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_carbon
         PROCEDURE :: equilibrate       => aed2_equilibrate_carbon
!        PROCEDURE :: mobility          => aed2_mobility_carbon
!        PROCEDURE :: light_extinction  => aed2_light_extinction_carbon
!        PROCEDURE :: delete            => aed2_delete_carbon

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_carbon(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS

   INTEGER  :: status

   AED_REAL          :: pH_initial=7.5
   AED_REAL          :: ionic = 0.0
   AED_REAL          :: dic_initial=4.5
   AED_REAL          :: Fsed_dic = 3.5
   AED_REAL          :: Ksed_dic = 30.0
   AED_REAL          :: theta_sed_dic = 1.0
   CHARACTER(len=64) :: Fsed_dic_variable=''
   AED_REAL          :: ch4_initial=4.5
   AED_REAL          :: Fsed_ch4 = 3.5
   AED_REAL          :: Ksed_ch4 = 30.0
   AED_REAL          :: theta_sed_ch4 = 1.0
   CHARACTER(len=64) :: Fsed_ch4_variable=''
   AED_REAL          :: Rch4ox = 0.01
   AED_REAL          :: Kch4ox = 0.01
   AED_REAL          :: vTch4ox= 1.05
!  AED_REAL          :: atmco2 = 367e-6
   AED_REAL          :: atmco2 = 367.
   CHARACTER(len=64) :: methane_reactant_variable=''


   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed2_carbon/ dic_initial,pH_initial,ionic,Fsed_dic,Ksed_dic,theta_sed_dic,Fsed_dic_variable, &
                         ch4_initial,Fsed_ch4,Ksed_ch4,theta_sed_ch4,Fsed_ch4_variable, &
                         atmco2,Rch4ox,Kch4ox,vTch4ox,methane_reactant_variable

!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_carbon,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist aed2_carbon'
      STOP
   ENDIF

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%Fsed_dic      = Fsed_dic/secs_pr_day
   data%Ksed_dic      = Ksed_dic
   data%theta_sed_dic = theta_sed_dic
   data%ionic         = ionic
   data%Fsed_ch4      = Fsed_ch4/secs_pr_day
   data%Ksed_ch4      = Ksed_ch4
   data%theta_sed_ch4 = theta_sed_ch4
   data%Rch4ox        = Rch4ox/secs_pr_day
   data%Kch4ox        = Kch4ox
   data%vTch4ox       = vTch4ox
   data%atmco2        = atmco2
   data%simDIC        = .false.
   data%simCH4        = .false.



   ! Register state variables
   IF (dic_initial>MISVAL) THEN
      data%id_dic = aed2_define_variable('dic','mmol/m**3','dissolved inorganic carbon',     &
                                       dic_initial,minimum=zero_)
      data%simDIC = .true.
      data%id_pH = aed2_define_variable('pH','-','pH',     &
                                       pH_initial,minimum=zero_)
   ENDIF

   IF (ch4_initial>MISVAL) THEN
      data%id_ch4 = aed2_define_variable('ch4','mmol/m**3','methane',    &
                                     ch4_initial,minimum=zero_)
      data%simCH4 = .true.
   ENDIF

   !# Register external state variable dependencies
   data%use_oxy = methane_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) THEN
      data%id_oxy = aed2_locate_variable(methane_reactant_variable)
   ENDIF

   data%use_sed_model_dic = Fsed_dic_variable .NE. ''
   IF (data%use_sed_model_dic) &
      data%id_Fsed_dic = aed2_locate_global_sheet(Fsed_dic_variable)

   data%use_sed_model_ch4 = Fsed_ch4_variable .NE. ''
   IF (data%use_sed_model_ch4) &
      data%id_Fsed_ch4 = aed2_locate_global_sheet(Fsed_ch4_variable)

   !# Register diagnostic variables
   data%id_ch4ox = aed2_define_diag_variable('ch4ox','/d', 'methane oxidation rate')
   data%id_sed_dic = aed2_define_sheet_diag_variable('sed_dic','mmol/m**2/d',        &
                                                      'Filterable reactive carbon')

   data%id_atm_co2_exch = aed2_define_sheet_diag_variable('atm_co2_exch',            &
                             'mmol/m**2/d', 'CO2 exchange across atm/water interface')
   data%id_atm_ch4_exch = aed2_define_sheet_diag_variable('atm_ch4_exch',            &
                             'mmol/m**2/d', 'CH4 exchange across atm/water interface')

   !# Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')
   data%id_wind = aed2_locate_global_sheet('wind_speed')
END SUBROUTINE aed2_define_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: dic,ch4,oxy,temp
   AED_REAL :: ch4oxidation
!
!-------------------------------------------------------------------------------
!BEGIN

   IF(data%simDIC .AND. data%simCH4) THEN
      ! Retrieve current (local) state variable values.
      dic = _STATE_VAR_(data%id_dic)! carbon
      ch4 = _STATE_VAR_(data%id_ch4)! carbon

      ! Retrieve current dependent state variable values.
      IF (data%use_oxy) THEN ! & use_oxy
         oxy = _STATE_VAR_(data%id_oxy)! oxygen
      ELSE
         oxy = 0.0
      ENDIF

      ! Retrieve current environmental conditions.
      temp = _STATE_VAR_(data%id_temp) ! temperature

      ! Define some intermediate quantities units mmol C/m3/day
      ch4oxidation = aed2_carbon_fch4ox(data%use_oxy,data%Rch4ox,data%Kch4ox,data%vTch4ox,oxy,temp)

      ! Set temporal derivatives
      _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (ch4*ch4oxidation)
      _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (-ch4*ch4oxidation)

      ! If an externally maintained oxygen pool is present, take nitrification from it
      IF (data%use_oxy) then ! & use_oxy
         _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-(32./12.)*ch4*ch4oxidation)
      ENDIF

      ! Export diagnostic variables
      _DIAG_VAR_(data%id_ch4ox) =  ch4oxidation
   ENDIF

END SUBROUTINE aed2_calculate_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind

   ! State
   AED_REAL :: dic,ph,ch4

   ! Temporary variables
   AED_REAL :: pCO2,FCO2,FCH4
   AED_REAL :: Ko,kCH4,KCO2, CH4solub
   AED_REAL :: Tabs,windHt,atm
   AED_REAL :: A1,A2,A3,A4,B1,B2,B3,logC

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simDIC .AND. .NOT.data%simCH4) RETURN


   !Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind = _STATE_VAR_S_(data%id_wind) ! Wind speed at 10 m above surface (m/s)
   windHt = 10.


   Tabs = temp + 273.15

   ! CO2 flux
   IF(data%simDIC) THEN

      ! Retrieve current (local) state variable values.
     dic = _STATE_VAR_(data%id_dic)! Concentration of carbon in surface layer
     ph = _STATE_VAR_(data%id_pH)  ! Concentration of carbon in surface layer

     kCO2 = aed2_gas_piston_velocity(windHt,wind,temp,salt,schmidt_model=2)

     ! Solubility, Ko (mol/L/atm)
     Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
          + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
     Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
     Ko = exp(Ko)

     ! pCO2 in surface water layer
     pCO2 = aed2_carbon_co2(data%ionic,temp,dic,ph) / Ko

     ! FCO2 = kCO2 * Ko * (pCO2 - PCO2a)
     ! pCO2a = 367e-6 (Keeling & Wharf, 1999)

     !------ Yanti correction (20/5/2013) ----------------------------------------
     ! pCO2 is actually in uatm (=ppm)
     ! mmol/m2/s = m/s * mmol/L/atm * atm
     FCO2 = kCO2 * Ko * (pCO2 - data%atmco2)

     ! FCO2 = - kCO2 * Ko*1e6 * ((pCO2 * 1e-6) - data%atmco2) ! dCO2/dt
     !----------------------------------------------------------------------------

     ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
     _FLUX_VAR_T_(data%id_dic) = -FCO2

     ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).
     _DIAG_VAR_S_(data%id_atm_co2_exch) = FCO2
   END IF

   ! CH4 flux
   IF(data%simCH4) THEN
     ! Algorithm from Arianto Santoso <abs11@students.waikato.ac.nz>

     ! Concentration of methane in surface layer
     ch4 = _STATE_VAR_(data%id_ch4)

     ! Piston velocity for CH4
     kCH4 = aed2_gas_piston_velocity(windHt,wind,temp,salt,schmidt_model=4)

     ! Solubility, Ko (mol/L/atm)

     atm = 1.76 * 1e-6 !## current atmospheric CH4 data from NOAA (in ppm)

     A1 = -415.2807
     A2 = 596.8104
     A3 = 379.2599
     A4 = -62.0757
     B1 = -0.05916
     B2 = 0.032174
     B3 = -0.0048198

     logC = (log(atm)) + A1 + (A2 * (100./Tabs)) + (A3 * log (Tabs/100.)) + (A4 * (Tabs/100.)) + &
                 salt * (B1 + (B2  * (Tabs/100.)) + (B3 * (Tabs/100.)*(Tabs/100.)))

     CH4solub = exp(logC) * 1e-3


     ! mmol/m2/s = m/s * mmol/m3
     FCH4 = kCH4 *  (ch4 - CH4solub)

     !----------------------------------------------------------------------------

     ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
     _FLUX_VAR_T_(data%id_ch4) = -FCH4

     ! Also store ch4 flux across the atm/water interface as diagnostic variable (mmmol/m2).
     _DIAG_VAR_S_(data%id_atm_ch4_exch) = FCH4
   END IF
END SUBROUTINE aed2_calculate_surface_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED carbon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: dic,oxy

   ! Temporary variables
   AED_REAL :: dic_flux, ch4_flux, Fsed_dic, Fsed_ch4

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simDIC) RETURN


   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   dic = _STATE_VAR_(data%id_dic)! carbon

   IF ( data%use_sed_model_dic ) THEN
      Fsed_dic = _STATE_VAR_S_(data%id_Fsed_dic)
   ELSE
       Fsed_dic = data%Fsed_dic
   ENDIF
   IF ( data%use_sed_model_ch4 ) THEN
      Fsed_ch4 = _STATE_VAR_S_(data%id_Fsed_ch4)
   ELSE
       Fsed_ch4 = data%Fsed_ch4
   ENDIF

   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      dic_flux = Fsed_dic * oxy/(data%Ksed_dic+oxy) * (data%theta_sed_dic**(temp-20.0))
      ch4_flux = Fsed_ch4 * data%Ksed_ch4/(data%Ksed_ch4+oxy) * (data%theta_sed_ch4**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      dic_flux = Fsed_dic * (data%theta_sed_dic**(temp-20.0))
      ch4_flux = Fsed_ch4 * (data%theta_sed_ch4**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   !_SET_BOTTOM_FLUX_(data%id_dic,dic_flux/secs_pr_day)
   !_SET_SED_FLUX_(data%id_dic,dic_flux)
   _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (dic_flux)
   _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (ch4_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_dic) = _FLUX_VAR_B_(data%id_ben_dic) + (-dic_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_dic) = dic_flux

END SUBROUTINE aed2_calculate_benthic_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_equilibrate_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update pH after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! State
   AED_REAL :: dic, pH

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT.data%simDIC) RETURN


    ! Retrieve current (local) state variable values.
!  dic = _STATE_VAR_(data%id_dic)! Concentration of carbon in surface layer
!   pH = _STATE_VAR_(data%id_pH)! Concentration of carbon in surface layer

!print*,"new pH = ",pH
!  _STATE_VAR_(data%id_pH) =  pH

END SUBROUTINE aed2_equilibrate_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_carbon_fch4ox(use_oxy,Rch4ox,Kch4ox,vTch4ox,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for methane oxidation
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rch4ox,Kch4ox,vTch4ox,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      aed2_carbon_fch4ox = Rch4ox * oxy/(Kch4ox+oxy) * (vTch4ox**(temp-20.0))
   ELSE
      aed2_carbon_fch4ox = Rch4ox * (vTch4ox**(temp-20.0))
   ENDIF

END FUNCTION aed2_carbon_fch4ox
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_carbon_co2(ionic, temp, dic, pH)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL, INTENT(IN) :: ionic, dic, temp, pH
!
!LOCALS
   ! Temporary variables
   AED_REAL :: K_h, Kw, Ka1, Ka2, i_f
   AED_REAL :: H, CO2, HCO3, CO3, TA
!-------------------------------------------------------------------------------
!BEGIN

   ! Acidity constants temperature dependence

   ! pKh =  -0.000075324675x2 + 0.016279653680x + 1.110424242424
   ! pKa1 = 0.000142121212x2 - 0.012648181818x + 6.577539393939
   ! pKa2 =  0.000113679654x2 - 0.014687186147x + 10.625769696970
   ! pKw =   0.000201991342x2 - 0.043419653680x + 14.949709090909

   K_h = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
   Ka1 = 0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
   Ka2 = 0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
   Kw = 0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909


   ! Ionic strength dependence

   ! 1st calculate function f
   i_f = (((SQRT(ionic)) / (1+SQRT(ionic))) -0.20*ionic) * &
                       (298.0/(temp+273.))**0.666667

   ! pKh = pKh(0) + bI
   ! b = 0.105 (Butler, 1982)
   K_h = K_h + 0.105*ionic

   ! pKw = pKw(0) - f
   Kw = Kw - i_f

   ! pKa1 = pKa1(0) - f - bI
   Ka1 = Ka1 - i_f - 0.105*ionic

   !pKa2 = pKa2(0) - 2f
   Ka2 = Ka2 + 2.0*i_f

   ! Convert from pK etc to Kh, Kw, Ka1, Ka2
   K_h  = 10.**(-K_h)
   Ka1 = 10.**(-Ka1)
   Ka2 = 10.**(-Ka2)
   Kw  = 10.**(-Kw)


   ! Calculate the speciation to know the molar mass of DIC                                                             !
   H    = 10.**(-pH)
   CO3  = (Ka1*Ka2)/(H*H + Ka1*H + Ka1*Ka2)
   HCO3 = (Ka1*H)/(H*H + Ka1*H + Ka1*Ka2)
   CO2  = (H*H)/(H*H + Ka1*H + Ka1*Ka2)


   ! and update speciation (mol C/L)
   CO3  = dic*CO3
   HCO3 = dic*HCO3
   CO2  = dic*CO2

   ! calculate TA for the previous timestep
   TA = dic * (Ka1*H + 2.0*Ka1*Ka2) / (H*H + Ka1*H + Ka1*Ka2)
   TA = TA + (Kw/H) - H

   aed2_carbon_co2 = CO2

END FUNCTION aed2_carbon_co2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed2_carbon
