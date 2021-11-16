!###############################################################################
!#                                                                             #
!# aed_carbon.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2021 -  The University of Western Australia               #
!#                                                                             #
!#   GLM is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   GLM is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2012                                                          #
!#  Track changes on GitHub @ https://github.com/AquaticEcoDynamics/libaed-water
!#                                                                             #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ______   | || |      __      | || |  _______     | |         !
!         | |   .' ___  |  | || |     /  \     | || | |_   __ \    | |         !
!         | |  / .'   \_|  | || |    / /\ \    | || |   | |__) |   | |         !
!         | |  | |         | || |   / ____ \   | || |   |  __ /    | |         !
!         | |  \ `.___.'\  | || | _/ /    \ \_ | || |  _| |  \ \_  | |         !
!         | |   `._____.'  | || ||____|  |____|| || | |____| |___| | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_carbon
!-------------------------------------------------------------------------------
! aed_carbon --- (inorganic) carbon biogeochemical model
!
! The AED module carbon contains equations that describe exchange of
! dissolved inorganic carbon across the air/water interface and sediment flux,
! and simulation of methane.
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,  ONLY: aed_gas_piston_velocity

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_carbon_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_carbon_data_t
      !# Variable identifiers
      INTEGER  :: id_dic, id_pH, id_ch4, id_oxy, id_talk, id_ch4_bub
      INTEGER  :: id_Fsed_dic, id_Fsed_ch4, id_Fsed_ch4_ebb
      INTEGER  :: id_temp, id_salt
      INTEGER  :: id_wind, id_vel, id_depth
      INTEGER  :: id_par, id_extc, id_dz, id_tau
      INTEGER  :: id_ch4ox, id_pco2
      INTEGER  :: id_sed_dic, id_sed_ch4, id_sed_ch4_ebb, id_sed_ch4_ebb_3d
      INTEGER  :: id_atm_co2, id_atm_ch4, id_atm_ch4_ebb, id_ch4_ebb_df

      !# Model parameters
      AED_REAL :: Fsed_dic, Ksed_dic, theta_sed_dic
      AED_REAL :: Fsed_ch4, Ksed_ch4, theta_sed_ch4, Fsed_ch4_ebb
      AED_REAL :: Rch4ox, Kch4ox, vTch4ox
      AED_REAL :: atm_co2, atm_ch4, ionic
      AED_REAL :: ch4_bub_aLL, ch4_bub_cLL, ch4_bub_kLL
      AED_REAL :: ch4_bub_disf1, ch4_bub_disf2, ch4_bub_disdp
      AED_REAL :: ch4_bub_tau0, ch4_bub_ws
      AED_REAL :: maxMPBProdn, IkMPB


      !# Model options
      LOGICAL  :: use_oxy, use_sed_model_dic, use_sed_model_ch4, use_sed_model_ebb
      LOGICAL  :: simDIC, simCH4, simCH4ebb
      INTEGER  :: alk_mode, co2_model, ebb_model
      INTEGER  :: co2_piston_model, ch4_piston_model

     CONTAINS
         PROCEDURE :: define            => aed_define_carbon
         PROCEDURE :: calculate_surface => aed_calculate_surface_carbon
         PROCEDURE :: calculate         => aed_calculate_carbon
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_carbon
         PROCEDURE :: equilibrate       => aed_equilibrate_carbon
   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_carbon(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_carbon_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER           :: status

!  %% NAMELIST   %%  /aed_carbon/
!  %% Last Checked 05/09/2021
   AED_REAL          :: dic_initial      = 1000.0
   AED_REAL          :: pH_initial       = 7.5
   AED_REAL          :: ch4_initial      = 4.5
   AED_REAL          :: ionic            = 0.0
   AED_REAL          :: Fsed_dic         = 0.0
   AED_REAL          :: Ksed_dic         = 30.0
   AED_REAL          :: theta_sed_dic    = 1.0
   CHARACTER(len=64) :: Fsed_dic_variable=''
   AED_REAL          :: Fsed_ch4         = 0.0
   AED_REAL          :: Ksed_ch4         = 30.0
   AED_REAL          :: theta_sed_ch4    = 1.0
   CHARACTER(len=64) :: Fsed_ch4_variable=''
   AED_REAL          :: atm_co2          = 367e-6
   AED_REAL          :: atm_ch4          = 1.76e-6
   AED_REAL          :: Rch4ox           = 0.01
   AED_REAL          :: Kch4ox           = 0.01
   AED_REAL          :: vTch4ox          = 1.05
   CHARACTER(len=64) :: methane_reactant_variable=''

   INTEGER           :: co2_model        = 1
   INTEGER           :: alk_mode         = 1
   INTEGER           :: ebb_model        = 0
   INTEGER           :: co2_piston_model = 1
   INTEGER           :: ch4_piston_model = 1

   LOGICAL           :: simCH4ebb
   AED_REAL          :: Fsed_ch4_ebb     = zero_
   CHARACTER(len=64) :: Fsed_ebb_variable=''
   AED_REAL          :: ch4_bub_aLL      = 42.9512677
   AED_REAL          :: ch4_bub_cLL      = 0.634
   AED_REAL          :: ch4_bub_kLL      = -0.8247
   AED_REAL          :: ch4_bub_disf1    = 0.07
   AED_REAL          :: ch4_bub_disf2    = 0.33
   AED_REAL          :: ch4_bub_disdp    = 20.0

   AED_REAL          :: ch4_bub_ws       = zero_
   AED_REAL          :: ch4_bub_tau0     = one_
   AED_REAL          :: maxMPBProdn      =  40.0   ! mmolC/m2/day
   AED_REAL          :: IkMPB            = 180.0   ! Light sensitivity of MPB

! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs

!  %% END NAMELIST   %%  /aed_carbon/

   NAMELIST /aed_carbon/ dic_initial,pH_initial,ch4_initial,                &
                         ionic,Rch4ox,Kch4ox,vTch4ox,                       &
                         methane_reactant_variable,                         &
                         Fsed_dic,Ksed_dic,theta_sed_dic,Fsed_dic_variable, &
                         Fsed_ch4,Ksed_ch4,theta_sed_ch4,Fsed_ch4_variable, &
                         atm_co2,atm_ch4,                                   &
                         co2_piston_model, ch4_piston_model,                &
                         co2_model, alk_mode, ebb_model,                    &
                         Fsed_ch4_ebb, Fsed_ebb_variable,                   &
                         ch4_bub_aLL,ch4_bub_cLL, ch4_bub_kLL,              &
                         ch4_bub_disf1, ch4_bub_disf2, ch4_bub_disdp,       &
                         ch4_bub_ws, ch4_bub_tau0,                          &
                         maxMPBProdn, IkMPB,                                &
                         diag_level


!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_carbon initialization"

   !# Set defaults
   data%simDIC        = .false.
   data%simCH4        = .false.
   data%simCH4ebb     = .false.
   IF (ebb_model>0) data%simCH4ebb = .true.

   !# Read the namelist
   read(namlst,nml=aed_carbon,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist aed_carbon'
      STOP
   ENDIF

   !# Store parameter values in modules own derived type
   !  Note: rates are provided in values per day, and
   !        are converted here to values per second.
   data%Fsed_dic         = Fsed_dic/secs_per_day
   data%Ksed_dic         = Ksed_dic
   data%theta_sed_dic    = theta_sed_dic
   data%ionic            = ionic
   data%co2_model        = co2_model
   data%alk_mode         = alk_mode
   data%atm_co2          = atm_co2
   data%co2_piston_model = co2_piston_model

   data%Fsed_ch4         = Fsed_ch4/secs_per_day
   data%Ksed_ch4         = Ksed_ch4
   data%theta_sed_ch4    = theta_sed_ch4
   data%Rch4ox           = Rch4ox/secs_per_day
   data%Kch4ox           = Kch4ox
   data%vTch4ox          = vTch4ox
   data%atm_ch4          = atm_ch4
   data%ch4_piston_model = ch4_piston_model
   data%ebb_model        = ebb_model
   data%Fsed_ch4_ebb     = Fsed_ch4_ebb
   data%ch4_bub_aLL      = ch4_bub_aLL
   data%ch4_bub_cLL      = ch4_bub_cLL
   data%ch4_bub_kLL      = ch4_bub_kLL
   data%ch4_bub_disf1    = ch4_bub_disf1
   data%ch4_bub_disf2    = ch4_bub_disf2
   data%ch4_bub_disdp    = ch4_bub_disdp

   data%ch4_bub_ws       = ch4_bub_ws    !unused
   data%ch4_bub_tau0     = ch4_bub_tau0  !unused
   data%maxMPBProdn      = maxMPBProdn   !unused
   data%IkMPB            = IkMPB         !unused


   !# Register state variables
   IF (dic_initial>MISVAL) THEN
      data%simDIC = .true.
      data%id_dic = aed_define_variable('dic','mmol/m**3','dissolved inorganic carbon', &
                                       dic_initial,minimum=zero_)
      data%id_pH = aed_define_variable('pH','-','pH',     &
                                       pH_initial,minimum=zero_)
   ENDIF

   IF (ch4_initial>MISVAL) THEN
      data%simCH4 = .true.
      data%id_ch4 = aed_define_variable('ch4','mmol/m**3','methane',    &
                                     ch4_initial,minimum=zero_)
      IF( data%simCH4ebb ) THEN
        data%id_ch4_bub = aed_define_variable('ch4_bub','mmol/m**3', &
                                     'methane bubbles',zero_,minimum=zero_)
      ENDIF
   ENDIF

   !# Register external state variable dependencies
   data%use_oxy = methane_reactant_variable .NE. '' !i.e., oxygen module engaged
   IF (data%use_oxy) THEN
      data%id_oxy = aed_locate_variable(methane_reactant_variable)
   ENDIF

   data%use_sed_model_dic = Fsed_dic_variable .NE. ''
   IF (data%use_sed_model_dic) &
      data%id_Fsed_dic = aed_locate_sheet_variable(Fsed_dic_variable)

   data%use_sed_model_ch4 = Fsed_ch4_variable .NE. ''
   IF (data%use_sed_model_ch4) &
      data%id_Fsed_ch4 = aed_locate_sheet_variable(Fsed_ch4_variable)

   data%use_sed_model_ebb = Fsed_ebb_variable .NE. ''
   IF (data%use_sed_model_ebb) &
      data%id_Fsed_ch4_ebb = aed_locate_sheet_variable(Fsed_ebb_variable)

   !# Register diagnostic variables
   data%id_pco2 = aed_define_diag_variable('pCO2','atm', 'pCO2')
   IF (diag_level>0) THEN
     data%id_sed_dic = aed_define_sheet_diag_variable('sed_dic','mmol/m**2/d', &
                            'CO2 exchange across sed/water interface')
     data%id_atm_co2 = aed_define_sheet_diag_variable('atm_co2_flux',          &
                            'mmol/m**2/d', 'CO2 exchange across atm/water interface')

     IF( data%simCH4 ) THEN
       data%id_ch4ox   = aed_define_diag_variable('ch4ox','mmol/m**3/d', 'methane oxidation rate')
       data%id_sed_ch4 = aed_define_sheet_diag_variable('sed_ch4','mmol/m**2/d', &
                            'CH4 exchange across sed/water interface')
       data%id_atm_ch4 = aed_define_sheet_diag_variable('atm_ch4_flux',        &
                            'mmol/m**2/d', 'CH4 exchange across atm/water interface')
       IF( data%simCH4ebb ) THEN
         data%id_sed_ch4_ebb_3d = aed_define_diag_variable('sed_ch4_ebb_3d','mmol/m**3/d', &
                            'CH4 ebullition release rate')
         data%id_ch4_ebb_df = aed_define_diag_variable('ch4_ebb_df','mmol/m**3/d', &
                            'CH4 bubble dissolution rate')
         data%id_sed_ch4_ebb = aed_define_sheet_diag_variable('sed_ch4_ebb','mmol/m**2/d', &
                            'CH4 ebullition across sed/water interface')
         data%id_atm_ch4_ebb = aed_define_sheet_diag_variable('atm_ch4_ebb_flux', &
                            'mmol/m**2/d', 'CH4 ebullition across atm/water interface')
       ENDIF
     ENDIF
   ENDIF

   !# Register environmental dependencies
   data%id_temp = aed_locate_global('temperature')
   data%id_salt = aed_locate_global('salinity')
   data%id_extc = aed_locate_global('extc_coef')
   data%id_par  = aed_locate_global('par')
   data%id_dz   = aed_locate_global('layer_ht')
   data%id_vel  = aed_locate_global('cell_vel')           ! needed for k600
   data%id_depth= aed_locate_global('depth')
!  data%id_depth= aed_locate_global('layer_ht')
   data%id_wind = aed_locate_sheet_global('wind_speed')
   IF( data%simCH4ebb ) data%id_tau  = aed_locate_sheet_global('taub')

END SUBROUTINE aed_define_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_carbon_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
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
      dic = _STATE_VAR_(data%id_dic)    ! DIC
      ch4 = _STATE_VAR_(data%id_ch4)    ! CH4

      !# Retrieve current dependent state variable values.
      IF (data%use_oxy) THEN
         oxy = _STATE_VAR_(data%id_oxy) ! O2
      ELSE
         oxy = zero_
      ENDIF

      !# Retrieve current environmental conditions.
      temp = _STATE_VAR_(data%id_temp)  ! temperature

      !# Compute rates of change (mmol C/m3/day)
      ch4oxidation = aed_carbon_fch4ox(data%use_oxy,data%Rch4ox,data%Kch4ox,data%vTch4ox,oxy,temp)

      !# Set temporal derivatives
      _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (ch4*ch4oxidation)
      _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (-ch4*ch4oxidation)

      !# If a linked oxygen pool is present, take oxidation from it assume 1:1 stoichometry
      IF (data%use_oxy) THEN
         _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - ch4*ch4oxidation  ! *(32./12.) not mass
      ENDIF

      !# Export diagnostic variables
      IF (diag_level>0) _DIAG_VAR_(data%id_ch4ox) =  ch4*ch4oxidation*secs_per_day
   ENDIF

END SUBROUTINE aed_calculate_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_carbon_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, S, T, depth, vel

   ! State
   AED_REAL :: dic,pHin,ch4,talk = 0.,TCO2 = 0.

   ! Temporary variables

   AED_REAL :: pCO2 = 0.,FCO2,FCH4,henry
   AED_REAL :: Ko, kCH4, KCO2, CH4solub
   AED_REAL :: Tabs,windHt,atm
   AED_REAL :: A1,A2,A3,A4,B1,B2,B3,logC
   AED_REAL :: a,b,c,dcf
   AED_REAL :: ca, bc, cb, carba, bicarb, carb, om_cal, om_arg
   AED_REAL :: p00,p10,p01,p20,p11,p02
   AED_REAL :: phdum,tadum,co2dum,pHout,pH,deltapH
   INTEGER  :: iter

   INTEGER  :: ii

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simDIC .AND. .NOT.data%simCH4) RETURN

   Ko = 0.

   !----------------------------------------------------------------------------
   !# Get dependent state variables from physical driver
   windHt = 10.
   wind   = _STATE_VAR_S_(data%id_wind) ! Wind speed at 10 m above surface (m/s)
   temp   = _STATE_VAR_(data%id_temp)   ! Temperature (degrees Celsius)
   salt   = _STATE_VAR_(data%id_salt)   ! Salinity (psu)
   depth  = MAX( _STATE_VAR_(data%id_depth), one_ )
   IF (data%id_vel > 0 ) THEN
     vel = _STATE_VAR_(data%id_vel)
   ELSE
    ! vel  = SQRT(_STATE_VAR_(data%id_E_tau)/_STATE_VAR_(data%id_E_dens))
    ! vel = vel/0.41 * log(depth/0.01)
     vel = 0.0001
   ENDIF

   Tabs = temp + 273.15

   !----------------------------------------------------------------------------
   !# CO2 concentration in the surface layer, depends on DIC, TA
   IF(data%simDIC) THEN

      !# Retrieve current (local) state variable values.
     dic  = _STATE_VAR_(data%id_dic)! Concentration of carbon in surface layer
     pHin = _STATE_VAR_(data%id_pH)  ! Concentration of carbon in surface layer

     IF ( data%co2_model == 1 ) THEN
       !# Use the CO2SYS code for computing pCO2 & pH
       S=salt; T=temp
       a    =  8.24493d-1 - 4.0899d-3*T + 7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
       b    = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2
       c    =  4.8314d-4
       dcf  = (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
                - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3
       TCO2 = dic / (1.0D6*dcf) ! change unit to mol/kgSW

       IF( data%alk_mode == 0 ) THEN
         ! Freshwater system - carbonate alkalinity assumed
         pH = pHin
         deltapH = 5.
         ! As carbonate alkalinity depends on H+, lets iterate
         iter = 0
         DO WHILE (abs(deltapH) > 1e-10)
           ! Use CO2SYS to estimate TA, from past pH
           CALL CO2SYS(1,T,S,zero_,tadum,TCO2,pH,pHdum,co2dum,talk)
           ! Use CO2SYS to re-estimate pH, from TA
           CALL CO2SYS(0,T,S,zero_,talk,TCO2,phdum,pHout,pCO2,tadum)
           deltapH = pH-pHout
           pH = pHout
          ! print *,'ss',talk
           iter = iter+1
           IF(iter>100) then
             print *, 'note pH-TA convergance failure',pH
             exit
           ENDIF
         END DO

       ELSEIF( data%alk_mode == 1 ) THEN
         ! talk = 520.1 + 51.24*S  ! Atlantic (Millero 1998) from fabm, not suitable for estuaries
         ! talk = 1136.1 + 1.2*S*S + 2.8*S !Chesapeake Bay (George et al., 2013)
         talk =  1627.4 + 22.176*S   !regression from Naomi's data on Caboolture
         talk  = talk / 1.0D6      ! change unit to mol/kgSW

       ELSEIF( data%alk_mode == 2 ) THEN
         p00  =       1063
         p10  =      1.751
         p01  =   -0.05369
         p20  =     0.2266
         p11  =  -0.001252
         p02  =  0.0002546
         talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
         talk = talk / 1.0D6      ! change unit to mol/kgSW

       ELSEIF( data%alk_mode == 3 ) THEN
         p00 =      -258.8
         p10 =       34.59
         p01 =      0.9923
         p20 =      0.8186
         p11 =    -0.03101
         p02 =   0.0001045
         talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
         talk = talk / 1.0D6      ! change unit to mol/kgSW

        ELSEIF( data%alk_mode == 4 ) THEN
          p00 =      -47.51
          p10 =      -17.21
          p01 =        1.32
          p20 =      0.1439
          p11 =     0.01224
          p02 =  -0.0002055
          talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
          talk = talk / 1.0D6      ! change unit to mol/kgSW

        ELSEIF( data%alk_mode == 5 ) THEN
          p00 =       157.7
          p10 =       4.298
          p01 =      0.6448
          p20 =      0.2107
          p11 =   -0.002072
          p02 =   0.0001239
          talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
          talk = talk / 1.0D6      ! change unit to mol/kgSW

       ENDIF

       !CALL CO2DYN ( TCO2, talk, T, S, pCO2, pH, HENRY, ca, bc, cb)
       !CALL CO2SYS(T,S,talk,TCO2,pCO2,pH)

       ! Return pH and pCO2 based on TA (alkalinity)
       CALL CO2SYS(0,T,S,zero_,talk,TCO2,phdum,pHout,pCO2,tadum)

       _DIAG_VAR_(data%id_pco2) = pCO2
       ! _STATE_VAR_(data%id_talk) = talk*(1.0D6)           ! total alkalinity (umol/kg)

    !   print *,'pHout',pHout, talk, pCO2

     ELSEIF ( data%co2_model == 2 ) THEN
       !# Use the Butler CO2 code for computing pCO2 & pH

       !# Solubility, Ko (mol/L/atm)
       Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
                                        + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
       Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
       Ko = exp(Ko)

       pCO2 = aed_carbon_co2(data%ionic,temp,dic,ph)*1e-6 / Ko  !(=atm), use Yanti's script for pCO2
       _DIAG_VAR_(data%id_pco2) = pCO2

     ELSEIF ( data%co2_model == 0 ) THEN
       !# Use the aed_geochemistry module for computing pCO2 & pH
       pCO2 = _DIAG_VAR_(data%id_pco2)  ! this diagnostic is getting set in aed_geocehmistry

     ENDIF

     !# Now compute piston velocity, k
     kCO2 = aed_gas_piston_velocity(windHt,wind,temp,salt,                     &
         vel=vel,depth=depth,schmidt_model=2,piston_model=data%co2_piston_model)

     !# Now compute the CO2 flux
     !  FCO2 = kCO2 * Ko * (pCO2 - PCO2a)
     !  pCO2a = 367e-6 atm (Keeling & Wharf, 1999)
     !  mmol/m2/s = m/s * mmol/m3/atm * atm
     !  FCO2 = - kCO2 * Ko*1e6 * ((pCO2 * 1e-6) - data%atm_co2) ! dCO2/dt
     FCO2 = kCO2 * (1e6*Ko) * (pCO2 - data%atm_co2)

     !--------------------------------------------------------------------------
     !# Transfer surface exchange value to AED2 (mmmol/m2/s) converted by driver
     _FLUX_VAR_T_(data%id_dic) = -FCO2

     !# Also store co2 flux across the atm/water interface as a
     !  diagnostic variable (mmmol/m2/d)
     IF (diag_level>0) _DIAG_VAR_S_(data%id_atm_co2) = FCO2*secs_per_day

   END IF


   !----------------------------------------------------------------------------
   !# CH4 flux
   IF(data%simCH4) THEN
     ! Algorithm from Arianto Santoso <abs11@students.waikato.ac.nz>

     ! Concentration of methane in surface layer
     ch4 = _STATE_VAR_(data%id_ch4)

     ! Piston velocity for CH4
     kCH4 = aed_gas_piston_velocity(windHt,wind,temp,salt, &
         vel=vel,depth=depth,schmidt_model=4,piston_model=data%ch4_piston_model)

     ! Solubility, Ko (mol/L/atm)
     atm = data%atm_ch4   ! 1.76 e-6 !## recent atmospheric CH4 from NOAA (in atm)

     A1 = -415.2807
     A2 =  596.8104
     A3 =  379.2599
     A4 =  -62.0757
     B1 =   -0.05916
     B2 =    0.032174
     B3 =   -0.0048198

     logC = (log(atm)) + A1                                                   &
          + (A2 * (100./Tabs)) + (A3 * log (Tabs/100.)) + (A4 * (Tabs/100.))  &
          + salt * (B1 + (B2  * (Tabs/100.)) + (B3 * (Tabs/100.)*(Tabs/100.)))

     CH4solub = exp(logC) * 1e-3

     !# Now compute methane atm flux  (mmol/m2/s = m/s * mmol/m3)
     FCH4 = kCH4 *  (ch4 - CH4solub)

     !----------------------------------------------------------------------------
     !# Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
     _FLUX_VAR_T_(data%id_ch4) = -FCH4

     !# Also store CH4 flux across the atm/water interface as
     !  diagnostic variable (mmmol/m2/d)
     IF (diag_level>0) _DIAG_VAR_S_(data%id_atm_ch4) = FCH4*secs_per_day

   END IF
   !----------------------------------------------------------------------------

END SUBROUTINE aed_calculate_surface_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED carbon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_carbon_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, par, extc, dz, depth

   ! State
   AED_REAL :: dic, oxy, mpb, ph

   ! Temporary variables
   AED_REAL :: dic_flux, ch4_flux, Fsed_dic, Fsed_ch4, ebb_flux, Fsed_ch4_ebb, ch4_bub_disf
   !AED_REAL, PARAMETER :: maxMPBProdn = 40.     ! mmolC/m2/day                     !
   !AED_REAL, PARAMETER :: IkMPB       = 180.0   ! Light sensitivity of MPB  !

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simDIC) RETURN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   par  = _STATE_VAR_(data%id_par)  ! local par
   dz   = _STATE_VAR_(data%id_dz)   ! local layer thickness
   depth= _STATE_VAR_(data%id_depth)! local layer depth
   extc = _STATE_VAR_(data%id_extc) ! local extinction

    ! Retrieve current (local) state variable values.
   dic  = _STATE_VAR_(data%id_dic)! carbon
   pH   = _STATE_VAR_(data%id_pH)! pH

   IF ( data%use_sed_model_dic ) THEN
      Fsed_dic = _STATE_VAR_S_(data%id_Fsed_dic)
   ELSE
      Fsed_dic = data%Fsed_dic
   ENDIF
   IF ( data%use_sed_model_ch4 ) THEN
      Fsed_ch4 = _STATE_VAR_S_(data%id_Fsed_ch4)
      IF( data%simCH4ebb ) Fsed_ch4_ebb = _STATE_VAR_S_(data%id_Fsed_ch4_ebb)
   ELSE
      Fsed_ch4 = data%Fsed_ch4
      IF( data%simCH4ebb ) Fsed_ch4_ebb = data%Fsed_ch4_ebb
   ENDIF

   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      dic_flux = Fsed_dic * oxy/(data%Ksed_dic+oxy) * (data%theta_sed_dic**(temp-20.0))
      ch4_flux = Fsed_ch4 * data%Ksed_ch4/(data%Ksed_ch4+oxy) * (data%theta_sed_ch4**(temp-20.0))
      IF( data%simCH4ebb ) ebb_flux = Fsed_ch4_ebb * (data%theta_sed_ch4**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      dic_flux = Fsed_dic * (data%theta_sed_dic**(temp-20.0))
      ch4_flux = Fsed_ch4 * (data%theta_sed_ch4**(temp-20.0))
      IF( data%simCH4ebb ) THEN
        ebb_flux = Fsed_ch4_ebb * (data%theta_sed_ch4**(temp-20.0))
        ! Kinneret special lake level equations
        ebb_flux = ebb_flux * data%ch4_bub_cLL * exp(data%ch4_bub_kLL*(data%ch4_bub_aLL-depth))
      ENDIF
   ENDIF


  !! Allow photosynthetic production of CO2 in the benthos due to MPB if light and suitable pH
  !par = par * (exp(-extc*dz))
  ! IF( par > 50. .AND. pH > 5.5 .AND. pH < 9.6 ) THEN
  !   mpb = (data%maxMPBProdn/secs_per_day)*(1.0-exp(-par/data%IkMPB)) * (data%theta_sed_dic**(temp-20.0))
  !   dic_flux = Fsed_dic - mpb
  !   IF (data%use_oxy) THEN
  !     _FLUX_VAR_(data%id_oxy) =  _FLUX_VAR_(data%id_oxy) + mpb
  !   ENDIF
  ! ENDIF


   ! Set bottom fluxes for the pelagic (flux per surface area, per second)
   ! Increment sediment flux value into derivative of water column variable
   _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (dic_flux)
   IF( data%simCH4 .and. diag_level>0) _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (ch4_flux)

   ! Store dissolved sediment fluxes as diagnostic variables (flux per surface area, per day)
   _DIAG_VAR_S_(data%id_sed_dic) = dic_flux * secs_per_day
   IF( data%simCH4) _DIAG_VAR_S_(data%id_sed_ch4) = ch4_flux * secs_per_day

   ! Re-distribute bubbles to the water or atmosphere, or dissolve
   IF( data%simCH4ebb ) THEN
     ! Add bubbles to layer
      !_FLUX_VAR_(data%id_ch4_bub) = _FLUX_VAR_(data%id_ch4_bub) + (ebb_flux)

      ! Dissolve bubbles in this bottom layer, depending on depth
      ch4_bub_disf = data%ch4_bub_disf1
      IF( depth > data%ch4_bub_disdp) ch4_bub_disf = data%ch4_bub_disf2
      IF( data%simCH4) _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + ebb_flux*ch4_bub_disf

      IF (diag_level>0) THEN
        _DIAG_VAR_(data%id_ch4_ebb_df) = ebb_flux*ch4_bub_disf * secs_per_day !/ dz. MH Something wrong with dz here?

       ! Release the remainder to the atmosphere (mmol/m2/day)
        IF (diag_level>0) _DIAG_VAR_S_(data%id_atm_ch4_ebb) = ebb_flux * (1-ch4_bub_disf) * secs_per_day

       ! Note the bubble flux, as the zone sees it  (mmol/m2/day)
        IF (diag_level>0) _DIAG_VAR_S_(data%id_sed_ch4_ebb) = ebb_flux * secs_per_day

       ! Note the bubble flux, as the water sees it  (mmol/m3/day)
        IF (diag_level>0) _DIAG_VAR_(data%id_sed_ch4_ebb_3d) = ebb_flux * secs_per_day / dz
      ENDIF
    ENDIF


   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_dic) = _FLUX_VAR_B_(data%id_ben_dic) + (-dic_flux/secs_per_day)


END SUBROUTINE aed_calculate_benthic_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update pH after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_carbon_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! State
   AED_REAL :: dic, pHin, pCO2, temp, salt
   AED_REAL :: S,T,a,b,c,dcf,talk = 0.,TCO2 = 0.,ca,bc,cb,HENRY
   AED_REAL :: p00,p10,p01,p20,p11,p02

   AED_REAL :: Tabs, Ko
   AED_REAL :: K_h, Kw, Ka1, Ka2, i_f
   AED_REAL :: H, CO2, HCO3, CO3, TA
   AED_REAL :: phdum,tadum,co2dum,pHout,pH,deltapH
   INTEGER  :: iter

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT.data%simDIC) RETURN

    pCO2 = zero_

    !# Retrieve current (local) state variable values.
    dic  = _STATE_VAR_(data%id_dic)  ! Concentration of DIC in the cell
    salt = _STATE_VAR_(data%id_salt) ! Salinity
    temp = _STATE_VAR_(data%id_temp) ! Temperature
    pHin = _STATE_VAR_(data%id_pH)   ! pH (from previous time-step)

    pH = pHin

    IF ( data%co2_model == 1 ) THEN

      S=salt; T=temp

      !# Use the CO2SYS code for computing pCO2 & pH
      a    =  8.24493d-1 - 4.0899d-3*T + 7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
      b    = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2
      c    =  4.8314d-4

      dcf  = (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
               - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3
      TCO2  = dic / (1.0D6*dcf) ! change unit to mol/kgSW


      IF( data%alk_mode == 0 ) THEN
        ! Freshwater system - carbonate alkalinity assumed
        pH = pHin
        deltapH = 5.
        ! As carbonate alkalinity depends on H+, lets iterate
        iter = 0
        DO WHILE (abs(deltapH) > 1e-10)
          ! Use CO2SYS to estimate TA, from past pH
          CALL CO2SYS(1,T,S,zero_,tadum,TCO2,pH,pHdum,co2dum,talk)
          ! Use CO2SYS to re-estimate pH, from TA
          CALL CO2SYS(0,T,S,zero_,talk,TCO2,phdum,pHout,pCO2,tadum)
          deltapH = pH-pHout
          pH = pHout
          iter = iter+1
          IF(iter>100) then
            print *, 'note pH-TA convergance failure',pH
            exit
          ENDIF
        END DO

      ELSEIF( data%alk_mode == 1 ) THEN
        ! talk = 520.1 + 51.24*S  ! Atlantic (Millero 1998) from fabm, not suitable for estuaries
        ! talk = 1136.1 + 1.2*S*S + 2.8*S !Chesapeake Bay (George et al., 2013)
        talk =  1627.4 + 22.176*S   !regression from Naomi's data on Caboolture
        talk  = talk / 1.0D6      ! change unit to mol/kgSW

      ELSEIF( data%alk_mode == 2 ) THEN
        p00  =       1063
        p10  =      1.751
        p01  =   -0.05369
        p20  =     0.2266
        p11  =  -0.001252
        p02  =  0.0002546
        talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
        talk = talk / 1.0D6      ! change unit to mol/kgSW

      ELSEIF( data%alk_mode == 3 ) THEN
        p00 =      -258.8
        p10 =       34.59
        p01 =      0.9923
        p20 =      0.8186
        p11 =    -0.03101
        p02 =   0.0001045
        talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
        talk = talk / 1.0D6      ! change unit to mol/kgSW

       ELSEIF( data%alk_mode == 4 ) THEN
         p00 =      -47.51
         p10 =      -17.21
         p01 =        1.32
         p20 =      0.1439
         p11 =     0.01224
         p02 =  -0.0002055
         talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
         talk = talk / 1.0D6      ! change unit to mol/kgSW

       ELSEIF( data%alk_mode == 5 ) THEN
         p00 =       157.7
         p10 =       4.298
         p01 =      0.6448
         p20 =      0.2107
         p11 =   -0.002072
         p02 =   0.0001239
         talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
         talk = talk / 1.0D6      ! change unit to mol/kgSW

      ENDIF

      !CALL CO2DYN ( TCO2, talk, T, S, pCO2, pH, HENRY, ca, bc, cb)
      !CALL CO2SYS(T,S,talk,TCO2,pCO2,pH)

      ! Return pH and pCO2 based on TA (alkalinity)
      CALL CO2SYS(0,T,S,zero_,talk,TCO2,phdum,pHout,pCO2,tadum)
      pH = pHout
  !    print *,'pH',pH,talk,pCO2

    ELSEIF ( data%co2_model == 2 ) THEN
      !# Use the Butler CO2 code for computing pCO2 & pH

      !# Solubility, Ko (mol/L/atm)
      Tabs = temp + 273.15
      Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
                                    + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
      Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
      Ko = exp(Ko)

      ! Acidity constants temperature dependence
      K_h = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
      Ka1 =  0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
      Ka2 =  0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
      Kw  =  0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909

      ! Ionic strength dependence, 1st calculate function f
      i_f = (((SQRT(data%ionic)) / (1+SQRT(data%ionic))) -0.20*data%ionic) * (298.0/(Tabs))**0.666667

      ! pKh = pKh(0) + bI;  b = 0.105 (Butler, 1982)
      K_h = K_h + 0.105*data%ionic

      ! pKw = pKw(0) - f
      Kw = Kw - i_f

      ! pKa1 = pKa1(0) - f - bI
      Ka1 = Ka1 - i_f - 0.105*data%ionic

      !pKa2 = pKa2(0) - 2f
      Ka2 = Ka2 + 2.0*i_f

      ! Convert from pK etc to Kh, Kw, Ka1, Ka2
      K_h = 10.**(-K_h)
      Ka1 = 10.**(-Ka1)
      Ka2 = 10.**(-Ka2)
      Kw  = 10.**(-Kw)

      ! Calculate H, from pH at the last timestep
      H    = 10.**(-pH)

      ! Calculate TA
      TA = dic * (Ka1*H + 2.0*Ka1*Ka2) / (H*H + Ka1*H + Ka1*Ka2)
      TA = TA + (Kw/H) - H

      ! Assuming TA is constant (CO2 has changed due to reactions, but this
      ! doesn't change the alkalinity directly), now calculate new H conc
      H =  calcPH(dic,TA,H,Ka1,Ka2,Kw)
      IF(H > zero_) pH = -log10(H)

      pCO2 = aed_carbon_co2(data%ionic,temp,dic,pH)*1e-6 / Ko  !(=atm)

    ENDIF


    !# Set pCO2 & pH as returned
    _DIAG_VAR_(data%id_pco2) = pCO2
    _STATE_VAR_(data%id_pH)  = pH

END SUBROUTINE aed_equilibrate_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed_carbon_fch4ox(use_oxy,Rch4ox,Kch4ox,vTch4ox,oxy,temp)
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
      aed_carbon_fch4ox = Rch4ox * oxy/(Kch4ox+oxy) * (vTch4ox**(temp-20.0))
   ELSE
      aed_carbon_fch4ox = Rch4ox * (vTch4ox**(temp-20.0))
   ENDIF

END FUNCTION aed_carbon_fch4ox
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed_carbon_co2(ionic, temp, dic, pH)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T, from Butler
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
   K_h = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
   Ka1 =  0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
   Ka2 =  0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
   Kw  =  0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909


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

   aed_carbon_co2 = CO2

END FUNCTION aed_carbon_co2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
! Function to calculate H+ conc given total alkalinity and DIC
!-----------------------------------------------------------------------
 FUNCTION calcPH(CT,TA,H,Ka1,Ka2,Kw) RESULT(rett)
   !-- Incoming
   AED_REAL,  INTENT(INOUT) :: CT      ! DIC
   AED_REAL,  INTENT(INOUT) :: TA      ! Total Alkalinity
   AED_REAL,  INTENT(INOUT) :: H       ! old H+ conc
   AED_REAL,  INTENT(INOUT) :: Ka1     ! 1st Acidity Const
   AED_REAL,  INTENT(INOUT) :: Ka2     ! 2nd Acidity Const
   AED_REAL,  INTENT(INOUT) :: Kw      ! ion product H2O
   !-- Returns the H+ conc at the forward time step
   AED_REAL :: rett
   !-- Local
   AED_REAL :: tmp

   ! 1st iterative sequence using high search bounds and coarse resolution
   tmp = 2.5
   rett = iteratepH(CT,TA,H,Ka1,Ka2,Kw,100,tmp)
   ! 2nd iterative sequence using low search bounds and fine resolution
   tmp = 0.5
   rett = iteratepH(CT,TA,rett,Ka1,Ka2,Kw,200,tmp)

 END FUNCTION calcPH
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
! Solves the Total Alkalinity equation for H+ using an iterative        !
! solution                                                              !
!-----------------------------------------------------------------------!
 FUNCTION iteratePH(CT,TA,H,Ka1,Ka2,Kw,numElements,bounds) RESULT(pH)
   !-- Incoming
   AED_REAL,  INTENT(IN) :: CT                      ! DIC
   AED_REAL,  INTENT(IN) :: TA                      ! Total Alkalinity
   AED_REAL,  INTENT(IN) :: H                       ! old H+ conc
   AED_REAL,  INTENT(IN) :: Ka1                     ! 1st Acidity Const
   AED_REAL,  INTENT(IN) :: Ka2                     ! 2nd Acidity Const
   AED_REAL,  INTENT(IN) :: Kw                      ! ion product H2O
   AED_REAL,  INTENT(IN) :: bounds
   INTEGER               :: numElements
   !-- Returns the H+ conc at the forward time step
   AED_REAL   :: pH
   !-- Local
   AED_REAL   :: startpH,endpH,increment,lowValue
   INTEGER    :: lowValueLoc1,lowValueLoc2,c
   AED_REAL, DIMENSION(:), ALLOCATABLE :: Hvec
   AED_REAL, DIMENSION(:), ALLOCATABLE :: estTA
   AED_REAL, DIMENSION(:), ALLOCATABLE :: residual
   LOGICAL :: pmode

   pmode = .true.

   ALLOCATE(Hvec(numElements))
   ALLOCATE(estTA(numElements))
   ALLOCATE(residual(numElements))

     startpH = -log10(H) - bounds
     endpH   = -log10(H) + bounds
     increment = (endpH - startpH)/REAL(numElements)

     DO c=1, numElements
       Hvec(c) = startpH + ((c-1)*increment)
     END DO

     ! Get into H+ conc
     Hvec = 10**(-Hvec)

     ! Now estimate TA based on the spectrum of H+ concs
     estTA = CT * (Ka1*Hvec + 2.0*Ka1*Ka2) / &
                   (Hvec*Hvec + Ka1*Hvec + Ka1*Ka2) + (Kw/Hvec) - Hvec

     ! now calc residual
     residual = estTA - TA

     lowValue = 100.0
     lowValueLoc1 = 0
     lowValueLoc2 = 0
     DO c = 1,numElements
       IF(ABS(residual(c)) < ABS(lowValue)) THEN
         lowValueLoc1 = c
         lowValue = residual(c)
       END IF
     END DO


     IF(lowValueLoc1 >= SIZE(residual)) THEN
       ! Solution hasn't converged and residuals are high
       lowValueLoc1 = SIZE(residual) - 1
       IF(pmode) THEN
         print *,'pH iteration convergence problems -'
         print *,'   consider reducing your timestep'
       END IF
     END IF
     IF(lowValueLoc1 <= 1) THEN
       ! Solution hasn't converged and residuals are high
       lowValueLoc1 = 2
       IF(pmode) THEN
         print *,'pH iteration convergence problems -'
         print *,'   consider reducing your timestep'
       END IF
     END IF


     IF (lowValue <= zero_ .AND. residual(lowValueLoc1+1)>zero_) THEN
        lowValueLoc2 = lowValueLoc1 + 1
     ELSE IF (lowValue <= zero_ .AND. residual(lowValueLoc1-1)>zero_) THEN
        lowValueLoc2 = lowValueLoc1 - 1
     ELSE IF (lowValue > zero_ .AND. residual(lowValueLoc1+1)<zero_) THEN
        lowValueLoc2 = lowValueLoc1 + 1
     ELSE IF (lowValue > zero_ .AND. residual(lowValueLoc1-1)<zero_) THEN
        lowValueLoc2 = lowValueLoc1 - 1
     ELSE
        ! emergency
        lowValueLoc2 = lowValueLoc1 + 1
     END IF

     pH = Hvec(lowValueLoc1) + &
     (residual(lowValueLoc1) * ABS(Hvec(lowValueLoc1) - Hvec(lowValueLoc2))  &
             / (ABS(residual(lowValueLoc1)) + ABS(residual(lowValueLoc2))))


   DEALLOCATE(Hvec)
   DEALLOCATE(estTA)
   DEALLOCATE(residual)

 END FUNCTION iteratePH
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE CO2SYS(mode,temp,salt,pres,TAin,TCin,pHin,pHout,fCO2out,TAout)

  INTEGER, INTENT(IN)   :: mode
  REAL,    INTENT(IN)   :: temp, salt, pres
  REAL,    INTENT(IN)   :: TCin, TAin, pHin
  REAL,    INTENT(OUT)  :: fCO2out, pHout, TAout
  ! LOCAL
  REAL                  :: PRE, K0, KS, kF, fH, KB, KW, KP1, KP2, KP3, KSi=0., &
                           K1, K2, TB, TP, TS, TF, TSi, TC, TA

  !===========Initialize the conditions =========================!

  TB  = 0.
  TS  = 0.
  TF  = 0.
  TP  = 0./1.e6
  TSi = 0./1.e6

  PRE = 0.  ! Fix up link to pressure, pres

  Call Cal_constants(temp, salt, PRE, K0, KS, kF, fH, KB, KW, KP1, KP2, KP3,   &
                         & KSi, K1, K2, TB, TP, TS, TF)

  IF (mode==0) THEN
    TC  = TCin !/1.e6
    TA  = TAin !/1.e6
    Call Cal_pHfromTATC(TA,TC,pHout,K1,K2,TB,KB,KW,KP1,KP2,KP3,TP,TSi,TS,KS,KSi,TF,KF)
    Call Cal_fCO2fromTCpH(TC,pHout,fCO2out,K1,K2,K0)

  ELSEIF (mode==1) THEN

    Call Cal_TAfromTCpH(TCin,pHin,TAout,K1,K2,KW,KP1,KP2,KP3,TP,TSi,KSi,TB,KB,TS,KS,TF,KF)

  ENDIF

END SUBROUTINE CO2SYS


SUBROUTINE Cal_constants(TempC, Sal, PRE, K0, KS, kF, fH, KB, KW, KP1, KP2,    &
                         & KP3, KSi, K1, K2, TB, TP, TS, TF)

  REAL,   INTENT(in)      :: TempC, Sal, PRE
  REAL,   INTENT(out)     :: K0, KS, kF, fH, KB, KW, KP1, KP2, KP3, KSi,  K1, K2
  REAL,   INTENT(inout)   :: TB, TP, TS, TF
  ! local
  REAL                    :: TempK, RT, logTempK, sqrSal, TempK100, Pbar
  REAL                    :: lnK0, IonS, lnKS, lnKF, lnKBtop, lnKB, lnKW, lnKP1, &
                           & lnKP2, lnKP3, lnKSi, lnK1, lnK2, pK1, pK2
  REAL                    :: SWStoTOT, FREEtoTOT
  REAL                    ::   deltaV, kappa, lnK1fac, lnK2fac, lnKWfac, lnKFfac,&
                           & lnKSfac, lnKP1fac, lnKP2fac, lnKP3fac, lnKSifac,    &
                           & K1fac, K2fac, KWfac, KFfac, KSfac, KP1fac, KP2fac,  &
                           & KP3fac, KSifac, pHfactor, Delta, b, P1atm, FugFac,  &
                           & VPWP, VPCorrWP, VPSWWP, VPFac, lnKBfac, KBfac

  TempK    = TempC + 273.15
  RT       = 83.1451*TempK
  logTempK = log(TempK)
  sqrSal   = sqrt(Sal)
  Pbar     = PRE/10.

  TB       = (0.000232/10.811)*(Sal/1.80655) ! Total Borate, mol/kg-sw
  TF       = (0.000067/18.998)*(Sal/1.80655) ! Total Fluoride, mol/kg-sw
  TS       = (0.14/96.062)*(Sal/1.80655)     ! Total Sulfate, mol/kg-sw

  !---------- CO2 solubility -------------------!
  TempK100 = TempK/100.
  lnK0     = -60.2409 + 93.4517 / TempK100 + 23.3585 * log(TempK100) + Sal*    &
           & (0.023517 - 0.023656*TempK100 + 0.0047036 * (TempK100**2.))
  K0       = exp(lnK0)   ! mol/kg-sw/atm

  !----------- Ion Strength --------------------!
  IonS     = 19.924 * Sal / (1000. - 1.005 * Sal)

  !-------- KS for the bisulfate ion -----------!
  lnKS     = -4276.1/TempK + 141.328 - 23.093*logTempK +                       &
           & (-13856./TempK + 324.57 - 47.986*logTempK)*sqrt(IonS) +           &
           & (35474./TempK - 771.54 + 114.723*logTempK)*IonS +                 &
           & (-2698./TempK)*sqrt(IonS)*IonS + (1776./TempK)*(IonS**2)
  KS       = exp(lnKS) * (1. - 0.001005*Sal)  ! mol/kg-sw

  !-------- KF for hydrogen fluoride -----------!
  lnKF     = 1590.2/TempK - 12.641 + 1.525*sqrt(IonS)
  KF       = exp(lnKS) * (1. - 0.001005*Sal)  ! mol/kg-sw
  SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF)
  FREEtoTOT =  1 + TS/KS

  !---fH, activity coefficient of the H+ ion ---!
  fH      = 1.2948 - 0.002036*TempK + (0.0004607 - 0.000001475*TempK)*(Sal**2)

  !--------KB for boric acid --------------------!
  lnKBtop = -8966.9 - 2890.53*sqrSal - 77.942*Sal + 1.728*sqrSal*Sal - 0.0996*(Sal**2)
  lnKB    = lnKBtop/TempK + 148.0248 + 137.1942*sqrSal + 1.62142*Sal + (-24.4344 - 25.085*sqrSal &
        & - 0.2474*Sal)*logTempK + 0.053105*sqrSal*TempK
  KB      = exp(lnKB)/SWStoTOT          ! convert to SWS pH scale

  !--------- KW for water -----------------------!
  lnKW = 148.9802 - 13847.26/TempK - 23.6521*logTempK + (-5.977 + &
      & 118.67/TempK + 1.0495*logTempK)*sqrSal - 0.01615*Sal
  KW   = exp(lnKW)

  !------KP1, KP2, KP3 for phosphoric acid-------!
  lnKP1 = -4576.752/TempK + 115.54 - 18.453*logTempK + (-106.736/TempK +  &
       & 0.69171)*sqrSal + (-0.65643/TempK - 0.01844)*Sal
  KP1   = exp(lnKP1)
  lnKP2 = -8814.715/TempK + 172.1033 - 27.927*logTempK + (-160.34/TempK + &
       & 1.3566)*sqrSal + (0.37335/TempK - 0.05778)*Sal
  KP2   = exp(lnKP2)
  lnKP3 = -3070.75/TempK - 18.126 + (17.27039/TempK + 2.81197)*sqrSal + &
       & (-44.99486/TempK - 0.09984)*Sal
  KP3   = exp(lnKP3)

  !-------KSi for silicic acid--------------------!
  lnKSi = -8904.2/TempK + 117.4 - 19.334*logTempK + (-458.79/TempK + 3.5913)*sqrt(IonS) + &
       & (188.74/TempK - 1.5998)*IonS + (-12.1652/TempK + 0.07871)*(IonS**2)

  !--------K1 and K2 for carbonic acid------------!
  pK1 = 3670.7/TempK-62.008+9.7944*logTempK-0.0118*Sal+0.000116*(Sal**2)
  K1  = 10.**(-pK1)  ! this is on the SWS pH scale in mol/kg-SW
  pK2 = 1394.7/TempK + 4.777 - 0.0184*Sal + 0.000118*(Sal**2)
  K2  = 10.**(-pK2)

  !============correct constants for pressure=================!
  !------correct K1, k2, kB for pressure----------!
  deltaV  = -25.5 + 0.1271*TempC
  Kappa   = (-3.08 + 0.0877*TempC)/1000.
  lnK1fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -15.82 - 0.0219*TempC
  Kappa   = (1.13 - 0.1475*TempC)/1000
  lnK2fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -29.48 + 0.1622*TempC - 0.002608*(TempC**2)
  Kappa   = -2.84/1000.
  lnKBfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  !----correct other constants for pressure--------!

  deltaV  = -20.02 + 0.1119*TempC - 0.001409*(TempC**2)
  Kappa   = (-5.13 + 0.0794*TempC)/1000
  lnKWfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -9.78 - 0.009*TempC - 0.000942*(TempC**2)
  Kappa   = (-3.91 + 0.054*TempC)/1000
  lnKFfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -18.03 + 0.0466*TempC + 0.000316*(TempC**2)
  Kappa   = (-4.53 + 0.09*TempC)/1000
  lnKSfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -14.51 + 0.1211*TempC - 0.000321*(TempC**2)
  Kappa   = (-2.67 + 0.0427*TempC)/1000
  lnKP1fac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -23.12 + 0.1758*TempC - 0.002647*(TempC**2)
  Kappa   = (-5.15 + 0.09  *TempC)/1000
  lnKP2fac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -26.57 + 0.202 *TempC - 0.003042*(TempC**2)
  Kappa   = (-4.08 + 0.0714*TempC)/1000
  lnKP3fac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  deltaV  = -29.48 + 0.1622*TempC - 0.002608*(TempC**2)
  Kappa   = -2.84/1000
  lnKSifac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

  K1fac  = exp(lnK1fac);  K1  = K1 *K1fac
  K2fac  = exp(lnK2fac);  K2  = K2 *K2fac
  KWfac  = exp(lnKWfac);  KW  = KW *KWfac
  KBfac  = exp(lnKBfac);  KB  = KB *KBfac
  KFfac  = exp(lnKFfac);  KF  = KF *KFfac
  KSfac  = exp(lnKSfac);  KS  = KS *KSfac
  KP1fac = exp(lnKP1fac); KP1 = KP1*KP1fac
  KP2fac = exp(lnKP2fac); KP2 = KP2*KP2fac
  KP3fac = exp(lnKP3fac); KP3 = KP3*KP3fac
  KSifac = exp(lnKSifac); KSi = KSi*KSifac

  SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF)
  FREEtoTOT =  1 + TS/KS

  pHfactor = SWStoTOT

  !----convert from SWS pH scale to chosen scale------!
  K1  = K1* pHfactor; K2  = K2* pHfactor;
  KW  = KW* pHfactor; KB  = KB* pHfactor;
  KP1 = KP1*pHfactor; KP2 = KP2*pHfactor;
  KP3 = KP3*pHfactor; KSi = KSi*pHfactor;

  !-----Calculate Fugacity constants------------------!
  Delta = (57.7 - 0.118*TempK)
  b = -1636.75 + 12.0408*TempK - 0.0327957*(TempK**2) + 3.16528*0.00001*(TempK**3)

  P1atm = 1.01325 ! in bar
  FugFac = exp((b + 2.*Delta)*P1atm/RT)

  !------Calculate VP faction
  VPWP = exp(24.4543 - 67.4509*(100./TempK) - 4.8489*log(TempK/100.))
  VPCorrWP = exp(-0.000544*Sal)
  VPSWWP = VPWP*VPCorrWP
  VPFac = 1. - VPSWWP

END SUBROUTINE Cal_constants
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE Cal_pHfromTATC(TAx, TCx, pHx, K1F, K2F, TBF, KBF, KWF, KP1F, KP2F, KP3F,&
                           & TPF, TSiF, TSF, KSF, KSiF, TFF, KFF)

  REAL,     INTENT(in) :: TAx, TCx
  REAL,     INTENT(in) :: K1F, K2F, TBF, KBF, KWF, KP1F, KP2F, KP3F, TPF, TSiF, TSF, KSF, TFF, KFF, KSiF
  REAL,     INTENT(out):: pHx
  !local
  REAL                 :: H, Denom, CAlk, BAlk, OH, PhosTop, PhosBot, PAlk, SiAlk, &
                        & FREEtoTOT, Hfree, HSO4, HF, Residual, Slope, deltapH, pHTol, ln10
  INTEGER              :: count
  INTEGER              :: max_count

  pHx       = 8.        !This is the first guess
  pHTol     = 0.0001    !tolerance for iterations end
  ln10      = log(10.)
  deltapH   = pHTol + 1

  count     = 0
  max_count = 100

  DO WHILE (abs(deltapH) > pHTol .AND. count < max_count)
    H         = 10.**(-pHx)
    Denom     = (H*H + K1F*H + K1F*K2F)
    CAlk      = TCx*K1F*(H + 2.*K2F)/Denom
    BAlk      = TBF*KBF/(KBF + H)
    OH        = KWF/H
    PhosTop   = KP1F*KP2F*H + 2.*KP1F*KP2F*KP3F - H*H*H
    PhosBot   = H*H*H + KP1F*H*H + KP1F*KP2F*H + KP1F*KP2F*KP3F
    PAlk      = TPF*PhosTop/PhosBot
    SiAlk     = TSiF*KSiF/(KSiF + H)
    FREEtoTOT = (1 + TSF/KSF) !% pH scale conversion factor
    Hfree     = H/FREEtoTOT !% for H on the total scale
    HSO4      = TSF/(1 + KSF/Hfree) !% since KS is on the free scale
    HF        = TFF/(1 + KFF/Hfree) !% since KF is on the free scale
    Residual  = TAx - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF
    Slope     = ln10*(TCx*K1F*H*(H*H + K1F*K2F + 4.*H*K2F)/Denom/Denom &
              & + BAlk*H/(KBF + H) + OH + H)
    deltapH   = Residual/Slope !% this is Newton's method
    ! to keep the jump from being too big;
    DO WHILE (abs(deltapH) > 1)
       deltapH = deltapH/2
    END DO

    pHx       = pHx + deltapH  !% Is on the same scale as K1 and K2 were calculated...
    count     = count + 1
  END DO

END SUBROUTINE Cal_pHfromTATC
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE Cal_fCO2fromTCpH(TCx,pHx,fCO2x,K1,K2,K0)

  REAL, INTENT(in) :: TCx, pHx, K1, K2, K0
  REAL, INTENT(out):: fCO2x
  REAL             :: H  ! local

  H           = 10.**(-pHx)
  fCO2x       = TCx*H*H/(H*H + K1*H + K1*K2)/K0

END SUBROUTINE Cal_fCO2fromTCpH
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE Cal_TAfromTCpH(TCx,pHx,TAx,K1,K2,KW,KP1,KP2,KP3,TP,TSi,KSi,TB,KB,TS,KS,TF,KF)
! SUB CalculateTAfromTCpH, version 02.02, 10-10-97, written by Ernie Lewis.
! Inputs: TC, pH, K(...), T(...)
! Output: TA
! This calculates TA from TC and pH.
! Though it is coded for H on the total pH scale, for the pH values occuring
! in seawater (pH > 6) it will be equally valid on any pH scale (H terms
! negligible) as long as the K Constants are on that scale.

  REAL, INTENT(in) :: TCx, pHx
  REAL, INTENT(in) :: K1, K2, KW, KP1, KP2, KP3, TP, TSi, KSi, TB, KB, TS, KS, TF, KF
  REAL, INTENT(out):: TAx
  ! local
  REAL             :: H,CAlk,BAlk,OH,PhosTop,PhosBot,PAlk,SiAlk,FREEtoTOT,Hfree,HSO4,HF

  H         = 10.**(-pHx)
  CAlk      = TCx*K1*(H + 2.*K2)/(H*H + K1*H + K1*K2)
  BAlk      = TB*KB/(KB + H)
  OH        = KW/H
  PhosTop   = KP1*KP2*H + 2.*KP1*KP2*KP3 - H*H*H
  PhosBot   = H*H*H + KP1*H*H + KP1*KP2*H + KP1*KP2*KP3
  PAlk      = TP*PhosTop/PhosBot
  SiAlk     = TSi*KSi/(KSi + H)
  FREEtoTOT = (1 + TS/KS)             ! pH scale conversion factor
  Hfree     = H/FREEtoTOT             ! for H on the total scale
  HSO4      = TS/(1 + KS/Hfree)       ! since KS is on the free scale
  HF        = TF/(1 + KF/Hfree)       ! since KF is on the free scale
  TAx       = CAlk + BAlk + OH + PAlk + SiAlk - Hfree - HSO4 - HF

END SUBROUTINE Cal_TAfromTCpH
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_carbon
