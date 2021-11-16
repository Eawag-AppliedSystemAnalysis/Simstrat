!###############################################################################
!#                                                                             #
!# aed_nitrogen.F90                                                            #
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
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created 9 May 2011                                                          #
!#                                                                             #
!###############################################################################
!        .-----------------. .----------------.  .----------------.            !
!        | .--------------. || .--------------. || .--------------. |          !
!        | | ____  _____  | || |     _____    | || |  _________   | |          !
!        | ||_   \|_   _| | || |    |_   _|   | || | |  _   _  |  | |          !
!        | |  |   \ | |   | || |      | |     | || | |_/ | | \_|  | |          !
!        | |  | |\ \| |   | || |      | |     | || |     | |      | |          !
!        | | _| |_\   |_  | || |     _| |_    | || |    _| |_     | |          !
!        | ||_____|\____| | || |    |_____|   | || |   |_____|    | |          !
!        | |              | || |              | || |              | |          !
!        | '--------------' || '--------------' || '--------------' |          !
!         '----------------'  '----------------'  '----------------'           !
!###############################################################################

#include "aed.h"

!
MODULE aed_nitrogen
!------------------------------------------------------------------------------+
! aed_nitrogen --- nitrogen biogeochemical model
!
! Nitrogen module contains equations for nitrification and denitrification
!------------------------------------------------------------------------------+
   USE aed_core
   USE aed_util,  ONLY: aed_gas_piston_velocity, aed_n2o_sat

   IMPLICIT NONE

   PRIVATE    ! By default make everything private
!
   PUBLIC aed_nitrogen_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_nitrogen_data_t
      !# Variable identifiers
      INTEGER  :: id_nox, id_amm, id_n2o, id_no2
      INTEGER  :: id_oxy, id_ph, id_temp, id_salt, id_denit_product
      INTEGER  :: id_wind, id_E_depth, id_E_tau, id_E_dens, id_E_rain
      INTEGER  :: id_Fsed_amm, id_Fsed_nit, id_Fsed_n2o, id_Fsed_no2
      INTEGER  :: id_nitrf, id_denit, id_n2op, id_anammox, id_dnra
      INTEGER  :: id_sed_amm, id_sed_nit, id_sed_n2o, id_sed_no2
      INTEGER  :: id_atm_n2o, id_atm_dep
      INTEGER  :: id_cell_vel

      !# Model parameters
      AED_REAL :: Rnitrif,Rdenit,Ranammox,Rn2o,Rdnra
      AED_REAL :: Knitrif,Kdenit,Kanmx_nit,Kanmx_amm,Kdnra_oxy
      AED_REAL :: Kpart_ammox, Kin_deamm, Rno2o2, Rnh4o2, Rnh4no2
      AED_REAL :: theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit
      AED_REAL :: Fsed_amm,Fsed_nit,Fsed_n2o,Fsed_no2,Ksed_amm,Ksed_nit,Ksed_n2o
      AED_REAL :: atm_din_dd, atm_din_conc, atm_pn_dd, f_dindep_nox
      AED_REAL :: atm_n2o

      LOGICAL  :: use_oxy, use_ph
      LOGICAL  :: use_sed_model_amm,use_sed_model_nit,use_sed_model_n2o,use_sed_model_no2
      LOGICAL  :: simNitrfpH,simNitrfLight,simDryDeposition,simWetDeposition
      INTEGER  :: oxy_lim, simN2O, n2o_piston_model, Fsed_nit_model

     CONTAINS
         PROCEDURE :: define            => aed_define_nitrogen
         PROCEDURE :: calculate         => aed_calculate_nitrogen
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_nitrogen
         PROCEDURE :: calculate_surface => aed_calculate_surface_nitrogen
!        PROCEDURE :: mobility          => aed_mobility_nitrogen
!        PROCEDURE :: light_extinction  => aed_light_extinction_nitrogen
!        PROCEDURE :: delete            => aed_delete_nitrogen

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
SUBROUTINE aed_define_nitrogen(data, namlst)
!------------------------------------------------------------------------------+
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!------------------------------------------------------------------------------+
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_nitrogen_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER           :: status

!  %% NAMELIST    %%  /aed_nitrogen/
!  %% Last Checked 20/08/2021
   AED_REAL          :: n_min         = zero_
   AED_REAL          :: n_max         =   1e6
   AED_REAL          :: nit_initial   =   4.5
   AED_REAL          :: amm_initial   =   4.5
   AED_REAL          :: n2o_initial   =   0.0
   AED_REAL          :: no2_initial   =   0.0
   AED_REAL          :: Rnitrif       =   0.01
   AED_REAL          :: Rdenit        =   0.01
   AED_REAL          :: Rn2o          =   0.0015
   AED_REAL          :: Ranammox      =   0.0
   AED_REAL          :: Rdnra         =   0.0
   AED_REAL          :: Knitrif       = 150.0
   AED_REAL          :: Kdenit        = 150.0
   AED_REAL          :: Kanmx_nit     = 150.0
   AED_REAL          :: Kanmx_amm     = 150.0
   AED_REAL          :: Kdnra_oxy     = 150.0
   AED_REAL          :: Kpart_ammox   =  20.0     !Ko2_4
   AED_REAL          :: Kin_deamm     =   0.886   !Ko2_5
   AED_REAL          :: Rno2o2
   AED_REAL          :: Rnh4o2
   AED_REAL          :: Rnh4no2
   AED_REAL          :: theta_nitrif  = 1.0
   AED_REAL          :: theta_denit   = 1.0
   AED_REAL          :: Fsed_amm      = zero_
   AED_REAL          :: Fsed_nit      = zero_
   AED_REAL          :: Fsed_n2o      = zero_
   AED_REAL          :: Fsed_no2      = zero_
   AED_REAL          :: Ksed_amm      = 30.0
   AED_REAL          :: Ksed_nit      = 30.0
   AED_REAL          :: Ksed_n2o      = 30.0
   AED_REAL          :: Ksed_no2      = 30.0
   AED_REAL          :: theta_sed_amm = 1.0
   AED_REAL          :: theta_sed_nit = 1.0
   AED_REAL          :: atm_n2o       = 0.32*1e-6 !## recent atmospheric N2O conc (in ppm)
   AED_REAL          :: atm_din_dd    = zero_
   AED_REAL          :: atm_din_conc  = zero_
   AED_REAL          :: atm_pn_dd     = zero_
   AED_REAL          :: f_dindep_nox  = 0.5

   INTEGER           :: simN2O  = 0
   INTEGER           :: oxy_lim = 1
   INTEGER           :: n2o_piston_model = 4
   INTEGER           :: Fsed_nit_model = 1

   LOGICAL           :: simNitrfpH = .false.
   LOGICAL           :: simNitrfLight = .false.
   LOGICAL           :: simDryDeposition = .false.
   LOGICAL           :: simWetDeposition = .false.

   CHARACTER(len=64) :: nitrif_reactant_variable=''
   CHARACTER(len=64) :: denit_product_variable=''
   CHARACTER(len=64) :: nitrif_ph_variable=''
   CHARACTER(len=64) :: Fsed_amm_variable=''
   CHARACTER(len=64) :: Fsed_nit_variable=''
   CHARACTER(len=64) :: Fsed_n2o_variable=''
   CHARACTER(len=64) :: Fsed_no2_variable=''
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST    %%  /aed_nitrogen/

   NAMELIST /aed_nitrogen/     n_min, n_max,                                  &
                nit_initial, amm_initial, n2o_initial, no2_initial,            &
                Rnitrif,Rdenit,Fsed_amm,Fsed_nit,                              &
                Knitrif,Kdenit,Ksed_amm,Ksed_nit,                              &
                Kpart_ammox, Kin_deamm, Rno2o2, Rnh4o2, Rnh4no2,               &
                theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit,          &
                nitrif_reactant_variable,denit_product_variable,nitrif_ph_variable,&
                Fsed_amm_variable,Fsed_nit_variable,Fsed_n2o_variable,Fsed_no2_variable,&
                simN2O, atm_n2o, oxy_lim, Rn2o, Fsed_n2o, Ksed_n2o, n2o_piston_model,&
                Ranammox, Rdnra, Kanmx_nit, Kanmx_amm, Kdnra_oxy,              &
                simNitrfpH, simNitrfLight, simDryDeposition, simWetDeposition, &
                atm_din_dd, atm_din_conc, atm_pn_dd, f_dindep_nox, Fsed_nit_model, &
                diag_level
!
!------------------------------------------------------------------------------+
!BEGIN
   print *,"        aed_nitrogen initialization"

   !---------------------------------------------------------------------------+
   ! Read the namelist
   read(namlst,nml=aed_nitrogen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_nitrogen'

   !---------------------------------------------------------------------------+
   ! Store config options and parameter values in module's own derived type
   ! Note: all rates must be provided in values per day in the nml,
   ! and are converted for internal use as values per second.
   data%simNitrfpH       = simNitrfpH
   data%simNitrfLight    = simNitrfLight
   data%simN2O           = simN2O
   data%oxy_lim          = oxy_lim
   data%Fsed_nit_model   = Fsed_nit_model

   data%Rnitrif          = Rnitrif/secs_per_day
   data%Rdenit           = Rdenit/secs_per_day
   data%Rn2o             = Rn2o/secs_per_day
   data%Ranammox         = Ranammox/secs_per_day
   data%Rdnra            = Rdnra/secs_per_day
   data%Knitrif          = Knitrif
   data%Kdenit           = Kdenit
   data%Kanmx_nit        = Kanmx_nit
   data%Kanmx_amm        = Kanmx_amm
   data%Kdnra_oxy        = Kdnra_oxy
   data%Kpart_ammox      = Kpart_ammox
   data%Kin_deamm        = Kin_deamm
   data%Rno2o2           = Rno2o2/secs_per_day
   data%Rnh4o2           = Rnh4o2/secs_per_day
   data%Rnh4no2          = Rnh4no2/secs_per_day
   data%theta_nitrif     = theta_nitrif
   data%theta_denit      = theta_denit

   data%atm_n2o          = atm_n2o
   data%n2o_piston_model = n2o_piston_model

   data%Fsed_amm         = Fsed_amm/secs_per_day
   data%Fsed_nit         = Fsed_nit/secs_per_day
   data%Fsed_n2o         = Fsed_n2o/secs_per_day
   data%Fsed_no2         = Fsed_no2/secs_per_day
   data%Ksed_amm         = Ksed_amm
   data%Ksed_nit         = Ksed_nit
   data%Ksed_n2o         = Ksed_n2o
   data%theta_sed_amm    = MIN(MAX(0.5,theta_sed_amm),1.5)
   data%theta_sed_nit    = MIN(MAX(0.5,theta_sed_nit),1.5)

   data%simDryDeposition = simDryDeposition
   data%simWetDeposition = simWetDeposition
   data%atm_din_dd       = MAX(zero_,atm_din_dd/secs_per_day)
   data%atm_pn_dd        = MAX(zero_,atm_pn_dd/secs_per_day)  ! Remains unused ?
   data%atm_din_conc     = MAX(zero_,atm_din_conc)
   data%f_dindep_nox     = MIN(MAX(zero_,f_dindep_nox),one_)

   !---------------------------------------------------------------------------+
   ! Register state variables
   data%id_amm = aed_define_variable('amm','mmol/m**3','ammonium',            &
                                   amm_initial,minimum=n_min, maximum=n_max)
   data%id_nox = aed_define_variable('nit','mmol/m**3','nitrate',             &
                                   nit_initial,minimum=n_min, maximum=n_max)

   IF( simN2O>0 ) THEN
      ! TODO check here to see if oxy is simulated if simN2O is also on
     IF (nitrif_reactant_variable .NE. '') THEN
       print *,'          advanced nitrogen redox model linking to ',TRIM(nitrif_reactant_variable)
       ! this allocation is done below with other state var dependancies
       !data%id_oxy = aed_locate_variable(nitrif_reactant_variable)
       print *,'          ... found'
      ELSE
        PRINT *,'  ERROR advanced nitrogen redox set (simN2O) but no oxygen target variable is set'
        STOP
      ENDIF

      data%id_n2o = aed_define_variable('n2o','mmol/m**3','nitrous oxide',   &
                                  n2o_initial, minimum=n_min, maximum=n_max)

     IF( simN2O>1 ) THEN ! Advanced model requires NO2
       data%id_no2 = aed_define_variable('no2','mmol/m**3','nitrate',         &
                                   no2_initial, minimum=n_min, maximum=n_max)
     ENDIF
   ENDIF

   !---------------------------------------------------------------------------+
   ! Register external state variable dependencies
   data%use_oxy = nitrif_reactant_variable .NE. '' !Check for link to OXY module
   IF (data%use_oxy) data%id_oxy = aed_locate_variable(nitrif_reactant_variable)
   data%use_ph = nitrif_ph_variable .NE. ''        !Check for link to CAR module
   IF (.NOT.data%use_ph) data%simNitrfpH = .false. !Disable if no link possible
   IF (data%simNitrfpH) data%id_ph = aed_locate_variable(nitrif_ph_variable)

   data%id_Fsed_amm = -1 ; data%id_Fsed_nit = -1
   data%id_Fsed_n2o = -1 ; data%id_Fsed_no2 = -1
   data%use_sed_model_amm = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_amm) &
     data%id_Fsed_amm = aed_locate_sheet_variable(Fsed_amm_variable)
   data%use_sed_model_nit = Fsed_nit_variable .NE. ''
   IF (data%use_sed_model_nit) &
     data%id_Fsed_nit = aed_locate_sheet_variable(Fsed_nit_variable)
   data%use_sed_model_n2o = Fsed_n2o_variable .NE. ''
   IF (data%use_sed_model_n2o .and. simN2O>0 ) &
     data%id_Fsed_n2o = aed_locate_sheet_variable(Fsed_n2o_variable)
   data%use_sed_model_no2 = Fsed_no2_variable .NE. ''
   IF (data%use_sed_model_no2 .and. simN2O>1 ) &
     data%id_Fsed_no2 = aed_locate_sheet_variable(Fsed_no2_variable)

   !---------------------------------------------------------------------------+
   ! Register diagnostic variables
   data%id_sed_amm = aed_define_sheet_diag_variable('sed_amm','mmol/m**2/d','ammonium sediment flux')
   data%id_sed_nit = aed_define_sheet_diag_variable('sed_nit','mmol/m**2/d','nitrate sediment flux')
   data%id_nitrf   = aed_define_diag_variable('nitrif','mmol/m**3/d','nitrification rate')
   data%id_denit   = aed_define_diag_variable('denit','mmol/m**3/d','de-nitrification rate')
   data%id_anammox = aed_define_diag_variable('anammox','mmol/m**3/d','anammox rate')
   data%id_dnra    = aed_define_diag_variable('dnra','mmol/m**3/d','dnra rate')

   IF( simN2O>0 ) THEN
    data%id_n2op    = aed_define_diag_variable('n2oprod','mmol/m**3/d','n2o prod rate')
    data%id_atm_n2o = aed_define_sheet_diag_variable('atm_n2o_flux','mmol/m**2/d','n2o atmospheric flux')
    data%id_sed_n2o = aed_define_sheet_diag_variable('sed_n2o','mmol/m**2/d','n2o sediment flux')
    IF( simN2O>1 ) data%id_sed_no2 = aed_define_sheet_diag_variable('sed_no2','mmol/m**2/d','no2 sediment flux')
   ENDIF

   IF( simWetDeposition .OR. simDryDeposition ) THEN
    data%id_atm_dep = aed_define_sheet_diag_variable('atm_din_flux','mmol/m**2/d','din atmospheric deposition flux')
   ENDIF

   !---------------------------------------------------------------------------+
   ! Register environmental dependencies
   data%id_temp    = aed_locate_global('temperature')
   data%id_salt    = aed_locate_global('salinity')
   IF( simWetDeposition ) data%id_E_rain  = aed_locate_sheet_global('rain')
   IF( simN2O>0 )         data%id_wind    = aed_locate_sheet_global('wind_speed')
   IF( simN2O>0 )         data%id_E_depth = aed_locate_global('layer_ht')
   IF( simN2O>0 )         data%id_cell_vel= aed_locate_global('cell_vel')! needed for k600
  !IF( simN2O>0 )         data%id_E_tau   = aed_locate_global('taub')    ! tau to be converted to velocity
  !IF( simN2O>0 )         data%id_E_dens  = aed_locate_global('density') ! density needed for tau-vel

END SUBROUTINE aed_define_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_nitrogen(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Right hand sides of aed_nitrogen model
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: amm,nit,n2o,no2   ! State variables
   AED_REAL           :: oxy,temp,pH       ! Dependant variables

   AED_REAL           :: nitrification,denitrification,anammox,dnra
   AED_REAL           :: denit_n2o_prod,denit_n2o_cons,nit_n2o_prod
   AED_REAL           :: nitratation, denitratation, &
                         nitritation, denitritation, &
                         deammonification,     & !dnra,
                         nitrousation, denitrousation, &
                         ammonium_oxidation, ammonium_release

   AED_REAL,PARAMETER :: Xon  =  3.0     !ratio of O2 to N utilised during nitrification
   AED_REAL,PARAMETER :: Knev =  3.0     !Nevison nitrification O2 threshold
   AED_REAL,PARAMETER :: aa   =  0.26    !Nevison nitrification parameter
   AED_REAL,PARAMETER :: bb   = -0.0006  !Nevison nitrification parameter
   AED_REAL,PARAMETER :: Xnc  = 16./106. !OM stoichiomtery
   AED_REAL,PARAMETER :: Kn2oc=  0.3     !N2O consumption O2 poisoning
   AED_REAL,PARAMETER :: Kno3 = 5.0      !Denit NO3 half-sat
   AED_REAL,PARAMETER :: Xon1 = 1.
   AED_REAL,PARAMETER :: Xon2 = 1.
   AED_REAL,PARAMETER :: Xon3 = 1.!
!------------------------------------------------------------------------------+
!BEGIN
   n2o = zero_ ; no2 = zero_

   !-----------------------------------------------
   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp)  ! temperature

   !-----------------------------------------------
   ! Set current (local) state variable values.
   amm = _STATE_VAR_(data%id_amm)                                      ! ammonium
   nit = _STATE_VAR_(data%id_nox)                                      ! nitrate
   oxy = 300.0 ;  IF( data%id_oxy>0 )   oxy = _STATE_VAR_(data%id_oxy) ! oxygen
   ph  = 7.0   ;  IF( data%simNitrfpH )  ph = _STATE_VAR_(data%id_ph)  ! pH

   !-----------------------------------------------
   !# Define process rates in units mmol N/m3/s, using either the fully resolved
   !  GHG model, or traditional simple nitrification/denitrification model
   IF ( data%simN2O == 2 ) THEN
     !# Full N GHG model with NO2 and N2O; hetertrophy is done in OGM

     no2 = _STATE_VAR_(data%id_no2)    ! nitrite
     n2o = _STATE_VAR_(data%id_n2o)    ! nitrous oxide

     !# Nitrification (autotrophic so done here)
     nitratation        = data%Rno2o2 * no2*oxy
     ammonium_oxidation = data%Rnh4o2 * amm*oxy
     nitritation        = data%Rnh4o2 * amm*oxy * oxy/(data%Kpart_ammox+oxy)
     nitrousation       = data%Rnh4o2 * amm*oxy * data%Kpart_ammox/(data%Kpart_ammox+oxy) * 0.5
     deammonification   = data%Rnh4no2* no2*amm * data%Kin_deamm/(data%Kin_deamm+oxy)

    !IF( data%simNitrfpH ) nitrousation = nitrousation * NitrfpHFunction(pH)  ! Ask Dan

     !# De-nitrification (for info only, they are hetertrophic so are done in OGM)
     !   denitratation         = kom*om * nit/(data%K_nit+nit) * data%Ko2_1/(data%Ko2_1+oxy)
     !   denitritation         = kom*om * nit/(data%K_nit+nit) * data%Ko2_3/(data%Ko2_3+oxy)
     !   denitrousation        = kom*om * n2o/(data%K_n2o+n2o) * data%Ko2_6/(data%Ko2_6+oxy)
     !   dnra                  = denitritation * data%K_no2/(data%K_no2+no2)
     !   nitrous_denitritation = denitritation * no2/(data%K_no2+no2) * data%K_n2o/(data%K_n2o+n2o)
     !   ammonium_release      =(denitratation + denitritation + denitrousation) * data%X_nc

     ! Update nitrogen pools with the flux calculations
     _FLUX_VAR_(data%id_nox) = _FLUX_VAR_(data%id_nox) + nitratation ! - denitratation
     _FLUX_VAR_(data%id_no2) = _FLUX_VAR_(data%id_no2) + &
                                   nitritation - nitratation - deammonification ! + denitritation - denitritation
     _FLUX_VAR_(data%id_n2o) = _FLUX_VAR_(data%id_n2o) + &
                                   nitrousation ! + nitrous_denitritation - denitrousation
     _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) - &
                                   ammonium_oxidation - deammonification ! + dnra + ammonium_release
    !_FLUX_VAR_(data%id_n2 ) = _FLUX_VAR_(data%id_n2 ) + deammonification

    !# Take nitrification consumption of oxygen
    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - &
                        Xon1*nitrousation - Xon2*nitritation - Xon3*nitratation

    !# Set diagnostics, based on this sub-model
    _DIAG_VAR_(data%id_nitrf)   = ammonium_oxidation * secs_per_day
    _DIAG_VAR_(data%id_anammox) = deammonification * secs_per_day
    _DIAG_VAR_(data%id_n2op)    = nitrousation* secs_per_day
    _DIAG_VAR_(data%id_denit)   = zero_ * secs_per_day     ! Heterotrophic done in OGM
    _DIAG_VAR_(data%id_dnra)    = zero_ * secs_per_day     ! Heterotrophic done in OGM

   ELSE
     !# Simple model

     !# Nitrification
     nitrification = amm * NitrfO2Function(data%use_oxy,data%Rnitrif,data%Knitrif,data%theta_nitrif,oxy,temp)
     IF( data%simNitrfpH ) nitrification = nitrification * NitrfpHFunction(pH)
    !IF( data%simNitrfLight ) nitrification = nitrification* NitrfLightFunction(I) ! Capone Pg 238 add one day

     !# De-nitrification
     denitrification = nit * DenitO2Function(data%Rdenit,                      &
                                    data%use_oxy,data%oxy_lim,                 &
                                    data%Kdenit,data%theta_denit,Kno3,         &
                                    oxy,temp,nit)

     !# Leakign of N2O, using the Babbin model
     IF( data%simN2O == 1 ) THEN
       ! Babbin style model to capture intermediate N2O pool
       n2o = _STATE_VAR_(data%id_n2o)            ! nitrous oxide
       denit_n2o_prod = 0.5 * denitrification    !* Xnc
       denit_n2o_cons = data%Rn2o * n2o * exp(-oxy/Kn2oc)
       nit_n2o_prod   = zero_
       IF(oxy>Knev) nit_n2o_prod = ((aa/oxy)+bb)*nitrification
       !IF(oxy>Knev) nit_n2o_prod = ((aa/oxy)+bb)*remin*Xnc
     ENDIF

     !# Anammox (NO2 + NH4 + CO2 -> N2); assuming nit ~ NO2 for now
     anammox = zero_
     IF( data%use_oxy .AND. oxy < 1e-1*(1e3/32.) ) THEN
       anammox = data%Ranammox * nit/(data%Kanmx_nit+nit) * amm/(data%Kanmx_amm+amm)
     ENDIF

     !## Dissasimilatory nitrate reduction to ammonia (DNRA)
     dnra = zero_
     IF( data%use_oxy ) THEN
       dnra = data%Rdnra * data%Kdnra_oxy/(data%Kdnra_oxy+oxy) * nit
     ENDIF

     !-----------------------------------------------
     ! Set temporal derivatives
     _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) - nitrification - anammox + dnra
     _FLUX_VAR_(data%id_nox) = _FLUX_VAR_(data%id_nox) &
                             + nitrification - denitrification - anammox - dnra
     IF( data%simN2O==1 ) &
       _FLUX_VAR_(data%id_n2o) = _FLUX_VAR_(data%id_n2o)  &
                             +  (denit_n2o_prod - denit_n2o_cons + nit_n2o_prod)

     ! If an externally maintained oxygen pool is linked, take nitrification from it
     IF (data%use_oxy) &
       _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-Xon*nitrification)

     !-----------------------------------------------
     ! Export diagnostic variables
     _DIAG_VAR_(data%id_nitrf)   = nitrification * secs_per_day
     _DIAG_VAR_(data%id_denit)   = denitrification * secs_per_day
     _DIAG_VAR_(data%id_anammox) = anammox * secs_per_day
     _DIAG_VAR_(data%id_dnra)    = dnra * secs_per_day
     IF( data%simN2O==1 ) THEN
       _DIAG_VAR_(data%id_n2op) = (denit_n2o_prod + nit_n2o_prod) * secs_per_day
     ENDIF

   ENDIF

END SUBROUTINE aed_calculate_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface_nitrogen(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Air-water exchange for the aed nitrogen model. Includes N2O atmos exchange
! and wet/dry deposition, depending on the configuration
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_nitrogen_data_t),INTENT(in) :: data
   TYPE  (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, depth, rain
   AED_REAL :: vel = 0.0001
   ! State
   AED_REAL :: n2o, nox, no2, nh4
   ! Temporary variables
   AED_REAL :: n2o_atm_flux = zero_ ! N2O atmos exchange with water
   AED_REAL :: Cn2o_air = zero_     ! N2O in the air phase
   AED_REAL :: kn2o_trans = zero_   ! N2O piston velocity
   AED_REAL :: wind_hgt = 10.       ! Height of wind data, for correction
   AED_REAL :: f_pres = 1.          ! Pressure correction function (currently
                                    !     unused but needed for high altitudes)
!
!------------------------------------------------------------------------------+
!BEGIN

  !----------------------------------------------------------------------------+
  !# Water-Atmosphere exchange of N2O
  IF( data%simN2O>0 ) THEN

     !-----------------------------------------------
     ! Get the necessary environmental variables (from physical driver)
     wind_hgt = 10.
     wind = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10m above surface (m/s)
     temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
     salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
     depth = MAX( _STATE_VAR_(data%id_E_depth), one_ )
     IF (data%id_cell_vel > 0 ) vel = _STATE_VAR_(data%id_cell_vel)

     !-----------------------------------------------
     ! Retrieve current (local) state variable values.
     n2o = _STATE_VAR_(data%id_n2o)      ! Concentration of N2O in surface layer

     !-----------------------------------------------
     ! Get the surface piston velocity  (Note : THIS IS BASED ON 2D FLOWS)
     ! transfer velocity k of Ho et al. 2016
     !kn2o_trans = 0.77*( (vel**0.5)*(depth**(-0.5)) + 0.266*wind**2. ) / 3.6e5
     kn2o_trans = aed_gas_piston_velocity(wind_hgt,wind,temp,salt,vel=vel,    &
                 depth=depth,schmidt_model=2,piston_model=data%n2o_piston_model)

     !-----------------------------------------------
     ! First get the N2O concentration in the air-phase at atm interface
     ! C_N2O = F x P ; Capone (2008) pg 56
     f_pres = 1.0 ! add function here for altitude correction.
     Cn2o_air = data%atm_n2o
     Cn2o_air = Cn2o_air * aed_n2o_sat(salt,temp)

     !-----------------------------------------------
     ! Get the N2O flux:    [ mmol/m2/s = m/s * mmol/m3 ]
     n2o_atm_flux = kn2o_trans * (n2o - Cn2o_air) * f_pres

     !-----------------------------------------------
     ! Set surface exchange value (mmmol/m2/s) for AED2 ODE solution.
     _FLUX_VAR_T_(data%id_n2o) = -n2o_atm_flux

     !-----------------------------------------------
     ! Also store N2O flux across the atm/water interface as a
     ! diagnostic variable (mmmol/m2/d).
     _DIAG_VAR_S_(data%id_atm_n2o) = n2o_atm_flux * secs_per_day
  ENDIF

  !----------------------------------------------------------------------------+
  !# Atmosphere loading of DIN to the water, due to dry or wet deposition
  IF( data%simDryDeposition ) THEN
    !-----------------------------------------------
    ! Set surface exchange value (mmmol/m2/s) for AED2 ODE solution.
    _FLUX_VAR_T_(data%id_nox) = data%atm_din_dd * data%f_dindep_nox
    _FLUX_VAR_T_(data%id_amm) = data%atm_din_dd * (1.-data%f_dindep_nox)
  ENDIF

  IF( data%simWetDeposition ) THEN
    !-----------------------------------------------
    ! Get the necessary environmental variables (from physical driver)
    rain = _STATE_VAR_S_(data%id_E_rain) / secs_per_day   ! Rain (m/s)

    !-----------------------------------------------
    ! Set surface exchange value (mmmol/m2/s) for AED2 ODE solution.
    _FLUX_VAR_T_(data%id_nox) = _FLUX_VAR_T_(data%id_nox) &
                              + rain * data%atm_din_conc * data%f_dindep_nox
    _FLUX_VAR_T_(data%id_amm) = _FLUX_VAR_T_(data%id_amm) &
                              + rain * data%atm_din_conc *(1.-data%f_dindep_nox)
  ENDIF

  IF( data%simDryDeposition .OR. data%simWetDeposition ) THEN
    !-----------------------------------------------
    ! Also store deposition across the atm/water interface as a
    ! diagnostic variable (mmmol/m2/day).
    _DIAG_VAR_S_(data%id_atm_dep) = _DIAG_VAR_S_(data%id_atm_dep) &
        + (_FLUX_VAR_T_(data%id_nox) + _FLUX_VAR_T_(data%id_amm)) * secs_per_day
  ENDIF

END SUBROUTINE aed_calculate_surface_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_nitrogen(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Calculate pelagic bottom fluxes and benthic sink / source terms dis nitrogen.
! Everything in units per surface area (not volume!) per time.
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp
   ! State
   AED_REAL :: amm,nit,n2o,oxy
   ! Temporary variables
   AED_REAL :: amm_flux, nit_flux, n2o_flux, no2_flux
   AED_REAL :: Fsed_amm, Fsed_nit, Fsed_n2o, Fsed_no2
   AED_REAL :: fTa, fTo, fNO3
   
   AED_REAL,PARAMETER :: Kno3 = 5.0      !Denit NO3 half-sat
!
!------------------------------------------------------------------------------+
!BEGIN

   no2_flux = zero_ ; n2o_flux = zero_ ; amm_flux = zero_ ; nit_flux = zero_

   !-----------------------------------------------
   ! Retrieve local environmental conditions for this bottom water layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature


   !-----------------------------------------------
   ! Set the maximum flux (@20C) to use in this cell, either constant or linked
   Fsed_amm = data%Fsed_amm
   IF (data%id_Fsed_amm>0) Fsed_amm = _STATE_VAR_S_(data%id_Fsed_amm)
   Fsed_nit = data%Fsed_nit ;
   IF (data%id_Fsed_nit>0) Fsed_nit = _STATE_VAR_S_(data%id_Fsed_nit)
   Fsed_n2o = data%Fsed_n2o ;
   IF (data%id_Fsed_n2o>0) Fsed_n2o = _STATE_VAR_S_(data%id_Fsed_n2o)
   Fsed_no2 = data%Fsed_no2 ;
   IF (data%id_Fsed_no2>0) Fsed_no2 = _STATE_VAR_S_(data%id_Fsed_no2)

   !-----------------------------------------------
   ! Compute temperature scaling
   fTa = data%theta_sed_amm**(temp-20.0)
   fTo = data%theta_sed_nit**(temp-20.0)

   !-----------------------------------------------
   ! Compute actual flux based on oxygen and temperature
   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      amm_flux = Fsed_amm * data%Ksed_amm/(data%Ksed_amm+oxy) * fTa
      if(data%Fsed_nit_model == 1) THEN
         nit_flux = Fsed_nit *           oxy/(data%Ksed_nit+oxy) * fTo
      ELSE
         nit = _STATE_VAR_(data%id_nox)  
         IF(Kno3==zero_)THEN
           fNO3 = one_
         ELSE
           fNO3 = nit/(Kno3+nit)
         ENDIF
         nit_flux = Fsed_nit * data%Ksed_nit/(data%Ksed_nit+oxy) * fTo * fNO3 
      ENDIF

      IF( data%simN2O>0 ) n2o_flux = Fsed_n2o * data%Ksed_n2o/(data%Ksed_n2o+oxy) * fTa
   ELSE
      ! Sediment flux dependent on temperature only
      amm_flux = Fsed_amm * fTa
      nit_flux = Fsed_nit * fTo
   ENDIF

   !-----------------------------------------------
   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + amm_flux
   _FLUX_VAR_(data%id_nox) = _FLUX_VAR_(data%id_nox) + nit_flux
   IF( data%simN2O>0 ) _FLUX_VAR_(data%id_n2o)=_FLUX_VAR_(data%id_n2o) + n2o_flux

   ! Needs NO2
   IF( data%simN2O>1 ) _FLUX_VAR_(data%id_no2)=_FLUX_VAR_(data%id_no2) + no2_flux

   !-----------------------------------------------
   ! Store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_amm) = amm_flux*secs_per_day
   _DIAG_VAR_S_(data%id_sed_nit) = nit_flux*secs_per_day
   IF( data%simN2O>0 ) _DIAG_VAR_S_(data%id_sed_n2o) = n2o_flux*secs_per_day

   !----------------------------------------------- NEEDS ADDING
   ! Set sink and source terms for the benthos (change per surface area per sec)
   ! Note that this must include the fluxes to AND from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_amm) = _FLUX_VAR_B_(data%id_ben_amm) &
   !                              + (-amm_flux/secs_per_day)


END SUBROUTINE aed_calculate_benthic_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION NitrfO2Function(use_oxy,Rnitrif,Knitrif,theta_nitrif,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rnitrif,Knitrif,theta_nitrif,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      NitrfO2Function = Rnitrif * oxy/(Knitrif+oxy) * (theta_nitrif**(temp-20.0))
   ELSE
      NitrfO2Function = Rnitrif * (theta_nitrif**(temp-20.0))
   ENDIF
END FUNCTION NitrfO2Function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION DenitO2Function(Rdenit,use_oxy,oxy_lim,                 &
                              Kdenit,theta_denit,Kdenitnit,oxy,temp,nit)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification sensitivity to O2, NO3, & T
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL, INTENT(in) :: use_oxy
   INTEGER, INTENT(in) :: oxy_lim
   AED_REAL,INTENT(in) :: Rdenit,Kdenit,theta_denit,oxy,temp,nit,Kdenitnit
   AED_REAL :: fT, fDO, fNO3
!
!-------------------------------------------------------------------------------
!BEGIN

   fT = (theta_denit**(temp-20.0))

   fDO = one_
   IF (use_oxy) THEN
     IF (oxy_lim == 1) fDO = Kdenit/(Kdenit+oxy)
     IF (oxy_lim == 2) fDO = exp(-oxy/Kdenit)
   ENDIF

   IF(Kdenitnit==zero_)THEN
     fNO3 = one_
   ELSE
     fNO3 = nit/(Kdenitnit+nit)
   ENDIF

   DenitO2Function = Rdenit * fDO * fT * fNO3

END FUNCTION DenitO2Function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION NitrfpHFunction(pH) RESULT(limitation)
!------------------------------------------------------------------------------+
! Dependence of nitrification rate on pH
!------------------------------------------------------------------------------+
   !-- Incoming
   AED_REAL, INTENT(IN) :: pH                      ! pH in the water column
   !-- Returns the salinity function
   AED_REAL :: limitation
   !-- Local
   AED_REAL, PARAMETER :: NITpHOptMin = 7.1        ! min pH of optimum range
   AED_REAL, PARAMETER :: NITpHOptMax = 7.9        ! max pH of optimum range
   AED_REAL, PARAMETER :: NITpHTolMax = 9.0        ! upper pH tolerance
   AED_REAL, PARAMETER :: NITpHTolMin = 5.5        ! lower pH tolerance
   AED_REAL :: tmp1,tmp2
   !____________________________________________________________________!
   !                                                                    !
   !  =1                    .!---------------!.                         !
   !                      !                     !                       !
   !                    !                         !                     !
   !                   !                           !                    !
   !                  !                             !                   !
   !                 !                               !                  !
   !                !                                 !                 !
   !_______________!___________________________________!________________!
   !               !         !               !         !                !
   !                    NITpHOptMin      NITpHOptMax                    !
   !           NITpHTolMax                         NITpHTolMax          !
   !____________________________________________________________________!

   !---------------------------------------------------------------------------+
   ! pH is within the tolerance; no limitation.
   IF(pH >= NITpHOptMin .AND. pH <= NITpHOptMax) THEN
     limitation = one_
   ENDIF

   ! pH is greater than the upper bound of the optimum region
   IF(pH > NITpHOptMax) THEN
     limitation = (-pH*pH+2.0*NITpHOptMax*pH -                                 &
                  2.0*NITpHOptMax*NITpHTolMax+NITpHTolMax*NITpHTolMax)/        &
                                         ((NITpHOptMax)*(NITpHOptMax))
   ENDIF

   ! pH is less than the lower bound of optimum region (NITpHOptMin)
   ! but greater than minimum tolerance NITpHTolMin
   tmp1 = zero_
   tmp2 = zero_
   IF(pH < NITpHOptMin .AND. pH > NITpHTolMin) THEN
     tmp1 = pH-NITpHTolMin
     tmp2 = NITpHOptMin-NITpHTolMin
     limitation = (2.0*tmp1/NITpHOptMin-(tmp1*tmp1/(NITpHOptMin*NITpHOptMin))) &
                / (2.0*tmp2/NITpHOptMin-(tmp2*tmp2/(NITpHOptMin*NITpHOptMin)))
   ENDIF

   ! pH is less than the NITpHTolMin
   IF(pH <= NITpHTolMin) THEN
     limitation = zero_
   ENDIF

   ! ensure we don't go negative
   IF(limitation <= zero_) THEN
     limitation = zero_
   ENDIF

END FUNCTION NitrfpHFunction
!------------------------------------------------------------------------------!


!###############################################################################
 FUNCTION NitrfLightFunction(light) RESULT(limitation)
   !----------------------------------------------------------------------------
   ! Dependence of nitrification rate on light
   !----------------------------------------------------------------------------
   !-- Incoming
   AED_REAL, INTENT(IN) :: light              ! pH in the water column
   !-- Returns the function between 0 and 1
   AED_REAL :: limitation
   !-- Local

     limitation = one_   ! TO BE COMPLETED - SEE CAPONE 2008

 END FUNCTION NitrfLightFunction
!------------------------------------------------------------------------------!


END MODULE aed_nitrogen
