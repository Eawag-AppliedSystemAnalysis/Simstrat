!###############################################################################
!#                                                                             #
!# aed_oxygen.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
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
!#   Created May 2011                                                          #
!#   Track changes on GitHub @ https://github.com/AquaticEcoDynamics/libaed2   #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ____     | || |  ____  ____  | || |  ____  ____  | |         !
!         | |   .'    `.   | || | |_  _||_  _| | || | |_  _||_  _| | |         !
!         | |  /  .--.  \  | || |   \ \  / /   | || |   \ \  / /   | |         !
!         | |  | |    | |  | || |    > `' <    | || |    \ \/ /    | |         !
!         | |  \  `--'  /  | || |  _/ /'`\ \_  | || |    _|  |_    | |         !
!         | |   `.____.'   | || | |____||____| | || |   |______|   | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

MODULE aed_oxygen
!-------------------------------------------------------------------------------
! aed_oxygen --- oxygen biogeochemical model
!
! The AED module oxygen contains equations that describe exchange of
! oxygen across the air/water interface and sediment flux. Other modules can
! also add or consume oxygen. A summary fo the module is provided online at:
! http://aquatic.science.uwa.edu.au/research/models/AED/aed_oxygen.html
!
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,  ONLY: aed_gas_piston_velocity, aed_oxygen_sat

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_oxygen_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_oxygen_data_t
      !# Variable identifiers
      INTEGER  :: id_oxy
      INTEGER  :: id_Fsed_oxy
      INTEGER  :: id_oxy_sat
      INTEGER  :: id_atm_oxy_exch, id_sed_oxy
      INTEGER  :: id_sed_oxy_pel, id_atm_oxy_exch3d
      INTEGER  :: id_temp, id_salt, id_wind
      INTEGER  :: id_larea, id_lht, id_cell_vel

      !# Model parameters
      AED_REAL :: Fsed_oxy,Ksed_oxy,theta_sed_oxy
      INTEGER  :: oxy_piston_model
      LOGICAL  :: use_sed_model
      AED_REAL :: altitude

     CONTAINS
         PROCEDURE :: define            => aed_define_oxygen
         PROCEDURE :: calculate_surface => aed_calculate_surface_oxygen
         PROCEDURE :: calculate         => aed_calculate_oxygen
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_oxygen

   END TYPE

! MODULE GLOBALS
   INTEGER :: diag_level = 10             ! 0 = no diagnostic outputs
                                          ! 1 = basic diagnostic outputs
                                          ! 2-10 = most diagnostic outputs
                                          ! >10 = debug/checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_oxygen(data, namlst)
!-------------------------------------------------------------------------------
! Setup and initialise the aed_oxygen model
!
!  Here, the oxygen namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_oxygen_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

!  %% NAMELIST   %%  /aed_oxygen/
!  %% Last Checked 20/08/2021
   AED_REAL          :: oxy_initial       = 300.  !% initial dissolved oxygen (DO) concentration
                                                  !% $$mmol\,m^{-3}$$
                                                  !% float
                                                  !%
                                                  !% 0 - 400
                                                  !% Note: will be overwritten by GLM or TFV IC

   AED_REAL          :: oxy_min           = 0.    !% minimum dissolved oxygen (DO) concentration
                                                  !% $$mmol\,m^{-3}$$
                                                  !% float
                                                  !%
                                                  !%
                                                  !% Optional variable to enforce negative number clipping

   AED_REAL          :: oxy_max           = nan_  !% maxmium dissolved oxygen (DO) concentration
                                                  !% $$mmol\,m^{-3}$$
                                                  !% float
                                                  !% -
                                                  !% 1000
                                                  !% Optional variable to enforce high number clipping

   AED_REAL          :: Fsed_oxy          = -20.0 !% sediment oxygen demand (SOD)
                                                  !% $$mmol\,m^{-2}\,day^{-1}$$
                                                  !% float
                                                  !%
                                                  !% -100
                                                  !% Note: unused if Fsed_oxy_variable is activated via aed_sedflux

   AED_REAL          :: Ksed_oxy          = 30.0  !% half-saturation concentration of oxygen sediment flux
                                                  !% $$mmol\,m^{-3}$$
                                                  !% float
                                                  !%   50
                                                  !% 10-100
                                                  !% Changes the sensitivity of the oxygen flux to the
                                                  !-     overlying oxygen concentration

   AED_REAL          :: theta_sed_oxy     = 1.0   !% Arrhenius temperature multiplier for sediment oxygen flux
                                                  !% -
                                                  !% float
                                                  !% 1e+00
                                                  !% 1.04 - 1.12
                                                  !% Changes the sensitivity of the oxygen flux to
                                                  !-     the overlying temperature

   CHARACTER(len=64) :: Fsed_oxy_variable = ''    !% oxygen sediment flux variable link
                                                  !% -
                                                  !% string
                                                  !% '
                                                  !% e.g.: SDF_Fsed_oxy
                                                  !% will use the value supplied by the aed_sedflux
                                                  !-     model for Fsed_oxy; use this option to allow for
                                                  !-     spatial or temperal variation

   INTEGER           :: oxy_piston_model  = 1     !% specifies the atm exchange piston model
                                                  !% -
                                                  !% integer
                                                  !% 1
                                                  !% 1 - X
                                                  !% Choice depends on waterbody type


   AED_REAL          :: altitude     = 0.0        !% Altitude of site above sea level
                                                  !% -
                                                  !% float
                                                  !% 1e+00
                                                  !% 0 - 4000
                                                  !% Changes oxygen solubility

! %% From Module Globals
!  INTEGER :: diag_level = 10             ! 0 = no diagnostic outputs
!                                         ! 1 = basic diagnostic outputs
!                                         ! 2-10 = most diagnostic outputs
!                                         ! >10 = debug/checking outputs
!  %% END NAMELIST   %%  /aed_oxygen/

   NAMELIST /aed_oxygen/ oxy_initial, oxy_min, oxy_max,            &
                          Fsed_oxy, Ksed_oxy, theta_sed_oxy,       &
                          Fsed_oxy_variable, oxy_piston_model,     &
                          altitude,                                &
                          diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_oxygen initialization"

   ! Read the namelist
   read(namlst,nml=aed_oxygen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_oxygen'

   ! Store parameter values in the modules own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%Ksed_oxy = Ksed_oxy
   data%Fsed_oxy = Fsed_oxy/secs_per_day
   data%theta_sed_oxy = theta_sed_oxy
   data%use_sed_model = Fsed_oxy_variable .NE. ''
   data%oxy_piston_model = oxy_piston_model
   data%altitude = altitude

   ! Register state variables
   data%id_oxy = aed_define_variable('oxy','mmol/m**3','oxygen',   &
                                    oxy_initial,minimum=oxy_min,maximum=oxy_max)

   ! Register the link to external variables
   IF (data%use_sed_model) data%id_Fsed_oxy = aed_locate_sheet_variable(Fsed_oxy_variable)

   ! Register diagnostic variables
   IF (diag_level>0) THEN
     data%id_oxy_sat = aed_define_diag_variable(                   &
                     'sat', '%', 'oxygen saturation')

     data%id_sed_oxy = aed_define_sheet_diag_variable(             &
                     'sed_oxy', 'mmol/m**2/d', 'O2 exchange across sed/water interface')

     data%id_atm_oxy_exch = aed_define_sheet_diag_variable(        &
                     'atm_oxy_flux', 'mmol/m**2/d', 'O2 exchange across atm/water interface')
    IF (diag_level>10) THEN
     data%id_sed_oxy_pel = aed_define_diag_variable(               &
                     'sed_oxy_pel', 'mmol/m**2/d', 'O2 exchange across sed/water interface')

     data%id_atm_oxy_exch3d = aed_define_diag_variable(      &
                     'atm_oxy_exch3d', 'mmol/m**3/d', 'Oxygen exchange across atm/water interface')
    ENDIF
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed_locate_global('temperature') ! Temperature (degrees Celsius)
   data%id_salt = aed_locate_global('salinity') ! Salinity (psu)
!  data%id_pres = aed_locate_sheet_global('pressure') ! Pressure (dbar = 10 kPa)
   data%id_wind = aed_locate_sheet_global('wind_speed') ! Wind speed at 10 m above surface (m/s)
   data%id_larea = aed_locate_sheet_global('layer_area')
   data%id_lht = aed_locate_global('layer_ht')
   data%id_cell_vel = -1
   IF( oxy_piston_model>3 )data%id_cell_vel= aed_locate_global('cell_vel')! needed for k600

END SUBROUTINE aed_define_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface_oxygen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-water exchange for the aed oxygen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_oxygen_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, depth
   AED_REAL :: vel = 0.0001

   ! State
   AED_REAL :: oxy

   ! Temporary variables
   AED_REAL :: oxy_atm_flux = zero_  ! Surface atm flux of O2
   AED_REAL :: Coxy_air = zero_      ! Dissolved oxygen in the air phase
   AED_REAL :: koxy_trans = zero_    ! k600 for O2
   AED_REAL :: windHt  = 10.0        ! Height of U10 sensor
   AED_REAL :: f_pres  = 1.0         ! Pressure correction, only applicable at high altitudes
!
!-------------------------------------------------------------------------------
!BEGIN
   !Get dependent state variables from physical driver
   temp  = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt  = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind  = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10 m above surface (m/s)
   windHt= 10.                          ! Assumed wind height of 10m
   depth = MAX( _STATE_VAR_(data%id_lht), one_ )
   IF (data%id_cell_vel > 0 )  vel = _STATE_VAR_(data%id_cell_vel)

   ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy) ! Concentration of oxygen in surface layer

  !koxy_trans = aed_gas_piston_velocity(windHt,wind,temp,salt)
   koxy_trans = aed_gas_piston_velocity(windHt,wind,temp,salt,               &
                                         vel=vel,                             &
                                         depth=depth,                         &
                                         schmidt_model=2,                     &
                                         piston_model=data%oxy_piston_model)

   ! First get the oxygen concentration in the air phase at the interface
   ! (taken from Riley and Skirrow, 1974)
   f_pres = 1.0    ! set pressure function here using data%id_pres
   IF(data%altitude>1) f_pres = aed_oxygen_fp(data%altitude,air_temp=10.)
   Coxy_air = f_pres * aed_oxygen_sat(salt,temp)

   ! Get the oxygen flux
   oxy_atm_flux = koxy_trans * (Coxy_air - oxy)

   ! Transfer surface exchange value to AED2 (mmmol/m2/s) converted by driver
   _FLUX_VAR_T_(data%id_oxy) = oxy_atm_flux

   ! Also store oxygen flux across the atm/water interface as a diagnostic (mmmol/m2/day)
   IF (diag_level>0) THEN
     _DIAG_VAR_S_(data%id_atm_oxy_exch) = oxy_atm_flux * secs_per_day
     _DIAG_VAR_(data%id_oxy_sat) =  Coxy_air
   ENDIF

   ! Also store oxygen flux across the atm/water interface as a diagnostic (mmmol/m2/day)
   IF (diag_level>10) &
     _DIAG_VAR_(data%id_atm_oxy_exch3d) = oxy_atm_flux * secs_per_day / _STATE_VAR_(data%id_lht)

END SUBROUTINE aed_calculate_surface_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_oxygen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_oxygen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS(aed_oxygen_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: oxy, temp, salt
   AED_REAL :: f_pres, coxy_sat

!-------------------------------------------------------------------------------
!BEGIN
   ! Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)

   ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen

   ! Compute the oxygen saturation for diagnostic output
   f_pres = 1.0 ! set pressure function here using data%id_pres
   IF(data%altitude>1) f_pres = aed_oxygen_fp(data%altitude,air_temp=10.)
   coxy_sat = f_pres * aed_oxygen_sat(salt,temp)

   ! Export diagnostic variables
   _DIAG_VAR_(data%id_oxy_sat) =  (oxy/coxy_sat)*100.

END SUBROUTINE aed_calculate_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_oxygen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED oxygen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_oxygen_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: oxy

   ! Temporary variables
   AED_REAL :: oxy_flux, Fsed_oxy
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen

   IF (data%use_sed_model) THEN
       Fsed_oxy = _DIAG_VAR_S_(data%id_Fsed_oxy)
   ELSE
       Fsed_oxy = data%Fsed_oxy
   ENDIF

   ! Compute the sediment flux dependent on overlying oxygen & temperature
   oxy_flux = Fsed_oxy * MIN(3.,oxy/(data%Ksed_oxy+oxy) * (data%theta_sed_oxy**(temp-20.0)))
!print*, "Oxy oxy ben = ", oxy, "oxy_flux ", oxy_flux

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (oxy_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_oxy) = _FLUX_VAR_B_(data%id_ben_oxy) + (-oxy_flux)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_oxy) = oxy_flux * secs_per_day
   IF (diag_level>10) _DIAG_VAR_(data%id_sed_oxy_pel) = oxy_flux * secs_per_day

END SUBROUTINE aed_calculate_benthic_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed_oxygen_fp(altitude,air_temp)
!-------------------------------------------------------------------------------
! An extra function is included to account for the effect of a
! non-standard atmosphere (i.e. high altitudes)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: altitude,air_temp
!LOCAL
   AED_REAL, PARAMETER  :: p_SL = 101.325      ! MSL Pressure
   AED_REAL             :: p_vap, p_H
!
!-------------------------------------------------------------------------------
!BEGIN

  ! Pressure at altitude
  p_H = p_SL * exp((9.81/(287.0*0.0065)) * log((288.0-0.0065*altitude)/288.0))

  p_vap = zero_
  IF( air_temp > zero_ ) THEN
     IF (((-216961.*(1./(air_temp))-3840.7)*(1./(air_temp))+16.4754) > &
                                              1+MINEXPONENT(air_temp)/2) THEN
       p_vap = exp((-216961.*(1./(air_temp))-3840.7) * (1./(air_temp))+16.4754)
     END IF
  END IF

  aed_oxygen_fp = (p_H/p_SL) * ( (1.0-(p_vap/p_H)) / (1.0-(p_vap/p_SL)) )

END FUNCTION aed_oxygen_fp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_oxygen
