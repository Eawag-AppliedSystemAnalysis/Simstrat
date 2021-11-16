!###############################################################################
!#                                                                             #
!# aed_tracer.F90                                                              #
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
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |  _________   | || |  _______     | || |     ______   | |         !
!         | | |  _   _  |  | || | |_   __ \    | || |   .' ___  |  | |         !
!         | | |_/ | | \_|  | || |   | |__) |   | || |  / .'   \_|  | |         !
!         | |     | |      | || |   |  __ /    | || |  | |         | |         !
!         | |    _| |_     | || |  _| |  \ \_  | || |  \ `.___.'\  | |         !
!         | |   |_____|    | || | |____| |___| | || |   `._____.'  | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_tracer
!-------------------------------------------------------------------------------
! aed_tracer --- tracer biogeochemical model
!
! The AED2 module tracer contains equations that describe a
! soluble or particle tracer, including decay, sediment interaction, and
! resupension and settling
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,ONLY : water_viscosity

   IMPLICIT NONE

   PRIVATE

   PUBLIC aed_tracer_data_t

   TYPE,extends(aed_model_data_t) :: aed_tracer_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_trc(:), id_sfss(:), id_trc_vvel(:)
      INTEGER :: id_age
      INTEGER :: id_l_bot, id_tau_0, id_epsilon, id_resus
      INTEGER :: id_temp, id_taub, id_salt, id_rho
      INTEGER :: id_d_taub
      INTEGER :: id_E_sedzone

      !# Module configuration
      INTEGER :: num_tracers
      INTEGER :: resuspension, settling

      !# Model parameters
      AED_REAL,ALLOCATABLE :: decay(:), Fsed(:), Ke_ss(:)
      AED_REAL,ALLOCATABLE :: w_ss(:), rho_ss(:), d_ss(:)
      AED_REAL,ALLOCATABLE :: fs(:), tau_0(:)
      AED_REAL             :: epsilon, kTau_0, tau_r

     CONTAINS
         PROCEDURE :: define            => aed_define_tracer
         PROCEDURE :: calculate         => aed_calculate_tracer
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_tracer
         PROCEDURE :: mobility          => aed_mobility_tracer
         PROCEDURE :: light_extinction  => aed_light_extinction_tracer
        !PROCEDURE :: delete            => aed_delete_tracer

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
SUBROUTINE aed_define_tracer(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read in and the variables simulated
!  by the model are registered with AED2 core.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_tracer_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status,i
   CHARACTER(4) :: trac_name

!  %% NAMELIST   %%  /aed_tracer/
!  %% Last Checked 20/08/2021
   INTEGER           :: num_tracers
   LOGICAL           :: retention_time = .FALSE.
   INTEGER           :: resuspension   = 0
   INTEGER           :: settling       = 0
   AED_REAL          :: trace_initial  = zero_
   AED_REAL          :: decay(100)     = zero_
   AED_REAL          :: Fsed(100)      = zero_
   AED_REAL          :: Ke_ss(100)     = 0.02
   AED_REAL          :: w_ss(100)      = zero_
   AED_REAL          :: rho_ss(100)    = 1.6e3
   AED_REAL          :: d_ss(100)      = 1e-6
   AED_REAL          :: epsilon        = 0.02
   AED_REAL          :: tau_r          = 1.0
   AED_REAL          :: kTau_0         = 1.0
   AED_REAL          :: tau_0(100)     = 0.04
   AED_REAL          :: fs(100)        = 1.0
   CHARACTER(len=64) :: macrophyte_link_var = ''
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_tracer/

   NAMELIST /aed_tracer/ num_tracers, decay, Fsed, Ke_ss, &
                          settling, w_ss, rho_ss, d_ss, &
                          resuspension, epsilon, tau_0, tau_r, Ktau_0, &
                          macrophyte_link_var, fs, &
                          trace_initial, retention_time
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_tracer initialization"

   ! Read the namelist
   read(namlst,nml=aed_tracer,iostat=status)
   IF (status /= 0) STOP 'ERROR reading namelist aed_tracer'

   ! Store parameter values in our own derived type
   data%num_tracers = num_tracers
   data%resuspension = resuspension
   data%settling = settling

   data%epsilon = epsilon
   data%tau_r = tau_r
   data%kTau_0 = kTau_0

   ! Setup tracers
   IF ( num_tracers > 0 ) THEN
      ALLOCATE(data%id_trc(num_tracers)) ; ALLOCATE(data%id_trc_vvel(num_tracers))
      ALLOCATE(data%decay(num_tracers)) ; data%decay(1:num_tracers) = decay(1:num_tracers)
      ALLOCATE(data%Fsed(num_tracers))  ; data%Fsed(1:num_tracers)  = Fsed(1:num_tracers)
      ALLOCATE(data%Ke_ss(num_tracers)) ; data%Ke_ss(1:num_tracers) = Ke_ss(1:num_tracers)

      ALLOCATE(data%w_ss(num_tracers))  ; data%w_ss(1:num_tracers)  = w_ss(1:num_tracers)/secs_per_day
      ALLOCATE(data%d_ss(num_tracers))  ; data%d_ss(1:num_tracers)  = d_ss(1:num_tracers)
      ALLOCATE(data%rho_ss(num_tracers)); data%rho_ss(1:num_tracers)= rho_ss(1:num_tracers)

      ALLOCATE(data%tau_0(num_tracers)) ; data%tau_0(1:num_tracers) = tau_0(1:num_tracers)
      ALLOCATE(data%fs(num_tracers))    ; data%fs(1:num_tracers)    = fs(1:num_tracers)

      trac_name = 'tr0'
      ! Register state variables
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
                                             ! divide settling by secs_per_day to convert m/d to m/s
         data%id_trc(i) = aed_define_variable(TRIM(trac_name),'mmol/m**3','tracer', &
                                              trace_initial,minimum=zero_,maximum=1e3,mobility=(w_ss(i)/secs_per_day))
         data%id_trc_vvel(i) = aed_define_diag_variable(TRIM(trac_name)//'_vvel','m/s','vertical velocity')
      ENDDO
   ENDIF

   ! Setup bottom arrays if spatially variable resuspension
   IF ( resuspension == 2 ) THEN
      data%id_tau_0 =  aed_define_sheet_diag_variable('tau_0','N/m**2', 'dynamic bottom drag')
      data%id_epsilon =  aed_define_sheet_diag_variable('epsilon','g/m**2/s', 'max resuspension rate')

      ALLOCATE(data%id_sfss(num_tracers))

      trac_name = 'fs0'
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_sfss(i) =  aed_define_sheet_diag_variable(TRIM(trac_name),'-', 'sediment fraction of sed size')
      ENDDO

      IF ( macrophyte_link_var .NE. '' ) THEN
         data%id_l_bot = aed_locate_sheet_variable(macrophyte_link_var)
         IF ( data%id_l_bot .LE. 0 ) THEN
            print *, "Macrophyte Link Variable ", TRIM(macrophyte_link_var), " is not defined."
            STOP
         ENDIF
      ELSE
         data%id_l_bot = 0
      ENDIF
   ENDIF

   ! Retention time
   IF (retention_time) THEN
      data%id_age = aed_define_variable("age",'sec','tracer',trace_initial,minimum=zero_)
   ELSE
      data%id_age = -1
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed_locate_global('temperature')
   data%id_salt = aed_locate_global('salinity')
   IF ( settling > 1 ) THEN
      data%id_rho = aed_locate_global('density')
   ENDIF
   IF ( resuspension > 0 ) THEN
      data%id_taub = aed_locate_sheet_global('taub')
      data%id_E_sedzone = aed_locate_sheet_global('material')
      data%id_d_taub = aed_define_sheet_diag_variable('d_taub','N/m**2',  'taub diagnostic')
      data%id_resus =  aed_define_sheet_diag_variable('resus','g/m**2/s', 'resuspension rate')
   ENDIF

END SUBROUTINE aed_define_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_tracer(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_tracer_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: trc
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i=1,data%num_tracers
      trc = _STATE_VAR_(data%id_trc(i))
      _FLUX_VAR_(data%id_trc(i)) = _FLUX_VAR_(data%id_trc(i)) + data%decay(i)*trc
   ENDDO

   IF (data%id_age < 1) RETURN
   _FLUX_VAR_(data%id_age) = _FLUX_VAR_(data%id_age) + 1.0
END SUBROUTINE aed_calculate_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_tracer(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED tracer.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_tracer_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss, bottom_stress, matz

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0, resus_flux = 0.
   AED_REAL :: dummy_eps, dummy_tau
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. ALLOCATED(data%id_trc) ) RETURN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   IF ( data%resuspension  > 0) THEN
      bottom_stress = _STATE_VAR_S_(data%id_taub)
      bottom_stress = MIN(bottom_stress, 1.)
      _DIAG_VAR_S_(data%id_d_taub) = bottom_stress
      _DIAG_VAR_S_(data%id_resus) = zero_
   ENDIF

   IF ( data%resuspension == 2 .AND. data%id_l_bot > 0 ) &
      _DIAG_VAR_S_(data%id_tau_0) = data%tau_0(1) + data%kTau_0 * _STATE_VAR_S_(data%id_l_bot)


   DO i=1,ubound(data%id_trc,1)
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_trc(i))

      ! Resuspension
      IF ( data%resuspension > 0 ) THEN
         !IF ( data%resuspension == 2 .AND. i == 1 ) THEN
         !   dummy_tau = _DIAG_VAR_S_(data%id_tau_0)
         !ELSE
         !   dummy_tau =  data%tau_0(i)
         !ENDIF

         IF ( data%resuspension == 2 ) THEN
            IF (data%id_l_bot > 0) THEN
               dummy_tau = data%tau_0(i) + data%kTau_0 * _STATE_VAR_S_(data%id_l_bot)
            ELSE
               dummy_tau = data%tau_0(i)
            ENDIF
            dummy_eps = data%epsilon * _DIAG_VAR_S_(data%id_sfss(i))
         ELSE
            dummy_tau = data%tau_0(i)
            dummy_eps = data%epsilon * data%fs(i)

            !MH CORRONG account for low clay conetent in more sandy MTAZ
            matz = _STATE_VAR_S_(data%id_E_sedzone)
            IF (matz >3) dummy_eps = dummy_eps* 0.3
         ENDIF

         IF ( bottom_stress > dummy_tau ) THEN
            resus_flux = dummy_eps * (bottom_stress - dummy_tau) / data%tau_r
         ELSE
            resus_flux = zero_
         ENDIF
         _DIAG_VAR_S_(data%id_resus) = _DIAG_VAR_S_(data%id_resus) + resus_flux
      ENDIF

      ! Sediment "flux" (not sedimentation) dependent on temperature only.
      ss_flux = data%Fsed(i) * (theta_sed_ss**(temp-20.0))

      ! Transfer sediment flux value to model.
      _FLUX_VAR_(data%id_trc(i)) = _FLUX_VAR_(data%id_trc(i)) + ss_flux + resus_flux
   ENDDO

END SUBROUTINE aed_calculate_benthic_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_tracer(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_tracer_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: ss
   INTEGER  :: ss_i
!
!-----------------------------------------------------------------------
!BEGIN
   DO ss_i=1,ubound(data%id_trc,1)
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_trc(ss_i))

      ! Self-shading with contribution from background tracer concentration.
      extinction = extinction + (data%Ke_ss(ss_i)*ss)
   ENDDO
END SUBROUTINE aed_light_extinction_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_mobility_tracer(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement velocities (+ve up; -ve down)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_tracer_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: vvel
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: temp, rho_s
!
!-------------------------------------------------------------------------------
!BEGIN
   ! settling = 0 : no settling
   ! settling = 1 : constant settling @ w_pom
   ! settling = 2 : constant settling @ w_pom, corrected for variable density
   ! settling = 3 : settling based on Stoke's Law (calculated below)

   DO i=1,data%num_tracers
      SELECT CASE (data%settling)

         CASE ( _MOB_OFF_ )
            ! disable settling by setting vertical velocity to 0
            vvel = zero_

         CASE ( _MOB_CONST_ )
            ! constant settling velocity using user provided value
            vvel = data%w_ss(i)

         CASE ( _MOB_TEMP_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_rho)
            temp = _STATE_VAR_(data%id_temp)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%w_ss(i)*mu20*pw / ( mu*pw20 )

         CASE ( _MOB_STOKES_ )
            ! settling velocity based on Stokes Law calculation and cell density
            pw = _STATE_VAR_(data%id_rho)              ! water density
            temp = _STATE_VAR_(data%id_temp)
            mu = water_viscosity(temp)                 ! water dynamic viscosity
            rho_s = data%rho_ss(i)
            vvel = -9.807*(data%d_ss(i)**2.)*( rho_s-pw ) / ( 18.*mu )

         CASE DEFAULT
            ! unknown settling/migration option selection
            vvel = data%w_ss(i)

      END SELECT
      ! set global mobility array
      mobility(data%id_trc(i)) = vvel
      _DIAG_VAR_(data%id_trc_vvel(i)) = vvel
   ENDDO

END SUBROUTINE aed_mobility_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_tracer
