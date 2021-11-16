!###############################################################################
!#                                                                             #
!# aed_noncohesive.F90                                                         #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2018 - 2021 - The University of Western Australia                #
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
!# Created Aug 2018                                                            #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .-----------------. .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | | ____  _____  | || |     ______   | || |    _______   | |         !
!         | ||_   \|_   _| | || |   .' ___  |  | || |   /  ___  |  | |         !
!         | |  |   \ | |   | || |  / .'   \_|  | || |  |  (__ \_|  | |         !
!         | |  | |\ \| |   | || |  | |         | || |   '.___`-.   | |         !
!         | | _| |_\   |_  | || |  \ `.___.'\  | || |  |`\____) |  | |         !
!         | ||_____|\____| | || |   `._____.'  | || |  |_______.'  | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!          '----------------' '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_noncohesive
!-------------------------------------------------------------------------------
! aed_noncohesive --- noncohesive sediment model
!
! The AED2 module noncohesive contains equations that describe a
! particle, noncohesive, sediment. It is subject to processes of
! resupension and settling, and other modules may link to these pools.
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,ONLY : water_viscosity

   IMPLICIT NONE

   PRIVATE

   PUBLIC aed_noncohesive_data_t

   TYPE,extends(aed_model_data_t) :: aed_noncohesive_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_ss(:), id_ss_vvel(:)
      INTEGER,ALLOCATABLE :: id_ss_sed(:), id_sfss(:)
      INTEGER :: id_l_bot, id_tau_0, id_epsilon, id_resus
      INTEGER :: id_e_temp, id_e_taub, id_e_salt, id_e_rho
      INTEGER :: id_e_sedzone , id_sed, id_swi_dz
      INTEGER :: id_d_taub

      !# Module configuration
      INTEGER :: num_ss
      INTEGER :: resuspension, settling
      LOGICAL :: simSedimentMass

      !# Model parameters
      AED_REAL,ALLOCATABLE :: decay(:), Ke_ss(:)
      AED_REAL,ALLOCATABLE :: w_ss(:), rho_ss(:), d_ss(:)
      AED_REAL,ALLOCATABLE :: fs(:), tau_0(:)
      AED_REAL,ALLOCATABLE :: Fsed(:)
      AED_REAL             :: epsilon, kTau_0, tau_r
      AED_REAL             :: sed_porosity

     CONTAINS
         PROCEDURE :: define             => aed_define_noncohesive
         PROCEDURE :: initialize_benthic => aed_initialize_benthic_noncohesive
         PROCEDURE :: calculate          => aed_calculate_noncohesive
         PROCEDURE :: calculate_benthic  => aed_calculate_benthic_noncohesive
         PROCEDURE :: mobility           => aed_mobility_noncohesive
         PROCEDURE :: light_extinction   => aed_light_extinction_noncohesive
        !PROCEDURE :: delete             => aed_delete_noncohesive

   END TYPE

! MODULE GLOBALS
   AED_REAL :: sed_depth = 1.0
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_noncohesive(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read in and the variables simulated
!  by the model are registered with AED2 core.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_noncohesive_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status,i
   CHARACTER(4) :: ncs_name

!  %% NAMELIST   %%  /aed_noncohesive/
!  %% Last Checked 20/08/2021
   ! Set default parameter values
   INTEGER           :: num_ss          = 0
   INTEGER           :: resuspension    = 0
   INTEGER           :: settling        = 0
   LOGICAL           :: simSedimentMass = .FALSE.
   AED_REAL          :: ss_initial(100) = zero_
   AED_REAL          :: decay(100)      = zero_
   AED_REAL          :: Ke_ss(100)      = 0.02
   AED_REAL          :: w_ss(100)       = 0.
   AED_REAL          :: rho_ss(100)     = 1.6e3
   AED_REAL          :: d_ss(100)       = 1e-6
   AED_REAL          :: Fsed(100)       = zero_
   AED_REAL          :: tau_0(100)      = 0.04
   AED_REAL          :: epsilon         = 0.02
   AED_REAL          :: tau_r           = 1.0
   AED_REAL          :: kTau_0          = 1.0
   AED_REAL          :: fs(100)         = 1.0
   AED_REAL          :: sed_porosity    = 0.3
   AED_REAL          :: sed_initial     = zero_
   CHARACTER(len=64) :: macrophyte_link_var = ''
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_noncohesive/

   NAMELIST /aed_noncohesive/ num_ss, decay, Ke_ss,                         &
                              settling, w_ss, rho_ss, d_ss,                 &
                              resuspension, epsilon, tau_0, tau_r, Ktau_0,  &
                              macrophyte_link_var, Fsed, fs,                &
                              simSedimentMass, ss_initial, sed_porosity, diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_noncohesive initialization"

   ! Read the namelist
   read(namlst,nml=aed_noncohesive,iostat=status)
   IF (status /= 0) STOP 'ERROR reading namelist aed_noncohesive'

   ! Store parameter values in our own derived type
   data%num_ss = num_ss
   data%settling = settling
   data%resuspension = resuspension
   data%simSedimentMass = simSedimentMass
   data%sed_porosity = sed_porosity
   data%epsilon = epsilon
   data%kTau_0 = kTau_0
   data%tau_r = tau_r

   ! Check if it is worth going further
   IF ( num_ss < 1 ) RETURN

   ! Setup non-cohesive particle groups
   ALLOCATE(data%id_ss(num_ss)) ; ALLOCATE(data%id_ss_vvel(num_ss))
   IF ( simSedimentMass ) ALLOCATE(data%id_ss_sed(num_ss))

   ALLOCATE(data%decay(num_ss)) ; data%decay(1:num_ss) = decay(1:num_ss)
   ALLOCATE(data%Ke_ss(num_ss)) ; data%Ke_ss(1:num_ss) = Ke_ss(1:num_ss)

   ALLOCATE(data%w_ss(num_ss))  ; data%w_ss(1:num_ss)  = w_ss(1:num_ss)/secs_per_day
   ALLOCATE(data%d_ss(num_ss))  ; data%d_ss(1:num_ss)  = d_ss(1:num_ss)
   ALLOCATE(data%rho_ss(num_ss)); data%rho_ss(1:num_ss)= rho_ss(1:num_ss)

   ALLOCATE(data%Fsed(num_ss))  ; data%Fsed(1:num_ss)  = Fsed(1:num_ss)
   ALLOCATE(data%tau_0(num_ss)) ; data%tau_0(1:num_ss) = tau_0(1:num_ss)
   ALLOCATE(data%fs(num_ss))    ; data%fs(1:num_ss)    = fs(1:num_ss)

   ! Register state variables
   ncs_name = 'ss0'
   DO i=1,num_ss
     ncs_name(3:3) = CHAR(ICHAR('0') + i)
                                             ! divide settling by secs_per_day to convert m/d to m/s
     data%id_ss(i) = aed_define_variable(TRIM(ncs_name),'g/m**3','noncohesive particle group', &
                         ss_initial(i),minimum=zero_,maximum=1e4,mobility=(w_ss(i)/secs_per_day))
     data%id_ss_vvel(i) = aed_define_diag_variable(TRIM(ncs_name)//'_vvel','m/s','vertical velocity')

     IF ( simSedimentMass ) THEN
       sed_initial = data%sed_porosity * fs(i) * sed_depth * data%rho_ss(i)
       data%id_ss_sed(i) = aed_define_sheet_variable(TRIM(ncs_name)//'_sed',&
                                'g/m**2','sedimented noncohesive particles', &
                                sed_initial,minimum=zero_)
     ENDIF
   ENDDO

   ! Setup bottom diag arrays for sediment and spatially variable resuspension
   IF ( simSedimentMass ) THEN
     data%id_sed = aed_define_sheet_diag_variable('ss_sed','g/m**2','total non-cohesive sediment mass')
   ENDIF
   IF ( resuspension == 2 ) THEN
      data%id_tau_0 = aed_define_sheet_diag_variable('tau_0','N/m**2','dynamic bottom drag')
      data%id_epsilon = aed_define_sheet_diag_variable('epsilon','g/m**2/s','max resuspension rate')

      ALLOCATE(data%id_sfss(num_ss))

      ncs_name = 'fs0'
      DO i=1,num_ss
         ncs_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_sfss(i) =  aed_define_sheet_diag_variable(TRIM(ncs_name),'-', 'sediment fraction of sed size')
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

   ! Register environmental dependencies
   data%id_e_temp = aed_locate_global('temperature')
   data%id_e_salt = aed_locate_global('salinity')
   IF ( settling > 1 ) THEN
      data%id_e_rho = aed_locate_global('density')
   ENDIF

   data%id_swi_dz =  aed_define_sheet_diag_variable('swi_dz','m/s','cum. swi position change')
   IF ( resuspension > 0 ) THEN
      data%id_resus = aed_define_sheet_diag_variable('resus','g/m**2/s','resuspension rate')
      data%id_d_taub = aed_define_sheet_diag_variable('d_taub','N/m**2','taub diagnostic')
      data%id_e_taub = aed_locate_sheet_global('taub')
      data%id_e_sedzone = aed_locate_sheet_global('sed_zone')
   ENDIF

END SUBROUTINE aed_define_noncohesive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic_noncohesive(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to set initial state of NCS variables (in the sediment)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_noncohesive_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
!-------------------------------------------------------------------------------
!BEGIN
   !---------------------------------------------------------------------------+
   ! If users update fs (sediment fraction) via benthic initialisaiton file,
   ! then we need to recompute the starting mass of each sediment particle type
   !---------------------------------------------------------------------------+
   IF ( data%resuspension == 2 ) THEN
    DO i=1,data%num_ss
     _STATE_VAR_S_(data%id_ss_sed(i))=_DIAG_VAR_S_(data%id_sfss(i)) * sed_depth &
                                      * data%sed_porosity * (data%rho_ss(i)*1e3)
    ENDDO
   ENDIF
   !---------------------------------------------------------------------------+
END SUBROUTINE aed_initialize_benthic_noncohesive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_noncohesive(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_noncohesive_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: ss
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i=1,data%num_ss
      ss = _STATE_VAR_(data%id_ss(i))
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + data%decay(i)*ss
   ENDDO
 END SUBROUTINE aed_calculate_noncohesive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_noncohesive(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate bottom fluxes and benthic sink and source terms of AED noncohesive
! Everything in units per surface area (not volume!) per time, eg. g/m2
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_noncohesive_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: bottom_stress, matz

   ! State
   AED_REAL :: ss

   ! Temporary variables
   AED_REAL :: resus_flux, ss_flux, dummy_eps, dummy_tau
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. ALLOCATED(data%id_ss) ) RETURN

   resus_flux = zero_

   ! Retrieve current environmental conditions for the bottom pelagic layer
   IF ( data%resuspension  > 0) THEN
      bottom_stress = MIN( _STATE_VAR_S_(data%id_e_taub), one_ )
      _DIAG_VAR_S_(data%id_d_taub) = bottom_stress
      _DIAG_VAR_S_(data%id_resus)  = zero_
   ENDIF

   ! If (spatially variable) resuspension has plant stabilisation (e.g.,
   ! seagrass-sediment-turbidity feedback), then update tau_0
   IF ( data%resuspension == 2 .AND. data%id_l_bot > 0 ) &
      _DIAG_VAR_S_(data%id_tau_0) = data%tau_0(1) + data%kTau_0 * _STATE_VAR_S_(data%id_l_bot)

   IF ( data%simSedimentMass ) _DIAG_VAR_S_(data%id_sed) = zero_

   ! Loop through each type of particle and lift it up via resuspension
   DO i=1,ubound(data%id_ss,1)
      ! Resuspension
      IF ( data%resuspension > 0 ) THEN
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

            !!MH COORONG account for low clay conetent in more sandy MTAZ
            !matz = _STATE_VAR_S_(data%id_E_sedzone)
            !IF (matz >3) dummy_eps = dummy_eps* 0.3
         ENDIF

         IF ( bottom_stress > dummy_tau ) THEN
            resus_flux = dummy_eps * (bottom_stress - dummy_tau) / data%tau_r
         ELSE
            resus_flux = zero_
         ENDIF
         _DIAG_VAR_S_(data%id_resus) = _DIAG_VAR_S_(data%id_resus) + resus_flux
      ENDIF

      ! Constant sediment "flux"  (not for general consumption)
      ss_flux = data%Fsed(i)

      ! Transfer sediment flux value to model for ODE
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + ss_flux + resus_flux

      ! Keep track of the cumulative deviation in SWI position due to
      ! resuspension of this particle class
      _DIAG_VAR_S_(data%id_swi_dz) = _DIAG_VAR_S_(data%id_swi_dz) &
                                   - (resus_flux+ss_flux) / (data%sed_porosity * (data%rho_ss(i)*1e3))

      IF ( data%simSedimentMass ) THEN
        ! Remove/add sediment fluxes value from the sediment vars
        _FLUX_VAR_B_(data%id_ss_sed(i)) = _FLUX_VAR_B_(data%id_ss_sed(i)) &
                                        - resus_flux - ss_flux

        ! Recompute the total sediment mass, adding this group
        _DIAG_VAR_S_(data%id_sed) = _DIAG_VAR_S_(data%id_sed) + &
                                  _STATE_VAR_S_(data%id_ss_sed(i))
      ENDIF
   ENDDO

END SUBROUTINE aed_calculate_benthic_noncohesive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_noncohesive(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_noncohesive_data_t),INTENT(in) :: data
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
   DO ss_i=1,ubound(data%id_ss,1)
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(ss_i))

      ! Self-shading with contribution from background noncohesive concentration.
      extinction = extinction + (data%Ke_ss(ss_i)*ss)
   ENDDO
END SUBROUTINE aed_light_extinction_noncohesive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_mobility_noncohesive(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement velocities (+ve up; -ve down)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_noncohesive_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: ss, vvel
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: temp, rho_s
!
!-------------------------------------------------------------------------------
!BEGIN
   ! settling = 0 : no settling
   ! settling = 1 : constant settling @ w_pom
   ! settling = 2 : constant settling @ w_pom, corrected for variable density
   ! settling = 3 : settling based on Stoke's Law (calculated below)
   DO i=1,data%num_ss
      SELECT CASE (data%settling)

         CASE ( _MOB_OFF_ )
            ! disable settling by setting vertical velocity to 0
            vvel = zero_

         CASE ( _MOB_CONST_ )
            ! constant settling velocity using user provided value
            vvel = data%w_ss(i)

         CASE ( _MOB_TEMP_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_e_rho)
            temp = _STATE_VAR_(data%id_e_temp)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%w_ss(i)*mu20*pw / ( mu*pw20 )

         CASE ( _MOB_STOKES_ )
            ! settling velocity based on Stokes Law calculation and cell density
            pw = _STATE_VAR_(data%id_e_rho)            ! water density
            temp = _STATE_VAR_(data%id_e_temp)         ! water temperature
            mu = water_viscosity(temp)                 ! water dynamic viscosity
            rho_s = data%rho_ss(i)
            vvel = -9.807*(data%d_ss(i)**2.)*( rho_s-pw ) / ( 18.*mu )

         CASE DEFAULT
            ! unknown settling/migration option selection
            vvel = data%w_ss(i)

      END SELECT

      !------------------------------------------------------------------------+
      ! Set global mobility array
      mobility(data%id_ss(i)) = vvel
      _DIAG_VAR_(data%id_ss_vvel(i)) = vvel

      !------------------------------------------------------------------------+
      ! EXPERIMENTAL : SEDIMENT CUMULATION
      ! Keep track of the cumulative deviation in SWI position due to
      ! sedimentation of this particle class
      ss = _STATE_VAR_(data%id_ss(i))
      _DIAG_VAR_S_(data%id_swi_dz) = _DIAG_VAR_S_(data%id_swi_dz) - (vvel*ss) &
                                   / (data%sed_porosity * (data%rho_ss(i)*1e3))

      IF ( data%simSedimentMass ) THEN
        ! Remove/add sediment fluxes value from the sediment vars
        _FLUX_VAR_B_(data%id_ss_sed(i)) = _FLUX_VAR_B_(data%id_ss_sed(i)) - vvel*ss
      ENDIF
      !------------------------------------------------------------------------+

   ENDDO

END SUBROUTINE aed_mobility_noncohesive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_noncohesive
