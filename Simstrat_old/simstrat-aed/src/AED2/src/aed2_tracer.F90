!###############################################################################
!#                                                                             #
!# aed2_tracer.F90                                                             #
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
MODULE aed2_tracer
!-------------------------------------------------------------------------------
! aed2_tracer --- tracer biogeochemical model
!
! The AED2 module tracer contains equations that describe exchange of
! soluable reactive tracer across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_tracer_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_tracer_data_t
      INTEGER :: num_tracers
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_ss(:), id_sfss(:)

      AED_REAL :: tau_0_min, kTau_0
      INTEGER  :: id_lnk_id, id_tau_0, id_epsilon

      INTEGER :: id_temp, id_retain, id_taub
      INTEGER :: id_d_taub, id_resus
      INTEGER :: resuspension

      !# Model parameters
      AED_REAL,ALLOCATABLE :: decay(:), settling(:), Fsed(:), Ke_ss(:)
      AED_REAL,ALLOCATABLE :: epsilon(:), tau_0(:), tau_r(:), fs(:)

     CONTAINS
         PROCEDURE :: define            => aed2_define_tracer
         PROCEDURE :: calculate         => aed2_calculate_tracer
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_tracer
!        PROCEDURE :: mobility          => aed2_mobility_tracer
         PROCEDURE :: light_extinction  => aed2_light_extinction_tracer
!        PROCEDURE :: delete            => aed2_delete_tracer

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_tracer(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_tracer_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: num_tracers
   AED_REAL :: decay(100)
   AED_REAL :: settling(100)
   AED_REAL :: Fsed(100)
   AED_REAL :: epsilon(100)
   AED_REAL :: tau_0(100)
   AED_REAL :: tau_r(100)
   AED_REAL :: Ke_ss(100)
   AED_REAL :: fs(100)
   AED_REAL :: trace_initial = zero_
   AED_REAL :: kTau_0
   CHARACTER(len=64) :: macrophyte_link_var = ''
   INTEGER  :: i
   LOGICAL  :: retention_time = .FALSE.
   INTEGER  :: resuspension = 0
   CHARACTER(4) :: trac_name

   NAMELIST /aed2_tracer/ num_tracers,decay,settling,Fsed,Ke_ss, &
                          resuspension, epsilon, tau_0, tau_r, fs, &
                          Ktau_0, macrophyte_link_var, &
                          trace_initial, retention_time

!
!-------------------------------------------------------------------------------
!BEGIN
   decay = zero_
   settling = zero_
   Fsed = zero_
   epsilon = 0.02
   tau_0 = 0.04
   tau_r = 1.0
   Ke_ss = 0.02
   kTau_0 = 1.0
   fs = 1.0

   ! Read the namelist
   read(namlst,nml=aed2_tracer,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_tracer'

   ! Store parameter values in our own derived type
   data%num_tracers = num_tracers

   ! Setup tracers
   IF ( num_tracers > 0 ) THEN
      ALLOCATE(data%id_ss(num_tracers))
      ALLOCATE(data%decay(num_tracers))    ; data%decay(1:num_tracers)    = decay(1:num_tracers)
      ALLOCATE(data%settling(num_tracers)) ; data%settling(1:num_tracers) = settling(1:num_tracers)
      ALLOCATE(data%Fsed(num_tracers))     ; data%Fsed(1:num_tracers)     = Fsed(1:num_tracers)
      ALLOCATE(data%Ke_ss(num_tracers))    ; data%Ke_ss(1:num_tracers)    = Ke_ss(1:num_tracers)

      ALLOCATE(data%epsilon(num_tracers))  ; data%epsilon(1:num_tracers)  = epsilon(1:num_tracers)
      ALLOCATE(data%tau_0(num_tracers))    ; data%tau_0(1:num_tracers)    = tau_0(1:num_tracers)
      ALLOCATE(data%tau_r(num_tracers))    ; data%tau_r(1:num_tracers)    = tau_r(1:num_tracers)
      ALLOCATE(data%fs(num_tracers))       ; data%fs(1:num_tracers)       = fs(1:num_tracers)

      trac_name = 'ss0'
      ! Register state variables
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
                                             ! divide settling by secs_per_day to convert m/d to m/s
         data%id_ss(i) = aed2_define_variable(TRIM(trac_name),'mmol/m**3','tracer', &
                                                  trace_initial,minimum=zero_,mobility=(settling(i)/secs_per_day))
      ENDDO
   ENDIF

   ! Resuspension stuff
   data%kTau_0    = kTau_0
   ! Setup bottom arrays if spatially variable resuspension
   IF ( resuspension == 2 ) THEN
      data%id_tau_0 =  aed2_define_sheet_diag_variable('tau_0','N/m**2', 'dynamic bottom drag')
      data%id_epsilon =  aed2_define_sheet_diag_variable('epsilon','g/m**2/s', 'max resuspension rate')

      ALLOCATE(data%id_sfss(num_tracers))
      trac_name = 'fs0'
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_sfss(i) =  aed2_define_sheet_diag_variable(TRIM(trac_name),'-', 'sediment fraction of sed size')
      ENDDO

      IF ( macrophyte_link_var .NE. '' ) THEN
         data%id_lnk_id = aed2_locate_sheet_variable(macrophyte_link_var)
         IF ( data%id_lnk_id .LE. 0 ) THEN
            print *, "Macrophyte Link Variable ", TRIM(macrophyte_link_var), " is not defined."
            STOP
         ENDIF
      ENDIF
   ENDIF


   ! Retention time
   IF (retention_time) THEN
      data%id_retain = aed2_define_variable("ret",'sec','tracer',trace_initial,minimum=zero_)
   ELSE
      data%id_retain = -1
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   IF ( resuspension > 0 ) THEN
      data%id_taub = aed2_locate_global_sheet('taub')
      data%id_d_taub = aed2_define_sheet_diag_variable('d_taub','N/m**2',  'taub diagnostic')
      data%id_resus =  aed2_define_sheet_diag_variable('resus','g/m**2/s', 'resuspension rate')
   ENDIF

   data%resuspension = resuspension
END SUBROUTINE aed2_define_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_tracer(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: trc
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i=1,data%num_tracers
      trc = _STATE_VAR_(data%id_ss(i))
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + data%decay(i)*trc
   ENDDO

   IF (data%id_retain < 1) RETURN
   _FLUX_VAR_(data%id_retain) = _FLUX_VAR_(data%id_retain) + 1.0
END SUBROUTINE aed2_calculate_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_tracer(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED tracer.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss, bottom_stress

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0, resus_flux = 0.
   AED_REAL :: dummy_eps, dummy_tau
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. ALLOCATED(data%id_ss) ) RETURN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   IF ( data%resuspension  > 0) THEN
      bottom_stress = _STATE_VAR_S_(data%id_taub)
      bottom_stress = MIN(bottom_stress, 100.)
      _DIAG_VAR_S_(data%id_d_taub) = bottom_stress
      _DIAG_VAR_S_(data%id_resus) = zero_
   ENDIF

   IF ( data%resuspension == 2 .AND. data%id_lnk_id > 0 ) &
      _DIAG_VAR_S_(data%id_tau_0) = data%tau_0(1) + data%kTau_0 * _STATE_VAR_S_(data%id_lnk_id)


   DO i=1,ubound(data%id_ss,1)

      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(i))

      ! Resuspension
      IF ( data%resuspension > 0 ) THEN

         !IF ( data%resuspension == 2 .AND. i == 1 ) THEN
         !   dummy_tau = _DIAG_VAR_S_(data%id_tau_0)
         !ELSE
         !   dummy_tau =  data%tau_0(i)
         !ENDIF

         IF ( data%resuspension == 2 ) THEN
            dummy_tau = data%tau_0(i) + data%kTau_0 * _STATE_VAR_S_(data%id_lnk_id)
            dummy_eps = data%epsilon(i) * _DIAG_VAR_S_(data%id_sfss(i))
         ELSE
            dummy_tau =  data%tau_0(i)
            dummy_eps = data%epsilon(i) * data%fs(i)
         ENDIF

         IF ( bottom_stress > dummy_tau ) THEN
            resus_flux = dummy_eps * (bottom_stress - dummy_tau) / data%tau_r(i)
         ELSE
            resus_flux = 0.
         ENDIF
         _DIAG_VAR_S_(data%id_resus) = _DIAG_VAR_S_(data%id_resus) + resus_flux
      ENDIF

      ! Sediment "flux" (not sedimentation) dependent on temperature only.
      ss_flux = data%Fsed(i) * (theta_sed_ss**(temp-20.0))

      ! Transfer sediment flux value to model.
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + ss_flux + resus_flux

   ENDDO

END SUBROUTINE aed2_calculate_benthic_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_tracer(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
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

      ! Self-shading with explicit contribution from background tracer concentration.
      extinction = extinction + (data%Ke_ss(ss_i)*ss)
   ENDDO
END SUBROUTINE aed2_light_extinction_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_tracer
