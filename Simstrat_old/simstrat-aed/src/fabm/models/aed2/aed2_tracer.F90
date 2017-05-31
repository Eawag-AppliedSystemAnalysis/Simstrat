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
      INTEGER,ALLOCATABLE :: id_ss(:)
      INTEGER :: id_temp, id_retain, id_taub
      INTEGER :: id_d_taub
      LOGICAL :: resuspension

      !# Model parameters
      AED_REAL,ALLOCATABLE :: decay(:), settling(:), Fsed(:)
      AED_REAL,ALLOCATABLE :: epsilon(:), tau_0(:), tau_r(:), Ke_ss(:)

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
   AED_REAL :: trace_initial = zero_
   INTEGER  :: i
   LOGICAL  :: retention_time = .FALSE.
   LOGICAL  :: resuspension = .FALSE.
   CHARACTER(4) :: trac_name

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed2_tracer/ num_tracers,decay,settling,Fsed,resuspension,epsilon,tau_0,tau_r,Ke_ss,retention_time
!
!-------------------------------------------------------------------------------
!BEGIN
   decay = zero_
   settling = -0.1
   Fsed = zero_
   epsilon = 0.02
   tau_0 = 0.01
   tau_r = 1.0
   Ke_ss = 0.02

   ! Read the namelist
   read(namlst,nml=aed2_tracer,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_tracer'

   ! Store parameter values in our own derived type

   data%num_tracers = num_tracers
   IF ( num_tracers > 0 ) THEN
      ALLOCATE(data%id_ss(num_tracers))
      ALLOCATE(data%decay(num_tracers))    ; data%decay(1:num_tracers)    = decay(1:num_tracers)
      ALLOCATE(data%settling(num_tracers)) ; data%settling(1:num_tracers) = settling(1:num_tracers)
      ALLOCATE(data%Fsed(num_tracers))     ; data%Fsed(1:num_tracers)     = Fsed(1:num_tracers)

      ALLOCATE(data%epsilon(num_tracers))  ; data%epsilon(1:num_tracers)  = epsilon(1:num_tracers)
      ALLOCATE(data%tau_0(num_tracers))    ; data%tau_0(1:num_tracers)    = tau_0(1:num_tracers)
      ALLOCATE(data%tau_r(num_tracers))    ; data%tau_r(1:num_tracers)    = tau_r(1:num_tracers)
      ALLOCATE(data%Ke_ss(num_tracers))    ; data%Ke_ss(1:num_tracers)    = Ke_ss(1:num_tracers)

      trac_name = 'ss0'
      ! Register state variables
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
                                             ! divide settling by secs_pr_day to convert m/d to m/s
         data%id_ss(i) = aed2_define_variable(TRIM(trac_name),'mmol/m**3','tracer', &
                                                  trace_initial,minimum=zero_,mobility=(settling(i)/secs_pr_day))
      ENDDO
   ENDIF
   IF (retention_time) THEN
      data%id_retain = aed2_define_variable("ret",'sec','tracer', &
                                   trace_initial,minimum=zero_)
   ELSE
      data%id_retain = -1
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   IF ( resuspension ) THEN
      data%id_taub = aed2_locate_global_sheet('taub')
      data%id_d_taub = aed2_define_sheet_diag_variable('d_taub','n/m**2',  'taub diagnostic')
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

!-------------------------------------------------------------------------------
!BEGIN
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
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. ALLOCATED(data%id_ss) ) RETURN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   IF ( data%resuspension ) THEN
      bottom_stress = _STATE_VAR_S_(data%id_taub)
      bottom_stress = MIN(bottom_stress, 100.)
      _DIAG_VAR_S_(data%id_d_taub) = bottom_stress
   ENDIF

   DO i=1,ubound(data%id_ss,1)
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(i))

      IF ( data%resuspension ) THEN
         IF (bottom_stress > data%tau_0(i)) THEN
            resus_flux = data%epsilon(i) * ( bottom_stress - data%tau_0(i)) / data%tau_r(i)
         ELSE
            resus_flux = 0.
         ENDIF
      ENDIF

      ! Sediment flux dependent on temperature only.
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
