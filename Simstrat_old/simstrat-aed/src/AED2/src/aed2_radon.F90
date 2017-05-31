!###############################################################################
!#                                                                             #
!# aed2_radon.F90                                                              #
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
!# Created February 2013                                                       #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_radon
!-------------------------------------------------------------------------------
! aed2_radon --- radon biogeochemical model
!
! The AED module radon contains equations that describe exchange of
! soluable reactive radon across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_radon_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_radon_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_ss(:)
      INTEGER  :: id_temp

      !# Model parameters
      AED_REAL,ALLOCATABLE :: decay(:),settling(:), Fsed(:)

     CONTAINS
         PROCEDURE :: define            => aed2_define_radon
!        PROCEDURE :: calculate         => aed2_calculate_radon
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_radon
!        PROCEDURE :: mobility          => aed2_mobility_radon
!        PROCEDURE :: light_extinction  => aed2_light_extinction_radon
!        PROCEDURE :: delete            => aed2_delete_radon

   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_radon(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_radon_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in)        :: namlst
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_radons
   AED_REAL :: decay(100)
   AED_REAL :: settling(100)
   AED_REAL :: Fsed(100)
   AED_REAL :: trace_initial = zero_
   INTEGER  :: i
   CHARACTER(4) :: trac_name

   NAMELIST /aed2_radon/ num_radons,decay,settling,Fsed
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_radon,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_radon'

   ! Store parameter values in our own derived type

   ALLOCATE(data%id_ss(num_radons))
   ALLOCATE(data%decay(num_radons))    ; data%decay(1:num_radons)    = decay(1:num_radons)
   ALLOCATE(data%settling(num_radons)) ; data%settling(1:num_radons) = settling(1:num_radons)
   ALLOCATE(data%Fsed(num_radons))     ; data%Fsed(1:num_radons)     = Fsed(1:num_radons)

   trac_name = 'ss0'
   ! Register state variables
   DO i=1,num_radons
      trac_name(3:3) = CHAR(ICHAR('0') + i)
      data%id_ss(i) = aed2_define_variable( TRIM(trac_name),'mmol/m**3','radon', &
                                   trace_initial,minimum=zero_)
   ENDDO

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global( 'temperature')

END SUBROUTINE aed2_define_radon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_radon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_radon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_radon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS

!-------------------------------------------------------------------------------
!BEGIN

END SUBROUTINE aed2_calculate_radon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_radon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED radon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_radon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN

RETURN
   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    DO i=1,ubound(data%id_ss,1)
    ! Retrieve current (local) state variable values.
       ss = _STATE_VAR_(data%id_ss(i))

      ! Sediment flux dependent on temperature only.
      ss_flux = data%Fsed(i) * (theta_sed_ss**(temp-20.0))

      ! Transfer sediment flux value to AED2.
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + (ss_flux)
   ENDDO
END SUBROUTINE aed2_calculate_benthic_radon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_radon
