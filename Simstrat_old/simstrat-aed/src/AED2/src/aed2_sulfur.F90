!###############################################################################
!#                                                                             #
!# aed2_sulfur.F90                                                             #
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
MODULE aed2_sulfur
!-------------------------------------------------------------------------------
! aed2_sulfur --- sulfur biogeochemical model
!
! The AED module sulfur contains equations that describe exchange of
! soluable reactive sulfur across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_sulfur_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_sulfur_data_t
      !# Variable identifiers
      INTEGER  :: id_so4
      INTEGER  :: id_temp
      INTEGER  :: id_sed_so4

      !# Model parameters
      AED_REAL :: Fsed_so4,Ksed_so4,theta_sed_so4
      LOGICAL  :: use_oxy,use_so4

     CONTAINS
         PROCEDURE :: define            => aed2_define_sulfur
         PROCEDURE :: calculate         => aed2_calculate_sulfur
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_sulfur
!        PROCEDURE :: mobility          => aed2_mobility_sulfur
!        PROCEDURE :: light_extinction  => aed2_light_extinction_sulfur
!        PROCEDURE :: delete            => aed2_delete_sulfur

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_sulfur(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_sulfur_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: num_sulfurs
   AED_REAL :: decay(100)
   AED_REAL :: settling(100)
   AED_REAL :: Fsed(100)

   NAMELIST /aed2_sulfur/ num_sulfurs,decay,settling,Fsed
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_sulfur,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_sulfur'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   PRINT *,'AED_SULFUR : Note this module has not been completed. Stopping.'


!  ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')

END SUBROUTINE aed2_define_sulfur
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_sulfur(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_sulfur model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_sulfur_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
!  AED_REAL           :: so4,diff_so4

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
!  so4 = _STATE_VAR_(data%id_so4)! sulfur

   ! Set temporal derivatives
!  diff_so4 = 0.

!  _FLUX_VAR_(data%id_so4) = _FLUX_VAR_(data%id_so4) + (diff_so4)


END SUBROUTINE aed2_calculate_sulfur
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_sulfur(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED sulfur.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_sulfur_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
!  AED_REAL :: so4,oxy

   ! Temporary variables
!  AED_REAL :: so4_flux
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
!  temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
!  so4 = _STATE_VAR_(data%id_so4)! sulfur

!  IF (data%use_oxy) THEN
!     ! Sediment flux dependent on oxygen and temperature
!     oxy = _STATE_VAR_(data%id_oxy)
!     so4_flux = data%Fsed_so4 * data%Ksed_so4/(data%Ksed_so4+oxy) * (data%theta_sed_so4**(temp-20.0))
!  ELSE
!     ! Sediment flux dependent on temperature only.
!     so4_flux = data%Fsed_so4 * (data%theta_sed_so4**(temp-20.0))
!  ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   !_SET_BOTTOM_FLUX_(data%id_so4,so4_flux/secs_per_day)
   !_SET_SED_FLUX_(data%id_so4,so4_flux)
!  _FLUX_VAR_(data%id_so4) = _FLUX_VAR_(data%id_so4) + (so4_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_so4) = _FLUX_VAR_B_(data%id_ben_so4) + (-so4_flux/secs_per_day)

   ! Also store sediment flux as diagnostic variable.
!  _DIAG_VAR_S_(data%id_sed_so4) = so4_flux


END SUBROUTINE aed2_calculate_benthic_sulfur
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_sulfur
