!###############################################################################
!#                                                                             #
!# aed2_iron.F90                                                               #
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
!# Created June 2012                                                           #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_iron
!-------------------------------------------------------------------------------
! aed2_iron --- iron biogeochemical model
!
! The AED module iron contains equations that describe exchange of
! soluable reactive iron across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_iron_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_iron_data_t
      !# Variable identifiers
      INTEGER  :: id_fe3
      INTEGER  :: id_temp
      INTEGER  :: id_sed_fe3

      !# Model parameters
      AED_REAL :: Fsed_fe3,Ksed_fe3,theta_sed_fe3
      LOGICAL  :: use_oxy,use_fe3

     CONTAINS
         PROCEDURE :: define            => aed2_define_iron
         PROCEDURE :: calculate         => aed2_calculate_iron
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_iron
!        PROCEDURE :: mobility          => aed2_mobility_iron
!        PROCEDURE :: light_extinction  => aed2_light_extinction_iron
!        PROCEDURE :: delete            => aed2_delete_iron

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_iron(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_iron_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status

   AED_REAL          :: fe3_initial=4.5
   AED_REAL          :: fe3_min = zero_
   AED_REAL          :: fe3_max = nan_
   AED_REAL          :: Fsed_fe3 = 3.5
   AED_REAL          :: Ksed_fe3 = 30.0
   AED_REAL          :: theta_sed_fe3 = 1.0
   CHARACTER(len=64) :: iron_reactant_variable=''

   INTEGER  :: num_irons
   AED_REAL          :: decay(100)
   AED_REAL          :: settling(100)
   AED_REAL          :: Fsed(100)

   NAMELIST /aed2_iron/ num_irons,decay,settling,Fsed
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_iron,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_iron'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   PRINT *,'AED_IRON : Note this module has not been completed. Stopping.'
   STOP

!  ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
END SUBROUTINE aed2_define_iron
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_iron(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_iron model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_iron_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: fe3,diff_fe3

!-------------------------------------------------------------------------------
!BEGIN
!  ! Retrieve current (local) state variable values.
!  fe3 = _STATE_VAR_(data%id_fe3)! iron

!  ! Set temporal derivatives
!  diff_fe3 = 0.

!  _FLUX_VAR_(data%id_fe3) = _FLUX_VAR_(data%id_fe3) + (diff_fe3)
END SUBROUTINE aed2_calculate_iron
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_benthic_iron(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED iron.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_iron_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
!  AED_REAL :: temp

   ! State
!  AED_REAL :: fe3,oxy

   ! Temporary variables
!  AED_REAL :: fe3_flux


!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
!  temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
!  fe3 = _STATE_VAR_(data%id_fe3)! iron

!  IF (data%use_oxy) THEN
!     ! Sediment flux dependent on oxygen and temperature
!     oxy = _STATE_VAR_(data%id_oxy)
!     fe3_flux = data%Fsed_fe3 * data%Ksed_fe3/(data%Ksed_fe3+oxy) * (data%theta_sed_fe3**(temp-20.0))
!  ELSE
!     ! Sediment flux dependent on temperature only.
!     fe3_flux = data%Fsed_fe3 * (data%theta_sed_fe3**(temp-20.0))
!  ENDIF

!  ! TODO:
!  ! (1) Get benthic sink and source terms (sccb?) for current environment
!  ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

!  ! Set bottom fluxes for the pelagic (change per surface area per second)
!  ! Transfer sediment flux value to AED2.
!  !_SET_BOTTOM_FLUX_(data%id_fe3,fe3_flux/secs_per_day)
!  !_SET_SED_FLUX_(data%id_fe3,fe3_flux)
!  _FLUX_VAR_(data%id_fe3) = _FLUX_VAR_(data%id_fe3) + (fe3_flux)

!  ! Set sink and source terms for the benthos (change per surface area per second)
!  ! Note that this must include the fluxes to and from the pelagic.
!  !_FLUX_VAR_B_(data%id_ben_fe3) = _FLUX_VAR_B_(data%id_ben_fe3) + (-fe3_flux/secs_per_day)

!  ! Also store sediment flux as diagnostic variable.
!  _DIAG_VAR_S_(data%id_sed_fe3) = fe3_flux


END SUBROUTINE aed2_calculate_benthic_iron
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_iron
