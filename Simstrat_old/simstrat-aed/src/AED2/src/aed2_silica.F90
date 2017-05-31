!###############################################################################
!#                                                                             #
!# aed2_silica.F90                                                             #
!#                                                                             #
!# Created 24 August 2011                                                      #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_silica
!-------------------------------------------------------------------------------
! Nitrogen module contains equations for nitrification and deitrification
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE   ! By default make everything private.
!
   PUBLIC aed2_silica_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_silica_data_t
      !# Variable identifiers
      INTEGER  :: id_rsi,id_oxy
      INTEGER  :: id_Fsed_rsi
      INTEGER  :: id_temp
      INTEGER  :: id_sed_rsi

      !# Model parameters
      AED_REAL :: Fsed_rsi,Ksed_rsi,theta_sed_rsi
      LOGICAL  :: use_oxy,use_sed_model

     CONTAINS
         PROCEDURE :: define            => aed2_define_silica
         PROCEDURE :: calculate         => aed2_calculate_silica
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_silica
!        PROCEDURE :: mobility          => aed2_mobility_silica
!        PROCEDURE :: light_extinction  => aed2_light_extinction_silica
!        PROCEDURE :: delete            => aed2_delete_silica

   END TYPE

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_silica(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: namlst
   CLASS (aed2_silica_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status

   AED_REAL          :: rsi_initial=4.5
   AED_REAL          :: rsi_min=zero_
   AED_REAL          :: rsi_max=nan_
   AED_REAL          :: Fsed_rsi = 3.5
   AED_REAL          :: Ksed_rsi = 30.0
   AED_REAL          :: theta_sed_rsi = 1.0
   CHARACTER(len=64) :: silica_reactant_variable=''
   CHARACTER(len=64) :: Fsed_rsi_variable=''

   NAMELIST /aed2_silica/ rsi_initial,rsi_min,rsi_max,Fsed_rsi,Ksed_rsi,theta_sed_rsi, &
                         silica_reactant_variable,Fsed_rsi_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_silica,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_silica'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%Fsed_rsi  = Fsed_rsi/secs_per_day
   data%Ksed_rsi  = Ksed_rsi
   data%theta_sed_rsi = theta_sed_rsi
   data%use_oxy = silica_reactant_variable .NE. '' !This means oxygen module switched on

   ! Register state variables
   data%id_rsi = aed2_define_variable('rsi','mmol/m**3', 'silica',     &
                                    rsi_initial,minimum=rsi_min,maximum=rsi_max)

   ! Register external state variable dependencies
   IF (data%use_oxy) &
      data%id_oxy = aed2_locate_variable(silica_reactant_variable)

   data%use_sed_model = Fsed_rsi_variable .NE. ''
   IF (data%use_sed_model) &
      data%id_Fsed_rsi = aed2_locate_global_sheet(Fsed_rsi_variable)

   ! Register diagnostic variables
   data%id_sed_rsi = aed2_define_sheet_diag_variable('sed_rsi','mmol/m**2/d', &
                     'Si exchange across sed/water interface')

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
END SUBROUTINE aed2_define_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_silica(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_silica model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_silica_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!ARGUMENT
   !AED_REAL                   :: rsi,oxy,temp,tss !State variables
!
!-------------------------------------------------------------------------------
!BEGIN


END SUBROUTINE aed2_calculate_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_silica(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED silica.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_silica_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: rsi,oxy

   ! Temporary variables
   AED_REAL :: rsi_flux, Fsed_rsi
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   rsi = _STATE_VAR_(data%id_rsi)! silica

   IF (data%use_sed_model) THEN
       Fsed_rsi = _STATE_VAR_S_(data%id_Fsed_rsi)
   ELSE
       Fsed_rsi = data%Fsed_rsi
   ENDIF

   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
       oxy = _STATE_VAR_(data%id_oxy)
       rsi_flux = Fsed_rsi * data%Ksed_rsi/(data%Ksed_rsi+oxy) * (data%theta_sed_rsi**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
       rsi_flux = Fsed_rsi * (data%theta_sed_rsi**(temp-20.0))
   ENDIF

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   _FLUX_VAR_(data%id_rsi) = _FLUX_VAR_(data%id_rsi) + (rsi_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this should include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_rsi) = _FLUX_VAR_B_(data%id_ben_rsi) + (-rsi_flux)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_rsi) = rsi_flux*secs_per_day


END SUBROUTINE aed2_calculate_benthic_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_silica
