!###############################################################################
!#                                                                             #
!# aed2_nitrogen.F90                                                           #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created 9 May 2011                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_nitrogen
!-------------------------------------------------------------------------------
! aed2_nitrogen --- nitrogen biogeochemical model
!
! Nitrogen module contains equations for nitrification and deitrification
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE    ! By default make everything private
!
   PUBLIC aed2_nitrogen_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_nitrogen_data_t
      !# Variable identifiers
      INTEGER  :: id_nit, id_amm !nitrate & ammonium
      INTEGER  :: id_oxy,id_denit_product
      INTEGER  :: id_temp
      INTEGER  :: id_Fsed_amm,id_Fsed_nit
      INTEGER  :: id_nitrif,id_denit
      INTEGER  :: id_sed_amm,id_sed_nit

      !# Model parameters
      AED_REAL :: Rnitrif,Rdenit,Fsed_amm,Fsed_nit,Knitrif,Kdenit,Ksed_amm,Ksed_nit, &
                          theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit
      LOGICAL  :: use_oxy,use_no2,use_sed_model_amm, use_sed_model_nit

     CONTAINS
         PROCEDURE :: define            => aed2_define_nitrogen
         PROCEDURE :: calculate         => aed2_calculate_nitrogen
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_nitrogen
!        PROCEDURE :: mobility          => aed2_mobility_nitrogen
!        PROCEDURE :: light_extinction  => aed2_light_extinction_nitrogen
!        PROCEDURE :: delete            => aed2_delete_nitrogen

   END TYPE

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_nitrogen(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_nitrogen_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status

   AED_REAL          :: nit_initial=4.5
   AED_REAL          :: nit_min=zero_
   AED_REAL          :: nit_max=nan_
   AED_REAL          :: amm_initial=4.5
   AED_REAL          :: amm_min=zero_
   AED_REAL          :: amm_max=nan_
   AED_REAL          :: Rnitrif = 0.01
   AED_REAL          :: Rdenit = 0.01
   AED_REAL          :: Fsed_amm = 3.5
   AED_REAL          :: Fsed_nit = 3.5
   AED_REAL          :: Knitrif = 150.0
   AED_REAL          :: Kdenit = 150.0
   AED_REAL          :: Ksed_amm = 30.0
   AED_REAL          :: Ksed_nit = 30.0
   AED_REAL          :: theta_nitrif = 1.0
   AED_REAL          :: theta_denit = 1.0
   AED_REAL          :: theta_sed_amm = 1.0
   AED_REAL          :: theta_sed_nit = 1.0
   CHARACTER(len=64) :: nitrif_reactant_variable=''
   CHARACTER(len=64) :: denit_product_variable=''
   CHARACTER(len=64) :: Fsed_amm_variable=''
   CHARACTER(len=64) :: Fsed_nit_variable=''


   NAMELIST /aed2_nitrogen/ nit_initial,nit_min, nit_max,                 &
                    amm_initial, amm_min, amm_max,                        &
                    Rnitrif,Rdenit,Fsed_amm,Fsed_nit,                     &
                    Knitrif,Kdenit,Ksed_amm,Ksed_nit,                     &
                    theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit, &
                    nitrif_reactant_variable,denit_product_variable,      &
                    Fsed_amm_variable, Fsed_nit_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_nitrogen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_nitrogen'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%Rnitrif  = Rnitrif/secs_per_day
   data%Rdenit   = Rdenit/secs_per_day
   data%Fsed_amm = Fsed_amm/secs_per_day
   data%Fsed_nit = Fsed_nit/secs_per_day
   data%Knitrif  = Knitrif
   data%Kdenit   = Kdenit
   data%Ksed_amm  = Ksed_amm
   data%Ksed_nit  = Ksed_nit
   data%theta_nitrif = theta_nitrif
   data%theta_denit  = theta_denit
   data%theta_sed_amm = theta_sed_amm
   data%theta_sed_nit = theta_sed_nit

   ! Register state variables
   data%id_amm = aed2_define_variable('amm','mmol/m**3','ammonium',            &
                                    amm_initial,minimum=amm_min, maximum=amm_max)
   data%id_nit = aed2_define_variable('nit','mmol/m**3','nitrate',             &
                                    nit_initial,minimum=nit_min, maximum=nit_max)
   ! Register external state variable dependencies
   data%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) THEN
     data%id_oxy = aed2_locate_variable(nitrif_reactant_variable)
   ENDIF
   data%use_no2 = denit_product_variable .NE. '' !This means n2 module switched on
   IF (data%use_no2) data%id_denit_product = aed2_locate_variable(denit_product_variable)

   data%use_sed_model_amm = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_amm) &
     data%id_Fsed_amm = aed2_locate_global_sheet(Fsed_amm_variable)
   data%use_sed_model_nit = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_nit) &
     data%id_Fsed_nit = aed2_locate_global_sheet(Fsed_nit_variable)

   ! Register diagnostic variables
   data%id_nitrif = aed2_define_diag_variable('nitrif','mmol/m**3/d',       &
                                                         'nitrification rate')
   data%id_denit = aed2_define_diag_variable('denit','mmol/m**3/d',         &
                                                         'de-nitrification rate')
   data%id_sed_amm = aed2_define_sheet_diag_variable('sed_amm','mmol/m**2/d',      &
                                                         'ammonium sediment flux')
   data%id_sed_nit = aed2_define_sheet_diag_variable('sed_nit','mmol/m**2/d',      &
                                                         'nitrate sediment flux')

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
END SUBROUTINE aed2_define_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_nitrogen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_nitrogen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: amm,nit,oxy,temp !State variables
   AED_REAL           :: nitrification,denitrification
   AED_REAL,PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to nitrogen utilised during nitrification
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Retrieve current (local) state variable values.
   amm = _STATE_VAR_(data%id_amm)! ammonium
   nit = _STATE_VAR_(data%id_nit)! nitrate
   IF (data%use_oxy) THEN ! & use_oxy
      oxy = _STATE_VAR_(data%id_oxy)! oxygen
   ELSE
      oxy = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp) ! temperature

   ! Define some intermediate quantities units mmol N/m3/day
   nitrification = fnitrif(data%use_oxy,data%Rnitrif,data%Knitrif,data%theta_nitrif,oxy,temp)
   denitrification = fdenit(data%use_oxy,data%Rdenit,data%Kdenit,data%theta_denit,oxy,temp)

   ! Set temporal derivatives
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (-amm*nitrification)
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + (amm*nitrification - nit*denitrification)

   ! If an externally maintained oxygen pool is present, take nitrification from it
   IF (data%use_oxy) then ! & use_oxy
      _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-Yoxy_nitrif*amm*nitrification)
   ENDIF
   !if (data%use_no2) then
   !   _FLUX_VAR_(data%id_denit_product) = _FLUX_VAR_(data%id_denit_product) + (denitrification)
   !end if

   ! Export diagnostic variables
   _DIAG_VAR_(data%id_nitrif) =  amm*nitrification*secs_per_day
   _DIAG_VAR_(data%id_denit) =  nit*denitrification*secs_per_day


END SUBROUTINE aed2_calculate_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_nitrogen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: amm,nit,oxy

   ! Temporary variables
   AED_REAL :: amm_flux,nit_flux
   AED_REAL :: Fsed_amm, Fsed_nit
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   amm = _STATE_VAR_(data%id_amm)! ammonium
   nit = _STATE_VAR_(data%id_nit)! nitrate

   IF (data%use_sed_model_amm) THEN
      Fsed_amm = _STATE_VAR_S_(data%id_Fsed_amm)
   ELSE
      Fsed_amm = data%Fsed_amm
   ENDIF
   IF (data%use_sed_model_nit) THEN
      Fsed_nit = _STATE_VAR_S_(data%id_Fsed_nit)
   ELSE
      Fsed_nit = data%Fsed_nit
   ENDIF

   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      amm_flux = Fsed_amm * data%Ksed_amm/(data%Ksed_amm+oxy) * (data%theta_sed_amm**(temp-20.0))
      nit_flux = Fsed_nit * oxy/(data%Ksed_nit+oxy) * (data%theta_sed_nit**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      oxy = 0.
      amm_flux = Fsed_amm * (data%theta_sed_amm**(temp-20.0))
      nit_flux = Fsed_nit * (data%theta_sed_nit**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   !_SET_BOTTOM_FLUX_(data%id_amm,amm_flux/secs_per_day)
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (amm_flux)
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + (nit_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_amm) = _FLUX_VAR_B_(data%id_ben_amm) + (-amm_flux/secs_per_day)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_amm) = amm_flux*secs_per_day
   _DIAG_VAR_S_(data%id_sed_nit) = nit_flux*secs_per_day
END SUBROUTINE aed2_calculate_benthic_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fnitrif(use_oxy,Rnitrif,Knitrif,theta_nitrif,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rnitrif,Knitrif,theta_nitrif,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fnitrif = Rnitrif * oxy/(Knitrif+oxy) * (theta_nitrif**(temp-20.0))
   ELSE
      fnitrif = Rnitrif * (theta_nitrif**(temp-20.0))
   ENDIF

END FUNCTION fnitrif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdenit(use_oxy,Rdenit,Kdenit,theta_denit,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rdenit,Kdenit,theta_denit,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdenit = Rdenit * Kdenit/(Kdenit+oxy) * (theta_denit**(temp-20.0))
   ELSE
      fdenit = Rdenit * (theta_denit**(temp-20.0))
   ENDIF

END FUNCTION fdenit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed2_nitrogen
