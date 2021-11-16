!###############################################################################
!#                                                                             #
!# aed_habitat.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2016 - 2021 -  The University of Western Australia               #
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
!# Created March 2016                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_habitat_water
!-------------------------------------------------------------------------------
! aed_habitat --- habitat model
!
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_habitat_water_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_habitat_water_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_mtox
      INTEGER, ALLOCATABLE :: id_fish(:)

      !# Dependencies
      INTEGER :: id_l_ph, id_l_hab, id_l_aass, id_l_rveg, id_l_bveg
      INTEGER :: id_l_oxys
      INTEGER :: id_l_otrc, id_l_oxy, id_l_sav
      INTEGER, ALLOCATABLE :: id_l_mtox(:)

      !# Environment variables
      INTEGER :: id_E_temp, id_E_salt, id_E_bathy, id_E_matz, id_E_depth
      INTEGER :: id_E_nearlevel, id_E_extc, id_E_Io, id_E_stress, id_E_airtemp

      !# Model switches
      LOGICAL :: simFishTolerance
      LOGICAL :: simMosquitoRisk,simCyanoRisk
      LOGICAL :: simMetalTox,simClearWater

      !# Model parameters
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      AED_REAL, ALLOCATABLE :: fish_alpha(:), fish_Tmax(:), fish_Taccl(:), fish_Ocrit(:), fish_KO(:)
      INTEGER :: num_mtox, num_fish

     CONTAINS
         PROCEDURE :: define             => aed_define_habitat_water
         PROCEDURE :: calculate          => aed_calculate_habitat_water
!        PROCEDURE :: calculate_benthic  => aed_calculate_benthic_habitat_water
!        PROCEDURE :: calculate_riparian => aed_calculate_riparian_habitat_water
!        PROCEDURE :: mobility           => aed_mobility_habitat_water
!        PROCEDURE :: light_extinction   => aed_light_extinction_habitat_water
!        PROCEDURE :: delete             => aed_delete_habitat_water

   END TYPE

!-------------------------------------------------------------------------------
!MODULE VARIABLES
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_habitat_water(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED HABITAT (WATER) module
!
!  Here, the aed namelist is read and the variables exported
!  are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_habitat_water_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, f, z, status

!  %% NAMELIST   %%  /aed_habitat_water/
!  %% Last Checked 20/08/2021
   LOGICAL           :: simFishTolerance
   INTEGER           :: num_mtox
   INTEGER           :: num_fish
   LOGICAL           :: simCyanoRisk
   LOGICAL           :: simClearWater
   LOGICAL           :: simMetalTox
   AED_REAL          :: mtox_lims(10)
   AED_REAL          :: fish_alpha(20)
   AED_REAL          :: fish_Tmax(20)
   AED_REAL          :: fish_Taccl(20)
   AED_REAL          :: fish_Ocrit(20)
   AED_REAL          :: fish_KO(20)
   CHARACTER(len=40) :: mtox_vars(10)
   CHARACTER(6)      :: fish_name

! From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_habitat_water/

   CHARACTER(len=40) :: mtox_acid_link
   CHARACTER(len=40) :: mtox_aass_link
   CHARACTER(len=40) :: fish_oxy_link
   LOGICAL           :: simMosquitoRisk

   NAMELIST /aed_habitat_water/                                                   &
                               simFishTolerance, num_fish,                        &
                               fish_alpha,fish_Tmax,fish_Taccl,fish_Ocrit,fish_KO,&
                               simClearWater,                                     &
                               simCyanoRisk,                                      &
                               simMetalTox, mtox_vars, mtox_lims,                 &
                               diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_habitat_water initialization"
   print *,"          WARNING! aed_habitat model is currently under development"

   ! Default
   simClearWater = .false.
   simCyanoRisk = .false.
   simMetalTox = .false.

   simFishTolerance = .false.
   fish_alpha(:) = 1.
   fish_Tmax(:) = 40.
   fish_Taccl(:) = 0.
   fish_Ocrit(:) = 10.
   fish_KO(:) = 10.

   ! Read the namelist
   read(namlst,nml=aed_habitat_water,iostat=status)
   IF (status /= 0) STOP 'ERROR reading namelist aed_habitat_water'

   ! Update module level switches
   data%num_habitats = 0
   data%simFishTolerance = simFishTolerance ; IF(simFishTolerance) data%num_habitats=data%num_habitats+1
   data%simMetalTox      = simMetalTox      ; IF(simMetalTox) data%num_habitats=data%num_habitats+1
   data%simCyanoRisk     = simCyanoRisk     ; IF(simCyanoRisk) data%num_habitats=data%num_habitats+1
   data%simClearWater    = simClearWater    ; IF(simClearWater) data%num_habitats=data%num_habitats+1

   print *,"          ... # habitat templates simulated: ",data%num_habitats

   !----------------------------------------------------------------------------
   ! Define variables and dependencies


   !-- CONTAMINATION
   IF( simMetalTox ) THEN
     data%id_mtox =  aed_define_sheet_diag_variable('toxicity','-', 'Suitability')

     mtox_acid_link = 'CAR_pH'
     mtox_aass_link = 'ASS_uzaass'

     mtox_vars = '' ;  mtox_lims = 1.0
     num_mtox = 0
     DO i=1,10 ; IF (mtox_vars(i)  .EQ. '' ) THEN ; num_mtox = i-1 ; EXIT ; ENDIF ; ENDDO
     ALLOCATE(data%id_l_mtox(num_mtox)); ALLOCATE(data%mtox_lims(num_mtox))
     data%num_mtox = num_mtox
     DO i=1,data%num_mtox
       data%id_l_mtox(i) =  aed_locate_variable(mtox_vars(i))
       data%mtox_lims(i) =  mtox_lims(i)
       !print*,'Tox : ', TRIM(tfe_vars(i)), ' * ', data%tfe_varscale(i)
     ENDDO
   ENDIF

   !-- FISH TOLERANCE
   IF( simFishTolerance ) THEN
     fish_name = 'fish0'
     data%num_fish = MAX(MIN(num_fish,20),0)
     ALLOCATE(data%id_fish(num_fish));
     DO f=1,num_fish
       fish_name(5:5) = CHAR(ICHAR('0') + f)
       data%id_fish(f) =  aed_define_diag_variable(TRIM(fish_name),'-', 'Fish Habittat Suitability')
     ENDDO

     fish_oxy_link = 'OXY_sat'
     data%id_l_oxys  = aed_locate_global(TRIM(fish_oxy_link))

     ALLOCATE(data%fish_alpha(num_fish))
     ALLOCATE(data%fish_Tmax(num_fish))
     ALLOCATE(data%fish_Taccl(num_fish))
     ALLOCATE(data%fish_Ocrit(num_fish))
     ALLOCATE(data%fish_KO(num_fish))
     DO f=1,num_fish
       data%fish_alpha(f) =  fish_alpha(f)
       data%fish_Tmax(f)  =  fish_Tmax(f)
       data%fish_Taccl(f)  = fish_Taccl(f)
       data%fish_Ocrit(f) =  fish_Ocrit(f)
       data%fish_KO(f)    =  fish_KO(f)
     ENDDO
     print *,'yo'
   ENDIF

!   !-- GENERAL
!   IF( simGalaxiidSpawning .OR. simCharaHabitat .OR. simRuppiaHabitat) THEN
!     data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
!     data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')
!   ENDIF

   ! Register environmental dependencies
   data%id_E_salt      = aed_locate_global('salinity')
   data%id_E_extc      = aed_locate_global('extc_coef')
   data%id_E_temp      = aed_locate_global('temperature')
   data%id_E_depth     = aed_locate_global('layer_ht')
!   data%id_E_bathy     = aed_locate_global_sheet('bathy')
!   data%id_E_matz      = aed_locate_global_sheet('material')
!   data%id_E_Io        = aed_locate_global_sheet('par_sf')
!   data%id_E_airtemp   = aed_locate_global_sheet('air_temp')
!   data%id_E_stress    = aed_locate_global_sheet('taub')
!   data%id_E_nearlevel = aed_locate_global_sheet('nearest_depth')

END SUBROUTINE aed_define_habitat_water
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_habitat_water(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic habitat quality metrics.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_water_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wlevel, extc, bathy, Io, vel, tau0

   ! State
   AED_REAL :: osat

   ! Temporary variables
   INTEGER  :: i, f
   AED_REAL :: euphotic, drytime, light
   AED_REAL :: OPCTmax, Ocrit, alpha, Tmax , HSI_f, Ko

   ! Parameters
   AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
   AED_REAL, PARAMETER :: crit_salinity = 50.0
   AED_REAL, PARAMETER :: crit_leg_depth = 0.12
   AED_REAL, PARAMETER :: crit_hab_conc = 500.

   AED_REAL :: fs_sdepth , fs_substr, fs_spntem, fs_stress, fs_dewatr, fs_mattem

   AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp = 0.,falg
   AED_REAL :: pshpl, pshfl, pshsd, pass, height
   AED_REAL :: crns = 0.,creg = 0.,crht = 0.,crml = 0.
   AED_REAL :: limitation(5,6)

!-------------------------------------------------------------------------------
!BEGIN
   salt = 0.0 ; euphotic = 0.0 ; bathy = 0.0


   !---------------------------------------------------------------------------+
   !-- HABITAT TEMPLATE : FISH TOLERANCE
   IF( data%simFishTolerance ) THEN

     salt  = _STATE_VAR_(data%id_E_salt)   ! salinity g/L
     temp  = _STATE_VAR_(data%id_E_temp)   ! degC
     osat  = _STATE_VAR_(data%id_l_oxys)   ! oxy %


     DO f=1,data%num_fish
       Ocrit = data%fish_Ocrit(f)
       alpha = data%fish_alpha(f)
       Tmax  = data%fish_Tmax(f)
       Ko    = data%fish_Ko(f)

       OPCTmax = Ocrit - MIN(alpha * (Tmax - temp),zero_)
       HSI_f = zero_
       IF(osat>OPCTmax) THEN
         HSI_f   = (osat-OPCTmax)/(Ko+osat-OPCTmax)
         HSI_f   = MIN(HSI_f / ((100.-OPCTmax)/(Ko+100.-OPCTmax)),one_)
       ENDIF

       _DIAG_VAR_(data%id_fish(f)) = HSI_f
     ENDDO

   ENDIF



END SUBROUTINE aed_calculate_habitat_water
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_habitat_water
