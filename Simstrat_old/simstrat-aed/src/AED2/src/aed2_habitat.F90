!###############################################################################
!#                                                                             #
!# aed2_habitat.F90                                                            #
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
!# Created March 2016                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_habitat
!-------------------------------------------------------------------------------
! aed2_habitat --- habitat model
!
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_habitat_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_habitat_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_bird, id_mtox
      !# Dependencies
      INTEGER :: id_l_ph, id_l_hab, id_l_aass, id_l_rveg, id_l_bveg
      INTEGER, ALLOCATABLE :: id_l_mtox(:)
      !# Environment variables
      INTEGER :: id_E_temp, id_E_salt, id_E_bathy, id_E_matz, id_E_depth, id_E_nearlevel, id_E_extc

      !# Model parameters
      LOGICAL :: simBirdForaging,simBenthicProd,simFishTolerance,simMetalTox,simCyanoRisk !,simMBORisk, simRuppiaGerm
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      INTEGER :: num_mtox

     CONTAINS
         PROCEDURE :: define             => aed2_define_habitat
         PROCEDURE :: calculate          => aed2_calculate_habitat
         PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_habitat
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_habitat
!        PROCEDURE :: mobility           => aed2_mobility_habitat
!        PROCEDURE :: light_extinction   => aed2_light_extinction_habitat
!        PROCEDURE :: delete             => aed2_delete_habitat

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_habitat(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_habitat_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, status, num_habitats, num_mtox
   LOGICAL :: simBirdForaging, simBenthicProd, simFishTolerance, simMetalTox, simCyanoRisk
   AED_REAL          :: mtox_lims(10)
   CHARACTER(len=64) :: bird_acid_link, bird_habs_link, bird_aass_link, bird_rveg_link, bird_bveg_link
   CHARACTER(len=40) :: mtox_acid_link, mtox_aass_link, mtox_vars(10)


   NAMELIST /aed2_habitat/ simBirdForaging, simBenthicProd, simFishTolerance, simMetalTox, simCyanoRisk, mtox_vars, mtox_lims
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Default
   simBirdForaging = .false.
   simBenthicProd = .false.
   simFishTolerance = .false.
   simMetalTox = .false.
   simCyanoRisk = .false.

   ! Read the namelist
   read(namlst,nml=aed2_habitat,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_habitat'

   ! Store parameter values in our own derived type
!   data%num_habitats = num_habitats
   data%simFishTolerance = simFishTolerance
   data%simBirdForaging = simBirdForaging
   data%simMetalTox = simMetalTox
   data%simBenthicProd = simBenthicProd
   data%simCyanoRisk = simCyanoRisk


   ! Define variables and dependencies
   IF( simBirdForaging ) THEN
     data%id_bird =  aed2_define_sheet_diag_variable('bird_wading','%', 'Suitability')

     bird_acid_link = 'CAR_pH'
     bird_habs_link = 'PHY_grn'
     bird_aass_link = 'ASS_uzaass'
     bird_rveg_link = '' ! 'VEG_xxx'
     bird_bveg_link = '' ! 'MAC_xxx'

     data%id_l_ph    = aed2_locate_global(TRIM(bird_acid_link))
     data%id_l_hab   = aed2_locate_global(TRIM(bird_habs_link))
     data%id_l_aass  = aed2_locate_global_sheet(TRIM(bird_aass_link))
     data%id_l_rveg  = 0 !aed2_locate_global_sheet(TRIM(bird_rveg_link))
     data%id_l_bveg  = 0 !aed2_locate_global_sheet(TRIM(bird_rveg_link))
   ENDIF
   IF( simMetalTox ) THEN
     data%id_mtox =  aed2_define_sheet_diag_variable('toxicity','-', 'Suitability')

     mtox_acid_link = 'CAR_pH'
     mtox_aass_link = 'ASS_uzaass'

     mtox_vars = '' ;  mtox_lims = 1.0
     DO i=1,10 ; IF (mtox_vars(i)  .EQ. '' ) THEN ; num_mtox = i-1 ; EXIT ; ENDIF ; ENDDO
     ALLOCATE(data%id_l_mtox(num_mtox)); ALLOCATE(data%mtox_lims(num_mtox))
     data%num_mtox = num_mtox
     DO i=1,data%num_mtox
       data%id_l_mtox(i) =  aed2_locate_variable(mtox_vars(i))
       data%mtox_lims(i) =  mtox_lims(i)
       !print*,'Tox : ', TRIM(tfe_vars(i)), ' * ', data%tfe_varscale(i)
     ENDDO
     IF ( .NOT.simBirdForaging ) data%id_l_ph   = aed2_locate_global(TRIM(mtox_acid_link))
     IF ( .NOT.simBirdForaging ) data%id_l_aass = aed2_locate_global_sheet(TRIM(mtox_aass_link))
   ENDIF


   ! Register environmental dependencies
   data%id_E_depth = aed2_locate_global('layer_ht')
   data%id_E_temp = aed2_locate_global('temperature')
   data%id_E_salt = aed2_locate_global('salinity')
   data%id_E_extc = aed2_locate_global('extc_coef')
   data%id_E_matz = aed2_locate_global_sheet('material')
   data%id_E_bathy = aed2_locate_global_sheet('bathy')
   data%id_E_nearlevel = aed2_locate_global_sheet('nearest_depth')


END SUBROUTINE aed2_define_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_habitat(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: trc
!
!-------------------------------------------------------------------------------
!BEGIN
!  DO i=1,data%num_habitats
!     trc = _STATE_VAR_(data%id_ss(i))
!     _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + data%decay(i)*trc
!  ENDDO

!  IF (data%id_retain < 1) RETURN
!  _FLUX_VAR_(data%id_retain) = _FLUX_VAR_(data%id_retain) + 1.0
END SUBROUTINE aed2_calculate_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_habitat(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED habitat.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss, bottom_stress

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0, resus_flux = 0., dummy_tau
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN


END SUBROUTINE aed2_calculate_benthic_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian_habitat(data,column,layer_idx, pc_wet)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED habitat.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wlevel, extc, bathy, matz

   ! State
   AED_REAL :: depth, ph, hab, sdepth, uzaass, aass, conc, mtox

   ! Temporary variables
   INTEGER  :: i
   AED_REAL :: bird_acid, bird_soil, bird_rveg, bird_bveg, bird_salt, bird_dept, bird_habs
   AED_REAL :: euphotic

   ! Parameters
   AED_REAL, PARAMETER :: crit_soil_acidity = 100.000
   AED_REAL, PARAMETER :: crit_water_ph = 6.0
   AED_REAL, PARAMETER :: crit_soil_type = 5.5 ! (matz 6&7 is sand)
   AED_REAL, PARAMETER :: crit_rveg_depth = 0.6
   AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
   AED_REAL, PARAMETER :: crit_salinity = 50.0
   AED_REAL, PARAMETER :: crit_leg_depth = 0.12
   AED_REAL, PARAMETER :: crit_hab_conc = 500.


!-------------------------------------------------------------------------------
!BEGIN
   !-- HABITAT TEMPLATE 1: Bird Foraging
   IF( data%simBirdForaging ) THEN

     bird_acid=one_; bird_soil=one_; bird_rveg=one_;
     bird_bveg=one_; bird_salt=one_; bird_dept=one_; bird_habs=one_

     IF( data%id_l_ph>0 .AND. pc_wet>0.1 ) THEN
       ph = _STATE_VAR_(data%id_l_pH)
     ELSE
       ph = 8.0
     ENDIF
     IF( data%id_E_depth>0 ) THEN
       depth = _STATE_VAR_(data%id_E_depth)
     ELSE
       depth = zero_
     ENDIF
     sdepth = zero_
     IF( pc_wet < 0.1 ) THEN
       bathy = _STATE_VAR_S_(data%id_E_bathy)
       wlevel = _STATE_VAR_S_(data%id_E_nearlevel)
       sdepth = bathy - wlevel
       depth = zero_
       aass = zero_
     ELSE
       extc = _STATE_VAR_(data%id_E_extc)
       euphotic = log(0.9) / (-1.*extc)
       aass = zero_
     ENDIF
     IF( data%id_l_hab>0 ) THEN
       hab = _STATE_VAR_(data%id_l_hab)
     ELSE
       hab = zero_
     ENDIF

     ! Effect of acid
     IF( ph<crit_water_ph .OR. aass>crit_soil_acidity) bird_acid = zero_

     ! Requirement for mud-flat
     IF( matz > crit_soil_type ) bird_soil = zero_

     ! Presence of riparian/emergent vegetation (veg is located >0.6mAHD)
     IF( bathy < crit_rveg_depth ) bird_rveg = zero_

     ! Presence or liklihood of benthic vegetation
     IF( (depth-sav_height) > euphotic .OR. pc_wet < 0.1 ) bird_bveg = zero_

     ! Check for tolerable salinity
     IF( salt > crit_salinity ) bird_salt = zero_

     ! Check for relevant leg/beak depth
     IF( depth > crit_leg_depth .OR. sdepth>crit_leg_depth ) bird_dept = zero_

     ! Limit suitability of Harmful Algal Bloom (HAB)
     IF( hab > crit_hab_conc ) bird_habs = zero_

     _DIAG_VAR_S_(data%id_bird) =  (bird_acid+bird_soil+bird_rveg+bird_bveg+bird_habs)/5. &
                                   * bird_salt * bird_dept

   ENDIF

   !-- HABITAT TEMPLATE 2: Metal Toxicity
   IF( data%simMetalTox ) THEN

     mtox = zero_
     DO  i=1,data%num_mtox
        conc = _STATE_VAR_( data%id_l_mtox(i) )
        IF ( conc > data%mtox_lims(i) ) mtox = one_
     ENDDO
     IF( pc_wet < 0.1 ) THEN
        aass = _STATE_VAR_S_( data%id_l_aass )
        IF ( aass < crit_soil_acidity ) mtox = one_
     ELSE
        aass = _STATE_VAR_( data%id_l_ph )
        IF ( aass < crit_water_ph ) mtox = one_
     END IF

     _DIAG_VAR_S_(data%id_mtox) = mtox

   ENDIF

   !-- HABITAT TEMPLATE 3: Fish Tolerance
   IF( data%simFishTolerance ) THEN

   !  _DIAG_VAR_S_(data%id_fish) = mtox

   ENDIF

END SUBROUTINE aed2_calculate_riparian_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_habitat(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
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
!  DO ss_i=1,ubound(data%id_ss,1)
!     ! Retrieve current (local) state variable values.
!     ss = _STATE_VAR_(data%id_ss(ss_i))

!     ! Self-shading with explicit contribution from background habitat concentration.
!     extinction = extinction + (data%Ke_ss(ss_i)*ss)
!  ENDDO
END SUBROUTINE aed2_light_extinction_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_habitat
