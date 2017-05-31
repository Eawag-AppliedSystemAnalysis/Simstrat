!###############################################################################
!#                                                                             #
!# aed2_common.F90                                                             #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
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
!# Created Aug 2013                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!#define CELL_VAR    1
!#define CELL_SV     2
!#define CELL_DIAG   3
!#define CELL_SD     4
!#define CELL_EXTERN 5
!#define CELL_SE     6


!###############################################################################
MODULE aed2_common
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_sedflux
   USE aed2_land
   USE aed2_ass
!  USE aed2_cladophora
   USE aed2_chlorophylla
   USE aed2_oxygen
   USE ufz_oxygen
   USE aed2_silica
   USE aed2_carbon
   USE aed2_nitrogen
   USE aed2_phosphorus
   USE aed2_macrophyte
   USE aed2_organic_matter
   USE aed2_phytoplankton
   USE aed2_pathogens
   USE aed2_iron
   USE aed2_isotope
   USE aed2_isotope_c
   USE aed2_radon
   USE aed2_sulfur
   USE aed2_zooplankton
   USE aed2_bivalve
   USE aed2_tracer
   USE aed2_geochemistry
   USE aed2_seddiagenesis
   USE aed2_vegetation
   USE aed2_totals
   USE csiro_optical
   USE aed2_test
   USE aed2_habitat

   IMPLICIT NONE

   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   PUBLIC aed2_new_model, aed2_define_model, aed2_build_model

   !#---------------------------------------------------------------------------

   PUBLIC aed2_calculate, aed2_calculate_surface, aed2_calculate_benthic
   PUBLIC aed2_light_extinction, aed2_delete, aed2_equilibrate
   PUBLIC aed2_initialize, aed2_calculate_riparian, aed2_calculate_dry
   PUBLIC aed2_rain_loss, aed2_light_shading, aed2_bio_drag, aed2_particle_bgc

   !# Re-export these from aed2_core.
   PUBLIC aed2_model_data_t, aed2_variable_t, aed2_column_t
   PUBLIC aed2_init_core, aed2_get_var, aed2_core_status

   PUBLIC zero_, one_, nan_, secs_per_day, misval_

   !#---------------------------------------------------------------------------

   CLASS(aed2_model_data_t), POINTER :: model_list => null()
   CLASS(aed2_model_data_t), POINTER :: last_model => null()

CONTAINS
!===============================================================================


!###############################################################################
FUNCTION aed2_new_model(modelname) RESULT(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
   CHARACTER(len=4) :: prefix
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)

   SELECT CASE (modelname)
      CASE ('aed2_sedflux');        prefix = 'SDF'; ALLOCATE(aed2_sedflux_data_t::model)
      CASE ('aed2_land');           prefix = 'LND'; ALLOCATE(aed2_land_data_t::model)
      CASE ('aed2_ass');            prefix = 'ASS'; ALLOCATE(aed2_ass_data_t::model)
!     CASE ('aed2_cladophora');     prefix = 'CLD'; ALLOCATE(aed2_cladophora_data_t::model)
      CASE ('aed2_chlorophylla');   prefix = 'CHL'; ALLOCATE(aed2_chla_data_t::model)
      CASE ('aed2_oxygen');         prefix = 'OXY'; ALLOCATE(aed2_oxygen_data_t::model)
      CASE ('ufz_oxygen');          prefix = 'UOX'; ALLOCATE(ufz_oxygen_data_t::model)
      CASE ('aed2_silica');         prefix = 'SIL'; ALLOCATE(aed2_silica_data_t::model)
      CASE ('aed2_carbon');         prefix = 'CAR'; ALLOCATE(aed2_carbon_data_t::model)
      CASE ('aed2_nitrogen');       prefix = 'NIT'; ALLOCATE(aed2_nitrogen_data_t::model)
      CASE ('aed2_phosphorus');     prefix = 'PHS'; ALLOCATE(aed2_phosphorus_data_t::model)
      CASE ('aed2_organic_matter'); prefix = 'OGM'; ALLOCATE(aed2_organic_matter_data_t::model)
      CASE ('aed2_phytoplankton');  prefix = 'PHY'; ALLOCATE(aed2_phytoplankton_data_t::model)
      CASE ('aed2_zooplankton');    prefix = 'ZOO'; ALLOCATE(aed2_zooplankton_data_t::model)
      CASE ('aed2_bivalve');        prefix = 'BIV'; ALLOCATE(aed2_bivalve_data_t::model)
      CASE ('aed2_macrophyte');     prefix = 'MAC'; ALLOCATE(aed2_macrophyte_data_t::model)
      CASE ('aed2_pathogens');      prefix = 'PAT'; ALLOCATE(aed2_pathogens_data_t::model)
      CASE ('aed2_iron');           prefix = 'IRN'; ALLOCATE(aed2_iron_data_t::model)
      CASE ('aed2_isotope');        prefix = 'ISO'; ALLOCATE(aed2_isotope_data_t::model)
      CASE ('aed2_isotope_c');      prefix = 'ISC'; ALLOCATE(aed2_isotope_c_data_t::model)
      CASE ('aed2_radon');          prefix = 'RAD'; ALLOCATE(aed2_radon_data_t::model)
      CASE ('aed2_sulfur');         prefix = 'SLF'; ALLOCATE(aed2_sulfur_data_t::model)
      CASE ('aed2_tracer');         prefix = 'TRC'; ALLOCATE(aed2_tracer_data_t::model)
      CASE ('aed2_geochemistry');   prefix = 'GEO'; ALLOCATE(aed2_geochemistry_data_t::model)
      CASE ('aed2_seddiagenesis');  prefix = 'SDD'; ALLOCATE(aed2_seddiagenesis_data_t::model)
      CASE ('aed2_vegetation');     prefix = 'VEG'; ALLOCATE(aed2_vegetation_data_t::model)
      CASE ('aed2_totals');         prefix = 'TOT'; ALLOCATE(aed2_totals_data_t::model)
      CASE ('aed2_test');           prefix = 'TST'; ALLOCATE(aed2_test_data_t::model)
      CASE ('csiro_optical');       prefix = 'OPT'; ALLOCATE(csiro_optical_data_t::model)
      CASE ('aed2_habitat');        prefix = 'HAB'; ALLOCATE(aed2_habitat_data_t::model)
      CASE DEFAULT;                 print *,'*** Unknown module ', modelname
   END SELECT

   model%aed2_model_name = modelname
   model%aed2_model_prefix = prefix

   IF ( .NOT. ASSOCIATED(model_list) ) model_list => model
   IF ( ASSOCIATED(last_model) ) last_model%next => model
   last_model => model
END FUNCTION aed2_new_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_build_model(model, namlst, do_prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_model_data_t),POINTER :: model
   INTEGER,INTENT(in)      :: namlst
   LOGICAL,INTENT(in)      :: do_prefix
!
!LOCALS
   CHARACTER(len=4),POINTER :: prefix_p
!
!-------------------------------------------------------------------------------
!BEGIN
   cur_model_name = model%aed2_model_name
   IF ( do_prefix ) THEN
      prefix_p => model%aed2_model_prefix
      CALL aed2_set_prefix(prefix_p)
   ENDIF
   CALL model%define(namlst)
   IF ( do_prefix ) THEN
      prefix_p => null()
      CALL aed2_set_prefix(prefix_p)
   ENDIF
   cur_model_name = ''
END SUBROUTINE aed2_build_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_model(modelname, namlst)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
   INTEGER,INTENT(in)      :: namlst
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
   model => aed2_new_model(modelname)
   CALL aed2_build_model(model, namlst, .TRUE.)
END SUBROUTINE aed2_define_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# These are wrappers for the individual models.                               #
!#                                                                             #
!###############################################################################


!###############################################################################
SUBROUTINE aed2_initialize(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%initialize(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_initialize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_surface(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_surface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_benthic(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian(column, layer_idx, pc_wet)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_riparian(column, layer_idx, pc_wet)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_dry(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_dry(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_dry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_equilibrate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%equilibrate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_equilibrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_validate(column, layer_idx) RESULT(valid)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
   LOGICAL :: valid
!-------------------------------------------------------------------------------
   valid = .TRUE.
   model => model_list
   DO WHILE (ASSOCIATED(model))
      IF (.NOT. model%validate(column, layer_idx)) valid = .FALSE.
      model => model%next
   ENDDO
END FUNCTION aed2_validate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction(column, layer_idx, extinction)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   extinction = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_extinction(column, layer_idx, extinction)
      model => model%next
   ENDDO
END SUBROUTINE aed2_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_rain_loss(column, layer_idx, infil)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   infil = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%rain_loss(column, layer_idx, infil)
      model => model%next
   ENDDO
END SUBROUTINE aed2_rain_loss
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_shading(column, layer_idx, shade_frac)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: shade_frac
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   shade_frac = one_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_shading(column, layer_idx, shade_frac)
      model => model%next
   ENDDO
END SUBROUTINE aed2_light_shading
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_bio_drag(column, layer_idx, drag)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   drag = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%bio_drag(column, layer_idx, drag)
      model => model%next
   ENDDO
END SUBROUTINE aed2_bio_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_particle_bgc(column, layer_idx, ppid)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   ppid = 0
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%particle_bgc(column, layer_idx, ppid)
      model => model%next
   ENDDO
END SUBROUTINE aed2_particle_bgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_delete
!-------------------------------------------------------------------------------
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%delete
      model => model%next
   ENDDO
END SUBROUTINE aed2_delete
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed2_common
