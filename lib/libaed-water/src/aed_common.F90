!###############################################################################
!#                                                                             #
!# aed_common.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2021 -  The University of Western Australia               #
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
!# Created Aug 2013                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"


!###############################################################################
MODULE aed_common
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_water

   IMPLICIT NONE

   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   PUBLIC aed_new_model, aed_define_model, aed_build_model

   !#---------------------------------------------------------------------------

   PUBLIC aed_initialize, aed_initialize_benthic
   PUBLIC aed_calculate, aed_calculate_surface, aed_calculate_benthic
   PUBLIC aed_calculate_riparian, aed_calculate_dry
   PUBLIC aed_light_extinction, aed_delete, aed_equilibrate
   PUBLIC aed_mobility, aed_rain_loss, aed_light_shading
   PUBLIC aed_bio_drag, aed_particle_bgc

   !# Re-export these from aed_core.
   PUBLIC aed_model_data_t, aed_variable_t, aed_column_t
   PUBLIC aed_init_core, aed_get_var, aed_core_status
   PUBLIC aed_provide_global, aed_provide_sheet_global

   PUBLIC zero_, one_, nan_, secs_per_day, misval_


CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed_build_model(model, namlst, do_prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_model_data_t),POINTER :: model
   INTEGER,INTENT(in)      :: namlst
   LOGICAL,INTENT(in)      :: do_prefix
!
!LOCALS
   CHARACTER(len=4),POINTER :: prefix_p
!
!-------------------------------------------------------------------------------
!BEGIN
   cur_model_name = model%aed_model_name
   IF ( do_prefix ) THEN
      prefix_p => model%aed_model_prefix
      CALL aed_set_prefix(prefix_p)
   ENDIF
   CALL model%define(namlst)
   IF ( do_prefix ) THEN
      prefix_p => null()
      CALL aed_set_prefix(prefix_p)
   ENDIF
   cur_model_name = ''
END SUBROUTINE aed_build_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_model(modelname, namlst)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
   INTEGER,INTENT(in)      :: namlst
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
   model => aed_new_model(modelname)
   IF ( ASSOCIATED(model) ) CALL aed_build_model(model, namlst, .TRUE.)
END SUBROUTINE aed_define_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# These are wrappers for the individual models.                               #
!#                                                                             #
!###############################################################################


!###############################################################################
SUBROUTINE aed_initialize(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%initialize(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_initialize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%initialize_benthic(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_initialize_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_surface(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate_surface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic(column, layer_idx, do_zones)
!-------------------------------------------------------------------------------
! The benthic routine may be grouped in zones by the global do_zone_averaging
! flag, however a new model level flag allows us to not average in zones.
! This routine takes the optional argument "do_zones" which should only be
! passed if do_zone_averaging is on.
! If do_zones is not present we can call every models calculate_benthic
! routine.
! If it IS present we pass true when called from inside the zone calculations,
! but false when called from the normal flux calculation section.
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   LOGICAL,OPTIONAL,INTENT(in) :: do_zones
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   IF ( PRESENT(do_zones) ) THEN
      DO WHILE (ASSOCIATED(model))
         IF ( (.NOT. model%aed_model_no_zones) .EQV. do_zones ) &
            CALL model%calculate_benthic(column, layer_idx)
         model => model%next
      ENDDO
   ELSE
      DO WHILE (ASSOCIATED(model))
         CALL model%calculate_benthic(column, layer_idx)
         model => model%next
      ENDDO
   ENDIF
END SUBROUTINE aed_calculate_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian(column, layer_idx, pc_wet)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_riparian(column, layer_idx, pc_wet)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate_riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_dry(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_dry(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate_dry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%equilibrate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_equilibrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_validate(column, layer_idx) RESULT(valid)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   LOGICAL :: valid
!-------------------------------------------------------------------------------
   valid = .TRUE.
   model => model_list
   DO WHILE (ASSOCIATED(model))
      IF (.NOT. model%validate(column, layer_idx)) valid = .FALSE.
      model => model%next
   ENDDO
END FUNCTION aed_validate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction(column, layer_idx, extinction)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   extinction = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_extinction(column, layer_idx, extinction)
      model => model%next
   ENDDO
END SUBROUTINE aed_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility(column, layer_idx, mobility)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   !mobility = zero_ !MH leave this as is in case default settling vals provided
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%mobility(column, layer_idx, mobility)
      model => model%next
   ENDDO
END SUBROUTINE aed_mobility
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_rain_loss(column, layer_idx, infil)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   infil = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%rain_loss(column, layer_idx, infil)
      model => model%next
   ENDDO
END SUBROUTINE aed_rain_loss
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_shading(column, layer_idx, shade_frac)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: shade_frac
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   shade_frac = one_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_shading(column, layer_idx, shade_frac)
      model => model%next
   ENDDO
END SUBROUTINE aed_light_shading
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_drag(column, layer_idx, drag)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   drag = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%bio_drag(column, layer_idx, drag)
      model => model%next
   ENDDO
END SUBROUTINE aed_bio_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc(column, layer_idx, ppid, partcl)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   AED_REAL,DIMENSION(:),INTENT(inout) :: partcl
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%particle_bgc(column, layer_idx, ppid, partcl)
      model => model%next
   ENDDO
END SUBROUTINE aed_particle_bgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_delete
!-------------------------------------------------------------------------------
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%delete
      model => model%next
   ENDDO
END SUBROUTINE aed_delete
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed_common
