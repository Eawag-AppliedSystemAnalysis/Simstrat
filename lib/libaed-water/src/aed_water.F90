!###############################################################################
!#                                                                             #
!# aed_water.F90                                                               #
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
!# Created May 2020                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

!###############################################################################
MODULE aed_water
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_sedflux
   USE aed_oxygen
   USE aed_silica
   USE aed_carbon
   USE aed_nitrogen
   USE aed_phosphorus
   USE aed_organic_matter
   USE aed_phytoplankton
   USE aed_zooplankton
   USE aed_tracer
   USE aed_noncohesive
   USE aed_totals
   USE aed_dummy
   USE aed_bio_particles
   USE aed_geochemistry
   USE aed_pathogens
   USE aed_pesticides
   USE aed_habitat_water

   USE aed_benthic
   USE aed_riparian
   USE aed_demo
   USE aed_dev

   IMPLICIT NONE


   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   PUBLIC aed_new_model, aed_print_version

   INTEGER,PARAMETER :: NO_ZONES = 1

CONTAINS
!===============================================================================


!###############################################################################
FUNCTION scan_name(modeldef, flags) RESULT(modelname)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in)  :: modeldef
   INTEGER(4), INTENT(out)  :: flags
!
!LOCALS
   INTEGER :: len, i
   CHARACTER(len=64) :: modelname
!
!-------------------------------------------------------------------------------
!BEGIN
   flags = 0
   modelname = ''
   len = LEN_TRIM(modeldef)

   DO i=1,len
      IF (modeldef(i:i) == ':') EXIT
      modelname(i:i) = modeldef(i:i)
   ENDDO

   IF ( i >= len ) RETURN

   DO WHILE ( i <= len )
      IF ( modeldef(i:i+1) == 'nz' ) flags = IOR(flags, NO_ZONES)

      DO WHILE ( i <= len .and. modeldef(i:i) /= ':' ) ; i = i + 1 ; ENDDO
      IF ( i <= len .and. modeldef(i:i) == ':' ) i = i + 1
   ENDDO

   modelname = TRIM(modelname)
END FUNCTION scan_name
!===============================================================================


!###############################################################################
FUNCTION aed_new_model(modeldef) RESULT(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modeldef
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   CHARACTER(len=4) :: prefix
   CHARACTER(len=64) :: modelname
   INTEGER(4) :: flags = 0
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
   modelname = scan_name(modeldef, flags)

   SELECT CASE (modelname)
      CASE ('aed_bio_particles');  prefix = 'PTM'; ALLOCATE(aed_bio_particles_data_t::model)
      CASE ('aed_sedflux');        prefix = 'SDF'; ALLOCATE(aed_sedflux_data_t::model)
      CASE ('aed_oxygen');         prefix = 'OXY'; ALLOCATE(aed_oxygen_data_t::model)
      CASE ('aed_silica');         prefix = 'SIL'; ALLOCATE(aed_silica_data_t::model)
      CASE ('aed_carbon');         prefix = 'CAR'; ALLOCATE(aed_carbon_data_t::model)
      CASE ('aed_nitrogen');       prefix = 'NIT'; ALLOCATE(aed_nitrogen_data_t::model)
      CASE ('aed_phosphorus');     prefix = 'PHS'; ALLOCATE(aed_phosphorus_data_t::model)
      CASE ('aed_organic_matter'); prefix = 'OGM'; ALLOCATE(aed_organic_matter_data_t::model)
      CASE ('aed_phytoplankton');  prefix = 'PHY'; ALLOCATE(aed_phytoplankton_data_t::model)
      CASE ('aed_zooplankton');    prefix = 'ZOO'; ALLOCATE(aed_zooplankton_data_t::model)
      CASE ('aed_tracer');         prefix = 'TRC'; ALLOCATE(aed_tracer_data_t::model)
      CASE ('aed_noncohesive');    prefix = 'NCS'; ALLOCATE(aed_noncohesive_data_t::model)
      CASE ('aed_geochemistry');   prefix = 'GEO'; ALLOCATE(aed_geochemistry_data_t::model)
      CASE ('aed_pathogens');      prefix = 'PAT'; ALLOCATE(aed_pathogens_data_t::model)
      CASE ('aed_pesticides');     prefix = 'PST'; ALLOCATE(aed_pesticides_data_t::model)
      CASE ('aed_totals');         prefix = 'TOT'; ALLOCATE(aed_totals_data_t::model)
      CASE ('aed_dummy');          prefix = 'DUM'; ALLOCATE(aed_dummy_data_t::model)
      CASE ('aed_habitat_water');  prefix = 'HAB'; ALLOCATE(aed_habitat_water_data_t::model)
!     CASE DEFAULT;                print *,'*** Unknown module ', TRIM(modelname)
   END SELECT

   IF (ASSOCIATED(model)) THEN
      model%aed_model_name = modelname
      model%aed_model_prefix = prefix
   ELSE
      model => aed_new_ben_model(modelname)
      IF (.NOT. ASSOCIATED(model)) model => aed_new_rip_model(modelname)
      IF (.NOT. ASSOCIATED(model)) model => aed_new_dmo_model(modelname)
      IF (.NOT. ASSOCIATED(model)) model => aed_new_dev_model(modelname)
   ENDIF

   IF (ASSOCIATED(model)) THEN
      model%aed_model_no_zones = ( IAND(flags, NO_ZONES) /= 0 )

      IF ( .NOT. ASSOCIATED(model_list) ) model_list => model
      IF ( ASSOCIATED(last_model) ) last_model%next => model
      last_model => model
   ELSE
       print *,'*** Unknown module ', TRIM(modelname)
   ENDIF
END FUNCTION aed_new_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_print_version
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-water version ", TRIM(AED_VERSION)
#ifdef __INTEL_COMPILER
   print*,"    libaed built using intel fortran version ", __INTEL_COMPILER
#else
# ifdef __PGI
   print*,"    libaed built using pgfortran version ", __PGIC__, '.', __PGIC_MINOR__, '.', __PGIC_PATCHLEVEL__
# else
#  ifdef __GNUC__
    print*,"    libaed built using gfortran version ", __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
#  else
#   ifdef __clang__
     print*,"    libaed built using flang version ", __clang_major__, '.', __clang_minor__, '.', __clang_patchlevel__
#   else
     print*,"    libaed built using unknown fortran version "
#   endif
#  endif
# endif
#endif

   CALL aed_print_version_ben
   CALL aed_print_version_rip
   CALL aed_print_version_dmo
   CALL aed_print_version_dev
END SUBROUTINE aed_print_version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_water
