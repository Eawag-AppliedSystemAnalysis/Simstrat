!###############################################################################
!#                                                                             #
!# aed2_core.F90                                                               #
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
!# Created June 2013                                                           #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_core

   IMPLICIT NONE

   PRIVATE  !# Everything defaults to private

   PUBLIC aed2_model_data_t, aed2_variable_t, aed2_column_t
   PUBLIC aed2_init_core, aed2_get_var, aed2_core_status, aed2_set_prefix
   PUBLIC aed2_define_variable,        aed2_define_sheet_variable,    &
          aed2_locate_sheet_variable,  aed2_locate_variable,          &
          aed2_define_diag_variable,   aed2_define_sheet_diag_variable
   PUBLIC aed2_locate_global, aed2_locate_global_sheet
   PUBLIC zero_, one_, nan_, misval_, secs_per_day
   PUBLIC cur_model_name


   !#---------------------------------------------------------------------------
   TYPE :: aed2_variable_t
      CHARACTER(len=64) :: name
      CHARACTER(len=64) :: model_name
      CHARACTER(len=128):: longname
      CHARACTER(len=24) :: units
      AED_REAL          :: initial
      AED_REAL          :: minimum
      AED_REAL          :: maximum
      AED_REAL          :: mobility
      AED_REAL          :: light_extinction
      LOGICAL           :: sheet, diag, extern, found
      LOGICAL           :: top, bot
   END TYPE aed2_variable_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !#---------------------------------------------------------------------------
   TYPE :: aed2_model_data_t
      INTEGER           :: aed2_model_id
      CHARACTER(len=64) :: aed2_model_name
      CHARACTER(len=4)  :: aed2_model_prefix
      CLASS(aed2_model_data_t),POINTER :: next => null()
      !# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONTAINS
         procedure :: define             => aed2_define
         procedure :: calculate_surface  => aed2_calculate_surface
         procedure :: calculate          => aed2_calculate
         procedure :: calculate_benthic  => aed2_calculate_benthic
         procedure :: calculate_riparian => aed2_calculate_riparian
         procedure :: calculate_dry      => aed2_calculate_dry
         procedure :: equilibrate        => aed2_equilibrate
         procedure :: light_extinction   => aed2_light_extinction
         procedure :: mobility           => aed2_mobility
         procedure :: validate           => aed2_validate
         procedure :: delete             => aed2_delete
      !# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   END TYPE aed2_model_data_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !#---------------------------------------------------------------------------
   TYPE :: aed2_column_t
      AED_REAL,DIMENSION(:),POINTER :: cell
      AED_REAL,             POINTER :: cell_sheet
      AED_REAL,             POINTER :: flux_atm
      AED_REAL,DIMENSION(:),POINTER :: flux_pel
      AED_REAL,             POINTER :: flux_ben
   END TYPE aed2_column_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!-------------------------------------------------------------------------------
!MODULE VARIABLES
   INTEGER :: cur_mod_base = 0
   INTEGER :: n_aed_vars, a_vars
   INTEGER :: n_vars, n_sheet_vars
   INTEGER :: n_diags, n_sheet_diags

   TYPE(aed2_variable_t),DIMENSION(:),ALLOCATABLE,TARGET :: all_vars
!  CLASS(aed2_model_data_t), POINTER :: aed2_cur_mod => null()
   CHARACTER(len=4),POINTER :: aed2_cur_prefix => null()
   CHARACTER(len=80) :: base_aed_directory = '.'
   CHARACTER(len=64) :: cur_model_name = ''

!CONSTANTS
   AED_REAL,PARAMETER :: zero_ = 0., one_ = 1.
   AED_REAL,PARAMETER :: secs_per_day = 86400.
   AED_REAL,PARAMETER :: misval_ = -9999.

!# We need a NaN for initialisation purposes, but gfortran wont compile (it sees
!# we are creating a NaN and stops) so for gfortran we use a regular variable
!# rather than constant and start with one, then in init divide it by 0.
!# [There may be a flag that turns this into a warning rather than an error, but
!#  so far I haven't found it]
!# The preprocessor symbols predefined by compilers are :
!#   __GFORTRAN__      for gfortran
!#   __INTEL_COMPILER  for intel fortran (ifort)
!#ifdef __GFORTRAN__
!   AED_REAL           :: nan_ = zero_
!#else
   AED_REAL,PARAMETER :: nan_ = zero_ / zero_
!#endif

!===============================================================================
CONTAINS


!###############################################################################
INTEGER FUNCTION aed2_init_core(dname)
!-------------------------------------------------------------------------------
! Initialise the aed model library core routines
!-------------------------------------------------------------------------------
!ARGUMENTS - none
   CHARACTER(len=*) :: dname
!
!LOCAL
!   AED_REAL :: tmpr = zero_
!
!-------------------------------------------------------------------------------
!BEGIN
   base_aed_directory = dname
   n_aed_vars = 0 ; a_vars = 0
   n_vars = 0;  n_sheet_vars = 0
   n_diags = 0; n_sheet_diags = 0
!#ifdef __GFORTRAN__
!   nan_ = nan_ / tmpr;
!#endif
   print*,"libaed2 version ", TRIM(AED2_VERSION)
   aed2_init_core = 0
END FUNCTION aed2_init_core
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION aed2_core_status(n_v, n_sv, n_d, n_sd)
!-------------------------------------------------------------------------------
! Status of the aed model library core routines
!-------------------------------------------------------------------------------
!ARGUMENTS
     INTEGER,INTENT(out) :: n_v, n_sv, n_d, n_sd
!
!-------------------------------------------------------------------------------
!BEGIN
   n_v = n_vars;  n_sv = n_sheet_vars
   n_d = n_diags; n_sd = n_sheet_diags
   aed2_core_status = n_aed_vars
END FUNCTION aed2_core_status
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_set_prefix(prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=4),POINTER :: prefix
!
!-------------------------------------------------------------------------------
!BEGIN
   aed2_cur_prefix => prefix
END SUBROUTINE aed2_set_prefix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!# These create space for the variables :

!###############################################################################
SUBROUTINE extend_allocated_variables(count)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: count
!
!LOCALS
   TYPE(aed2_variable_t),DIMENSION(:),ALLOCATABLE :: tmp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( ALLOCATED(all_vars) ) THEN
      ALLOCATE(tmp(1:a_vars))
      tmp = all_vars
      DEALLOCATE(all_vars)
      ALLOCATE(all_vars(1:a_vars+count))
      all_vars(1:a_vars) = tmp(1:a_vars)
      DEALLOCATE(tmp);
   ELSE
      ALLOCATE(all_vars(1:count))
   ENDIF

   all_vars(a_vars+1:a_vars+count)%initial = nan_
   all_vars(a_vars+1:a_vars+count)%minimum = nan_
   all_vars(a_vars+1:a_vars+count)%maximum = nan_
   all_vars(a_vars+1:a_vars+count)%mobility = nan_
   all_vars(a_vars+1:a_vars+count)%light_extinction = nan_

   all_vars(a_vars+1:a_vars+count)%sheet = .false.
   all_vars(a_vars+1:a_vars+count)%diag = .false.
   all_vars(a_vars+1:a_vars+count)%extern = .false.
   all_vars(a_vars+1:a_vars+count)%found = .false.

   all_vars(a_vars+1:a_vars+count)%top = .false.
   all_vars(a_vars+1:a_vars+count)%bot = .false.

!  IF ( ALLOCATED(column) ) THEN
!     ALLOCATE(tmpc(1:a_vars))
!     tmpc = column
!     DEALLOCATE(column)
!     ALLOCATE(column(1:a_vars+count))
!     column(1:a_vars) = tmpc(1:a_vars)
!     DEALLOCATE(tmpc);
!  ELSE
!     ALLOCATE(column(1:count))
!  ENDIF

   a_vars = a_vars + count

END SUBROUTINE extend_allocated_variables
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!# These define names and defaults for variables, returning their index:

!###############################################################################
SUBROUTINE aed2_create_variable(name, longname, units)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
!
!-------------------------------------------------------------------------------
!BEGIN
!print*,"Variable '",trim(name),"' long name '",trim(longname),"' with units '",trim(units),"'"
   IF (n_aed_vars .GE. a_vars) CALL extend_allocated_variables(n_aed_vars+1-a_vars)
   n_aed_vars = n_aed_vars + 1
   IF ( ASSOCIATED(aed2_cur_prefix) .AND. units /= '' ) THEN
      all_vars(n_aed_vars)%name = TRIM(aed2_cur_prefix)//"_"//name
   ELSE
      all_vars(n_aed_vars)%name = name
   ENDIF
   all_vars(n_aed_vars)%model_name = cur_model_name
   all_vars(n_aed_vars)%longname = longname
   all_vars(n_aed_vars)%units = units

   all_vars(n_aed_vars)%sheet = .false.
   all_vars(n_aed_vars)%diag = .false.
   all_vars(n_aed_vars)%extern = .false.
   all_vars(n_aed_vars)%found = .false.
END SUBROUTINE aed2_create_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_define_variable(name, units, longname, initial, minimum, maximum, mobility) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   AED_REAL,INTENT(in),OPTIONAL :: initial, minimum, maximum
   AED_REAL,INTENT(in),OPTIONAL :: mobility
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed2_create_variable(name, longname, units)
   n_vars = n_vars + 1

   if ( present(initial) ) all_vars(n_aed_vars)%initial = initial
   if ( present(minimum) ) all_vars(n_aed_vars)%minimum = minimum
   if ( present(maximum) ) all_vars(n_aed_vars)%maximum = maximum
   if ( present(mobility) ) all_vars(n_aed_vars)%mobility = mobility

!  all_vars(n_aed_vars)%sheet = .false.
!  all_vars(n_aed_vars)%diag = .false.
!  all_vars(n_aed_vars)%extern = .false.
   all_vars(n_aed_vars)%found = .true.

   ret = n_aed_vars
END FUNCTION aed2_define_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_define_sheet_variable(name, units, longname, initial, minimum, maximum, surf) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   AED_REAL,INTENT(in),OPTIONAL :: initial, minimum, maximum
   LOGICAL,INTENT(in),OPTIONAL :: surf
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed2_create_variable(name, longname, units)
   n_sheet_vars = n_sheet_vars + 1

   if ( present(initial) ) all_vars(n_aed_vars)%initial = initial
   if ( present(minimum) ) all_vars(n_aed_vars)%minimum = minimum
   if ( present(maximum) ) all_vars(n_aed_vars)%maximum = maximum

   all_vars(n_aed_vars)%sheet = .true.
!  all_vars(n_aed_vars)%diag = .false.
!  all_vars(n_aed_vars)%extern = .false.
   all_vars(n_aed_vars)%found = .true.

   IF ( PRESENT(surf) ) THEN
      all_vars(n_aed_vars)%bot = .NOT. surf
      all_vars(n_aed_vars)%top = surf
   ELSE
      all_vars(n_aed_vars)%bot = .TRUE.
   ENDIF

   ret = n_aed_vars
END FUNCTION aed2_define_sheet_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_define_diag_variable(name, units, longname) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed2_create_variable(name, longname, units)
   n_diags = n_diags + 1

!  all_vars(n_aed_vars)%sheet = .false.
   all_vars(n_aed_vars)%diag = .true.
!  all_vars(n_aed_vars)%extern = .false.
   all_vars(n_aed_vars)%found = .true.

   ret = n_aed_vars
END FUNCTION aed2_define_diag_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_define_sheet_diag_variable(name, units, longname, surf) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   LOGICAL,INTENT(in),OPTIONAL :: surf
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed2_create_variable(name, longname, units)
   n_sheet_diags = n_sheet_diags + 1

   all_vars(n_aed_vars)%sheet = .true.
   all_vars(n_aed_vars)%diag = .true.
!  all_vars(n_aed_vars)%extern = .false.
   all_vars(n_aed_vars)%found = .true.

   IF ( PRESENT(surf) ) THEN
      all_vars(n_aed_vars)%bot = .NOT. surf
      all_vars(n_aed_vars)%top = surf
   ELSE
      all_vars(n_aed_vars)%bot = .TRUE.
   ENDIF

   ret = n_aed_vars
END FUNCTION aed2_define_sheet_diag_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!# These return the index for the named variable :

!###############################################################################
FUNCTION aed2_find_variable(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   DO ret=1,n_aed_vars
      IF ( all_vars(ret)%name == name ) RETURN
      IF ( trim(all_vars(ret)%model_name)//'_'//trim(all_vars(ret)%name) == name ) RETURN
!print*,trim(name), ' is neither ', trim(all_vars(ret)%name), ' nor ',  trim(all_vars(ret)%model_name)//'_'//trim(all_vars(ret)%name)
   ENDDO
   ret = 0
END FUNCTION aed2_find_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_locate_variable(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed2_find_variable(name)
   IF ( ret == 0 ) THEN
!print*,'Variable ',trim(name),' not found'
!stop
      CALL aed2_create_variable(name, '', '')
!     all_vars(n_aed_vars)%sheet = .false.
!     all_vars(n_aed_vars)%diag = .false.
!     all_vars(n_aed_vars)%extern = .true.

      ret = n_aed_vars
      RETURN
   ENDIF
   IF ( all_vars(ret)%sheet .OR. all_vars(ret)%diag ) ret = 0
END FUNCTION aed2_locate_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_locate_sheet_variable(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed2_find_variable(name)
   IF ( ret == 0 ) THEN
      CALL aed2_create_variable(name, '', '')
      all_vars(n_aed_vars)%name = name
      all_vars(n_aed_vars)%sheet = .true.
!     all_vars(n_aed_vars)%diag = .false.
!     all_vars(n_aed_vars)%extern = .true.

      ret = n_aed_vars
      RETURN
   ENDIF
   IF ( (.NOT. all_vars(ret)%sheet) .OR. all_vars(ret)%diag ) ret = 0
END FUNCTION aed2_locate_sheet_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_locate_global(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed2_find_variable(name)
   IF ( ret == 0 ) THEN
!     print *,"variable ",trim(name)," not found"
      CALL aed2_create_variable(name, '', '')
!     all_vars(n_aed_vars)%sheet = .false.
!     all_vars(n_aed_vars)%diag = .false.
      all_vars(n_aed_vars)%extern = .true.

      ret = n_aed_vars
      RETURN
   ENDIF
   IF ( all_vars(ret)%sheet ) ret = 0
END FUNCTION aed2_locate_global
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_locate_global_sheet(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed2_find_variable(name)
   IF ( ret == 0 ) THEN
!     print *,"sheet var ",trim(name)," not found"
      CALL aed2_create_variable(name, '', '')
      all_vars(n_aed_vars)%sheet = .true.
!     all_vars(n_aed_vars)%diag = .false.
      all_vars(n_aed_vars)%extern = .true.

      ret = n_aed_vars
      RETURN
   ENDIF
   IF ( .NOT. all_vars(ret)%sheet ) ret = 0
END FUNCTION aed2_locate_global_sheet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed2_get_var(which, tvar)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: which
   TYPE(aed2_variable_t),POINTER,INTENT(out) :: tvar
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   aed2_get_var = .TRUE. ; NULLIFY(tvar)
   IF (which > 0 .AND. which <= n_aed_vars) THEN ; tvar => all_vars(which)
   ELSE ; aed2_get_var = .FALSE. ; ENDIF
END FUNCTION aed2_get_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# These are dummies for the declarations.                                     #
!#                                                                             #
!###############################################################################


!###############################################################################
SUBROUTINE aed2_define(data, namlst)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!-------------------------------------------------------------------------------
!print*,"Default aed2_define ", trim(data%aed2_model_name)
END SUBROUTINE aed2_define
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed2_calculate ", trim(data%aed2_model_name)
END SUBROUTINE aed2_calculate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed2_calculate_surface ", trim(data%aed2_model_name)
END SUBROUTINE aed2_calculate_surface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed2_calculate_benthic ", trim(data%aed2_model_name)
END SUBROUTINE aed2_calculate_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian(data,column,layer_idx, pc_wet)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!-------------------------------------------------------------------------------
!print*,"Default aed2_calculate_riparian ", trim(data%aed2_model_name)
END SUBROUTINE aed2_calculate_riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_dry(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed2_calculate_dry ", trim(data%aed2_model_name)
END SUBROUTINE aed2_calculate_dry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_equilibrate(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed2_equilibrate ", trim(data%aed2_model_name)
END SUBROUTINE aed2_equilibrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed2_validate(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed2_validate ", trim(data%aed2_model_name)
   aed2_validate = .TRUE.
END FUNCTION aed2_validate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_mobility(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!-------------------------------------------------------------------------------
!print*,"Default aed2_mobility ", trim(data%aed2_model_name)
END SUBROUTINE aed2_mobility
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!-------------------------------------------------------------------------------
!print*,"Default aed2_light_extinction ", trim(data%aed2_model_name)
END SUBROUTINE aed2_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_delete(data)
!-------------------------------------------------------------------------------
   CLASS (aed2_model_data_t),INTENT(inout) :: data
!-------------------------------------------------------------------------------
!print*,"Default aed2_delete ", trim(data%aed2_model_name)
END SUBROUTINE aed2_delete
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_core
