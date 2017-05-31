!###############################################################################
!#                                                                             #
!# csiro_optical.F90                                                           #
!#                                                                             #
!# Developed by :                                                              #
!#                                                                             #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created Feb 2016                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE csiro_optical
!-------------------------------------------------------------------------------
! csiro_optical --- optical model
!
! The module .......
!
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC csiro_optical_data_t
!
   TYPE,extends(aed2_model_data_t) :: csiro_optical_data_t
      !# Variable identifiers
      INTEGER  :: id_opt_tn
      INTEGER  :: num_tn
      INTEGER  :: id_par, id_extc
      INTEGER,ALLOCATABLE :: id_dep_tn(:)
      AED_REAL,ALLOCATABLE :: tn_varscale(:)

      !# Model parameters
      INTEGER  :: num_tracers
      AED_REAL,ALLOCATABLE :: decay(:)
      INTEGER,ALLOCATABLE :: id_ss(:)

      CONTAINS
         PROCEDURE :: define            => aed2_define_optical
         PROCEDURE :: calculate         => aed2_calculate_optical
!        PROCEDURE :: calculate_benthic => aed2_calculate_benthic_optical
!        PROCEDURE :: mobility          => aed2_mobility_optical
!        PROCEDURE :: light_extinction  => aed2_light_extinction_optical
!        PROCEDURE :: delete            => aed2_delete_optical

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_optical(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the model
!
!  Here, the namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (csiro_optical_data_t),INTENT(inout) :: data
!
!LOCALS
!
   INTEGER           :: status
   INTEGER           :: i, num_tn, num_tracers
   CHARACTER(len=40) :: tn_vars(100)
   AED_REAL          :: tn_varscale(100), decay(100)
   AED_REAL          :: trace_initial = zero_
   CHARACTER(4)      :: trac_name

   NAMELIST /csiro_optical/ tn_vars,tn_varscale, &
                            num_tracers,decay
!
!-------------------------------------------------------------------------------
!BEGIN
   tn_vars = ''
   tn_varscale = 1.0


   ! TOTAL VARIABLE

   ! Read the namelist
   read(namlst,nml=csiro_optical,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist csiro_optical'


   DO i=1,100 ; IF (tn_vars(i)  .EQ. '' ) THEN ; num_tn  = i-1 ; EXIT ; ENDIF ; ENDDO

   ALLOCATE(data%id_dep_tn(num_tn))     ; ALLOCATE(data%tn_varscale(num_tn))

   data%num_tn   = num_tn

   ! Register external state variable dependencies
   DO i=1,data%num_tn
      data%id_dep_tn(i) =  aed2_locate_variable(tn_vars(i))
      data%tn_varscale(i) =  tn_varscale(i)
      print*,'TN : ', TRIM(tn_vars(i)), ' * ', data%tn_varscale(i)
   ENDDO

   ! Register environmental dependencies
   data%id_par = aed2_locate_global('par')
   data%id_extc = aed2_locate_global('extc_coef')

   ! Register diagnostic variables
   IF (data%num_tn>0) THEN
      data%id_opt_tn = aed2_define_diag_variable('tn',               &
                     'mmol/m**3', 'Total Nitrogen')
   ENDIF



   ! TRACER STATE VARIABLE

   data%num_tracers = num_tracers
   IF ( num_tracers > 0 ) THEN
      ALLOCATE(data%id_ss(num_tracers))
      ALLOCATE(data%decay(num_tracers))    ; data%decay(1:num_tracers)    = decay(1:num_tracers)

      trac_name = 'tr0'
      ! Register state variables
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_ss(i) = aed2_define_variable(TRIM(trac_name),'mmol/m**3','tracer', &
                                                  trace_initial,minimum=zero_,mobility=zero_)
      ENDDO
   ENDIF
END SUBROUTINE aed2_define_optical
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_optical(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of csiro_optical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (csiro_optical_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: i,count
   AED_REAL :: val, tot,trc

!-------------------------------------------------------------------------------
!BEGIN

   ! Sum over TN variables provided in the NML and export them to a diagnostic
   IF (data%num_tn>0) THEN
      tot = 0.
      count = ubound(data%id_dep_tn,1)
      DO i=1,count
         val = _STATE_VAR_(data%id_dep_tn(i))
         tot = tot + (val*data%tn_varscale(i))
      ENDDO

      _DIAG_VAR_(data%id_opt_tn) =  tot
   ENDIF

   ! Now compute decau of tarcer
   DO i=1,data%num_tracers
      trc = _STATE_VAR_(data%id_ss(i))
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + data%decay(i)*trc
   ENDDO
END SUBROUTINE aed2_calculate_optical
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE csiro_optical
