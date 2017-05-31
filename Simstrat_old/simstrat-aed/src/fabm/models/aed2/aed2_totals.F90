!###############################################################################
!#                                                                             #
!# aed2_totals.F90                                                             #
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
!# Created July 2012                                                           #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_totals
!-------------------------------------------------------------------------------
! aed2_totals --- totals model
!
! The AED module "totals" contains only diagnostic variables to provide
! totals of other variables (eg tss) or derived metrics of interest (eg turbidity)
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_totals_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_totals_data_t
      !# Variable identifiers
      INTEGER  :: id_totals_tn, id_totals_tkn, id_totals_tp,                   &
                  id_totals_toc, id_totals_tss, id_totals_turbidity,           &
                  id_totals_light, id_totals_par, id_totals_uv, id_totals_extc
      INTEGER  :: num_tn, num_tkn, num_tp, num_toc, num_tss, num_turb
      INTEGER  :: id_par, id_extc
      INTEGER,ALLOCATABLE :: id_dep_tn(:), id_dep_tkn(:), id_dep_tp(:),        &
                             id_dep_toc(:), id_dep_tss(:), id_dep_turb(:)
      AED_REAL,ALLOCATABLE :: tn_varscale(:), tkn_varscale(:), tp_varscale(:), &
                              toc_varscale(:), tss_varscale(:), turb_varscale(:)


      !# Model parameters
      LOGICAL  :: outputLight

     CONTAINS
         PROCEDURE :: define            => aed2_define_totals
         PROCEDURE :: calculate         => aed2_calculate_totals
!        PROCEDURE :: mobility          => aed2_mobility_totals
!        PROCEDURE :: light_extinction  => aed2_light_extinction_totals
!        PROCEDURE :: delete            => aed2_delete_totals

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_totals(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_totals_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status

   INTEGER           :: i, num_tn,num_tkn,num_tp,num_toc,num_tss,num_turb
   CHARACTER(len=40) :: tn_vars(100), tkn_vars(100), tp_vars(100)
   CHARACTER(len=40) :: toc_vars(100), tss_vars(100), turb_vars(100)
   AED_REAL          :: tn_varscale(100), tkn_varscale(100), tp_varscale(100)
   AED_REAL          :: toc_varscale(100), tss_varscale(100), turb_varscale(100)
   LOGICAL           :: outputLight = .FALSE.

   NAMELIST /aed2_totals/ tn_vars,  tn_varscale,  tkn_vars,  tkn_varscale,  &
                          tp_vars,  tp_varscale,  toc_vars, toc_varscale,   &
                          tss_vars, tss_varscale, turb_vars, turb_varscale, &
                          outputLight
!
!-------------------------------------------------------------------------------
!BEGIN
   tn_vars = '' ;       tkn_vars = '' ;      tp_vars = ''
   toc_vars = '' ;      tss_vars = '' ;      turb_vars = ''
   tn_varscale = 1.0 ;  tkn_varscale = 1.0;  tp_varscale = 1.0
   toc_varscale = 1.0 ; tss_varscale = 1.0;  turb_varscale = 1.0

   ! Read the namelist
   read(namlst,nml=aed2_totals,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_totals'

   data%outputLight = outputLight

   DO i=1,100 ; IF (tn_vars(i)  .EQ. '' ) THEN ; num_tn  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tkn_vars(i) .EQ. '' ) THEN ; num_tkn = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tp_vars(i)  .EQ. '' ) THEN ; num_tp  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (toc_vars(i) .EQ. '' ) THEN ; num_toc = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tss_vars(i) .EQ. '' ) THEN ; num_tss = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (turb_vars(i).EQ. '' ) THEN ; num_turb= i-1 ; EXIT ; ENDIF ; ENDDO

   ALLOCATE(data%id_dep_tn(num_tn))     ; ALLOCATE(data%tn_varscale(num_tn))
   ALLOCATE(data%id_dep_tkn(num_tkn))   ; ALLOCATE(data%tkn_varscale(num_tkn))
   ALLOCATE(data%id_dep_tp(num_tp))     ; ALLOCATE(data%tp_varscale(num_tp))
   ALLOCATE(data%id_dep_toc(num_toc))   ; ALLOCATE(data%toc_varscale(num_toc))
   ALLOCATE(data%id_dep_tss(num_tss))   ; ALLOCATE(data%tss_varscale(num_tss))
   ALLOCATE(data%id_dep_turb(num_turb)) ; ALLOCATE(data%turb_varscale(num_turb))

   data%num_tn   = num_tn
   data%num_tkn  = num_tkn
   data%num_tp   = num_tp
   data%num_toc  = num_toc
   data%num_tss  = num_tss
   data%num_turb = num_turb

   ! Register external state variable dependencies
   DO i=1,data%num_tn
      data%id_dep_tn(i) =  aed2_locate_variable(tn_vars(i))
      data%tn_varscale(i) =  tn_varscale(i)
     print*,'TN : ', TRIM(tn_vars(i)), ' * ', data%tn_varscale(i)
   ENDDO
   DO i=1,data%num_tkn
      data%id_dep_tkn(i) =  aed2_locate_variable(tkn_vars(i))
      data%tkn_varscale(i) =  tkn_varscale(i)
     print*,'TKN : ', TRIM(tkn_vars(i)), ' * ', data%tkn_varscale(i)
   ENDDO
   DO i=1,data%num_tp
      data%id_dep_tp(i) =  aed2_locate_variable(tp_vars(i))
      data%tp_varscale(i) =  tp_varscale(i)
     print*,'TP : ', TRIM(tp_vars(i)), ' * ', data%tp_varscale(i)
   ENDDO
   DO i=1,data%num_toc
      data%id_dep_toc(i) = aed2_locate_variable(toc_vars(i))
      data%toc_varscale(i) = toc_varscale(i)
     print*,'TOC : ', TRIM(toc_vars(i)), ' * ', data%toc_varscale(i)
   ENDDO
   DO i=1,data%num_tss
      data%id_dep_tss(i) = aed2_locate_variable(tss_vars(i))
      data%tss_varscale(i) = tss_varscale(i)
     print*,'TSS : ', TRIM(tss_vars(i)), ' * ', data%tss_varscale(i)
   ENDDO
   DO i=1,data%num_turb
      data%id_dep_turb(i) = aed2_locate_variable(turb_vars(i))
      data%turb_varscale(i) = turb_varscale(i)
     print*,'TURB : ', TRIM(turb_vars(i)), ' * ', data%turb_varscale(i)
   ENDDO
   ! Register environmental dependencies
   IF (data%outputLight) THEN
     data%id_par = aed2_locate_global('par')
     data%id_extc = aed2_locate_global('extc_coef')
   END IF

   ! Register diagnostic variables
   IF (data%num_tn>0) THEN
     data%id_totals_tn = aed2_define_diag_variable('tn',               &
                     'mmol/m**3', 'Total Nitrogen')
   ENDIF

   IF (data%num_tkn>0) THEN
     data%id_totals_tkn = aed2_define_diag_variable('tkn',             &
                     'mmol/m**3', 'Total Kjedhal Nitrogen')
   ENDIF

   IF (data%num_tp>0) THEN
     data%id_totals_tp = aed2_define_diag_variable('tp',               &
                     'mmol/m**3', 'Total Phosphorus')
   ENDIF

   IF (data%num_toc>0) THEN
     data%id_totals_toc = aed2_define_diag_variable('toc',             &
                     'mmol/m**3', 'Total Organic Carbon')
   ENDIF

   IF (data%num_tss>0) THEN
     data%id_totals_tss = aed2_define_diag_variable('tss',             &
                     'mg/L', 'Total Suspended Solids')
   ENDIF

   IF (data%num_turb>0) THEN
     data%id_totals_turbidity = aed2_define_diag_variable('turbidity', &
                     'NTU', 'Turbidity')
   ENDIF

   IF (data%outputLight) THEN
     data%id_totals_light = aed2_define_diag_variable('light',         &
                     'W/m2', 'Shortwave light flux')

     data%id_totals_par = aed2_define_diag_variable('par',             &
                     'W/m2', 'Photosynthetically active light flux')

     data%id_totals_uv = aed2_define_diag_variable('uv',               &
                     'W/m2', 'Ultraviolet light flux')

     data%id_totals_extc = aed2_define_diag_variable('extc',           &
                     'W/m2', 'Light extinction coefficient')
   ENDIF

END SUBROUTINE aed2_define_totals
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_totals(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_totals model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_totals_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: i,count
   AED_REAL :: val, tot

!-------------------------------------------------------------------------------
!BEGIN

   ! Sum over TN variables
   IF (data%num_tn>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tn,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tn(i)); tot = tot + (val*data%tn_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tn) =  tot
   ENDIF

   ! Sum over TKN variables
   IF (data%num_tkn>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tkn,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tkn(i)); tot = tot + (val*data%tkn_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tkn) =  tot
   ENDIF

   ! Sum over TP variables
   IF (data%num_tp>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tp,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tp(i)); tot = tot + (val*data%tp_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tp) =  tot
   ENDIF

   ! Sum over TOC variables
   IF (data%num_toc>0) THEN
     tot = 0.
     count = ubound(data%id_dep_toc,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_toc(i)); tot = tot + (val*data%toc_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_toc) =  tot
   ENDIF

   ! Sum over TSS variables
   IF (data%num_tss>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tss,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tss(i)); tot = tot + (val*data%tss_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tss) =  tot
   ENDIF

   ! Sum over Turbidity variables
   IF (data%num_turb>0) THEN
     tot = 0.
     count = ubound(data%id_dep_turb,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_turb(i)); tot = tot + (val*data%turb_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_turbidity) =  tot
   ENDIF

   ! Light and Extinction Coefficient
   IF (data%outputLight) THEN
     val = _STATE_VAR_(data%id_par)
     _DIAG_VAR_(data%id_totals_light) =  val/0.45
     _DIAG_VAR_(data%id_totals_par) =    val
     _DIAG_VAR_(data%id_totals_uv) =     (val/0.45)*0.03
     val = _STATE_VAR_(data%id_extc)
     _DIAG_VAR_(data%id_totals_extc) =   val
   END IF

END SUBROUTINE aed2_calculate_totals
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_totals
