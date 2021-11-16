!###############################################################################
!#                                                                             #
!# aed_totals.F90                                                              #
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
!# Created July 2012                                                           #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |  _________   | || |     ____     | || |  _________   | |         !
!         | | |  _   _  |  | || |   .'    `.   | || | |  _   _  |  | |         !
!         | | |_/ | | \_|  | || |  /  .--.  \  | || | |_/ | | \_|  | |         !
!         | |     | |      | || |  | |    | |  | || |     | |      | |         !
!         | |    _| |_     | || |  \  `--'  /  | || |    _| |_     | |         !
!         | |   |_____|    | || |   `.____.'   | || |   |_____|    | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_totals
!-------------------------------------------------------------------------------
! aed_totals --- totals model
!
! The AED module "totals" contains only diagnostic variables to provide totals
! of other variables (eg tn, tss) or derived metrics of interest (eg turbidity)
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_totals_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_totals_data_t
      !# Variable identifiers
      INTEGER  :: id_totals_tn, id_totals_tkn, id_totals_tp, id_tfe, id_tal,      &
                  id_totals_toc, id_totals_tss, id_totals_turbidity,              &
                  id_totals_light, id_totals_par, id_totals_uv, id_totals_extc
      INTEGER  :: num_tn, num_tkn, num_tp, num_toc, num_tss, num_turb, num_tfe, num_tal
      INTEGER  :: id_par, id_extc
      INTEGER,ALLOCATABLE :: id_dep_tn(:), id_dep_tkn(:), id_dep_tp(:),           &
                             id_dep_toc(:), id_dep_tss(:), id_dep_turb(:),        &
                             id_dep_tfe(:), id_dep_tal(:)
      AED_REAL,ALLOCATABLE :: tn_varscale(:), tkn_varscale(:), tp_varscale(:),    &
                              toc_varscale(:), tss_varscale(:), turb_varscale(:), &
                              tfe_varscale(:), tal_varscale(:)


      !# Model parameters
      LOGICAL  :: outputLight, outputMetals

     CONTAINS
         PROCEDURE :: define            => aed_define_totals
         PROCEDURE :: calculate         => aed_calculate_totals
!        PROCEDURE :: mobility          => aed_mobility_totals
!        PROCEDURE :: light_extinction  => aed_light_extinction_totals
!        PROCEDURE :: delete            => aed_delete_totals

   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_totals(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_totals_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER :: status
   INTEGER :: i, num_tn,num_tkn,num_tp,num_toc,num_tss,num_turb,num_tfe,num_tal

!  %% NAMELIST   %%  /aed_totals/
!  %% Last Checked 20/08/2021
   CHARACTER(len=40) :: tn_vars(100)       = ''
   AED_REAL          :: tn_varscale(100)   = 1.0
   CHARACTER(len=40) :: tkn_vars(100)      = ''
   AED_REAL          :: tkn_varscale(100)  = 1.0
   CHARACTER(len=40) :: tp_vars(100)       = ''
   AED_REAL          :: tp_varscale(100)   = 1.0
   CHARACTER(len=40) :: toc_vars(100)      = ''
   AED_REAL          :: toc_varscale(100)  = 1.0
   CHARACTER(len=40) :: tss_vars(100)      = ''
   AED_REAL          :: tss_varscale(100)  = 1.0
   CHARACTER(len=40) :: turb_vars(100)     = ''
   AED_REAL          :: turb_varscale(100) = 1.0
   CHARACTER(len=40) :: tfe_vars(10)       = ''
   AED_REAL          :: tfe_varscale(10)   = 1.0
   CHARACTER(len=40) :: tal_vars(10)       = ''
   AED_REAL          :: tal_varscale(10)   = 1.0
   LOGICAL           :: outputLight        = .FALSE.
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_totals/

   NAMELIST /aed_totals/ tn_vars,  tn_varscale,  tkn_vars,  tkn_varscale,  &
                         tp_vars,  tp_varscale,  toc_vars, toc_varscale,   &
                         tss_vars, tss_varscale, turb_vars, turb_varscale, &
                         tfe_vars, tfe_varscale, tal_vars, tal_varscale,   &
                         outputLight
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_totals initialization"

   ! Read the namelist
   read(namlst,nml=aed_totals,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_totals'

   data%outputLight = outputLight

   num_tn  = 0 ; num_tkn  = 0 ; num_tp  = 0 ; num_toc = 0
   num_tss = 0 ; num_turb = 0 ; num_tfe = 0 ; num_tal = 0

   DO i=1,100 ; IF (tn_vars(i)  .EQ. '' ) THEN ; num_tn  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tkn_vars(i) .EQ. '' ) THEN ; num_tkn = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tp_vars(i)  .EQ. '' ) THEN ; num_tp  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (toc_vars(i) .EQ. '' ) THEN ; num_toc = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tss_vars(i) .EQ. '' ) THEN ; num_tss = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (turb_vars(i).EQ. '' ) THEN ; num_turb= i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tfe_vars(i) .EQ. '' ) THEN ; num_tfe = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tal_vars(i) .EQ. '' ) THEN ; num_tal = i-1 ; EXIT ; ENDIF ; ENDDO

   ALLOCATE(data%id_dep_tn(num_tn))     ; ALLOCATE(data%tn_varscale(num_tn))
   ALLOCATE(data%id_dep_tkn(num_tkn))   ; ALLOCATE(data%tkn_varscale(num_tkn))
   ALLOCATE(data%id_dep_tp(num_tp))     ; ALLOCATE(data%tp_varscale(num_tp))
   ALLOCATE(data%id_dep_toc(num_toc))   ; ALLOCATE(data%toc_varscale(num_toc))
   ALLOCATE(data%id_dep_tss(num_tss))   ; ALLOCATE(data%tss_varscale(num_tss))
   ALLOCATE(data%id_dep_turb(num_turb)) ; ALLOCATE(data%turb_varscale(num_turb))
   ALLOCATE(data%id_dep_tfe(num_tfe))   ; ALLOCATE(data%tfe_varscale(num_tfe))
   ALLOCATE(data%id_dep_tal(num_tal))   ; ALLOCATE(data%tal_varscale(num_tal))

   data%num_tn   = num_tn
   data%num_tkn  = num_tkn
   data%num_tp   = num_tp
   data%num_toc  = num_toc
   data%num_tss  = num_tss
   data%num_turb = num_turb
   data%num_tfe  = num_tfe
   data%num_tal  = num_tal

   ! Register external state variable dependencies
   DO i=1,data%num_tn
      data%id_dep_tn(i) =  aed_locate_variable(tn_vars(i))
      data%tn_varscale(i) =  tn_varscale(i)
      !print*,'          TN : ', TRIM(tn_vars(i)), ' * ', data%tn_varscale(i)
   ENDDO
   DO i=1,data%num_tkn
      data%id_dep_tkn(i) =  aed_locate_variable(tkn_vars(i))
      data%tkn_varscale(i) =  tkn_varscale(i)
     !print*,'          TKN : ', TRIM(tkn_vars(i)), ' * ', data%tkn_varscale(i)
   ENDDO
   DO i=1,data%num_tp
      data%id_dep_tp(i) =  aed_locate_variable(tp_vars(i))
      data%tp_varscale(i) =  tp_varscale(i)
      !print*,'          TP : ', TRIM(tp_vars(i)), ' * ', data%tp_varscale(i)
   ENDDO
   DO i=1,data%num_toc
      data%id_dep_toc(i) = aed_locate_variable(toc_vars(i))
      data%toc_varscale(i) = toc_varscale(i)
      !print*,'          TOC : ', TRIM(toc_vars(i)), ' * ', data%toc_varscale(i)
   ENDDO
   DO i=1,data%num_tss
      data%id_dep_tss(i) = aed_locate_variable(tss_vars(i))
      data%tss_varscale(i) = tss_varscale(i)
      !print*,'          TSS : ', TRIM(tss_vars(i)), ' * ', data%tss_varscale(i)
   ENDDO
   DO i=1,data%num_turb
      data%id_dep_turb(i) = aed_locate_variable(turb_vars(i))
      data%turb_varscale(i) = turb_varscale(i)
      !print*,'          TURB : ', TRIM(turb_vars(i)), ' * ', data%turb_varscale(i)
   ENDDO
   DO i=1,data%num_tfe
      data%id_dep_tfe(i) =  aed_locate_variable(tfe_vars(i))
      data%tfe_varscale(i) =  tfe_varscale(i)
      !print*,'          TFE : ', TRIM(tfe_vars(i)), ' * ', data%tfe_varscale(i)
   ENDDO
   DO i=1,data%num_tal
      data%id_dep_tal(i) =  aed_locate_variable(tal_vars(i))
      data%tal_varscale(i) =  tal_varscale(i)
      !print*,'          TAL : ', TRIM(tal_vars(i)), ' * ', data%tal_varscale(i)
   ENDDO

   ! Register environmental dependencies
   IF (data%outputLight) THEN
     data%id_par = aed_locate_global('par')
     data%id_extc = aed_locate_global('extc_coef')
   END IF

   ! Register diagnostic variables
   IF (data%num_tn>0) THEN
     data%id_totals_tn = aed_define_diag_variable('tn',               &
                     'mmol/m**3', 'Total Nitrogen')
   ENDIF

   IF (data%num_tkn>0) THEN
     data%id_totals_tkn = aed_define_diag_variable('tkn',             &
                     'mmol/m**3', 'Total Kjedhal Nitrogen')
   ENDIF

   IF (data%num_tp>0) THEN
     data%id_totals_tp = aed_define_diag_variable('tp',               &
                     'mmol/m**3', 'Total Phosphorus')
   ENDIF

   IF (data%num_toc>0) THEN
     data%id_totals_toc = aed_define_diag_variable('toc',             &
                     'mmol/m**3', 'Total Organic Carbon')
   ENDIF

   IF (data%num_tss>0) THEN
     data%id_totals_tss = aed_define_diag_variable('tss',             &
                     'mg/L', 'Total Suspended Solids')
   ENDIF

   IF (data%num_turb>0) THEN
     data%id_totals_turbidity = aed_define_diag_variable('turbidity', &
                     'NTU', 'Turbidity')
   ENDIF

   IF (data%num_tfe>0) THEN
     data%id_tfe = aed_define_diag_variable('tfe',                    &
                     'mmol/m**3', 'Total Iron')
   ENDIF

   IF (data%num_tal>0) THEN
     data%id_tal = aed_define_diag_variable('tal',                    &
                     'mmol/m**3', 'Total Aluminium')
   ENDIF

   IF (data%outputLight) THEN
     data%id_totals_light = aed_define_diag_variable('light',         &
                     'W/m2', 'Shortwave light flux')

     data%id_totals_par = aed_define_diag_variable('par',            &
                     'W/m2', 'Photosynthetically active light flux')

     data%id_totals_uv = aed_define_diag_variable('uv',               &
                     'W/m2', 'Ultraviolet light flux')

     data%id_totals_extc = aed_define_diag_variable('extc',           &
                     'W/m2', 'Light extinction coefficient')
   ENDIF

END SUBROUTINE aed_define_totals
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_totals(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_totals model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_totals_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: i, count
   AED_REAL :: val, tot

!-------------------------------------------------------------------------------
!BEGIN
   ! Sum over TN variables
   IF (data%num_tn>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tn,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tn(i)); tot = tot + (val*data%tn_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tn) = tot
   ENDIF

   ! Sum over TKN variables
   IF (data%num_tkn>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tkn,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tkn(i)); tot = tot + (val*data%tkn_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tkn) = tot
   ENDIF

   ! Sum over TP variables
   IF (data%num_tp>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tp,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tp(i)); tot = tot + (val*data%tp_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_totals_tp) = tot
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

   ! Sum over TFe variables
   IF (data%num_tfe>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tfe,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tfe(i)); tot = tot + (val*data%tfe_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_tfe) =  tot
   ENDIF

   ! Sum over TAl variables
   IF (data%num_tal>0) THEN
     tot = 0.
     count = ubound(data%id_dep_tal,1)
     DO i=1,count ; val = _STATE_VAR_(data%id_dep_tal(i)); tot = tot + (val*data%tal_varscale(i)) ; ENDDO
     _DIAG_VAR_(data%id_tal) =  tot
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
END SUBROUTINE aed_calculate_totals
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_totals
