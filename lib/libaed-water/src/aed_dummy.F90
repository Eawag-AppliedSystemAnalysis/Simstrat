!###############################################################################
!#                                                                             #
!# aed_dummy.F90                                                               #
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
!# Created July 2019                                                           #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |  ________    | || | ____    ____ | || | ____    ____ | |         !
!         | | |_   ___ `.  | || ||_  |    |  _|| || ||_   \  /   _|| |         !
!         | |   | |   `. \ | || |  | |    | |  | || |  |   \/   |  | |         !
!         | |   | |    | | | || |  | |    , |  | || |  | |\  /| |  | |         !
!         | |  _| |___.' / | || |  | `.__.  ,  | || | _| |_\/_| |_ | |         !
!         | | |________.'  | || |  |_______.   | || ||_____||_____|| |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_dummy
!-------------------------------------------------------------------------------
! aed_dummy --- dummy model
!
! The AED module "dummy" contains only variables to provide vars required in
! other modules that we didnt provide (usually for debugging purposes)
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_dummy_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_dummy_data_t
      !# Variable identifiers
      INTEGER :: num_v, num_dv, num_sv, num_dsv
      INTEGER :: id_sine, id_vsine
      INTEGER,ALLOCATABLE :: id_dummy_v(:), id_dummy_dv(:),           &
                             id_dummy_sv(:), id_dummy_dsv(:)
      AED_REAL,ALLOCATABLE :: dm_max(:), dm_min(:)
      AED_REAL,ALLOCATABLE :: dm_smax(:), dm_smin(:)

     CONTAINS
         PROCEDURE :: define            => aed_define_dummy
         PROCEDURE :: calculate         => aed_calculate_dummy
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_dummy
   END TYPE

! MODULE GLOBALS
   AED_REAL :: today = 1.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_dummy(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_dummy_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i, num_v, num_dv, num_sv, num_dsv

   CHARACTER(len=4),POINTER :: prefix => null()

!  %% NAMELIST   %%  /aed_dummy/
!  %% Last Checked 20/08/2021
   CHARACTER(len=40) :: dm_vars(100) = ''
   AED_REAL          :: dm_max(100) = NaN_
   AED_REAL          :: dm_min(100) = NaN_
   AED_REAL          :: dm_init(100) = 0.
   CHARACTER(len=40) :: dm_dvars(100) = ''
   CHARACTER(len=40) :: dm_svars(100) = ''
   AED_REAL          :: dm_smax(100) = NaN_
   AED_REAL          :: dm_smin(100) = NaN_
   AED_REAL          :: dm_sinit(100) = 0.
   CHARACTER(len=40) :: dm_dsvars(100) = ''
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_dummy/

   NAMELIST /aed_dummy/ dm_vars, dm_max, dm_min, dm_init,             &
                        dm_dvars,                                     &
                        dm_svars, dm_smax, dm_smin, dm_sinit,         &
                        dm_dsvars, diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_dummy initialization"

   num_v = 0 ; num_dv = 0 ; num_sv = 0 ; num_dsv = 0

   ! Read the namelist
   read(namlst,nml=aed_dummy,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_dummy'

   DO i=1,100 ; IF (dm_vars(i)  .EQ. '' ) THEN  ; num_v  = i-1  ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_dvars(i) .EQ. '' ) THEN  ; num_dv = i-1  ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_svars(i)  .EQ. '' ) THEN ; num_sv  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_dsvars(i) .EQ. '' ) THEN ; num_dsv = i-1 ; EXIT ; ENDIF ; ENDDO

   ALLOCATE(data%id_dummy_v(num_v))
   ALLOCATE(data%id_dummy_dv(num_dv))
   ALLOCATE(data%id_dummy_sv(num_sv))
   ALLOCATE(data%id_dummy_dsv(num_dsv))

   ALLOCATE(data%dm_min(num_v))   ; ALLOCATE(data%dm_max(num_v))
   ALLOCATE(data%dm_smin(num_sv)) ; ALLOCATE(data%dm_smax(num_sv))

   data%num_v   = num_v
   data%num_dv  = num_dv
   data%num_sv  = num_sv
   data%num_dsv = num_dsv

   CALL aed_set_prefix(prefix)

   ! Register state variables
   DO i=1,data%num_v
      data%id_dummy_v(i) = aed_define_variable(dm_vars(i), '', '', dm_init(i), dm_min(i), dm_max(i), 0.)
      data%dm_min(i) = dm_min(i)
      data%dm_max(i) = dm_max(i)
   ENDDO

   DO i=1,data%num_sv
      data%id_dummy_sv(i) = aed_define_sheet_variable(dm_svars(i), '', '', dm_sinit(i), dm_smin(i), dm_smax(i), .FALSE.)
      data%dm_smin(i) = dm_smin(i)
      data%dm_smax(i) = dm_smax(i)
   ENDDO

   DO i=1,data%num_dv
      data%id_dummy_dv(i) = aed_define_diag_variable(dm_dvars(i), '', '')
   ENDDO

   DO i=1,data%num_dsv
      data%id_dummy_dsv(i) = aed_define_sheet_diag_variable(dm_dsvars(i), '', '', .FALSE.)
   ENDDO

   data%id_vsine = aed_define_diag_variable('DUM_vol_sine', 'no units', 'DBG volume sine between 0.0 and 1.0')
   data%id_sine = aed_define_sheet_diag_variable('DUM_sine', 'no units', 'DBG sine wave between 0.0 and 1.0', .FALSE.)
END SUBROUTINE aed_define_dummy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_dummy(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_dummy_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
!  INTEGER  :: i, count
!  AED_REAL :: val, tot
   INTEGER :: i
   AED_REAL :: scale, offs

!-------------------------------------------------------------------------------
!BEGIN
   _DIAG_VAR_(data%id_vsine) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * 0.5) + 0.5

   DO i=1,data%num_v
      scale = (data%dm_max(i) - data%dm_min(i)) / 2.
      offs = data%dm_min(i) + scale
      _STATE_VAR_(data%id_dummy_v(i)) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * scale) + offs
   ENDDO
END SUBROUTINE aed_calculate_dummy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_dummy(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_dummy_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: scale, offs

!-------------------------------------------------------------------------------
!BEGIN
   IF (layer_idx .EQ. 1) today = today + 1.0/36.5

   _DIAG_VAR_S_(data%id_sine) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * 0.5) + 0.5

   DO i=1,data%num_sv
      scale = (data%dm_smax(i) - data%dm_smin(i)) / 2.
      offs = data%dm_smin(i) + scale
      _STATE_VAR_S_(data%id_dummy_sv(i)) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * scale) + offs
   ENDDO
END SUBROUTINE aed_calculate_benthic_dummy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_dummy
