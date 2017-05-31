!###############################################################################
!#                                                                             #
!# aed2_riptypes.F90                                                           #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created March 2016                                                          #
!###############################################################################

#include "aed2.h"


MODULE aed2_riptypes
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  PUBLIC  ! all

  integer, parameter :: dble_prec = SELECTED_REAL_KIND(15,307)
  integer, parameter :: GCHP = dble_prec
  integer, parameter :: GCHM_STR_LEN = 32

  INTEGER,   PARAMETER :: IVOID    = -9999
  REAL,      PARAMETER :: VOID     = -9999.0
  CHARACTER(*), PARAMETER :: VOID_STR = "VOID_STRING"

  INTEGER, PARAMETER :: HORTON = 2
  INTEGER, PARAMETER :: GRNAMPT = 3



!-------------------------------------------------------------------------------
! Infiltration type - based on SWMM model approach
!-------------------------------------------------------------------------------
 TYPE Infiltration ! scalars

! infiltration parameters THorton
   AED_REAL     ::        f0                   !  initial infil. rate (m/d)
   AED_REAL     ::        fmin                 !  minimum infil. rate (m/d)
   AED_REAL     ::        Fmax                 !  maximum total infiltration (m/d)
   AED_REAL     ::        decay                !  decay coeff. of infil. rate (1/d)
   AED_REAL     ::        regen                !  regeneration coeff. of infil. rate (1/d)
  !-----------------------------
   AED_REAL     ::        tp                   !  present time on infiltration curve (d)
   AED_REAL     ::        Fe                   !  cumulative infiltration (m)

! infiltration parameters TGrnAmpt
   AED_REAL     ::        S               ! avg. capillary suction (m)
   AED_REAL     ::        Ks              ! saturated conductivity (m/d)
   AED_REAL     ::        IMDmax          ! max. soil moisture deficit (m/m)
   !-----------------------------
   AED_REAL     ::        IMD             ! current initial soil moisture deficit
   AED_REAL     ::        F               ! current cumulative infiltrated volume (m)
   AED_REAL     ::        Fu              ! current upper zone infiltrated volume (m)
   AED_REAL     ::        Lu              ! depth of upper soil zone (m)
   AED_REAL     ::        T               ! time until start of next rain event (d)
  CHARACTER(len=3)  ::    Sat = "Sat"      ! saturation flag

! infiltration parameters TCurveNum
   AED_REAL     ::        Smax         ! max. infiltration capacity (m)
!   AED_REAL     ::        regen        ! infil. capacity regeneration constant (1/d)
   AED_REAL     ::        Tmax         ! maximum inter-event time (d)
   !-----------------------------
!   AED_REAL     ::        S            ! current infiltration capacity (m)
!   AED_REAL     ::        F            ! current cumulative infiltration (m)
   AED_REAL     ::        P            ! current cumulative precipitation (m)
 !  AED_REAL     ::        T            ! current inter-event time (d)
   AED_REAL     ::        Se           ! current event infiltration capacity (m)
  ! AED_REAL     ::        f            ! previous infiltration rate (m/d)

END TYPE Infiltration
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! Main type storing hydrolgogic and soil parameters needs for land modules
!-------------------------------------------------------------------------------
TYPE SoilUnit ! scalars

    INTEGER   :: UnitNum = 0, Substrate = 0
    AED_REAL  :: Depth, Area, Bathy, Sb         ! Physical properties
    AED_REAL  :: Ucap, St, Ssat, theta, PhreaticHgt,PhreaticDepth, qss, qse, qie, recharge, qsuc ! Hydrologic dynamics
    AED_REAL  :: aep, ass, Bss, Debs, Dper, Porosity, Density, phi, qinf ! Hydrologic parameters
    AED_REAL  :: Dcap, Zcap, Dtrn, Ztrn, fc, Sus, S_top, S_trn, S_cap, pastMaxLevel, rn, et, es, CapHgt
    ! Moisture disaggregation model
    AED_REAL, DIMENSION(:), ALLOCATABLE :: Moisture
    AED_REAL, DIMENSION(:), ALLOCATABLE :: Temp
    INTEGER   :: nlay
    ! Infiltration sub-model
    TYPE ( Infiltration ) :: inf

END TYPE SoilUnit
!-------------------------------------------------------------------------------




END MODULE aed2_riptypes
!------------------------------------------------------------------------------!
