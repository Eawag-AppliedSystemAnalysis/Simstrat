!###############################################################################
!#                                                                             #
!# aed2_ass.F90                                                                #
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
!# Created Oct 2015                                                            #
!#                                                                             #
!###############################################################################

! /AED/StudySites/CoorongLowerLakes/simulations/CoorongLowerLakes_i7/SeawaterEIS/010_anc_896_ASS01_yAB/caedym
! ! -----------------------------------------------------------------------------!
! ! Depth to perching (m)                                                        !
!  1.50000E+00                 : ASS SOIL 1 - Sand with Low OM content           !
!  1.50000E+00                 : ASS SOIL 2 - Clay with Low OM content           !
!  1.50000E+00                 : ASS SOIL 3 - Sand with High OM content          !
!  1.50000E+00                 : ASS SOIL 4 - Clay with High OM content          !
! ! Sediment density (kg/m3)                                                     !
!  1530.000000                 : ASS SOIL 1 - Sand with Low OM content           !
!  1230.000000                 : ASS SOIL 2 - Clay with Low OM content           !
!  1530.000000                 : ASS SOIL 3 - Sand with High OM content          !
!  1230.000000                 : ASS SOIL 4 - Clay with High OM content          !
! ! Potential Acidity (mol H+/kg)                                                !
!  0.041000E+00                : ASS SOIL 1 - Sand with Low OM content           !
!  0.250000E+00                : ASS SOIL 2 - Clay with Low OM content           !
!  0.041000E+00                : ASS SOIL 3 - Sand with High OM content          !
!  0.250000E+00                : ASS SOIL 4 - Clay with High OM content          !
! ! Acidity production rate (- scaling factor for Taylor OxCOn algorithm)        !
!  1.000000000                 : ASS SOIL 1 - Sand with Low OM content           !
!  1.000000000                 : ASS SOIL 2 - Clay with Low OM content           !
!  1.000000000                 : ASS SOIL 3 - Sand with High OM content          !
!  1.000000000                 : ASS SOIL 4 - Clay with High OM content          !
! ! ASS flux parameter - baseflow (/day)                                         !
!  0.000000000                 : ASS SOIL 1 - Sand with Low OM content           !
!  0.000000000                 : ASS SOIL 2 - Clay with Low OM content           !
!  0.000000000                 : ASS SOIL 3 - Sand with High OM content          !
!  0.000000000                 : ASS SOIL 4 - Clay with High OM content          !
! ! ASS flux parameter - rain pulsing (/day) & (mm)                              !
!  0.000 0.000                 : ASS SOIL 1 - Sand with Low OM content           !
!  0.000 0.000                 : ASS SOIL 2 - Clay with Low OM content           !
!  0.000 0.000                 : ASS SOIL 3 - Sand with High OM content          !
!  0.000 0.000                 : ASS SOIL 4 - Clay with High OM content          !
! ! ASS flux parameter - rewetting (Day1: mol H+/m2/day & Day2-90:mol H+/m2/day) !
! 0.1380 0.0070                : ASS SOIL 1 - Sand with Low OM content           !
! 0.1610 0.0100                : ASS SOIL 2 - Clay with Low OM content           !
! 0.1380 0.0070                : ASS SOIL 3 - Sand with High OM content          !
! 0.1610 0.0100                : ASS SOIL 4 - Clay with High OM content          !
!#
!# Translates to :
!#
!     AED_REAL :: aep(MAX_ASS_PARAMS)
!     AED_REAL :: ass(MAX_ASS_PARAMS)
!     AED_REAL :: Bss(MAX_ASS_PARAMS)
!     AED_REAL :: Debs(MAX_ASS_PARAMS)
!     AED_REAL :: Dper(MAX_ASS_PARAMS)
!     AED_REAL :: Porosity(MAX_ASS_PARAMS)
!     AED_REAL :: Density(MAX_ASS_PARAMS)
!     AED_REAL :: zASS
!
! &aed2_ass
!    nPars   = 4
!    dPer    =  1.50000E+00,  1.50000E+00,  1.50000E+00,  1.50000E+00  ! Depth to perching (m)
!    Density =  1530.000000,  1230.000000,  1530.000000,  1230.000000  ! Sediment density (kg/m3)

!    pass_0 = 0.041000E+00, 0.250000E+00, 0.041000E+00, 0.250000E+00  ! Potential Acidity (mol H+/kg)
!    Racid    =  1.000000000,  1.000000000,  1.000000000,  1.000000000  ! Acidity production rate
!    flux_bf =  0.000000000,  0.000000000,  0.000000000,  0.000000000  ! ASS flux parameter - baseflow (/day)
!    flux_rn=        0.000,        0.000,        0.000,        0.000  ! ASS flux parameter - rain pulsing (/day)
!    flux_rn_max   =        0.000,        0.000,        0.000,        0.000  ! ASS flux parameter - rain pulsing (mm)
!    flux_rw_a =       0.1380,       0.1610,       0.1380,       0.1610  ! ASS flux parameter - rewetting (Day1
!    flux_rw_d =       0.0070,       0.0100,       0.0070,       0.0100  ! ASS flux parameter - rewetting (Day2
! /


#include "aed2.h"

!
MODULE aed2_ass
!------------------------------------------------------------------------------+
! AED2 module for Acid Sulfate Soils                                           |
!------------------------------------------------------------------------------+
   USE aed2_core
   USE aed2_util
   USE aed2_riptypes

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_ass_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_ass_data_t
      !# Variable identifiers
      INTEGER :: id_uzaass, id_szaass, id_taa, id_passt
      INTEGER,ALLOCATABLE :: id_anc0(:), id_pass0(:)
      INTEGER,ALLOCATABLE :: id_anc(:), id_pass(:)
      INTEGER :: id_wettime, id_drytime
      INTEGER :: id_asstracer

      !# Environmental variables
      INTEGER :: id_E_rain, id_E_area, id_E_matz, id_E_bath, id_E_salt, id_E_nearlevel

      !# Diagnostic variables
      INTEGER :: id_rwet, id_rchg, id_bflw, id_pyrox, id_so4r, id_pml, id_sflw

      !# Dependant variable IDs
      INTEGER :: id_g_ubalchg, id_g_so4
      INTEGER :: id_l_depth, id_l_Sb, id_l_St, id_l_Ssat, id_l_theta, id_l_Stop, id_l_capz
      INTEGER :: id_l_phreatic, id_l_qss, id_l_qse, id_l_qcap, id_l_qper, id_l_wt
      INTEGER :: id_l_soiltemp

      LOGICAL :: simProfiles

      !# ASS params
      INTEGER  :: n_zones, nlay
      AED_REAL,ALLOCATABLE :: active_zones(:)
      AED_REAL :: pass_0(MAX_ASS_PARAMS)
      AED_REAL :: Racid(MAX_ASS_PARAMS)
      AED_REAL :: flux_bf(MAX_ASS_PARAMS)
      AED_REAL :: flux_rn(MAX_ASS_PARAMS)
      AED_REAL :: flux_rn_max(MAX_ASS_PARAMS)
      AED_REAL :: flux_rw_a(MAX_ASS_PARAMS)
      AED_REAL :: flux_rw_d(MAX_ASS_PARAMS)
      AED_REAL :: zASS(MAX_ASS_PARAMS)
      AED_REAL :: Porosity(MAX_ASS_PARAMS),Density(MAX_ASS_PARAMS)
      AED_REAL :: X_hso4

     CONTAINS
         PROCEDURE :: define             => aed2_define_ass
         PROCEDURE :: initialize         => aed2_initialize_ass
!        PROCEDURE :: calculate_surface  => aed2_calculate_surface_ass
!        PROCEDURE :: calculate          => aed2_calculate_ass
!        PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_ass
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_ass
!        PROCEDURE :: calculate_dry      => aed2_calculate_dry_ass
!        PROCEDURE :: equilibrate        => aed2_equilibrate_ass
!        PROCEDURE :: mobility           => aed2_mobility_ass
!        PROCEDURE :: light_extinction   => aed2_light_extinction_ass
!        PROCEDURE :: delete             => aed2_delete_ass
   END TYPE


!-------------------------------------------------------------------------------
!MODULE VARIABLES

   TYPE(SoilUnit) :: SoilCol
   INTEGER, PARAMETER :: ASSCLAYL = 4
   INTEGER, PARAMETER :: ASSCLAYH = 5
   INTEGER, PARAMETER :: ASSSANDL = 6
   INTEGER, PARAMETER :: ASSSANDH = 7
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_ass(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_ass_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER :: status, i

   !# ASS params
   INTEGER  :: n_zones, nlay
   INTEGER  :: active_zones(MAX_ZONES)
   AED_REAL :: pass_0(MAX_ASS_PARAMS)
   AED_REAL :: Racid(MAX_ASS_PARAMS)
   AED_REAL :: flux_bf(MAX_ASS_PARAMS)
   AED_REAL :: flux_rn(MAX_ASS_PARAMS)
   AED_REAL :: flux_rn_max(MAX_ASS_PARAMS)
   AED_REAL :: flux_rw_a(MAX_ASS_PARAMS)
   AED_REAL :: flux_rw_d(MAX_ASS_PARAMS)
   AED_REAL :: zASS(MAX_ASS_PARAMS)
   AED_REAL :: Porosity(MAX_ASS_PARAMS),Density(MAX_ASS_PARAMS)
   AED_REAL :: X_hso4

   LOGICAL  :: simProfiles

   CHARACTER(len=64) :: acid_link = ''
   CHARACTER(len=64) :: so4_link = ''
   CHARACTER(4) :: trac_name

   NAMELIST /aed2_ass/ n_zones, active_zones, nlay, acid_link, simProfiles, &
                       pass_0, Racid, flux_bf, flux_rn, flux_rn_max, flux_rw_a, flux_rw_d, zASS, &
                       Porosity, Density, so4_link, X_hso4
!
!-------------------------------------------------------------------------------
!BEGIN

   simProfiles = .FALSE.

   ! Read the namelist
   read(namlst,nml=aed2_ass,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_ass'

   data%simProfiles = simProfiles
   data%n_zones = n_zones         ! THIS NEEDS TO BE 4 AT THE MOMENT!
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   data%X_hso4 = X_hso4
   DO i=1,n_zones
      data%pass_0(i)  = pass_0(i)
      data%Racid(i)  = Racid(i)
      data%flux_bf(i)  = flux_bf(i)
      data%flux_rn(i) = flux_rn(i)
      data%flux_rn_max(i) = flux_rn_max(i)
      data%flux_rw_a(i) = flux_rw_a(i)
      data%flux_rw_d(i)  = flux_rw_d(i)
      data%zASS(i)  = zASS(i)
      data%nlay  = nlay
      data%Porosity(i) = Porosity(i)
      data%Density(i) = Density(i)
   ENDDO

   ALLOCATE(data%id_pass0(nlay));  ALLOCATE(data%id_anc0(nlay))
   ALLOCATE(data%id_pass(nlay));  ALLOCATE(data%id_anc(nlay))

   !# Register state variables
   data%id_asstracer = aed2_define_variable('tracer', 'mmol/m**3', 'ass tracer', &
                                                                zero_,minimum=zero_)
   IF ( acid_link .EQ. '' ) &
   data%id_g_ubalchg = aed2_define_variable('acidity', 'mmol/m**3', 'acidity', &
                                                                zero_,minimum=zero_)
     ! Register diagnostic variables
   data%id_passt = aed2_define_sheet_diag_variable('pass_total','molH+/kg','soil total potential acidity (depth averaged)')
   data%id_uzaass = aed2_define_sheet_diag_variable('uzaass','molH+/m2','unsat zone available acidity')
   data%id_szaass = aed2_define_sheet_diag_variable('szaass','molH+/m2','sat zone available acidity')
   data%id_taa    = aed2_define_sheet_diag_variable('taa','molH+/L','total actual acidity (porewater)')
   trac_name = 'lay0'
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_pass0(i) = aed2_define_sheet_diag_variable('pass0_'//TRIM(trac_name),'molH+/kg','starting potential acidity (sulfides)')
   ENDDO
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_pass(i)  = aed2_define_sheet_diag_variable('pass_'//TRIM(trac_name),'molH+/m2','potential acidity (sulfides)')
   ENDDO
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_anc0(i) = aed2_define_sheet_diag_variable('anc0_'//TRIM(trac_name),'molH+/kg','starting acid neutralising capacity')
   ENDDO
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_anc(i)  = aed2_define_sheet_diag_variable('anc_'//TRIM(trac_name),'molH+/m2','acid neutralising capacity')
   ENDDO
   data%id_rwet  = aed2_define_sheet_diag_variable('rwet','molH+/m2/day','acid neutralising capacity')
   data%id_rchg  = aed2_define_sheet_diag_variable('rchg','molH+/m2/day','acid neutralising capacity')
   data%id_bflw  = aed2_define_sheet_diag_variable('bflw','molH+/m2/day','acid neutralising capacity')
   data%id_sflw  = aed2_define_sheet_diag_variable('sflw','molH+/m2/day','acid neutralising capacity')
   data%id_pyrox = aed2_define_sheet_diag_variable('pyrox','/day','acid neutralising capacity')
   data%id_so4r  = aed2_define_sheet_diag_variable('so4r','molH+/m2/day','acid neutralising capacity')
   data%id_pml   = aed2_define_sheet_diag_variable('pml','m','past maximum groundwater level')
   data%id_wettime = aed2_define_sheet_diag_variable('wettime','day','time cell has been innundated')
   data%id_drytime = aed2_define_sheet_diag_variable('drytime','day','time cell has been exposed')

   ! Register module dependencies
   IF ( .NOT. acid_link .EQ. '' ) &
   data%id_g_ubalchg  = aed2_locate_global(TRIM(acid_link))     ! ('GEO_ubalchg')       !,'meq/L','chg balance variable')
   IF ( .NOT. so4_link .EQ. '' ) &
   data%id_g_so4      = aed2_locate_global(TRIM(so4_link))      ! ('GEO_SO4')       !,'meq/L','chg balance variable')
   data%id_l_depth    = aed2_locate_global_sheet('LND_depth')   !,'m','soil depth (to datum)')
   data%id_l_phreatic = aed2_locate_global_sheet('LND_phreatic')!,'m','depth of phreatic surface below surface')
   data%id_l_wt       = aed2_locate_global_sheet('LND_wt')!,'m','depth of phreatic surface below surface')
   data%id_l_Sb       = aed2_locate_global_sheet('LND_Sb')      !,'mm','total capacity for water storage')
   data%id_l_St       = aed2_locate_global_sheet('LND_St')      ! 'mm', 'total soil water storage at time t')
   data%id_l_Ssat     = aed2_locate_global_sheet('LND_Ssat')    !,'mm','saturated zone soil water storage')
   data%id_l_Stop     = aed2_locate_global_sheet('LND_Stop')    !,'mm','unsaturated storage capacity')
   data%id_l_capz     = aed2_locate_global_sheet('LND_capz')    !,'mm','unsaturated storage capacity')
   data%id_l_theta    = aed2_locate_global_sheet('LND_theta')   !,'-','unsaturated moisture content')
   data%id_l_qss      = aed2_locate_global_sheet('LND_qss')     !,'mm','sat zone seepage')
   data%id_l_qse      = aed2_locate_global_sheet('LND_qs')     !,'mm','surface runoff')
   data%id_l_qcap     = aed2_locate_global_sheet('LND_qcap')    !,'mm','capillarity')
   data%id_l_qper     = aed2_locate_global_sheet('LND_qper')    !,'mm','recharge')
   data%id_l_soiltemp = aed2_locate_global_sheet('LND_soiltemp_10cm') !,'C','soiltemp_10cm')

   !# Register environmental dependencies
   data%id_E_rain = aed2_locate_global_sheet('rain')        ! daily rainfall
   data%id_E_area = aed2_locate_global_sheet('layer_area')  ! cell area
   data%id_E_matz = aed2_locate_global_sheet('material')    ! material index
   data%id_E_bath = aed2_locate_global_sheet('bathy')       ! cell bathy
   data%id_E_salt = aed2_locate_global('salinity')          ! salinity of overlying water
   data%id_E_nearlevel= aed2_locate_global_sheet('nearest_depth')

   ! Initialisation occurs in first call of calculate_riparian

END SUBROUTINE aed2_define_ass
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_initialize_ass(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_ass_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: PASSt, ANCt
!-------------------------------------------------------------------------------
!BEGIN
   !---------------------------------------------------------------------------
   ! Prime local cell hydrology object with data from global arrays
   CALL SetSoilHydrology(data, column, SoilCol)

   !---------------------------------------------------------------------------
   ! Initialise vertical profiles of ANC and PASS
   _DIAG_VAR_S_(data%id_pml) = SoilCol%PhreaticHgt  !Reset
   ANCt = _DIAG_VAR_S_(data%id_anc0(1))
   PASSt = _DIAG_VAR_S_(data%id_pass0(1))
   CALL SetSoilASSProfile(data, column, SoilCol, PASSt, ANCt)
   _DIAG_VAR_S_(data%id_uzaass)= zero_


END SUBROUTINE aed2_initialize_ass
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian_ass(data, column, layer_idx, pc_wet)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_ass_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   AED_REAL :: rain, salt, bathy, area, matz, soiltemp
   INTEGER  :: var, i, sub, assz

   AED_REAL :: newAcidity, acidityFlux, acid_anc, newANC, avgpyrox, pass_vol
   AED_REAL :: avgLevel, acid_gen, acid_potl, OxdnRate, SZDepth, UZDepth, NeutRate
   AED_REAL :: acid_flux0, acid_flux1, acid_flux2, acid_flux3, acid_flux4, acid_flux5, acid_flux6, acid_flux7
   AED_REAL :: Kass, sedDensity, pass_depth, moist, maxqse, acid_neut
   AED_REAL :: oldlevel, newlevel

   AED_REAL :: ASSRWET, ASSRCHG, ASSBFLW, ASSSFLW, PYROX, SO4REDN
   AED_REAL :: PASSt, TAA, ANCt, SZAASS, UZAASS, UZpH, SZpH, SO4, PH
!
!-------------------------------------------------------------------------------
!BEGIN

   !---------------------------------------------------------------------------
   ! Set local cell properties
   matz = _STATE_VAR_S_(data%id_E_matz)
   IF(.NOT.in_zone_set(matz,data%active_zones)) RETURN
   assz = 1
   DO i =1,data%n_zones
     IF( INT(matz) == data%active_zones(i) )THEN
       assz = i
     ENDIF
   ENDDO
   avgLevel = _STATE_VAR_S_(data%id_E_nearlevel)
   bathy = _STATE_VAR_S_(data%id_E_bath)
   IF( ABS(bathy) > 1e9 ) RETURN
   rain = _STATE_VAR_S_(data%id_E_rain)
   area = _STATE_VAR_S_(data%id_E_area)
   soiltemp = _STATE_VAR_S_(data%id_l_soiltemp)


  ! WRITE(*,"('Start acid sulfate soils calculations...',I4)"),INT(matz)

   !---------------------------------------------------------------------------
   ! Zero dummy (working) variables
   acid_anc   = zero_
   acid_gen   = zero_
   acid_potl  = zero_
   acid_flux0 = zero_
   acid_flux1 = zero_
   acid_flux2 = zero_
   acid_flux3 = zero_
   acid_flux4 = zero_
   acid_flux5 = zero_
   acid_flux6 = zero_
   acid_flux7 = zero_
   acid_neut  = zero_
   pass_vol   = zero_
   avgpyrox   = zero_
   acidityFlux = zero_
   newAcidity  = zero_
   avgpyrox    = zero_
   maxqse = zero_

   ! Zero diagnostic variables
   ASSRWET = zero_
   ASSRCHG = zero_
   ASSBFLW = zero_
   PYROX   = zero_
   SO4REDN = zero_
   UZpH = zero_
   SZpH = zero_

   !---------------------------------------------------------------------------
   ! Get local ASS properties for this cell
   TAA = _DIAG_VAR_S_(data%id_taa)
   ANCt = _DIAG_VAR_S_(data%id_anc0(1))
   PASSt = _DIAG_VAR_S_(data%id_pass0(1))
   UZAASS = _DIAG_VAR_S_(data%id_uzaass)
   SZAASS = _DIAG_VAR_S_(data%id_szaass)

   !---------------------------------------------------------------------------
   ! Prime local cell hydrology object with data from global arrays
   CALL SetSoilHydrology(data, column, SoilCol)

   IF(SoilCol%qse>maxqse) THEN
      maxqse = SoilCol%qse
   END IF

   !---------------------------------------------------------------------------
   ! Get PASS depth & volume
   pass_depth = MAX(SoilCol%Bathy - SoilCol%pastMaxLevel, zero_)
   pass_vol   = pass_depth * SoilCol%Area

   SZDepth    = SoilCol%Depth - SoilCol%PhreaticDepth
   UZDepth    = SoilCol%PhreaticDepth


   !---------------------------------------------------------------------------
   ! Dry sediment : Check for exposure of new PASS and do oxidation, etc
   IF(pc_wet<0.1) THEN

       ! Get WT level from SoilHydrology
       newLevel = SoilCol%PhreaticHgt

       !-----------------------------------------------------------------------!
       !-- Increase in potential ASS due to increase in depth & then oxidise

       CALL UpdatePASSProfile( SoilCol,          &
                               column,           &
                               newAcidity,       &
                               avgpyrox,         &
                               data%RAcid(assz), &
                               temp=soiltemp     )

       PASSt = zero_
       DO i=1,data%nlay
         PASSt = PASSt + _DIAG_VAR_S_(data%id_pass(i))
       ENDDO
       IF( pass_depth > 0.05 ) THEN
         _DIAG_VAR_S_(data%id_passt) = PASSt / (pass_depth * SoilCol%Density * 1e-3)
       ENDIF

       ! Increase available acidity
       UZAASS =  UZAASS + newAcidity
       acid_gen = acid_gen + newAcidity
       PYROX = avgpyrox

       newAcidity = TAA
       CALL UpdateANCProfile( SoilCol, &
                               column, &
                              newAcidity )

       ANCt = zero_
       DO i=1,data%nlay
         ANCt = ANCt + _DIAG_VAR_S_(data%id_anc(i))
       ENDDO

       ! Consume available acidity
       UZAASS = UZAASS - newAcidity
       acid_anc = acid_anc - newAcidity

       ! Reset pastMaxLevel if necessary
       oldLevel = newLevel
       IF(SoilCol%PhreaticHgt < SoilCol%pastMaxLevel) THEN
           ! Newly exposed sediment is yet to contibute
           SoilCol%pastMaxLevel = SoilCol%PhreaticHgt
       END IF

       !-----------------------------------------------------------------------
       !-- Decrease in available acidity in SZ due to ANC/SO4 reduction

       !newANC = c%SED%ASS%flux_rn_max(sIndex) * WQ%sed(i,SZAASS) * DDT
       !newANC = data%flux_rn_max(assz) * MIN(SZDepth,0.5) * area * DDT * 1e-3
       newANC = data%flux_rn_max(assz) * MIN(SZDepth,0.5) * DDT * 1e-3
       IF(newANC > SZAASS) newANC=SZAASS

       ! Consume available acidty
       SZAASS = SZAASS - newANC

       acid_neut = acid_neut - newANC

!print *,'newAcidity',newAcidity,newANC,PYROX


       !-----------------------------------------------------------------------
       !-- Decrease in available acidity in due to acidity recharge/discharge/etc

       ! R-RCG : recharge/percolation   [UNITS=> mol/L/timestep = mol/L * ((m/ts)/m) ]
       IF( SoilCol%recharge>zero_ .AND. UZDepth>0.1 ) THEN

         acidityFlux = UZAASS * 0.5 * ( MIN(SoilCol%recharge/UZDepth,1.0) **data%flux_rn(assz) )

         UZAASS = UZAASS - acidityFlux
         SZAASS = SZAASS + acidityFlux

         acid_flux0 = acid_flux0 + acidityFlux
         ASSRCHG = acidityFlux/DDT  ! molH+/L/day
       END IF


       ! R-BFLW : baseflow  [ UNITS>= mol/timestep = mol/L * m/day /m *day/ts ]
       IF( SoilCol%qss>zero_ .AND. SZDepth>0.01 ) THEN

        IF(SZAASS>zero_)THEN
           acidityFlux = data%flux_bf(assz) * SZAASS * ( SoilCol%qss / MIN(SZDepth,0.5) )
        ELSE
           acidityFlux = zero_
        ENDIF

         IF ( acidityFlux > SZAASS ) THEN
           acidityFlux = SZAASS
         END IF

         SZAASS = SZAASS - acidityFlux
         acid_flux1 = acid_flux1 + acidityFlux

        ! ! Store acidity flux in DissFluxRates  [ UNITS = mol /m2 /day ]
        ! ! Used in UpdateBoundaryWaterConcs where non-vdo2D cells are fluxed
        ! DissFluxRates(i,DICHM(chargeBalCol)) = DissFluxRates(i,DICHM(chargeBalCol)) &
        !                                      + acidityFlux / ( ivec(1) * DDT)

        ! set diagnostic flux and update seepage tracer into the water: mol/m2/day
         ASSBFLW = acidityFlux/DDT
        _FLUX_VAR_R_(data%id_asstracer) = _FLUX_VAR_R_(data%id_asstracer) + ASSBFLW/secs_per_day

        ! ! Store acidity flux in DissFluxRates  [ UNITS = mol /m2 /day ]
        ! ! Used in UpdateBoundaryWaterConcs where non-vdo2D cells are fluxed
        ! DissFluxRates(i,DICHM(chargeBalCol)) = DissFluxRates(i,DICHM(chargeBalCol)) &
        !                                      + acidityFlux / ( ivec(1) * DDT)

!        _FLUX_VAR_R_(data%id_g_ubalchg) = _FLUX_VAR_R_(data%id_g_ubalchg) + 1e-3*ASSBFLW/secs_per_day
        _FLUX_VAR_R_(data%id_g_ubalchg) = _FLUX_VAR_R_(data%id_g_ubalchg) + ASSBFLW/secs_per_day
        IF ( data%id_g_so4>0 ) THEN
         _FLUX_VAR_R_(data%id_g_so4) = _FLUX_VAR_R_(data%id_g_so4) + (ASSBFLW/secs_per_day) * data%X_hso4
        ENDIF

       END IF


       ! R-DIS : If water table moves up, move some UZAASS -> SZAASS
       IF(newLevel > oldLevel .AND. UZDepth > 0.01 ) THEN

         acidityFlux = MAX((( newLevel - oldLevel ) / UZDepth),1.0) * UZAASS

         UZAASS = UZAASS - acidityFlux
         SZAASS = SZAASS + acidityFlux

         acid_flux6 = acid_flux6 + acidityFlux
       END IF


       ! R-POND : ponding and sat excess  [ UNITS>= mol/timestep = mol/L * m/day /m *day/ts ]
       IF( SoilCol%qse>zero_ .AND. UZAASS>zero_ .AND. UZDepth>0.1 ) THEN

         !  [UNITS=> mol/L/timestep = mol/L *m/ts *day/tstep *m2 * ??]
         acidityFlux = UZAASS * ( MIN(SoilCol%qse/UZDepth,1.0) **data%flux_rn(assz) )

         UZAASS = UZAASS - acidityFlux
         acid_flux5 = acid_flux5 + acidityFlux

        ! set diagnostic flux and update seepage tracer into the water: mol/m2/day
         ASSSFLW = acidityFlux/DDT

        _FLUX_VAR_R_(data%id_asstracer) = _FLUX_VAR_R_(data%id_asstracer) + ASSSFLW/secs_per_day

        ! ! Store acidity flux in DissFluxRates  [ UNITS = mol /m2 /day ]
        ! ! Used in UpdateBoundaryWaterConcs where non-vdo2D cells are fluxed
        ! DissFluxRates(i,DICHM(chargeBalCol)) = DissFluxRates(i,DICHM(chargeBalCol)) &
        !                                      + acidityFlux / ( ivec(1) * DDT )
        _FLUX_VAR_R_(data%id_g_ubalchg) = _FLUX_VAR_R_(data%id_g_ubalchg) + ASSSFLW/secs_per_day
!        _FLUX_VAR_R_(data%id_g_ubalchg) = _FLUX_VAR_R_(data%id_g_ubalchg) + 1e-3*ASSSFLW/secs_per_day
        IF ( data%id_g_so4>0 ) THEN
         ! SO4 flux associated with pyrite oxidation (SHOULD THIS BE MULTIPLIED BY 1000 as SO4 in mmol/m3?)
         _FLUX_VAR_R_(data%id_g_so4) = _FLUX_VAR_R_(data%id_g_so4) + (ASSSFLW/secs_per_day) * data%X_hso4
        ENDIF



         ! WQ%sed(i,SZAASS) = WQ%sed(i,SZAASS) + SedimentData(i,1,SED(DICHM(chargeBalCol)))
         ! SedimentData(i,1,SED(DICHM(chargeBalCol))) = 0.0
         !
         !
         ! acidityFlux = WQ%sed(i,SZAASS) * ( DDT * SoilCol(i)%qse / MIN(SZDepth,0.5) )
         !
         ! IF (acidityFlux > WQ%sed(i,SZAASS) ) THEN
         !   acidityFlux = 0.9* WQ%sed(i,SZAASS)
         ! END IF
         !
         ! acid_flux5 = acid_flux5 + acidityFlux
         !
         ! ! Store acidity flux in DissFluxRates  [ UNITS = mol /m2 /day ]
         ! ! Used in UpdateBoundaryWaterConcs where non-vdo2D cells are fluxed
         ! DissFluxRates(i,DICHM(chargeBalCol)) = DissFluxRates(i,DICHM(chargeBalCol)) &
         !                                      + acidityFlux / ( ivec(1) * DDT )
         !
         ! WQ%sed(i,SZAASS) = WQ%sed(i,SZAASS) - acidityFlux
         !
         ! !WQ%sed(i,ASSBFLW)= acidityFlux/DDT


       END IF




       ASSRWET = 0.0
       _DIAG_VAR_S_(data%id_drytime) = _DIAG_VAR_S_(data%id_drytime) + DDT
       _DIAG_VAR_S_(data%id_wettime) = zero_

       TAA = zero_
       IF( SoilCol%Sus>0.01 ) TAA = UZAASS  / ( MAX(SoilCol%Sus,0.001) )
       IF(TAA>1e-12) THEN
         UZpH  = MAX (-LOG10(TAA), 2.5)
       ELSE
         UZpH  = 9.0
       END IF
       TAA = UZAASS * 1e3 / ( MAX(UZDepth,0.01)  * SoilCol%Density )

       SZpH = SZAASS / ( MAX(SoilCol%Ssat,0.001) )
       IF(SZpH>1e-9) THEN
         SZpH  = -LOG10(SZpH)
       ELSE
         SZpH  = 9.0
       END IF

     !-------------------------------------------------------------------------
     ! R-WET : Newly re-wetted sediment (ii is in vdo2D but wasn't last step)
     ! Flux rate is initially high
     ELSE IF(pc_wet>0.1 .AND. _DIAG_VAR_S_(data%id_wettime)<=0.25) THEN   ! .AND. .NOT. ANY(Prevvdo2D == i)) THEN

       salt = _STATE_VAR_(data%id_E_salt)
       acidityFlux = GetAcidityFluxRate(data%flux_rw_a(assz),salt,assz,1)

       ! Store acidity flux in DissFluxRates for update within routine
       ! Update BottomWaterConcs [ UNITS = mol /m2 /day ]

       ASSRWET = acidityFlux !* 1e-3

   !    DissFluxRates(i,DICHM(chargeBalCol)) = acidityFlux * 1e-3
       _FLUX_VAR_B_(data%id_g_ubalchg) = _FLUX_VAR_B_(data%id_g_ubalchg) + ASSRWET/secs_per_day
       IF ( data%id_g_so4>0 ) THEN
         _FLUX_VAR_B_(data%id_g_so4) = _FLUX_VAR_B_(data%id_g_so4) + (ASSRWET/secs_per_day) * data%X_hso4
       ENDIF

       acid_flux3 = acid_flux3 + acidityFlux*DDT

       _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT



     !  ANCt = zero_
     !  DO i=1,data%nlay
     !o    ANCt = ANCt + _DIAG_VAR_S_(data%id_anc(i))
     !  ENDDO
     !  IF(ANCt>UZAASS) THEN
     !    SZAASS = zero_
     !
     !

     ! ANCt    =  _DIAG_VAR_S_(data%id_anc(i))
     !  IF(ANCt>SZAASS) THEN
     !    SZAASS = zero_
     !  END IF
     SZAASS = SZAASS + UZAASS
     UZAASS = zero_


     !-------------------------------------------------------------------------
     ! R-WET : Old re-wetted sediment - keep fluxing remainder (via diffusion)
     ELSE IF(pc_wet>0.1 .AND. _DIAG_VAR_S_(data%id_wettime)>0.25 ) THEN

       SO4 = _STATE_VAR_(data%id_g_so4)
       PH = 8.
       salt = _STATE_VAR_(data%id_E_salt)
       IF(salt < 1.)  salt = 1.0
       IF(salt > 40.) salt = 40.0

       IF( _DIAG_VAR_S_(data%id_wettime)>90. .OR.  _DIAG_VAR_S_(data%id_drytime)<5.0) THEN
         ! Old innundation .or. soil not dry for very long before innundation

         ! SO4 reduction of 5mmol H+/day from Koschorreck M, Tittel J.
         acidityFlux = -0.005

         ! Scale SO4 reduction down if SO4 is limiting (>280mg/L SO4 is non-limiting)
         acidityFlux = acidityFlux * ( SO4 ) / ( SO4+(153.*(1e3/96.)) )

         SO4REDN = acidityFlux
         acid_flux7 = acid_flux7 + acidityFlux*DDT

         ! Past dryness and PASS production is reset by now
         _DIAG_VAR_S_(data%id_drytime) = zero_

       ELSEIF(_DIAG_VAR_S_(data%id_wettime) < 1.) THEN
         ! Innundated within last 1 days
         acidityFlux = GetAcidityFluxRate(data%flux_rw_a(assz),salt,assz,1)

         acid_flux3 = acid_flux3 + acidityFlux*DDT

       ELSE
         ! Innundated within last 1-90 days
         acidityFlux = GetAcidityFluxRate(data%flux_rw_d(assz),salt,assz,2)

         acid_flux4 = acid_flux4 + acidityFlux*DDT
       END IF

       UZpH = PH
       SZpH = PH

       ! Store acidity flux in ASSRWET for update within routine
       ! UpdateBottomWaterConcs [ UNITS = mol /m2 /day ]
       ASSRWET = acidityFlux * 1e-3

  !     !IF (acidityFlux*DDT > SedimentData(i,1,SED(DICHM(chargeBalCol)))) THEN
  !     !  acidityFlux = 0.9*SedimentData(i,1,SED(DICHM(chargeBalCol)))/DDT
  !     !END IF
  !     DissFluxRates(i,DICHM(chargeBalCol)) = acidityFlux * 1e-3
       _FLUX_VAR_B_(data%id_g_ubalchg) = _FLUX_VAR_B_(data%id_g_ubalchg) + ASSRWET/secs_per_day
       IF ( data%id_g_so4>0 ) THEN
         _FLUX_VAR_B_(data%id_g_so4) = _FLUX_VAR_B_(data%id_g_so4) + (ASSRWET/secs_per_day) * data%X_hso4
       ENDIF

       ! Update the time-counter for inundation
       _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT

     END IF
     !------------------------------




   ! Update the main arrays
   _DIAG_VAR_S_(data%id_uzaass) = UZAASS
   _DIAG_VAR_S_(data%id_szaass) = SZAASS
   _DIAG_VAR_S_(data%id_taa)    = TAA
   _DIAG_VAR_S_(data%id_rwet)   = ASSRWET
   _DIAG_VAR_S_(data%id_rchg)   = ASSRCHG
   _DIAG_VAR_S_(data%id_bflw)   = ASSBFLW
   _DIAG_VAR_S_(data%id_sflw)   = ASSSFLW
   _DIAG_VAR_S_(data%id_pyrox)  = PYROX
   _DIAG_VAR_S_(data%id_so4r)   = SO4REDN
   _DIAG_VAR_S_(data%id_pml)    = SoilCol%pastMaxLevel  ! _DIAG_VAR_S_(data%id_l_wt)



   !---------------------------------------------------------------------------!
   ! Write Acid Sulfate Soils ITS file                                         !
   !
   !   (1)  Time
   !   (2)  Total moles potential H+ remaining (PASS)
   !   (3)  Total soil volume containing potential H+
   !   (4)  Change in potential ASS (PASS) in this timestep
   !   (5)  Total moles available H+ in sediment (UZ)
   !   (6)  Total moles available H+ in sediment (SZ)
   !   (7)  Change in available acidity in this timestep
   !   (8)  Change in available ANC in this timestep
   !   (9)  ASS flux to SatZone 0  - recharge
   !   (10)  ASS flux to water 1    - baseflow
   !   (11) ASS flux to water 2    - rain pulsing
   !   (12) ASS flux to water 3a   - rewetting (initial advective pulse)
   !   (13) ASS flux to water 3b   - rewetting (subsequent diffusive flux)
   !   (14) ASS flux to water 7    - SO4 reduction
   !   (15) ASS flux to water 5    - ponding and runoff (subsequent diffusive flux)
   !   (16) ASS flux to water 6    - UZ->SZ from rising water table
   !   (17) CHGBAL change in water column due to boundary flux
   !   (18) XSoilFluxi change in water column due to boundary flux
   !   (19) Area of wet cells
   !   (20) Area of dry cells
   !   (21) Area of dry cells not topgraphically connected to water
   !   (22) Average water level of domain
   !   (23) Rainfall
   !   (24) Boundary
   !   (25) maxqse (ponding)
   !   (26) SZ acid consumption
   !                      SUM(WQ%sed(:,PASS)),                                  &
   !                      pass_vol,                                             &
   !                      acid_potl,                                            &
   !                      SUM(WQ%sed(:,UZAASS)),                                &
   !                      SUM(WQ%sed(:,SZAASS)),                                &
   !                      acid_gen,                                             &
   !                      acid_anc,                                             &
   !                      acid_flux0,                                           &
   !                      acid_flux1,                                           &
   !                      acid_flux2,                                           &
   !                      acid_flux3,                                           &
   !                      acid_flux4,                                           &
   !                      acid_flux7,                                           &
   !                      acid_flux5,                                           &
   !                      acid_flux6,                                           &
   !                      ass_chgbal,                                           &
   !                      XSoilFluxi(DICHM(chargeBalCol)),                      &
   !                      wetarea, dryarea, lostarea,                           &
   !                      avgLevel, rain,                                       &
   !                      boundary_flux/DDT,                                    &
   !                      maxqse,                                               &
   !                      acid_neut
   !
   !---------------------------------------------------------------------------!

CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
! SUBROUTINE UpdatePASSProfile(theSoil,thePASS,PASS0,newAASS, avgpyrox, pyrox)
 SUBROUTINE UpdatePASSProfile(theSoil,column,newAASS, avgpyrox, pyrox, temp)
   !-- Incoming
   TYPE(SoilUnit)       :: theSoil
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   !AED_REAL, DIMENSION(20)   :: thePASS
   !AED_REAL, DIMENSION(20)   :: PASS0
   AED_REAL                  :: newAASS, avgpyrox, pyrox
   AED_REAL, OPTIONAL        :: temp
   !-- Local
   INTEGER  :: lay
   AED_REAL :: dep, depm1, middep, newAcidity, OxdnRate

   !---------------------------------------------------------------------------!

   newAASS    = zero_

   DO lay = 1,theSoil%nlay

       newAcidity = zero_

!       dep   = theSoil%Bathy  + ( lay    * ( theSoil%Depth / theSoil%nlay ) )
!       depm1 = theSoil%Bathy  + ((lay-1) * ( theSoil%Depth / theSoil%nlay ) )
       dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
       depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )

       middep = (dep+depm1)/2.

       !-----------------------------------------------------------------------!
       ! Check for newly exposed layers and update PASScontent based on PASScontent0
       IF(middep > theSoil%pastMaxLevel .AND. _DIAG_VAR_S_(data%id_pass0(lay)) > zero_ ) THEN

         _DIAG_VAR_S_(data%id_pass(lay)) = (dep - depm1)  & !* theSoil%Area &
                      * _DIAG_VAR_S_(data%id_pass0(lay))  * theSoil%Density * 1e-3
         ! Once a layer is added to the PASScontent array it cannot be repeated.
         _DIAG_VAR_S_(data%id_pass0(lay)) = zero_

         IF( middep > 0.5 )  _DIAG_VAR_S_(data%id_pass(lay)) = zero_
       END IF


       !-----------------------------------------------------------------------!
       !-- Increase in actual acidity due to oxidation process
       OxdnRate = GetPyriteOxidnRate( pyrox, &
                                      theSoil%Moisture(lay),  theSoil%Substrate )
       IF( PRESENT(temp) )  OxdnRate = OxdnRate * 1.05**(temp-20.)


   !    IF(theTAA*1e3 < 2500) THEN
         ! [ UNITS = mol/L * /day * day/timestep = mol /L /timestep ]
         newAcidity = _DIAG_VAR_S_(data%id_pass(lay)) * OxdnRate * DDT
   !    ELSE
   !      newAcidity = wq_zero
   !    END IF


       ! Reduce PASS based on used amount
       _DIAG_VAR_S_(data%id_pass(lay)) = _DIAG_VAR_S_(data%id_pass(lay)) - newAcidity

       ! Increase available acidity
       newAASS = newAASS + newAcidity

       avgpyrox  =  avgpyrox + OxdnRate
   END DO


   ! Average PYROX for plots
   avgpyrox =  avgpyrox/theSoil%nlay

 END SUBROUTINE UpdatePASSProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
 SUBROUTINE UpdateANCProfile(theSoil,column,conAASS)
   !-- Incoming
   TYPE(SoilUnit)       :: theSoil
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   !AED_REAL, DIMENSION(20)   :: theANC
   !AED_REAL, DIMENSION(20)   :: ANC0
   AED_REAL                  :: conAASS
   !-- Local
   INTEGER :: lay
   AED_REAL :: dep, depm1, middep, newAcidity, NeutRate, theTAA

   !---------------------------------------------------------------------------!

   theTAA = conAASS
   conAASS = zero_

   DO lay = 1,theSoil%nlay

       newAcidity = zero_

      !dep   = theSoil%Bathy + ( lay * ( theSoil%Depth/theSoil%nlay ))
      !depm1 = theSoil%Bathy + ((lay-1) * ( theSoil%Depth/theSoil%nlay ))
       dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
       depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )

       middep = (dep+depm1)/2.

       !-----------------------------------------------------------------------!
       ! Check for newly exposed layers and update ANCcontent
       IF(middep > theSoil%pastMaxLevel .AND. _DIAG_VAR_S_(data%id_anc0(lay)) > zero_ ) THEN

         _DIAG_VAR_S_(data%id_anc(lay)) = (dep - depm1) & !* theSoil%Area &
                      * _DIAG_VAR_S_(data%id_anc0(lay))  * theSoil%Density * 1e-3
         ! Once a layer is added to the ANCcontent array it cannot be repeated.
         _DIAG_VAR_S_(data%id_anc0(lay)) = zero_

       END IF

       !-----------------------------------------------------------------------!
       !-- Decrease in actual acidity due to nuetralisation process
       NeutRate = 0.00027  * 2.0   !20% /year

       IF(theTAA*1e3 > 1.0) THEN
         ! [ UNITS = mol/L * /day * day/timestep = mol /L /timestep ]
         newAcidity = _DIAG_VAR_S_(data%id_anc(lay)) * NeutRate * DDT
       ELSE
         newAcidity = zero_
       END IF

       ! Reduce ANC based on used amount
        _DIAG_VAR_S_(data%id_anc(lay)) =  _DIAG_VAR_S_(data%id_anc(lay)) - newAcidity

       ! Increase consumed acidty
       conAASS = conAASS + newAcidity


   END DO


 END SUBROUTINE UpdateANCProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END SUBROUTINE aed2_calculate_riparian_ass


!###############################################################################
SUBROUTINE SetSoilHydrology(data, column, theSoil)
!-------------------------------------------------------------------------------
! Primes the local cell Bucket model "object" with relevant parameters based on!
! its material zone id number                                                  !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_ass_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   TYPE(SoilUnit),INTENT(inout) :: theSoil
   !INTEGER,INTENT(in) :: MatZoneID
!LOCALS
   INTEGER :: lay, zindex, i
   AED_REAL :: dep, depm1, middep, botmoist, topmoist

!-------------------------------------------------------------------------------
!BEGIN


   theSoil%UnitNum = -1
   theSoil%Substrate = INT(_STATE_VAR_S_(data%id_E_matz))
   zindex = 1
   DO i =1,data%n_zones
     IF( theSoil%Substrate == INT(data%active_zones(i)) )THEN
       zindex = i
     ENDIF
   ENDDO

   ! Physical properties
   theSoil%Depth = _DIAG_VAR_S_(data%id_l_depth)
   theSoil%Area = _STATE_VAR_S_(data%id_E_area)
   theSoil%Bathy = _STATE_VAR_S_(data%id_E_bath)
   theSoil%Sb = _STATE_VAR_S_(data%id_l_Sb)
   theSoil%Porosity = data%Porosity(zindex)
   theSoil%Density = data%Density(zindex)

   ! Hydrologic dynamics (these have been previously sent by aed2_land)
   theSoil%St = _DIAG_VAR_S_(data%id_l_St)
   theSoil%Ssat = _DIAG_VAR_S_(data%id_l_Ssat)
   theSoil%Sus = MAX( theSoil%St-theSoil%Ssat , zero_ )
   theSoil%S_top = MAX( _DIAG_VAR_S_(data%id_l_Stop) , zero_ )
   theSoil%theta = _DIAG_VAR_S_(data%id_l_theta)
   !theSoil%Ucap = _DIAG_VAR_S_(data%id_l_Ucap)

   theSoil%PhreaticDepth =  _DIAG_VAR_S_(data%id_l_phreatic)
   theSoil%PhreaticHgt = _DIAG_VAR_S_(data%id_l_wt) !theSoil%Bathy - theSoil%PhreaticDepth
   theSoil%CapHgt = _DIAG_VAR_S_(data%id_l_capz)
   theSoil%Dcap = 0.5*  (theSoil%CapHgt - theSoil%PhreaticHgt)
   theSoil%Dtrn = 0.5*  (theSoil%CapHgt - theSoil%PhreaticHgt)
   theSoil%pastMaxLevel = _DIAG_VAR_S_(data%id_pml)

   theSoil%qss = _DIAG_VAR_S_(data%id_l_qss)
   theSoil%qse = _DIAG_VAR_S_(data%id_l_qse)
   theSoil%recharge = _DIAG_VAR_S_(data%id_l_qper)
   theSoil%qsuc = _DIAG_VAR_S_(data%id_l_qcap)

   ! Moisture disaggregation model
   theSoil%nlay = data%nlay

   IF(.NOT.ALLOCATED(theSoil%Moisture)) ALLOCATE( theSoil%Moisture(theSoil%nlay) )

   theSoil%Moisture(:) = theSoil%theta


   ! Update Layer Moisture Fractions (wt% water)
   botmoist = 0.9
   topmoist = (theSoil%S_top) /  MAX( (theSoil%Bathy-theSoil%CapHgt) * theSoil%Porosity,0.001)

   DO lay = 1,theSoil%nlay

         dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
         depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )

         middep = (dep+depm1)/2.

         theSoil%Moisture(lay) = 0.0

         IF( middep  <= theSoil%PhreaticHgt ) THEN
            ! Saturated Zone: % VWC/theta
            theSoil%Moisture(lay) = 1.000
         ELSE IF(middep > theSoil%PhreaticHgt .AND. middep <=  theSoil%PhreaticHgt+theSoil%Dcap) THEN
            ! Near saturated region above the water table: % VWC/theta
            theSoil%Moisture(lay) = botmoist
         ELSE IF(middep >  theSoil%PhreaticHgt+theSoil%Dcap .AND. middep <= theSoil%CapHgt) THEN
            ! Transition region: % VWC/theta
            !print *,'middep',topmoist,theSoil%Ztrn, (botmoist-topmoist) , ((middep)-theSoil%Bathy),(theSoil%Dcap - middep)
            theSoil%Moisture(lay) = botmoist - (botmoist-topmoist) *  (middep-(theSoil%PhreaticHgt+theSoil%Dcap)) / theSoil%Dtrn  ! linear gradient
         ELSE
            ! Top of unsaturated zone: % VWC/theta
            theSoil%Moisture(lay) = topmoist !&
            ! MIN(theSoil%S_top / MAX(theSoil%Sb - theSoil%Ssat - theSoil%S_trn - theSoil%S_cap,1e-2),1.0)
         ENDIF

         !print *,'moist',lay,dep,middep,theSoil%Moisture(lay)

         ! convert to % wtWC
         theSoil%Moisture(lay) = theSoil%Moisture(lay)*theSoil%Porosity * 1e3 /        &
            ( theSoil%Moisture(lay)*theSoil%Porosity*1e3 + (1.-theSoil%Porosity)*theSoil%Density )

         !theSoil%Moisture(lay) =theSoil%theta  !TEMPP
   ENDDO

END SUBROUTINE SetSoilHydrology
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 FUNCTION GetPyriteOxidnRate(MaxRate, Moisture, SOIL) RESULT (OxdnRate)
   !-- Incoming
   AED_REAL :: MaxRate,Moisture   ! Max rate of FeS2 oxidation for this soil
   AED_REAL :: OxdnRate           ! The actual oxidation rate calculated
   INTEGER  :: SOIL               ! Integer soil type
   !-- Local
   AED_REAL,PARAMETER   :: X=0.0333 ! Parameter - Jeff Turner ASS study
   AED_REAL             :: MC

   !---------------------------------------------------------------------------!

   MC = Moisture

   OxdnRate = 0.0

   ! SANDY ASS
   IF(SOIL == ASSSANDL .OR. SOIL == ASSSANDH) THEN

     OxdnRate = -9.7011*MC**3. + 2.1949*MC**2. + 0.0025*MC + 0.0006

     OxdnRate = MaxRate* OxdnRate

   ! CLAYEY ASS
   ELSE IF(SOIL == ASSCLAYL .OR. SOIL == ASSCLAYH) THEN

     IF(MC > 0.225 .AND. MC < 0.48) THEN
       OxdnRate = -0.0142*MC + 0.0068
     ELSE IF(MC <= 0.225) THEN
       OxdnRate = 0.0142*MC
     ELSE
       OxdnRate = 0.0
     END IF

     OxdnRate = MaxRate* OxdnRate

   ! UNSURE
   ELSE

     OxdnRate = ( MaxRate - X*(Moisture*0.5) )/2.0

   END IF

   ! CHECK
   IF(OxdnRate<0.0) OxdnRate = 0.0
   IF(OxdnRate>0.5) OxdnRate = 0.5

 END FUNCTION GetPyriteOxidnRate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
 FUNCTION GetAcidityFluxRate(FWRate,Salinity,SOIL,theTime) RESULT (FluxRate)
   !-- Incoming
   AED_REAL   :: FWRate,Salinity    ! Rate of FW release of H for this soil
   AED_REAL   :: FluxRate           ! The actual flux rate calculated
   INTEGER    :: SOIL, theTime      ! Integer soil type & Time
   !-- Local
   AED_REAL,PARAMETER   :: BC1=0.436 ! Parameter - CSIRO ASS study CLAY
   AED_REAL,PARAMETER   :: BS1=0.152 ! Parameter - CSIRO ASS study SAND
   AED_REAL,PARAMETER   :: BCt=0.033 ! Parameter - CSIRO ASS study CLAY
   AED_REAL,PARAMETER   :: BSt=0.006 ! Parameter - CSIRO ASS study SAND
   AED_REAL             :: kk

   !---------------------------------------------------------------------------!
   IF(theTime == 1) THEN
     IF(SOIL == ASSCLAYL .OR. SOIL == ASSCLAYH) THEN
        kk = BC1
     ELSE
        kk = BS1
     END IF
   ELSE
     IF(SOIL == ASSCLAYL .OR. SOIL == ASSCLAYH) THEN
        kk = BCt
     ELSE
        kk = BSt
     END IF
   END IF

   IF(FWRate < 1e-10) THEN
     FluxRate = 0.00
   ELSE
     FluxRate =  kk * Salinity/35. + FWRate
   END IF

 END FUNCTION GetAcidityFluxRate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE SetSoilASSProfile(data, column, theSoil, PASSt, ANCt)
!-------------------------------------------------------------------------------
! Primes the local cell Bucket model "object" with relevant parameters based on!
! its material zone id number                                                  !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_ass_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   TYPE(SoilUnit),INTENT(inout) :: theSoil ! passing SoilUnit into theSoil
   AED_REAL,INTENT(in) :: PASSt, ANCt
!LOCALS
   INTEGER :: lay, matz, zindex, i
   AED_REAL :: dep(theSoil%nlay), depm1(theSoil%nlay), middep(theSoil%nlay)
!-------------------------------------------------------------------------------
!BEGIN

   matz = INT(_STATE_VAR_S_(data%id_E_matz))
   IF (.NOT.in_zone_set(REAL(matz),data%active_zones) ) RETURN
   zindex = 1
   DO i =1,data%n_zones
     IF( matz == INT(data%active_zones(i)) )THEN
       zindex = i
     ENDIF
   ENDDO


   ! Set the grid depths of layers below surface
   DO lay = 1,theSoil%nlay
     dep(lay)   =  ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
     depm1(lay) =  (  lay *    ( theSoil%Depth / theSoil%nlay ) )
!     dep(lay)   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
!     depm1(lay) = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )
!     dep(lay)    = theSoil%Bathy  + (  lay *    ( theSoil%Depth / theSoil%nlay ) )
!     depm1(lay)  = theSoil%Bathy  + ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
     middep(lay) = ( dep(lay)+depm1(lay) )/2.
   END DO

   ! Set active arrays pass and anc to 0 (no oxidation as yet) & pass0 and anc0 to top value
   DO lay = 1,theSoil%nlay
     _DIAG_VAR_S_(data%id_pass(lay)) = zero_
     _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt
     _DIAG_VAR_S_(data%id_anc(lay)) = zero_
     _DIAG_VAR_S_(data%id_anc0(lay)) = ANCt
   ENDDO
   _DIAG_VAR_S_(data%id_passt) = zero_

     !PASSContent(i,:)     = 0.0
     !PASSContent0(i,1)    = WQ%sed(i,PASS)
     !PASSContent0(i,2)    = WQ%sed(i,PASS)
     !ANCContent(i,:)     = 0.0
     !ANCContent0(i,1)    = WQ%sed(i,ANC)
     !ANCContent0(i,2)    = WQ%sed(i,ANC)



   IF (.NOT. data%simProfiles)  RETURN

!IF(theSoil%bathy > 0.5 .and. theSoil%PhreaticDepth>0.) THEN
!   print *,'data%zASS(matz)',theSoil%bathy,data%zASS(matz),theSoil%Depth,theSoil%PhreaticHgt,middep
!!PAUSE
!ENDIF

   ! Lake Alex == 1
   IF( data%zASS(zindex)>0.5 .AND. data%zASS(zindex)<1.5) THEN

       IF(theSoil%Substrate == ASSCLAYL .OR. theSoil%Substrate == ASSCLAYH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt       ! Surface concs by CSIRO
              _DIAG_VAR_S_(data%id_anc0(lay))  = ANCt        ! Surface concs by CSIRO
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 1.009705089
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.187140007
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.83412953
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.487488851
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.667303624
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.761085667
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.500477718
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.034682483
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.333651812
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.308279299
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.166825906
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.581876115
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.215260078
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.855472932
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.205900944
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.01769224
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           END IF
         END DO

       ELSE IF(theSoil%Substrate == ASSSANDL .OR. theSoil%Substrate == ASSSANDH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt       ! Surface concs by CSIRO
              _DIAG_VAR_S_(data%id_anc0(lay))  = ANCt        ! Surface concs by CSIRO
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.025127245
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.128028612
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.162224986
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.040779625
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.109189895
              _DIAG_VAR_S_(data%id_anc0(lay))  = 2.741532621
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.056154803
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.44407992
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.087351916
              _DIAG_VAR_S_(data%id_anc0(lay))  = 4.750826297
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 4.057572674
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.092343454
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.874113298
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.690556051
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_anc0(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_anc0(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.945188106
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.791124722
           END IF
         END DO

       END IF

     ! Currency == 2
!CAB was KASS
     ELSEIF( data%zASS(zindex)>1.5 .AND. data%zASS(zindex) <2.5) THEN

       IF(theSoil%Substrate == ASSCLAYL .OR. theSoil%Substrate == ASSCLAYH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt      ! Surface concs by CSIRO
              _DIAG_VAR_S_(data%id_anc0(lay))  = ANCt        ! Surface concs by CSIRO
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.582346104
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.094914577
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.808629162
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.092250308
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.579018412
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.108235921
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.290133148
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.275751823
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.288573293
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.274895451
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.106070183
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.517534219
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.172624024
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.716688314
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.205900944
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.927165551
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.287013437
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.119892097
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.287013437
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.119892097
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) =  0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  =  0.089919073
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) =  0.155985564
              _DIAG_VAR_S_(data%id_anc0(lay))  =  0.089919073
           END IF
         END DO

       ELSE IF(theSoil%Substrate == ASSSANDL .OR. theSoil%Substrate == ASSSANDH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt       ! Surface concs by CSIRO
              _DIAG_VAR_S_(data%id_anc0(lay))  = ANCt        ! Surface concs by CSIRO
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.026740382
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.188401867
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.03327692
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.254371066
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.043675958
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.365670896
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.047419611
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.463183135
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.054074995
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.744663137
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.147666334
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.731341792
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.2214995
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.006094515
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.2214995
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.006094515
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.215260078
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.670729677
           END IF
         END DO

       END IF


     ! Lake Albert == 3
!CAB was KASS
     ELSEIF( data%zASS(zindex)>2.5 .AND. data%zASS(zindex) <3.5) THEN

       IF(theSoil%Substrate == ASSCLAYL .OR. theSoil%Substrate == ASSCLAYH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt       ! Surface concs by CSIRO
              _DIAG_VAR_S_(data%id_anc0(lay))  = ANCt        ! Surface concs by CSIRO
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 1.009705089
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.187140007
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.83412953
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.487488851
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.667303624
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.761085667
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.500477718
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.034682483
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.333651812
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.308279299
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.166825906
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.581876115
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.215260078
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.855472932
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.205900944
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.01769224
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.342130858
           END IF
         END DO

       ELSE IF(theSoil%Substrate == ASSSANDL .OR. theSoil%Substrate == ASSSANDH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = PASSt       ! Surface concs by CSIRO
              _DIAG_VAR_S_(data%id_anc0(lay))  = ANCt        ! Surface concs by CSIRO
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.025127245
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.128028612
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.162224986
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.040779625
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.109189895
              _DIAG_VAR_S_(data%id_anc0(lay))  = 2.741532621
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.056154803
              _DIAG_VAR_S_(data%id_anc0(lay))  = 5.44407992
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.087351916
              _DIAG_VAR_S_(data%id_anc0(lay))  = 4.750826297
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_anc0(lay))  = 4.057572674
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.092343454
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.874113298
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_anc0(lay))  = 3.690556051
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_anc0(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_anc0(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_anc0(lay))  = 1.945188106
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pass0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_anc0(lay))  = 0.791124722
           END IF
         END DO

       END IF

     ! Constant
     ELSE

      ! PASSContent0(i,3:20) = PASSContent0(i,1) !+ c%SED%ASS%KASS(sIndex) * (middep(3:20)-dep(1))
      ! ANCContent0(i,3:20)  = ANCContent0(i,1)  !+  3.0 * (middep(3:20)-dep(1))

     END IF


END SUBROUTINE SetSoilASSProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




END MODULE aed2_ass
