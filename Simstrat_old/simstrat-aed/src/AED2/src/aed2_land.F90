!###############################################################################
!#                                                                             #
!# aed2_land.F90                                                               #
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
!# Created Feb 2016                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"
#define _MAX_ZONES_ 100

MODULE aed2_land
!------------------------------------------------------------------------------+
! AED2 module for land/soil dynamics                                           |
!------------------------------------------------------------------------------+
   USE  aed2_core
   USE  aed2_util
   USE  aed2_riptypes

   IMPLICIT NONE

   PRIVATE

   PUBLIC aed2_land_data_t

   TYPE,extends(aed2_model_data_t) :: aed2_land_data_t
      !# Variable identifiers
      INTEGER :: id_tracer, id_erosion
      INTEGER :: id_St, id_Ssat, id_Scap, id_Strn, id_Stop, id_Theta, id_smd
      INTEGER :: id_Phreatic, id_wt, id_Ucap, id_dh
      INTEGER :: id_soiltemp_10cm

      !# Environmental variables
      INTEGER :: id_E_rain,id_E_area,id_E_height,id_E_density,id_E_material,id_E_temp, &
                 id_E_bathy,id_E_airtemp, id_E_nearest, id_E_rainloss, id_E_nearlevel

      !# Diagnostic variables
      INTEGER :: id_depth,id_area,id_bathy,id_Sb,id_qss,id_qs,id_qcap,id_qper,id_qinf,id_tp
      INTEGER :: id_rain, id_matz, id_rainloss, id_evap

      !# Configuration options and soil parameters
      LOGICAL  :: simTemp, simSalt, simErosion
      INTEGER  :: infil_model
      INTEGER  :: nlay                     ! Moisture disaggregation model
      INTEGER  :: n_zones
      AED_REAL,ALLOCATABLE :: active_zones(:)
      AED_REAL, DIMENSION(:),ALLOCATABLE :: aep, ass, Bss, Dper, minbathy
      AED_REAL, DIMENSION(:),ALLOCATABLE :: Zcap, Ztrn, zASS, Zebs
      AED_REAL, DIMENSION(:),ALLOCATABLE :: fc, phi, Porosity, Density
      AED_REAL, DIMENSION(:),ALLOCATABLE :: f0, fmin, Fmax, decay, regen

      !# Active procedures
      CONTAINS
          PROCEDURE :: define            => aed2_define_land
          PROCEDURE :: initialize        => aed2_initialize_land
          PROCEDURE :: calculate_dry     => aed2_calculate_dry_land
          PROCEDURE :: rain_loss         => aed2_rain_loss_land
          PROCEDURE :: calculate_benthic => aed2_calculate_benthic_land
   END TYPE

!------------------------------------------------------------------------------+
!MODULE VARIABLES
   TYPE(SoilUnit) :: SoilCol
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep
!  AED_REAL, PARAMETER :: DRYDZ = 0.04      ! Host model considers a cell below this as "dry"
!   INTEGER, PARAMETER  :: Hydrology = 1, BUCKET = 1
!   AED_REAL :: day_rain, angfreq, epdayav, epamp, epphase
!------------------------------------------------------------------------------+


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_land(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model-aed namelist is read and the variables exported
! by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!LOCALS
   LOGICAL  :: simTemp,simSalt,simErosion
   INTEGER  :: status,i
   INTEGER  :: infil_model
   INTEGER  :: nlay, n_zones
   INTEGER  :: active_zones(MAX_ZONES)
   AED_REAL :: Dper(MAX_ASS_PARAMS),minbathy(MAX_ASS_PARAMS)
   AED_REAL :: aep(MAX_ASS_PARAMS),Zebs(MAX_ASS_PARAMS)
   AED_REAL :: ass(MAX_ASS_PARAMS),Bss(MAX_ASS_PARAMS)
   AED_REAL :: Zcap(MAX_ASS_PARAMS),Ztrn(MAX_ASS_PARAMS)
   AED_REAL :: Porosity(MAX_ASS_PARAMS),Density(MAX_ASS_PARAMS)
   AED_REAL :: fc(MAX_ASS_PARAMS),phi(MAX_ASS_PARAMS)
   AED_REAL :: f0(MAX_ASS_PARAMS), fmin(MAX_ASS_PARAMS), Fmax(MAX_ASS_PARAMS), decay(MAX_ASS_PARAMS), regen(MAX_ASS_PARAMS)
   CHARACTER(len=64) :: turbid_link = ''

   NAMELIST /aed2_land/ n_zones, active_zones, minbathy, Dper, &
                        aep, Zebs, ass, Bss,   &
                        Porosity, Density, Zcap, Ztrn, &
                        infil_model, fc, phi, f0, fmin, Fmax, decay, regen, &
                        simTemp, simSalt, simErosion, turbid_link


!-------------------------------------------------------------------------------
!BEGIN
   read(namlst,nml=aed2_land,iostat=status) ! Read the namelist
   IF (status /= 0) STOP 'Error reading namelist aed2_land'

   ALLOCATE(data%Zebs(n_zones),data%aep(n_zones))
   ALLOCATE(data%ass(n_zones),data%Bss(n_zones))
   ALLOCATE(data%Dper(n_zones),data%minbathy(n_zones),data%Ztrn(n_zones),data%Zcap(n_zones))
   ALLOCATE(data%Porosity(n_zones),data%Density(n_zones),data%fc(n_zones),data%phi(n_zones))
   ALLOCATE(data%f0(n_zones),data%fmin(n_zones))
   ALLOCATE(data%Fmax(n_zones),data%decay(n_zones),data%regen(n_zones))

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   ! Assign the locally read in values to the module data structure
   data%simTemp = simTemp
   data%simSalt = simSalt
   data%simErosion = simErosion
   data%minbathy(1:n_zones) = minbathy(1:n_zones)  !
   data%Dper(1:n_zones) = Dper(1:n_zones)
   data%Zebs(1:n_zones) = Zebs(1:n_zones)
   data%aep(1:n_zones) = aep(1:n_zones)
   data%ass(1:n_zones) = ass(1:n_zones)
   data%Bss(1:n_zones) = Bss(1:n_zones)
   data%Porosity(1:n_zones) = Porosity(1:n_zones)
   data%Density(1:n_zones) = Density(1:n_zones)
   data%fc(1:n_zones) = fc(1:n_zones)
   data%phi(1:n_zones) = phi(1:n_zones)
   data%Zcap(1:n_zones) = Zcap(1:n_zones)
   data%Ztrn(1:n_zones) = Ztrn(1:n_zones)
   data%f0(1:n_zones) = f0(1:n_zones)
   data%fmin(1:n_zones) = fmin(1:n_zones)
   data%Fmax(1:n_zones) = Fmax(1:n_zones)
   data%decay(1:n_zones) = decay(1:n_zones)
   data%regen(1:n_zones) = regen(1:n_zones)


   ! Register state variables (mass consevred quantities, subject to ODE)
   data%id_tracer = aed2_define_variable('tracer', 'mmol/m**3', 'local runoff fraction', &
                                                                zero_,minimum=zero_)
!   data%id_erosion = aed2_define_variable(turbid_link, 'g/m3', '...', &
!                                                                zero_,minimum=zero_)

   ! Register diagnostic variables
   !-- soil bucket dimensions and capacity
   data%id_area =     aed2_define_sheet_diag_variable('area','m2','soil area')
   data%id_bathy =    aed2_define_sheet_diag_variable('bathy','m','soil surface height')
   data%id_depth =    aed2_define_sheet_diag_variable('depth','m','soil depth (to datum)')
   data%id_Sb =       aed2_define_sheet_diag_variable('Sb','m','total capacity for water storage')
   !-- soil water storages
   data%id_St =       aed2_define_sheet_diag_variable('St', 'm', 'total soil water storage at time t')
   data%id_Ssat =     aed2_define_sheet_diag_variable('Ssat','m','saturated zone soil water storage')
   data%id_Scap =     aed2_define_sheet_diag_variable('Scap','m','unsaturated zone soil water storage')
   data%id_Strn =     aed2_define_sheet_diag_variable('Strn','m','unsaturated zone soil water storage')
   data%id_Stop =     aed2_define_sheet_diag_variable('Stop','m','unsaturated zone soil water storage')
   data%id_theta =    aed2_define_sheet_diag_variable('theta','-','unsaturated moisture content')
   data%id_smd =      aed2_define_sheet_diag_variable('smd','m','soil moisture deficit')
   !-- soil water heads
   data%id_Ucap =     aed2_define_sheet_diag_variable('capz','m','capillary fringe height')
   data%id_wt =       aed2_define_sheet_diag_variable('wt','m','water table height')
   data%id_phreatic = aed2_define_sheet_diag_variable('phreatic','m','depth of water table below soil')
   data%id_dh =       aed2_define_sheet_diag_variable('dh','m','head difference between soil and water')
   !-- soil water fluxes
   data%id_qs =       aed2_define_sheet_diag_variable('qs','m/ts','surface runoff')
   data%id_qinf =     aed2_define_sheet_diag_variable('qinf','m/ts','water infiltrating into the soil')
   data%id_qper =     aed2_define_sheet_diag_variable('qper','m/ts','recharge')
   data%id_qcap =     aed2_define_sheet_diag_variable('qcap','m/ts','capillarity')
   data%id_qss =      aed2_define_sheet_diag_variable('qss','m/ts','sat zone seepage')
   data%id_tp =       aed2_define_sheet_diag_variable('tp','sec','equivalent time on the Horton curve')
   !-- soil temperature and salinity
   data%id_soiltemp_10cm = aed2_define_sheet_diag_variable('soiltemp_10cm','degC','soiltemp_10cm')
   !-- outputs for debugging
   data%id_rain =     aed2_define_sheet_diag_variable('rain','m','input rain')
   data%id_evap =     aed2_define_sheet_diag_variable('evap','m','input rain')
   data%id_matz =     aed2_define_sheet_diag_variable('matz','-','material zone')
   data%id_rainloss = aed2_define_sheet_diag_variable('rainloss','m','rainloss')

   ! Register environmental dependencies (vars from host model we need)
   data%id_E_temp =     aed2_locate_global('temperature')
   data%id_E_rain =     aed2_locate_global_sheet('rain')
   data%id_E_area =     aed2_locate_global_sheet('layer_area')
   data%id_E_bathy =    aed2_locate_global_sheet('bathy')
   data%id_E_height =   aed2_locate_global('layer_ht')
   data%id_E_density =  aed2_locate_global('density')
   data%id_E_airtemp =  aed2_locate_global_sheet('air_temp')
   data%id_E_material = aed2_locate_global_sheet('material')
   data%id_E_rainloss = aed2_locate_global_sheet('rainloss')
   data%id_E_nearlevel =aed2_locate_global_sheet('nearest_depth')

END SUBROUTINE aed2_define_land
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_initialize_land(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
     AED_REAL :: rain, matz, bathy, avgLevel, airtemp
!-------------------------------------------------------------------------------
!BEGIN

   ! Set local cell properties
   rain = _STATE_VAR_S_(data%id_E_rain)
   matz = _STATE_VAR_S_(data%id_E_material)
   bathy = _STATE_VAR_S_(data%id_E_bathy)
   airtemp = _STATE_VAR_S_(data%id_E_airtemp)
   avgLevel = _STATE_VAR_S_(data%id_E_nearlevel)      !avgLevel = 0.6

   ! Update bucket settings
   CALL InitialHydrology(data, column, layer_idx, avgLevel)

   ! Link local updates back to global array
   _DIAG_VAR_S_(data%id_depth) = SoilCol%Depth
   _DIAG_VAR_S_(data%id_Sb) = SoilCol%Sb
   _DIAG_VAR_S_(data%id_St) = SoilCol%St
   _DIAG_VAR_S_(data%id_Ssat) = SoilCol%Ssat
   _DIAG_VAR_S_(data%id_Scap) = SoilCol%S_cap
   _DIAG_VAR_S_(data%id_Strn) = SoilCol%S_trn
   _DIAG_VAR_S_(data%id_Stop) = SoilCol%S_top
   _DIAG_VAR_S_(data%id_theta) = SoilCol%theta
   _DIAG_VAR_S_(data%id_Ucap) = SoilCol%CapHgt
   _DIAG_VAR_S_(data%id_phreatic) = SoilCol%PhreaticDepth
   _DIAG_VAR_S_(data%id_wt) = SoilCol%PhreaticHgt
   _DIAG_VAR_S_(data%id_dh) = SoilCol%PhreaticHgt - avgLevel
   _DIAG_VAR_S_(data%id_smd) = MAX((SoilCol%PhreaticDepth*SoilCol%porosity)-SoilCol%Sus,zero_)
   _DIAG_VAR_S_(data%id_qss) = SoilCol%qss
   _DIAG_VAR_S_(data%id_qs) = SoilCol%qse+SoilCol%qie
   _DIAG_VAR_S_(data%id_qper) = SoilCol%recharge
   _DIAG_VAR_S_(data%id_qinf) = zero_
   _DIAG_VAR_S_(data%id_tp) = zero_
   _DIAG_VAR_S_(data%id_qcap) = SoilCol%qsuc
   _DIAG_VAR_S_(data%id_area) = SoilCol%area
   _DIAG_VAR_S_(data%id_bathy) = bathy
   _DIAG_VAR_S_(data%id_rain) = rain
   _DIAG_VAR_S_(data%id_evap) = zero_
   _DIAG_VAR_S_(data%id_matz) = matz
   _DIAG_VAR_S_(data%id_rainloss) = MAX((SoilCol%PhreaticDepth*SoilCol%porosity)-SoilCol%Sus,zero_)
   _DIAG_VAR_S_(data%id_soiltemp_10cm) = airtemp


END SUBROUTINE aed2_initialize_land
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_dry_land(data, column, layer_idx)!
!-------------------------------------------------------------------------------
! Main interface routine for the hydrology/ass module. Assumes the module has
! already been configured (ConfigureHydrology) and initialised (InitialHydrology)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in)   :: data
   TYPE  (aed2_column_t),  INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!LOCALS
   AED_REAL :: bathy, matz, rainloss, airtemp, area
   AED_REAL :: rain, avgLevel
!-------------------------------------------------------------------------------

   ! Set local cell properties
   rain = _STATE_VAR_S_(data%id_E_rain)
   area = _STATE_VAR_S_(data%id_E_area)
   matz = _STATE_VAR_S_(data%id_E_material)
   bathy = _STATE_VAR_S_(data%id_E_bathy)
   airtemp = _STATE_VAR_S_(data%id_E_airtemp)
   avgLevel = _STATE_VAR_S_(data%id_E_nearlevel)

   ! Prime local cell object with data from global arrays
   CALL SetBucketParameters(data, SoilCol, INT(matz))

   IF(.NOT.in_zone_set(matz,data%active_zones)) RETURN

   ! Set bucket object form the global array
   SoilCol%Bathy = bathy
   SoilCol%Depth = _DIAG_VAR_S_(data%id_depth)
   SoilCol%Sb = _DIAG_VAR_S_(data%id_Sb)
   SoilCol%St = _DIAG_VAR_S_(data%id_St)
   SoilCol%Ssat = _DIAG_VAR_S_(data%id_Ssat)
   SoilCol%S_top = _DIAG_VAR_S_(data%id_Stop)
   SoilCol%CapHgt = _DIAG_VAR_S_(data%id_Ucap)
   SoilCol%PhreaticDepth = _DIAG_VAR_S_(data%id_phreatic)
   IF ( data%infil_model==HORTON )  SoilCol%inf%tp = _DIAG_VAR_S_(data%id_tp)

   ! Update the Soil Water Balance based on hydrological processes
   CALL UpdateSoilHydrology(data, SoilCol, rain, avgLevel, airtemp)

   ! runoff tracer into the water: (m rain)/m2/s
   _FLUX_VAR_R_(data%id_tracer) = (SoilCol%qse + SoilCol%qss + SoilCol%qie) / (0.25*3600)
   !IF(simErosion) &
   !_FLUX_VAR_R_(data%id_erosion) = ??

   ! Put the updated cell hydro variables back into the global arrays
   _DIAG_VAR_S_(data%id_Sb) = SoilCol%Sb
   _DIAG_VAR_S_(data%id_St) = SoilCol%St
   _DIAG_VAR_S_(data%id_Ssat) = SoilCol%Ssat
   _DIAG_VAR_S_(data%id_theta) = SoilCol%theta
   _DIAG_VAR_S_(data%id_Ucap) = SoilCol%CapHgt
   _DIAG_VAR_S_(data%id_Stop) = SoilCol%S_top
   _DIAG_VAR_S_(data%id_Scap) = SoilCol%S_cap
   _DIAG_VAR_S_(data%id_Strn) = SoilCol%S_trn
   _DIAG_VAR_S_(data%id_phreatic) = SoilCol%PhreaticDepth
   _DIAG_VAR_S_(data%id_depth) = SoilCol%Depth
   _DIAG_VAR_S_(data%id_wt) = SoilCol%PhreaticHgt
   _DIAG_VAR_S_(data%id_dh) = SoilCol%PhreaticHgt - avgLevel
   _DIAG_VAR_S_(data%id_smd) = MAX((SoilCol%PhreaticDepth*SoilCol%porosity)-SoilCol%Sus,zero_)
   _DIAG_VAR_S_(data%id_qss) = SoilCol%qss
   _DIAG_VAR_S_(data%id_qper) = SoilCol%recharge
   _DIAG_VAR_S_(data%id_qcap) = SoilCol%qsuc
   _DIAG_VAR_S_(data%id_qinf) = SoilCol%qinf
   _DIAG_VAR_S_(data%id_qs) = SoilCol%qse+SoilCol%qie
   _DIAG_VAR_S_(data%id_evap) = SoilCol%et
   _DIAG_VAR_S_(data%id_soiltemp_10cm) = airtemp
   rainloss = MAX(rain - (SoilCol%qse+SoilCol%qss+SoilCol%qie),zero_)
   _DIAG_VAR_S_(data%id_rainloss) = rainloss
   IF ( data%infil_model==HORTON )  _DIAG_VAR_S_(data%id_tp) = SoilCol%inf%tp
   ! Duplicating Env variables just for plotting
   _DIAG_VAR_S_(data%id_area) = avgLevel
   _DIAG_VAR_S_(data%id_bathy) = bathy
   _DIAG_VAR_S_(data%id_rain) = rain
   _DIAG_VAR_S_(data%id_matz) = matz


   ! If soil temperature them, update soil heating/cooling
   IF (data%simTemp) THEN
      airtemp = _STATE_VAR_S_(data%id_E_airtemp)
      !# CALL UpdateSoilTemperature(SoilCol, airtemp)
      _DIAG_VAR_S_(data%id_soiltemp_10cm) = 15. !SoilCol%Temp(10)
   END IF

   ! If soil salinity them, update soil solution salinities
   IF (data%simSalt) THEN
      !
   END IF

END SUBROUTINE aed2_calculate_dry_land
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_rain_loss_land(data,column,layer_idx,infil)
!-------------------------------------------------------------------------------
! Get the soil moisture deficit, so host model can infiltrate rain accordingly
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!LOCALS
   AED_REAL :: smd

!-------------------------------------------------------------------------------
!BEGIN
   infil = _DIAG_VAR_S_(data%id_rainloss)

END SUBROUTINE aed2_rain_loss_land
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE InitialHydrology(data, column, layer_idx, avgLevel)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in)    :: data
   TYPE  (aed2_column_t),  INTENT(inout) :: column(:)
   INTEGER,INTENT(in)  :: layer_idx
   AED_REAL,INTENT(in) :: avgLevel
  !LOGICAL,INTENT(in) :: CellIsDry
!LOCALS
   INTEGER :: i, sub,  n, chnk, mincell, sIndex,MatZoneID
   CHARACTER(LEN=4) :: filetype
   AED_REAL :: minbathy,area,bathy,initWTDepth
   AED_REAL :: matz
!-------------------------------------------------------------------------------
!BEGIN

   ! Assign global vals for this cell to the local variables

   ! -> Environment
   MatZoneID = INT(_STATE_VAR_S_(data%id_E_material))
   matz = _STATE_VAR_S_(data%id_E_material)
   bathy = _STATE_VAR_S_(data%id_E_bathy)
   area = _STATE_VAR_S_(data%id_E_area)

   ! -> Initial conditions (initialised diagnostics from file)
   SoilCol%theta = 0.2   !_DIAG_VAR_S_(data%id_theta)
   initWTDepth = 0.6 !MAX( _DIAG_VAR_S_(data%id_phreatic), zero_ )

   !---------------------------------------------------------------------------
   ! Setup jobs
   minbathy = data%minbathy(1)   !MH set by user, or computed as 2m below lowest cell in mesh
   SoilCol%nlay = 20

   IF(.NOT.ALLOCATED(SoilCol%Moisture)) ALLOCATE( SoilCol%Moisture(SoilCol%nlay) )
   IF(data%simTemp) THEN
      IF(.NOT.ALLOCATED(SoilCol%Temp)) ALLOCATE( SoilCol%Temp(SoilCol%nlay) )
    ENDIF

   ! Assign the soil parameters to this cell based on its MatZone
   CALL SetBucketParameters(data, SoilCol, MatZoneID)

   ! Cell geometry
   SoilCol%Area = area
   SoilCol%Bathy = bathy

   ! Initialise processes to 0
   SoilCol%qss            = zero_
   SoilCol%qse            = zero_
   SoilCol%recharge       = zero_
   SoilCol%qsuc           = zero_
   SoilCol%rn             = zero_
   SoilCol%et             = zero_
   SoilCol%es             = zero_

   !---------------------------------------------------------------------------
   ! Deep cells or cells in inactive zones - forget them

   IF ( bathy<minbathy .OR.  .NOT.in_zone_set(matz,data%active_zones) ) THEN
     SoilCol%Depth = zero_
     SoilCol%Sb = zero_
     SoilCol%Ssat  = zero_
     SoilCol%theta = SoilCol%porosity
     SoilCol%PhreaticDepth = zero_
     SoilCol%PhreaticHgt = bathy
     SoilCol%Sus = zero_
     RETURN
   ENDIF

   !---------------------------------------------------------------------------
   ! Potentially dry cells

   ! Depth
   SoilCol%Depth = MIN( bathy - minbathy , 10.)  !MH zASS needs to be read in
   SoilCol%pastMaxLevel = SoilCol%Bathy
   SoilCol%Debs = SoilCol%Depth - SoilCol%Debs
   IF(SoilCol%Debs > SoilCol%Depth) THEN
      SoilCol%Debs = SoilCol%Depth
   ENDIF

   ! Soil capacity and water table
   SoilCol%Sb = SoilCol%Depth * SoilCol%Porosity
   SoilCol%PhreaticHgt = SoilCol%Bathy - initWTDepth

   IF(bathy<avgLevel) THEN
     ! Cell is underwater
     SoilCol%PhreaticHgt = avgLevel
     SoilCol%PhreaticDepth = zero_
   ELSE
     IF(bathy > 1.0) THEN
       SoilCol%PhreaticHgt = 1.0
       SoilCol%PhreaticDepth = bathy - SoilCol%PhreaticHgt
     ELSE
       SoilCol%PhreaticHgt = MIN( bathy, 0.5)! avgLevel
       SoilCol%PhreaticDepth = bathy - SoilCol%PhreaticHgt
     ENDIF

!     print *,'yuppo'
!     ! Cell is out of water, check the water table
!     IF( bathy-initWTDepth > avgLevel ) THEN
!       SoilCol%PhreaticHgt = bathy - initWTDepth
!       SoilCol%PhreaticDepth = bathy - SoilCol%PhreaticHgt
!     ELSE
!       SoilCol%PhreaticDepth = bathy - avgLevel
!       SoilCol%PhreaticHgt = bathy - SoilCol%PhreaticDepth
!     ENDIF
   ENDIF


   SoilCol%Ssat  = (SoilCol%Depth-SoilCol%PhreaticDepth) &
                                                * SoilCol%Porosity

   !SoilCol%PhreaticHgt = SoilCol%Bathy+SoilCol%PhreaticDepth


    IF(SoilCol%PhreaticHgt+SoilCol%Zcap > SoilCol%Bathy) THEN
       ! water table is at or just below soil surface
       SoilCol%Dcap = SoilCol%Bathy - SoilCol%PhreaticHgt
       SoilCol%Dtrn = zero_
       SoilCol%CapHgt = SoilCol%Bathy
    ELSE
       ! water table is below surface by more than capillary rise height
       SoilCol%Dcap = SoilCol%Zcap
       IF( (SoilCol%PhreaticHgt + SoilCol%Dcap + SoilCol%Ztrn) > SoilCol%Bathy) THEN
         SoilCol%Dtrn = SoilCol%Bathy - (SoilCol%PhreaticHgt + SoilCol%Dcap)
       ELSE
         SoilCol%Dtrn = SoilCol%Ztrn
       ENDIF
       SoilCol%CapHgt = SoilCol%PhreaticHgt + SoilCol%Dtrn + SoilCol%Dcap
    ENDIF

    ! set an assumed unsaturated zone profile
    SoilCol%S_cap = 0.9 * SoilCol%Dcap * SoilCol%Porosity
    SoilCol%S_trn = 0.5 * SoilCol%Dtrn * SoilCol%Porosity
    SoilCol%S_top = 0.2 * (SoilCol%Bathy - SoilCol%CapHgt) * SoilCol%Porosity

    SoilCol%Sus = SoilCol%S_top + SoilCol%S_trn + SoilCol%S_cap
    IF( SoilCol%PhreaticDepth*SoilCol%Porosity > 1e-5 ) THEN
      SoilCol%theta = SoilCol%Sus / (SoilCol%PhreaticDepth*SoilCol%Porosity)
    ELSE
      SoilCol%theta = SoilCol%Porosity
    ENDIF
    !---------------------------------------------------------------------------

   SoilCol%St = SoilCol%Ssat + SoilCol%Sus
   IF( SoilCol%St  >  SoilCol%Sb )THEN
     print *,'St>Sb'
     SoilCol%St = SoilCol%Sb
     PAUSE
   ENDIF

   SoilCol%Moisture = SoilCol%theta

END SUBROUTINE InitialHydrology
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE UpdateSoilHydrology( data, theSoil, rain, avgLevel, air_temp )
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(IN)    :: data
   TYPE(SoilUnit), INTENT(INOUT) :: theSoil
   AED_REAL, INTENT(IN) :: rain, avgLevel, air_temp
!LOCALS
   INTEGER  :: lay
   AED_REAL :: flow, evap, evapsat, rainp, evapp, infil,  ponding
   AED_REAL :: Sttm1, Ssattm1, Susfc, Susp
   AED_REAL :: dep, depm1, middep, botmoist, topmoist, top_depth
!-------------------------------------------------------------------------------
!BEGIN

   flow = zero_
   evap = zero_
   evapsat = zero_
   infil = zero_
   ponding = 0.0

   rainp = rain * DDT                           ! m of rain / timestep
   evapp = GetPotlEvap( rain,air_temp ) * DDT   !* GetPotlEvap(jm,jd)


!print *,'rainp',rainp,evapp

   ! Set prev timestep vars
   Sttm1   =   theSoil%St
   Ssattm1 =   theSoil%Ssat

   theSoil%Sus = theSoil%St - theSoil%Ssat
   theSoil%PhreaticHgt = theSoil%Bathy - theSoil%PhreaticDepth

!print *,'Sus',theSoil%St, theSoil%Ssat

   ! Infiltration amount
   infil = GetInfiltration( data, theSoil, rain , ponding ) * DDT
  !  infil = MIN( rain, theSoil%phi ) * DDT

   !theSoil%qie = MAX( rainp - infil, zero_ )
   IF( rainp>infil ) THEN ! both in m/ts now
       theSoil%qie  = rainp - infil
    ELSE
       theSoil%qie = zero_
   ENDIF


!print *,'infil',infil, theSoil%qie

   ! Unsaturated Zone capacity
   !   Susfc = (theSoil%Sb - Ssattm1 - theSoil%S_trn - theSoil%S_cap ) * theSoil%fc
   !   Susp  = theSoil%S_top + rain
   top_depth =  ( theSoil%Bathy - theSoil%CapHgt )
   Susfc = top_depth * theSoil%fc
   Susp  = theSoil%Sus + infil


!print *,'Susfc',Susp, top_depth, Susfc

   ! Check for saturation excess runoff
   !IF(Susp > top_depth * theSoil%Porosity) THEN
   IF((theSoil%Ssat+Susp) > theSoil%Sb) THEN
      theSoil%qse = (theSoil%Ssat+Susp) - theSoil%Sb ! - (0.9 * (top_depth * theSoil%Porosity))
      infil = infil - theSoil%qse
      Susp  = theSoil%Sus + infil !theSoil%S_top + infil
      theSoil%qinf = infil
   ELSE
      theSoil%qinf = infil !zero_
      theSoil%qse = 0.0
   ENDIF
!print *,'Susp',Susp, theSoil%qse, infil

   ! Recharge/percolation
   IF(Susp > Susfc) THEN
      theSoil%recharge = 0.5 * (Susp - Susfc)
   ELSE
      theSoil%recharge = 0.0
   ENDIF
   !print *,'ss: ',Sttm1, Ssattm1, Susfc, Susp, theSoil%recharge
!print *,'recharge',theSoil%recharge

   ! Evaporation
   theSoil%Debs = MAX( theSoil%Depth-theSoil%Debs, zero_ )

!print *,'Debs',theSoil%Debs,theSoil%Depth

   IF(Ssattm1 <=  theSoil%Porosity*theSoil%Debs) THEN
      ! Water table below reach, only take from unsaturated volume

      !evap = theSoil%aep * theSoil%Sus / ( theSoil%Sb - Ssattm1 )
      evap = theSoil%aep * MIN( theSoil%Sus/MAX(theSoil%Sb - Ssattm1,1e-2),1.0 )
      evap = evap * evapp
      evapsat = 0.0

   ELSE
      ! Take from unsaturated and saturated regions
      evap = theSoil%aep * (theSoil%Sus + (Ssattm1 - theSoil%Porosity*theSoil%Debs)) &
                / MIN( theSoil%Sb - theSoil%Porosity*theSoil%Debs, 1e-2 )
      evap = MIN(evap,1.0) * evapp
      evapsat = (Ssattm1 - theSoil%Porosity*theSoil%Debs) &
                   / (theSoil%Sus + (Ssattm1-theSoil%Porosity*theSoil%Debs) )
      evapsat = MIN(evapsat,1.0) * evap

  ENDIF

!print *,'evap',evap,evapsat

!evap = 0.5*evapp
!evapsat = 0.5*evapp

   ! intermediate stage estimate of top region of unsat zone
   Susp = Susp - theSoil%recharge - (evap - evapsat) !* top_depth/theSoil%PhreaticDepth

   ! Seepage
   IF ( (theSoil%PhreaticHgt) > avgLevel ) THEN
      IF(avgLevel >= theSoil%Bathy) THEN
         theSoil%qss = 0.
      ELSE
         theSoil%qss = ( theSoil%PhreaticHgt - avgLevel ) / ( theSoil%Bathy - avgLevel  )
         theSoil%qss = theSoil%ass * theSoil%qss**theSoil%Bss
      ENDIF
   ELSE
      ! Backflow
      !theSoil%qss = 0.0
      theSoil%qss = ( avgLevel - theSoil%PhreaticHgt ) / ( avgLevel - (theSoil%Bathy-theSoil%Depth)  )
      theSoil%qss = - theSoil%ass * theSoil%qss**theSoil%Bss
   ENDIF

   flow =  theSoil%qss
!print *,'flow',flow

   !print *,'ee: ',evap, evapsat, theSoil%recharge, theSoil%qss

   ! Update bucket water balance for timestep
   theSoil%St = Sttm1 + infil - flow - evap


   ! Suction
   !IF( Susp < zero_ .AND. top_depth > 0.01 ) THEN
   !   theSoil%qsuc = -1.05*Susp
   !ELSE
   !   theSoil%qsuc = zero_
   !   evapsat = evapsat - Susp
   !END IF
      IF( theSoil%Sus <= zero_ .AND. theSoil%St < theSoil%Sb ) THEN
        theSoil%qsuc = MIN (theSoil%Depth-(theSoil%St/theSoil%Porosity),theSoil%Zcap) *0.9*theSoil%Porosity
      END IF




!print *,'St',theSoil%St , Sttm1

   IF ( theSoil%St <= zero_ ) THEN
      ! Empty
      theSoil%St = 0.0
      theSoil%Ssat = 0.0
      theSoil%recharge = 0.0
      evapsat = 0.0
      flow = 0.0
      theSoil%qse = 0.0
      theSoil%Moisture(1:theSoil%nlay) = 0.0
      theSoil%PhreaticDepth   = theSoil%Depth
      theSoil%PhreaticHgt     = theSoil%Bathy -  theSoil%Depth
      theSoil%Dcap = 0.0
      theSoil%Dtrn = 0.0
      theSoil%S_cap  = 0.0
      theSoil%S_trn  = 0.0
      theSoil%S_top  = 0.0
      theSoil%theta  = 0.0
      theSoil%Sus = 0.0
   ELSEIF( theSoil%St >= theSoil%Sb ) THEN
      ! Full
      theSoil%qse = theSoil%St - theSoil%Sb
      theSoil%St  = theSoil%Sb
      theSoil%Sus = 0.0
      theSoil%Ssat = theSoil%St
      theSoil%Moisture(1:theSoil%nlay) = 1.0 * theSoil%Porosity * 1000. &
                                          /(theSoil%Porosity * 1000. + (1.0-theSoil%Porosity)*theSoil%Density)
      theSoil%recharge = theSoil%Ssat - Ssattm1
      theSoil%PhreaticDepth = zero_
      theSoil%PhreaticHgt =  theSoil%Bathy - theSoil%PhreaticDepth
      theSoil%Dcap = 0.0
      theSoil%Dtrn = 0.0
      theSoil%S_cap = 0.0
      theSoil%S_trn = 0.0
      theSoil%S_top = 0.0
      theSoil%theta = theSoil%Porosity
   ELSE
      ! Normal
      theSoil%Ssat = Ssattm1 + theSoil%recharge - evapsat - flow - theSoil%qsuc
      theSoil%Sus  = theSoil%St - theSoil%Ssat !+ theSoil%qsuc


      theSoil%PhreaticDepth = (theSoil%Sb - theSoil%Ssat) / theSoil%Porosity
      theSoil%PhreaticHgt =  theSoil%Bathy - theSoil%PhreaticDepth

!print *,'Ssat',theSoil%Ssat, theSoil%Sus

      ! Update unsat zone levels
      IF(theSoil%PhreaticHgt+theSoil%Zcap > SoilCol%Bathy) THEN
       theSoil%Dcap = theSoil%Bathy - theSoil%PhreaticHgt
       theSoil%Dtrn = zero_
       theSoil%CapHgt = theSoil%Bathy
      ELSE
       theSoil%Dcap = theSoil%Zcap
       IF( (theSoil%PhreaticHgt + theSoil%Dcap + theSoil%Ztrn) > theSoil%Bathy) THEN
         theSoil%Dtrn = theSoil%Bathy - (theSoil%PhreaticHgt + theSoil%Dcap)
       ELSE
         theSoil%Dtrn = theSoil%Ztrn
       ENDIF
       theSoil%CapHgt = theSoil%PhreaticHgt + theSoil%Dtrn + theSoil%Dcap
      ENDIF


      !MAKE THIS A FUNCTION SO ASS CAN USE IT
      ! Calculate unsat zone storage values (disaggregation of bulk value)
      botmoist = 0.9
      IF( ABS(theSoil%Bathy-theSoil%CapHgt) <1e-3) THEN
           topmoist = 0.9
      ELSE
           topmoist = MAX(Susp,0.001) /  (( theSoil%Bathy-theSoil%CapHgt ) * theSoil%Porosity)
      END IF

      theSoil%S_cap = 0.9 * theSoil%Dcap * theSoil%Porosity
      IF( theSoil%S_cap > theSoil%Sus ) THEN
        theSoil%S_top = zero_
        theSoil%S_trn = zero_
        theSoil%S_cap = theSoil%Sus
        theSoil%Dcap = theSoil%S_cap * theSoil%Porosity * 1.05
        theSoil%Dtrn = zero_
      ELSE
        theSoil%S_trn = ((botmoist-topmoist)/2.) * (theSoil%Dtrn) * theSoil%Porosity
        IF( theSoil%S_trn+theSoil%S_cap > theSoil%Sus ) THEN
          theSoil%S_trn = theSoil%Sus - theSoil%S_cap
          theSoil%S_top = zero_
          theSoil%Dtrn = theSoil%S_trn * theSoil%Porosity * 2.
        ELSE
        theSoil%S_top = theSoil%Sus - theSoil%S_trn - theSoil%S_cap
        ENDIF
      ENDIF
      topmoist = (theSoil%S_top) /  MAX( (theSoil%Bathy-theSoil%CapHgt) * theSoil%Porosity,0.001)
      theSoil%CapHgt = theSoil%PhreaticHgt + theSoil%Dtrn + theSoil%Dcap
      IF(theSoil%PhreaticDepth*theSoil%Porosity >1e-10) THEN
        theSoil%theta = theSoil%Sus / (theSoil%PhreaticDepth*theSoil%Porosity)
      ELSE
        theSoil%theta = theSoil%Porosity
      ENDIF
   ENDIF

   theSoil%rn = rainp
   theSoil%et = evap
   theSoil%es = evapsat

END SUBROUTINE UpdateSoilHydrology
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_land(data, column, layer_idx)
!-------------------------------------------------------------------------------
!Force values of hydrology variables to 0 when too much water over the cell
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!LOCALS
   AED_REAL :: bathy, matz, rainloss, airtemp, area
   AED_REAL :: rain, avgLevel, temp, rainp
!-------------------------------------------------------------------------------

   ! Set local cell properties
   temp = _STATE_VAR_(data%id_E_temp)
   rain = _STATE_VAR_S_(data%id_E_rain)
   area = _STATE_VAR_S_(data%id_E_area)
   matz = _STATE_VAR_S_(data%id_E_material)
   bathy = _STATE_VAR_S_(data%id_E_bathy)
   avgLevel = _STATE_VAR_S_(data%id_E_nearlevel)

   ! Prime local cell object with data from global arrays
   CALL SetBucketParameters(data, SoilCol, INT(matz))

!! IF(.NOT.in_zone_set(matz,data%active_zones)) RETURN

   SoilCol%Depth = MIN( bathy - data%minbathy(1) , 10.)  !MH zASS needs to be read in

      ! Put the updated cell hydro variables back into the global arrays
   _DIAG_VAR_S_(data%id_Sb) = SoilCol%Depth * SoilCol%Porosity
   _DIAG_VAR_S_(data%id_St) = SoilCol%Depth * SoilCol%Porosity
   _DIAG_VAR_S_(data%id_Ssat) = SoilCol%Depth * SoilCol%Porosity
   _DIAG_VAR_S_(data%id_theta) = SoilCol%Porosity
   _DIAG_VAR_S_(data%id_Ucap) = bathy
   _DIAG_VAR_S_(data%id_Stop) = zero_
   _DIAG_VAR_S_(data%id_Scap) = zero_
   _DIAG_VAR_S_(data%id_Strn) = zero_
   _DIAG_VAR_S_(data%id_phreatic) = zero_
   _DIAG_VAR_S_(data%id_wt) = avgLevel
   _DIAG_VAR_S_(data%id_dh) = zero_
   _DIAG_VAR_S_(data%id_depth) = SoilCol%Depth
   _DIAG_VAR_S_(data%id_smd) = zero_
   _DIAG_VAR_S_(data%id_qss) = zero_
   _DIAG_VAR_S_(data%id_qs) = zero_
   _DIAG_VAR_S_(data%id_qper) = zero_
   _DIAG_VAR_S_(data%id_qcap) = zero_
   _DIAG_VAR_S_(data%id_qinf) = zero_
   _DIAG_VAR_S_(data%id_rainloss) = zero_
   _DIAG_VAR_S_(data%id_soiltemp_10cm) = temp
   ! Duplicating Env variables just for plotting
   _DIAG_VAR_S_(data%id_area) = avgLevel
   _DIAG_VAR_S_(data%id_bathy) = bathy
   _DIAG_VAR_S_(data%id_rain) = rain
   _DIAG_VAR_S_(data%id_evap) = zero_
   _DIAG_VAR_S_(data%id_matz) = zero_  !matz - ! temp changed to detrmine status

END SUBROUTINE aed2_calculate_benthic_land
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE SetBucketParameters(data, theSoil, MatZoneID)
!-------------------------------------------------------------------------------
! Primes the local cell Bucket model "object" with relevant parameters based on!
! its material zone id number                                                  !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in) :: data
   TYPE(SoilUnit),INTENT(inout) :: theSoil ! passing SoilUnit into theSoil
   INTEGER,INTENT(in) :: MatZoneID
!LOCALS
   INTEGER :: i, zindex
!-------------------------------------------------------------------------------
!BEGIN

   theSoil%UnitNum = 0
   theSoil%Substrate = MatZoneID

   zindex = 1
   DO i =1,data%n_zones
     IF( MatZoneID == INT(data%active_zones(i)) )THEN
       zindex = i
     ENDIF
   ENDDO

   ! Hydrologic dynamics (initialise here, they'll be updated with state vars later)
   theSoil%St = ZERO_
   theSoil%Ssat = ZERO_
   theSoil%PhreaticHgt = ZERO_
   theSoil%PhreaticDepth = ZERO_
   theSoil%qss = ZERO_
   theSoil%qse = ZERO_
   theSoil%qie = ZERO_
   theSoil%recharge = ZERO_
   theSoil%qsuc = ZERO_
   theSoil%Sus = ZERO_

   ! Physical properties
   theSoil%Depth = ZERO_
   theSoil%Area = ZERO_
   theSoil%Bathy = ZERO_
   theSoil%Sb = ZERO_
   theSoil%Debs = data%Zebs(zindex)
   theSoil%Dper = data%Dper(zindex)

   ! Hydrologic parameters
   theSoil%aep = data%aep(zindex)
   theSoil%ass = data%ass(zindex)
   theSoil%Bss = data%Bss(zindex)

   ! Soil parameters
   theSoil%Porosity = data%Porosity(zindex)
   theSoil%Density = data%Density(zindex)
   theSoil%fc = data%fc(zindex)

   ! Moisture disaggregation model
   theSoil%nlay = data%nlay
   theSoil%Dcap = data%Zcap(zindex)
   theSoil%Zcap = data%Zcap(zindex)
   theSoil%Dtrn = data%Ztrn(zindex)
   theSoil%Ztrn = data%Ztrn(zindex)

   ! Infiltration sub-model
   IF( data%infil_model==1 ) THEN
     theSoil%phi = data%phi(zindex)
   ELSEIF( data%infil_model==HORTON ) THEN
     theSoil%inf%f0 = data%f0(zindex)
     theSoil%inf%fmin = data%fmin(zindex)
     theSoil%inf%Fmax = data%Fmax(zindex)
     theSoil%inf%decay = data%decay(zindex)
     theSoil%inf%regen = data%regen(zindex)
   ELSEIF( data%infil_model==GRNAMPT ) THEN

   ENDIF

END SUBROUTINE SetBucketParameters
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION GetInfiltration( data, theSoil, RainIntensity , PondingDepth ) RESULT (infiltration)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_land_data_t),INTENT(in) :: data
   TYPE(SoilUnit),INTENT(inout) :: theSoil ! passing SoilUnit into theSoil
   AED_REAL,INTENT(in) ::  RainIntensity, PondingDepth
!LOCALS
   INTEGER, PARAMETER :: infil_model = 1
   AED_REAL :: infiltration
!BEGIN
!-------------------------------------------------------------------------------
  infiltration = zero_

  IF( data%infil_model==1 ) THEN

     ! PHI index model - all rain upto phi m/d can infiltrate
     infiltration = MIN( RainIntensity, theSoil%phi )

  ELSEIF( data%infil_model==HORTON ) THEN

     ! HORTON model
     infiltration = horton_getInfil( theSoil%inf, RainIntensity, PondingDepth, theSoil%theta, theSoil%et )

  ELSEIF( data%infil_model==GRNAMPT ) THEN

     ! GREEN-AMPT model


  ENDIF

END FUNCTION GetInfiltration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION horton_getInfil( infil, RainIntensity, PondingDepth, moist , et) RESULT (CumInf)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(Infiltration),INTENT(inout) :: infil !  TYPE Infiltration
   AED_REAL,INTENT(in) :: RainIntensity ! the depth falling over the cell in the current time step [m/time step]
   AED_REAL,INTENT(in) :: PondingDepth  ! the depth stading over the cell available within the current time step [m/time step]
!LOCALS
    INTEGER  ::  iter         ! counter
    AED_REAL ::  fa, fp = 0.0;
    AED_REAL ::  FFp, F1, t1, tlim, ex, kt,kreg
    AED_REAL ::  FF, FF1, r
    AED_REAL ::  f0           !  initial infil. rate [m/d]
    AED_REAL ::  fmin         !  minimum infil. rate [m/d]
    AED_REAL ::  Fmax         !  maximum total infiltration [m];
    AED_REAL ::  tp           !  present time on infiltration curve [m]
    AED_REAL ::  df           !  f0 - fmin
    AED_REAL ::  kd           !  decay coeff. of infil. rate [1/d]
    AED_REAL ::  kr
    AED_REAL :: CumInf
    AED_REAL :: moist, et
	AED_REAL :: accumul_depth
!BEGIN
!-------------------------------------------------------------------------------

        kr   = infil%regen  !* Evap.recoveryFactor; [1/d] depends on evaporation, which is calculate above
        kd   = infil%decay  !  decay coeff. of infil. rate [1/d]
		moist=max(0.1,min(moist, 0.9))
        f0   = (infil%f0 * moist) ** (0.2)   !
        fmin = infil%fmin   !* Adjust.hydconFactor;
        Fmax = infil%Fmax
        CumInf=infil%Fe
        kreg= kr* et;       ! kreg(i)=kr*evaporation(i)*evaporation_factor; to be linked to evaporation
        tp   = infil%tp     ! tp= present time on infiltration curve [d] passed from previous timestep
        df   = f0 - fmin
        fa = RainIntensity + PondingDepth     ! water available for infiltration [m]
        tlim = 1.0 / kd      ! for tp >= tlim, fp = fmin -  kd and time make the units consistent tlim in days!
		t1 = tp + DDT        ! future cumul. time


       IF  (CumInf <=0.0) THEN
	  accumul_depth=0.0
      ELSE
	  accumul_depth=fa+ accumul_depth
      ENDIF

        IF (df < 0.0 .or. kd < 0.0 .or. kr < 0.0)  THEN  ! no infil. or constant infil
               CumInf=0.0
         return
        ENDIF

        IF ( df == 0.0 .or. kd == 0.0 )    THEN
             fp = min(f0 ,fa)      ! actual infiltration rate is equal to the initial, but constrained by availability
             CumInf=MAX(0.0, fp)   ! fp=average infil. rate over time step
             return
        ENDIF

        IF ( fa > 0.0 )	THEN      ! if there is water to infiltrate, compute average infil. rate over time step

                IF ( tp >= tlim ) THEN  !   if tp, the present time on infiltration curve, is greater than the time limit, tlim
                 FFp = max(0.0, fmin * tp + df / kd) ![m] and df=f0- fmin
                 F1 = max(0.0, FFp + fmin * DDT) ! [m]=[m]+[m/d]*[d/time step]
                ELSE !% if tp, the present time on infiltration curve, is lesser than the time limit, tlim, exponential part of the graph
                 FFp = max(0.0,fmin * tp + df / kd * (1.0 - exp(-kd * tp)))
                 F1 = max(0.0,fmin * t1 + df / kd * (1.0 - exp(-kd * t1)))
                ENDIF

              fp = max(0.0,min((F1 - FFp) / DDT,fa)) ! average infil. rate over time step limited by fa

             IF (t1 > tlim .or.  fp < fa )   THEN         !  if fp on flat portion of curve
                  tp = t1                       !  then increase tp by the time step
             ELSE                             !  if infil limited by available capacity then
                 tp = tp + DDT / 2.0
                 F1 = max(0.0, FFp + fp * DDT)


                   DO  iter=1,20                   !the  solve F(tp) - F1 = 0 using Newton-Raphson method
                    kt = MIN ( 60.0, kd*tp )
                    ex = exp(-kt)
                    FF = fmin * tp + df / kd * (1.0 - ex) - F1
                    FF1 = fmin + df * ex
                    r = FF / FF1
                    tp = tp - r
                    IF ( abs(r) <= 0.001 * DDT ) EXIT
                   END DO
            ENDIF

         ELSE !% if fa <=0

		         fp=0.0;
                 CumInf=CumInf;
                 tp=max(0.0,1/kd*log(1-exp(-kreg*DDT)*(1-exp(-kd*tp))));
                 F1=0.0;
                 FFp=0.0;
         ENDIF



                 IF ( CumInf + fp * DDT > Fmax )  THEN! --- limit cumulative infiltration to Fmax
                   fp =0.0
				   CumInf=CumInf

				  ELSE
				  fp = max(0.0, min(fa, fmin) )  !(Fmax - CumInf) / DDT
                  CumInf =CumInf + fp * DDT      !FE=Fe +fp*DDT
                 ENDIF



        IF (kr > 0.0) THEN ! --- case where infil. capacity is regenerating;
           r = exp(-kreg * DDT)
           tp = 1.0 - exp(-kd * tp)
           tp = -log(1.0 - r*tp) / kd    !update tp.
        ENDIF

         CumInf =  max(CumInf, (df/kd)*(1.0 - exp(-kd*tp)))!!! reduction in cumulative infiltration

		infil%tp = tp

END FUNCTION horton_getInfil
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION GetPotlEvap( rain,air_temp ) RESULT (evaporation)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: rain
   AED_REAL,INTENT(in) :: air_temp
!LOCALS
   AED_REAL :: evaporation  ! m/day
!BEGIN
!-------------------------------------------------------------------------------
   !evaporation = 1. + epamp*sin(angfreq*((30.44*REAL(jm-1)+jd)+epphase))
   !evaporation = epdayav * evaporation

   evaporation = zero_

   ! if rain is less than 5mm/day then compute potential evap based on simple
   ! correlation with air temp (hipsey rough calc)
   if ( rain < 0.005 ) evaporation = air_temp / (3.5*1e3)

END FUNCTION GetPotlEvap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed2_land
