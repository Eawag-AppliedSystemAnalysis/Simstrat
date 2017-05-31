!###############################################################################
!#                                                                             #
!# aed2_macrophyte.F90                                                         #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!#                                                                             #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created July 2015                                                           #
!#                                                                             #
!###############################################################################

#include "aed2.h"


MODULE aed2_macrophyte
!-------------------------------------------------------------------------------
!  aed2_macrophyte --- multi-group macrophyte (/seagrass) model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util
   USE aed2_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_macrophyte_data_t

   TYPE :: macrophyte_data
       CHARACTER(64) :: m_name
       AED_REAL      :: m0
       AED_REAL      :: R_growth
       INTEGER       :: fT_Method
       AED_REAL      :: theta_growth
       AED_REAL      :: T_std
       AED_REAL      :: T_opt
       AED_REAL      :: T_max
       INTEGER       :: lightModel
       AED_REAL      :: I_K
       AED_REAL      :: I_S
       AED_REAL      :: KeMAC
       AED_REAL      :: f_pr
       AED_REAL      :: R_resp
       AED_REAL      :: theta_resp
       INTEGER       :: salTol
       AED_REAL      :: S_bep
       AED_REAL      :: S_maxsp
       AED_REAL      :: S_opt
       AED_REAL      :: K_CD
       AED_REAL      :: f_bg
       AED_REAL      :: k_omega
       AED_REAL      :: Xcc
       AED_REAL      :: K_N
       AED_REAL      :: X_ncon
       AED_REAL      :: K_P
       AED_REAL      :: X_pcon
   END TYPE

   TYPE,extends(aed2_model_data_t) :: aed2_macrophyte_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_mphy(:)
      INTEGER :: id_par, id_tem, id_sal, id_dz, id_extc, id_I_0
      INTEGER :: id_diag_par, id_gpp, id_p2r, id_mac, id_sed_zone
      INTEGER :: id_MAC_ag, id_MAC_bg, id_lai
      INTEGER :: n_zones
      AED_REAL,ALLOCATABLE :: active_zones(:)

      !# Model parameters
      INTEGER  :: num_mphy
      TYPE(macrophyte_data),DIMENSION(:),ALLOCATABLE :: mphydata
      LOGICAL  :: simMacFeedback, simStaticBiomass

     CONTAINS
         PROCEDURE :: define             => aed2_define_macrophyte
!        PROCEDURE :: calculate_riparian => aed2_calculate_riparian_macrophyte
         PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_macrophyte
!        PROCEDURE :: mobility           => aed2_mobility_macrophyte
!        PROCEDURE :: light_extinction   => aed2_light_extinction_macrophyte
!        PROCEDURE :: delete             => aed2_delete_macrophyte
         PROCEDURE :: bio_drag           => aed2_bio_drag_macrophyte
   END TYPE


   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================




!###############################################################################
SUBROUTINE aed2_macrophyte_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macrophyte_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(macrophyte_data) :: md(MAX_PHYTO_TYPES)
   NAMELIST /macrophyte_data/ md
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD', iostat=status)
    IF (status /= 0) STOP 'Cannot open macrophyte_data namelist file for macrophytes'
    read(tfil,nml=macrophyte_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist macrophyte_data for macrophytes'

    data%num_mphy = count
    ALLOCATE(data%mphydata(count))

    ALLOCATE(data%id_mphy(count))

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%mphydata(i)%m_name       = md(list(i))%m_name
       data%mphydata(i)%m0           = md(list(i))%m0
       data%mphydata(i)%R_growth     = md(list(i))%R_growth/secs_per_day
       data%mphydata(i)%fT_Method    = md(list(i))%fT_Method
       data%mphydata(i)%theta_growth = md(list(i))%theta_growth
       data%mphydata(i)%T_std        = md(list(i))%T_std
       data%mphydata(i)%T_opt        = md(list(i))%T_opt
       data%mphydata(i)%T_max        = md(list(i))%T_max
       data%mphydata(i)%lightModel   = md(list(i))%lightModel
       data%mphydata(i)%I_K          = md(list(i))%I_K
       data%mphydata(i)%I_S          = md(list(i))%I_S
       data%mphydata(i)%KeMAC        = md(list(i))%KeMAC
       data%mphydata(i)%f_pr         = md(list(i))%f_pr
       data%mphydata(i)%R_resp       = md(list(i))%R_resp/secs_per_day
       data%mphydata(i)%theta_resp   = md(list(i))%theta_resp
       data%mphydata(i)%salTol       = md(list(i))%salTol
       data%mphydata(i)%S_bep        = md(list(i))%S_bep
       data%mphydata(i)%S_maxsp      = md(list(i))%S_maxsp
       data%mphydata(i)%S_opt        = md(list(i))%S_opt
       data%mphydata(i)%K_CD         = md(list(i))%K_CD
       data%mphydata(i)%f_bg         = md(list(i))%f_bg
       data%mphydata(i)%k_omega      = md(list(i))%k_omega
       data%mphydata(i)%Xcc          = md(list(i))%Xcc
       data%mphydata(i)%K_N          = md(list(i))%K_N
       data%mphydata(i)%X_ncon       = md(list(i))%X_ncon
       data%mphydata(i)%K_P          = md(list(i))%K_P
       data%mphydata(i)%X_pcon       = md(list(i))%X_pcon

       ! Register group as a state variable
       data%id_mphy(i) = aed2_define_sheet_variable(                   &
                              md(list(i))%m_name,                      &
                              'mmolC/m**2', 'macrophyte',              &
                              md(list(i))%m0,                          &
                              minimum=zero_)

    ENDDO
END SUBROUTINE aed2_macrophyte_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_macrophyte(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the macrophyte/seagrass  model
!
!  Here, the aed2_ namelist is read and the variables exported
!  by the model are registered
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macrophyte_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_mphy
   INTEGER  :: the_mphy(MAX_PHYTO_TYPES)
   CHARACTER(len=128) :: dbase='aed2_macrophyte_pars.nml'
   INTEGER  :: n_zones = 0, active_zones(MAX_ZONES), i
   LOGICAL  :: simMacFeedback, simStaticBiomass

   NAMELIST /aed2_macrophyte/ num_mphy, the_mphy, dbase, n_zones, active_zones, &
                              simMacFeedback, simStaticBiomass

!-----------------------------------------------------------------------
!BEGIN

   simMacFeedback = .FALSE.
   simStaticBiomass = .FALSE.

   ! Read the namelist
   read(namlst,nml=aed2_macrophyte,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_macrophyte'

   data%simMacFeedback = simMacFeedback
   data%simStaticBiomass = simStaticBiomass

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   ! Macrophyte state variable allocated in here
   CALL aed2_macrophyte_load_params(data, dbase, num_mphy, the_mphy)

  ! CALL aed2_bio_temp_function(data%num_phytos,              &
  !                             data%phytos%theta_growth,     &
  !                             data%phytos%T_std,            &
  !                             data%phytos%T_opt,            &
  !                             data%phytos%T_max,            &
  !                             data%phytos%aTn,              &
  !                             data%phytos%bTn,              &
  !                             data%phytos%kTn,              &
  !                             data%phytos%p_name)


   ! Register diagnostic variables
   data%id_diag_PAR = aed2_define_sheet_diag_variable('par','W/m**2','benthic light intensity')
   data%id_GPP = aed2_define_sheet_diag_variable('gpp','/d',  'benthic plant productivity')
   data%id_P2R = aed2_define_sheet_diag_variable('p_r','-',  'macrophyte p:r ratio')
   data%id_MAC = aed2_define_sheet_diag_variable('mac','mmolC/m2',  'total macrophyte biomass')
   data%id_LAI = aed2_define_sheet_diag_variable('lai','m2/m2',  'macrophyte leaf area density')
   data%id_MAC_ag = aed2_define_sheet_diag_variable('mac_ag','mmolC/m2',  'total above ground macrophyte biomass')
   data%id_MAC_bg = aed2_define_sheet_diag_variable('mac_bg','mmolC/m2',  'total below ground macrophyte biomass')

   ! Register environmental dependencies
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_par = aed2_locate_global('par')
   data%id_I_0 = aed2_locate_global_sheet('par_sf')
   data%id_dz = aed2_locate_global('layer_ht')
   data%id_extc = aed2_locate_global('extc_coef')
   data%id_sed_zone = aed2_locate_global_sheet('sed_zone')

END SUBROUTINE aed2_define_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if 0
!# Already in aed2_util
!###############################################################################
LOGICAL FUNCTION in_zone_set(matz, active_zones)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: matz
   AED_REAL,INTENT(in) :: active_zones(:)
!
!LOCALS
   INTEGER :: i, l
   LOGICAL :: res
!BEGIN
!-------------------------------------------------------------------------------
   res = .FALSE.
   l = size(active_zones)
   DO i=1,l
      IF ( active_zones(i) == matz ) THEN
         res = .TRUE.
         EXIT
      ENDIF
   ENDDO

   in_zone_set = res
END FUNCTION in_zone_set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

!###############################################################################
SUBROUTINE aed2_calculate_benthic_macrophyte(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of phytoplankton.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: mphy        ! State
   INTEGER  :: mphy_i
   AED_REAL :: mphy_flux
   AED_REAL :: fT, fI, fSal, fDO
   AED_REAL :: extc, dz, par, Io, temp, salinity
   AED_REAL :: primprod(data%num_mphy)
   AED_REAL :: respiration(data%num_mphy)
   AED_REAL :: matz
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Retrieve current environmental conditions
   temp = _STATE_VAR_(data%id_tem)      ! local temperature
   salinity = _STATE_VAR_(data%id_sal)  ! local salinity
   par = _STATE_VAR_(data%id_par)       ! local photosynthetically active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation

   ! Initialise cumulative biomass diagnostics
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_
   _DIAG_VAR_S_(data%id_lai) = zero_

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i))! macrophyte group i

      ! LIGHT
      extc = _STATE_VAR_(data%id_extc)
      dz   = _STATE_VAR_(data%id_dz)     ! dz = 0.5
      fI   = photosynthesis_irradiance(data%mphydata(mphy_i)%lightModel, &
                     data%mphydata(mphy_i)%I_K, data%mphydata(mphy_i)%I_S, par, extc, Io, dz)
      fT   = 1.

      primprod(mphy_i) = data%mphydata(mphy_i)%R_growth * fI * fT

      ! Respiration and general metabolic loss
      respiration(mphy_i) = bio_respiration(data%mphydata(mphy_i)%R_resp, data%mphydata(mphy_i)%theta_resp, temp)

      ! Salinity stress effect on respiration
      fSal = 1.0 ! phyto_salinity(data%mphydata,mphy_i,salinity)
      fDO  = 1.0 ! phyto_oxygen(data%mphydata,mphy_i,oxygen)

      respiration(mphy_i) = respiration(mphy_i) * fSal * fDO

      IF( .NOT.data%simStaticBiomass ) THEN
        mphy_flux = (primprod(mphy_i) - respiration(mphy_i)) *  mphy

        ! Set bottom fluxes for the pelagic (change per surface area per second)
        _FLUX_VAR_B_(data%id_mphy(mphy_i)) = _FLUX_VAR_B_(data%id_mphy(mphy_i)) + mphy_flux
      ENDIF
      IF( data%simMacFeedback ) THEN
    !    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + mphy_flux
    !    _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - mphy_flux
      ENDIF

     _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + mphy
     _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + mphy*(one_-data%mphydata(mphy_i)%f_bg)
     _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + mphy*(data%mphydata(mphy_i)%f_bg)
     _DIAG_VAR_S_(data%id_lai) = _DIAG_VAR_S_(data%id_lai) + (one_ - exp(-data%mphydata(mphy_i)%k_omega * mphy*(one_-data%mphydata(mphy_i)%f_bg)))
   ENDDO

   ! Export diagnostic variables
   _DIAG_VAR_S_(data%id_diag_par)= par
   _DIAG_VAR_S_(data%id_gpp) = SUM(primprod)*secs_per_day
   IF( SUM(respiration(:)) > 1e-5 ) THEN
     _DIAG_VAR_S_(data%id_p2r)  = (SUM(primprod(:))/data%num_mphy) / (SUM(respiration(:))/data%num_mphy)
   ELSE
     _DIAG_VAR_S_(data%id_p2r)  = 9999.
   ENDIF


END SUBROUTINE aed2_calculate_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_macrophyte(data, column, layer_idx, extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to macrophyte biomass
!
!  WARNING - THIS IS ADDDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: dz,mphy,matz
   INTEGER  :: mphy_i
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i)) * (one_-mphy*data%mphydata(mphy_i)%f_bg) ! above ground density of macrophyte group i

      ! Self-shading depending on amount of carbon in water volume
      extinction = extinction + (data%mphydata(mphy_i)%KeMAC * (mphy/dz) )
   ENDDO
END SUBROUTINE aed2_light_extinction_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_bio_drag_macrophyte(data, column, layer_idx, drag)
!-------------------------------------------------------------------------------
! Get the effect of macrophyte biomass on benthic drag
!
!  WARNING - THIS IS ADDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
   CLASS (aed2_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   AED_REAL :: dz,mphy,matz
   INTEGER  :: mphy_i
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i)) * (one_-mphy*data%mphydata(mphy_i)%f_bg) ! above ground density of macrophyte group i

      ! additional drag due to biomass
      drag = drag + (data%mphydata(mphy_i)%K_CD * (mphy/dz) )
   ENDDO
END SUBROUTINE aed2_bio_drag_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_macrophyte
