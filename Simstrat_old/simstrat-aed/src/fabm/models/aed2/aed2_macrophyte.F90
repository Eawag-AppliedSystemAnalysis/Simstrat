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
   USE aed2_phyto_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_macrophyte_data_t

   TYPE,extends(aed2_model_data_t) :: aed2_macrophyte_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_mphy(:)
      INTEGER :: id_par, id_tem, id_sal, id_dz, id_extc, id_I_0
      INTEGER :: id_diag_par, id_gpp, id_p2r, id_sed_zone
      INTEGER :: n_zones
      AED_REAL,ALLOCATABLE :: active_zones(:)

      !# Model parameters
      INTEGER  :: num_mphy
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: mphydata  ! Using phyto data for now

     CONTAINS
         PROCEDURE :: define             => aed2_define_macrophyte
!        PROCEDURE :: calculate_riparian => aed2_calculate_riparian_macrophyte
         PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_macrophyte
!        PROCEDURE :: mobility           => aed2_mobility_macrophyte
!        PROCEDURE :: light_extinction   => aed2_light_extinction_macrophyte
!        PROCEDURE :: delete             => aed2_delete_macrophyte

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

   TYPE(phyto_nml_data) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /phyto_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD', iostat=status)
    IF (status /= 0) STOP 'Cannot open phyto_data namelist file for macrophytes'
    read(tfil,nml=phyto_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist phyto_data for macrophytes'

    data%num_mphy = count
    ALLOCATE(data%mphydata(count))

    ALLOCATE(data%id_mphy(count))

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%mphydata(i)%p_name       = pd(list(i))%p_name
       data%mphydata(i)%p0           = pd(list(i))%p0
     ! data%mphydata(i)%w_p          = pd(list(i))%w_p/secs_per_day
       data%mphydata(i)%Xcc          = pd(list(i))%Xcc
       data%mphydata(i)%R_growth     = pd(list(i))%R_growth/secs_per_day
       data%mphydata(i)%fT_Method    = pd(list(i))%fT_Method
       data%mphydata(i)%theta_growth = pd(list(i))%theta_growth
       data%mphydata(i)%T_std        = pd(list(i))%T_std
       data%mphydata(i)%T_opt        = pd(list(i))%T_opt
       data%mphydata(i)%T_max        = pd(list(i))%T_max
       data%mphydata(i)%lightModel   = pd(list(i))%lightModel
       data%mphydata(i)%I_K          = pd(list(i))%I_K
       data%mphydata(i)%I_S          = pd(list(i))%I_S
       data%mphydata(i)%KePHY        = pd(list(i))%KePHY
       data%mphydata(i)%f_pr         = pd(list(i))%f_pr
       data%mphydata(i)%R_resp       = pd(list(i))%R_resp/secs_per_day
       data%mphydata(i)%theta_resp   = pd(list(i))%theta_resp
     ! data%mphydata(i)%k_fres       = pd(list(i))%k_fres
     ! data%mphydata(i)%k_fdom       = pd(list(i))%k_fdom
       data%mphydata(i)%salTol       = pd(list(i))%salTol
       data%mphydata(i)%S_bep        = pd(list(i))%S_bep
       data%mphydata(i)%S_maxsp      = pd(list(i))%S_maxsp
       data%mphydata(i)%S_opt        = pd(list(i))%S_opt
     ! data%mphydata(i)%simDINUptake = pd(list(i))%simDINUptake
     ! data%mphydata(i)%simDONUptake = pd(list(i))%simDONUptake
     ! data%mphydata(i)%simNFixation = pd(list(i))%simNFixation
     ! data%mphydata(i)%simINDynamics= pd(list(i))%simINDynamics
     ! data%mphydata(i)%N_o          = pd(list(i))%N_o
     ! data%mphydata(i)%K_N          = pd(list(i))%K_N
     ! data%mphydata(i)%X_ncon       = pd(list(i))%X_ncon
     ! data%mphydata(i)%X_nmin       = pd(list(i))%X_nmin
     ! data%mphydata(i)%X_nmax       = pd(list(i))%X_nmax
     ! data%mphydata(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_per_day
     ! data%mphydata(i)%k_nfix       = pd(list(i))%k_nfix
     ! data%mphydata(i)%R_nfix       = pd(list(i))%R_nfix/secs_per_day
     ! data%mphydata(i)%simDIPUptake = pd(list(i))%simDIPUptake
     ! data%mphydata(i)%simIPDynamics= pd(list(i))%simIPDynamics
     ! data%mphydata(i)%P_0          = pd(list(i))%P_0
     ! data%mphydata(i)%K_P          = pd(list(i))%K_P
     ! data%mphydata(i)%X_pcon       = pd(list(i))%X_pcon
     ! data%mphydata(i)%X_pmin       = pd(list(i))%X_pmin
     ! data%mphydata(i)%X_pmax       = pd(list(i))%X_pmax
     ! data%mphydata(i)%R_puptake    = pd(list(i))%R_puptake/secs_per_day
     ! data%mphydata(i)%simSiUptake  = pd(list(i))%simSiUptake
     ! data%mphydata(i)%Si_0         = pd(list(i))%Si_0
     ! data%mphydata(i)%K_Si         = pd(list(i))%K_Si
     ! data%mphydata(i)%X_sicon      = pd(list(i))%X_sicon

       ! Register group as a state variable
       data%id_mphy(i) = aed2_define_sheet_variable(                   &
                              pd(list(i))%p_name,                      &
                              'mmolC/m**2', 'macrophyte',              &
                              pd(list(i))%p0,                          &
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
!  by the model are registered with FABM.
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


   NAMELIST /aed2_macrophyte/ num_mphy, the_mphy, dbase, n_zones, active_zones

!-----------------------------------------------------------------------
!BEGIN

   ! Read the namelist
   read(namlst,nml=aed2_macrophyte,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_macrophyte'

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

  ! CALL aed2_bio_temp_function(data%num_phytos,             &
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
   data%id_GPP = aed2_define_sheet_diag_variable('gpp','mmolC/m**3/day',  'benthic plant productivity')
   data%id_P2R = aed2_define_sheet_diag_variable('p_r','-',  'macrophyte P/R ratio')

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
   do i=1,l
      if ( active_zones(i) == matz ) then
         res = .TRUE.
         exit
      endif
   enddo

   in_zone_set = res
END FUNCTION in_zone_set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
   AED_REAL :: fT, fI, fSal
   AED_REAL :: extc, dz, par, Io, temp, salinity
   AED_REAL :: primprod(data%num_mphy)
   AED_REAL :: respiration(data%num_mphy)
   AED_REAL :: matz
!
!-------------------------------------------------------------------------------
!BEGIN
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)      ! local temperature
   salinity = _STATE_VAR_(data%id_sal)  ! local salinity
   par = _STATE_VAR_(data%id_par)       ! local photosynthetically active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values.
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i))! macrophyte group i

      ! LIGHT
      extc = _STATE_VAR_(data%id_extc)
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5
      fI   = phyto_light(data%mphydata, mphy_i, par, extc, Io, dz)

      primprod(mphy_i) = data%mphydata(mphy_i)%R_growth * fT * fI

      ! Respiration and general metabolic loss
      respiration(mphy_i) = phyto_respiration(data%mphydata,mphy_i,temp)

      ! Salinity stress effect on respiration
      fSal = 1.0 ! phyto_salinity(data%mphydata,mphy_i,salinity)
      respiration(mphy_i) = respiration(mphy_i) * fSal

      mphy_flux = (primprod(mphy_i) - respiration(mphy_i)) *  mphy

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     _FLUX_VAR_B_(data%id_mphy(mphy_i)) = _FLUX_VAR_B_(data%id_mphy(mphy_i)) + mphy_flux
   ENDDO

   ! Export diagnostic variables
   _DIAG_VAR_S_(data%id_gpp  ) = SUM(primprod)
   _DIAG_VAR_S_(data%id_diag_par )  = par
   _DIAG_VAR_S_(data%id_p2r )  = SUM(primprod(:)) / SUM(respiration(:))
END SUBROUTINE aed2_calculate_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_macrophyte
