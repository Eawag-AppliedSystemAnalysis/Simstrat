!###############################################################################
!#                                                                             #
!# aed2_cladophora.F90                                                         #
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


MODULE aed2_cladophora
!-------------------------------------------------------------------------------
!  aed2_cladophora --- multi-group cladophora (/macroalgae) model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util
   USE aed2_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_cladophora_data_t

   TYPE,extends(aed2_model_data_t) :: aed2_cladophora_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_malg(:),id_malg_in(:), id_malg_ip(:)
      INTEGER,ALLOCATABLE :: id_slough_trig(:)
      INTEGER :: id_slough, id_par, id_tem, id_sal, id_dz, id_extc, id_I_0
      INTEGER :: id_diag_par, id_gpp, id_p2r, id_mac, id_sed_zone, id_taub
      INTEGER :: n_zones
      AED_REAL,ALLOCATABLE :: active_zones(:)
      LOGICAL :: simSloughing

      !# Model parameters
      INTEGER  :: num_malg
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: malgdata  ! Using phyto data for now

     CONTAINS
         PROCEDURE :: define             => aed2_define_cladophora
!        PROCEDURE :: calculate_riparian => aed2_calculate_riparian_cladophora
         PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_cladophora
         PROCEDURE :: equilibrate        => aed2_slough_cladophora
!        PROCEDURE :: mobility           => aed2_mobility_cladophora
         PROCEDURE :: light_extinction   => aed2_light_extinction_cladophora
!        PROCEDURE :: delete             => aed2_delete_cladophora

   END TYPE

   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================




!###############################################################################
SUBROUTINE aed2_cladophora_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_cladophora_data_t),INTENT(inout) :: data
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
    IF (status /= 0) STOP 'Cannot open phyto_data namelist file for cladophora'
    read(tfil,nml=phyto_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist phyto_data for cladophora'

    data%num_malg = count
    ALLOCATE(data%malgdata(count))

    ALLOCATE(data%id_malg(count))
    ALLOCATE(data%id_slough_trig(count))

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%malgdata(i)%p_name       = pd(list(i))%p_name
       data%malgdata(i)%p0           = pd(list(i))%p0
     ! data%malgdata(i)%w_p          = pd(list(i))%w_p/secs_per_day
       data%malgdata(i)%Xcc          = pd(list(i))%Xcc
       data%malgdata(i)%R_growth     = pd(list(i))%R_growth/secs_per_day
       data%malgdata(i)%fT_Method    = pd(list(i))%fT_Method
       data%malgdata(i)%theta_growth = pd(list(i))%theta_growth
       data%malgdata(i)%T_std        = pd(list(i))%T_std
       data%malgdata(i)%T_opt        = pd(list(i))%T_opt
       data%malgdata(i)%T_max        = pd(list(i))%T_max
       data%malgdata(i)%lightModel   = pd(list(i))%lightModel
       data%malgdata(i)%I_K          = pd(list(i))%I_K
       data%malgdata(i)%I_S          = pd(list(i))%I_S
       data%malgdata(i)%KePHY        = pd(list(i))%KePHY
       data%malgdata(i)%f_pr         = pd(list(i))%f_pr
       data%malgdata(i)%R_resp       = pd(list(i))%R_resp/secs_per_day
       data%malgdata(i)%theta_resp   = pd(list(i))%theta_resp
     ! data%malgdata(i)%k_fres       = pd(list(i))%k_fres
     ! data%malgdata(i)%k_fdom       = pd(list(i))%k_fdom
       data%malgdata(i)%salTol       = pd(list(i))%salTol
       data%malgdata(i)%S_bep        = pd(list(i))%S_bep
       data%malgdata(i)%S_maxsp      = pd(list(i))%S_maxsp
       data%malgdata(i)%S_opt        = pd(list(i))%S_opt
     ! data%malgdata(i)%simDINUptake = pd(list(i))%simDINUptake
     ! data%malgdata(i)%simDONUptake = pd(list(i))%simDONUptake
     ! data%malgdata(i)%simNFixation = pd(list(i))%simNFixation
     ! data%malgdata(i)%simINDynamics= pd(list(i))%simINDynamics
     ! data%malgdata(i)%N_o          = pd(list(i))%N_o
     ! data%malgdata(i)%K_N          = pd(list(i))%K_N
       data%malgdata(i)%X_ncon       = pd(list(i))%X_ncon
       data%malgdata(i)%X_nmin       = pd(list(i))%X_nmin
       data%malgdata(i)%X_nmax       = pd(list(i))%X_nmax
     ! data%malgdata(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_per_day
     ! data%malgdata(i)%k_nfix       = pd(list(i))%k_nfix
     ! data%malgdata(i)%R_nfix       = pd(list(i))%R_nfix/secs_per_day
     ! data%malgdata(i)%simDIPUptake = pd(list(i))%simDIPUptake
     ! data%malgdata(i)%simIPDynamics= pd(list(i))%simIPDynamics
     ! data%malgdata(i)%P_0          = pd(list(i))%P_0
     ! data%malgdata(i)%K_P          = pd(list(i))%K_P
      data%malgdata(i)%X_pcon       = pd(list(i))%X_pcon
      data%malgdata(i)%X_pmin       = pd(list(i))%X_pmin
      data%malgdata(i)%X_pmax       = pd(list(i))%X_pmax
     ! data%malgdata(i)%R_puptake    = pd(list(i))%R_puptake/secs_per_day
     ! data%malgdata(i)%simSiUptake  = pd(list(i))%simSiUptake
     ! data%malgdata(i)%Si_0         = pd(list(i))%Si_0
     ! data%malgdata(i)%K_Si         = pd(list(i))%K_Si
     ! data%malgdata(i)%X_sicon      = pd(list(i))%X_sicon

       ! Register group as a state variable
       data%id_malg(i) = aed2_define_sheet_variable(                             &
                              data%malgdata(i)%p_name,                           &
                              'mmolC/m**2', 'cladophora',                        &
                              pd(list(i))%p0,                                    &
                              minimum=zero_)

       ! Register IN group as a state variable
       minNut = data%malgdata(i)%p0*data%malgdata(i)%X_nmin
       data%id_malg_in(i) = aed2_define_sheet_variable(                          &
                              TRIM(data%malgdata(i)%p_name)//'_IN',              &
                              'mmolN/m**3',                                      &
                              'cladophora '//TRIM(data%malgdata(i)%p_name)//'_IN', &
                              pd(list(i))%p_initial*data%malgdata(i)%X_ncon,     &
                              minimum=minNut)

       ! Register IP group as a state variable
       minNut = data%malgdata(i)%p0*data%malgdata(i)%X_pmin
       data%id_malg_in(i) = aed2_define_sheet_variable(                          &
                              TRIM(data%malgdata(i)%p_name)//'_IP',              &
                              'mmolP/m**3',                                      &
                              'cladophora '//TRIM(data%malgdata(i)%p_name)//'_IP', &
                              pd(list(i))%p_initial*data%malgdata(i)%X_pcon,     &
                              minimum=minNut)
       IF(data%simSloughing) THEN
         data%id_slough_trig(i) = aed2_define_sheet_variable(                          &
                              TRIM(data%malgdata(i)%p_name)//'_ST',              &
                              '-',                                      &
                              'cladophora '//TRIM(data%malgdata(i)%p_name)//' slough trigger', &
                              zero_,     &
                              minimum=-100.)
       ENDIF
    ENDDO
END SUBROUTINE aed2_cladophora_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_cladophora(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the cladophora/seagrass  model
!
!  Here, the aed2_ namelist is read and the variables exported
!  by the model are registered
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_cladophora_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_malg
   INTEGER  :: the_malg(MAX_PHYTO_TYPES)
   CHARACTER(len=128) :: dbase='aed2_cladophora_pars.nml'
   INTEGER  :: n_zones = 0, active_zones(MAX_ZONES), i
   LOGICAL  :: simSloughing = .FALSE.

   NAMELIST /aed2_cladophora/ num_malg, the_malg, dbase, n_zones, active_zones, simSloughing

!-----------------------------------------------------------------------
!BEGIN

   ! Read the namelist
   read(namlst,nml=aed2_cladophora,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_cladophora'

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   data%simSloughing = simSloughing

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   ! Main macroalgae state variables are allocated in here
   CALL aed2_cladophora_load_params(data, dbase, num_malg, the_malg)

   ! Add an extra 3D variable to store "slough" - broken off macroalgal material
   IF(data%simSloughing) THEN
     data%id_slough = aed2_define_sheet_variable(             &
                              'slough',                       &
                              'mmolC/m**3', 'slough',         &
                              zero_,                          &
                              minimum=zero_)
   END IF


   ! Register diagnostic variables
   data%id_diag_PAR = aed2_define_sheet_diag_variable('par','W/m**2','benthic light intensity')
   data%id_GPP = aed2_define_sheet_diag_variable('gpp','/d','benthic plant productivity')
   data%id_P2R = aed2_define_sheet_diag_variable('p_r','-','cladophora p/r ratio')
   data%id_MAC = aed2_define_sheet_diag_variable('mac','mmolC/m2','total cladophora biomass')
   data%id_tem_avg = aed2_define_sheet_diag_variable('avg_temp','degC','temperature moving average')
   data%id_par_avg = aed2_define_sheet_diag_variable('avg_light','W/m2','light moving average')
   data%id_tau_avg = aed2_define_sheet_diag_variable('avg_tau','N/m2','stress moving average')

   ! Register environmental dependencies
   data%id_dz = aed2_locate_global('layer_ht')
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_par = aed2_locate_global('par')
   data%id_I_0 = aed2_locate_global_sheet('par_sf')
   data%id_extc = aed2_locate_global('extc_coef')
   data%id_taub = aed2_locate_global_sheet('taub')
   data%id_sed_zone = aed2_locate_global_sheet('sed_zone')

END SUBROUTINE aed2_define_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!CAB This exists in a utils module
! !###############################################################################
! LOGICAL FUNCTION in_zone_set(matz, active_zones)
! !-------------------------------------------------------------------------------
! !ARGUMENTS
!    AED_REAL,INTENT(in) :: matz
!    AED_REAL,INTENT(in) :: active_zones(:)
! !
! !LOCALS
!    INTEGER :: i, l
!    LOGICAL :: res
! !BEGIN
! !-------------------------------------------------------------------------------
!    res = .FALSE.
!    l = size(active_zones)
!    do i=1,l
!       if ( active_zones(i) == matz ) then
!          res = .TRUE.
!          exit
!       endif
!    enddo
!
!    in_zone_set = res
! END FUNCTION in_zone_set
! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_cladophora(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of phytoplankton.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_cladophora_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL, PARAMETER :: DDT = (15./60.)/24.  ! 15 minutes

   AED_REAL, PARAMETER :: a01 = -0.7171732e-01 !
   AED_REAL, PARAMETER :: a02 =  0.6610416     !
   AED_REAL, PARAMETER :: a03 =  0.4797417     !
   AED_REAL, PARAMETER :: a04 = -0.2358428e+01 !
   AED_REAL, PARAMETER :: a05 =  0.4662903e+01 !
   AED_REAL, PARAMETER :: a06 = -0.2759110e+01 !
   AED_REAL, PARAMETER :: a07 =  0.3343594e+01 !
   AED_REAL, PARAMETER :: a08 = -0.5609280e+01 !
   AED_REAL, PARAMETER :: a09 = -0.3914500e+01 !
   AED_REAL, PARAMETER :: a10 =  0.3372881e+01 !
   AED_REAL, PARAMETER :: a11 = -0.1742290e+01 !
   AED_REAL, PARAMETER :: a12 =  0.1304549e+01 !
   AED_REAL, PARAMETER :: a13 =  0.2618247e+01 !
   AED_REAL, PARAMETER :: a14 =  0.8705460     !
   AED_REAL, PARAMETER :: a15 = -0.1207573e+01 !
   AED_REAL, PARAMETER :: b01 =  0.8204769E-02 !
   AED_REAL, PARAMETER :: b02 =  0.1622813E+00 !
   AED_REAL, PARAMETER :: b03 =  0.5599809E+00 !
   AED_REAL, PARAMETER :: b04 = -0.9230173E-01 !
   AED_REAL, PARAMETER :: b05 = -0.3701472E+00 !
   AED_REAL, PARAMETER :: b06 = -0.1466722E+01 !
   AED_REAL, PARAMETER :: b07 = -0.1021670E+00 !
   AED_REAL, PARAMETER :: b08 =  0.6366414E-01 !
   AED_REAL, PARAMETER :: b09 =  0.9596720E+00 !
   AED_REAL, PARAMETER :: b10 =  0.1632652E+01 !
   AED_REAL, PARAMETER :: b11 =  0.1431378E+00 !
   AED_REAL, PARAMETER :: b12 =  0.6283502E+00 !
   AED_REAL, PARAMETER :: b13 = -0.1520305E+01 !
   AED_REAL, PARAMETER :: b14 =  0.3120532E+00 !
   AED_REAL, PARAMETER :: b15 = -0.7267257E+00 !

   AED_REAL :: malg
   INTEGER  :: malg_i
   AED_REAL :: extc, dz, par, Io, temp, salinity
   AED_REAL :: matz

   AED_REAL :: X_maxp
   AED_REAL :: lght, pplt, prlt, pf_MB
   AED_REAL :: Kq

!
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for cladophoras
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Retrieve current environmental conditions
   salinity = _STATE_VAR_(data%id_sal)  ! local salinity
   temp = _STATE_VAR_(data%id_tem)      ! local temperature
   extc = _STATE_VAR_(data%id_extc)     ! extinction coefficent
   par = _STATE_VAR_(data%id_par)       ! local photosynthetically active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation
   dz = _STATE_VAR_(data%id_dz)         ! dz = 0.5

   ! Initialise cumulative biomass diagnostic
   _DIAG_VAR_S_(data%id_mac) = 0

   DO malg_i=1,data%num_malg
      ! Retrieve current (local) state variable values
      malg = _STATE_VAR_S_(data%id_malg(malg_i))! cladophora group i

      ! BED DEPTH AND EXTINCTION
      ! Emprical forumaltion from Higgins CGM model
      macroHgt = (1/100.0) * 1.467 * ( (1./0.25)*malg )**0.425
      macroExt = macroExt +  7.830 * ( (1./0.25)*malg )**0.250
   ENDDO

   ! LIGHT
   macroPAR_Top = par * exp(-extc*(dz-macroHgt))
   macroPAR_Bot = par * exp(-(extc+macroExt)*(dz-macroHgt))


    onm = one_; opm = zero_

    !--------------------------------------------------------------------------!
    !-- Update moving average for daily temp
    !   (averaged over the past 1 day)

    !_DIAG_VAR_S_(data%id_tem_avg) =  _DIAG_VAR_S_(data%id_tem_avg)*(1-(DDT/TempAvgTime)) + temp*(DDT/TempAvgTime)
    AvgTemp = temp ! _DIAG_VAR_S_(data%id_tem_avg)

    TempAvgTime = MIN(TempAvgTime+DDT,1.0)

    !--------------------------------------------------------------------------!
    !-- Now check light, and perform specific day or night activities
    IF(Io > 10.0) THEN

      !------------------------------------------------------------------------!
      !-- Calculate the nutrient limitation (phosporus and nitrogen) and
      !   find the most limiting
      sf = zero_

      IF (malg > zero_) THEN
        sf = malg_ip/malg
      ENDIF
      IF (sf > zero_) THEN
        sf = 1.0 - (data%IPmin/sf)
      ENDIF
      print *,'PHOS',sf


      !------------------------------------------------------------------------!
      !-- Update moving average for daily light
      !   (averaged over the photo-period, PP)

      !_DIAG_VAR_S_(data%id_par_avg) = _DIAG_VAR_S_(data%id_par_avg) * (1-(DDT/LgtAvgTime)) + par*(DDT/LgtAvgTime)
      AvgLight = par     ! _DIAG_VAR_S_(data%id_par_avg)

      LgtAvgTime = MIN(LgtAvgTime+DDT,0.5)

      !------------------------------------------------------------------------!
      !-- Self-shading

      ! empirical thing with *0.25 to get from g DM to g C
      X_maxp = ( 1.18 * AvgLight - 58.7 ) * 0.25

      lf = 1.0 - malg/X_maxp

      IF(lf < zero_) lf = zero_
      IF(lf > one_) lf = one_

      print *,'SHAD',lf
      print *,'LGHT',AvgLight
      print *,'TEMP',AvgTemp

      !--------------------------------------------------------------------------!
      !-- Now light and temperature function for photosynthesis
      temp = AvgTemp/35.0

      lght = MIN(600.0,AvgLight)/1235.0

      pplt = a01                             &
           + a02 * temp                      &
           + a03 * lght                      &
           + a04 * temp * temp               &
           + a05 * temp * lght               &
           + a06 * lght * lght               &
           + a07 * temp * temp * temp        &
           + a08 * temp * temp * lght        &
           + a09 * temp * lght * lght        &
           + a10 * lght * lght * lght        &
           + a11 * temp * temp * temp * temp &
           + a12 * temp * temp * temp * lght &
           + a13 * temp * temp * lght * lght &
           + a14 * temp * lght * lght * lght &
           + a15 * lght * lght * lght * lght

      pf = (data%Vmax(malg_i) * pplt * 5.43) * sf * lf     ! + (0.44 * prlt * 4.52)

      print *,'PROD',pf

      !--------------------------------------------------------------------------!
      !-- Now light and temperature function for daytime respiration

      prlt = b01                             &
           + b02 * temp                      &
           + b03 * lght                      &
           + b04 * temp * temp               &
           + b05 * temp * lght               &
           + b06 * lght * lght               &
           + b07 * temp * temp * temp        &
           + b08 * temp * temp * lght        &
           + b09 * temp * lght * lght        &
           + b10 * lght * lght * lght        &
           + b11 * temp * temp * temp * temp &
           + b12 * temp * temp * temp * lght &
           + b13 * temp * temp * lght * lght &
           + b14 * temp * lght * lght * lght &
           + b15 * lght * lght * lght * lght


      !rf = (0.44 * pplt * 4.52) * sf * lf

      rf= zero_

      !rf = 0.151 * (0.0025 * AvgTemp + 0.1)

      print *,'RDAY',rf

    ELSE
      !------------------------------------------------------------------------!
      !-- Now calculate basal respiration....
      rf = 0.151 * (0.0025 * AvgTemp + 0.1)

      pf = zero_

      print *,'PROD',pf
      print *,'RBAS',rf

    END IF


    !--------------------------------------------------------------------------!
    !-- Get the INTERNAL PHOSPHORUS stores for the macroalgae groups.          !
    !   Recall that internal nutrient is in mol and must be converted to mol P/mol C !
    !   by division by macroalgae biomass for the nitrogen limitation          !

    !-- Get the internal phosphorus ratio
    IF(malg>zero_) THEN
      pu = malg_ip/malg
    ELSEWHERE
      pu = zero_
    END WHERE

    print *,'IP/M',pu

    temp = AvgTemp
    WHERE(temp < 18.0)
      tf = exp((temp-18.0)/39.00)
    ELSEWHERE
      tf = exp((18.0-temp)/18.75)
    END WHERE


    Kq = 0.0028 ! 0.07% = 0.0028gP/gC     !!! FACTOR IN STOICHIOMETRY

    !-- IPmax(bb) = Kq; IPmin(bb) = Qo; KP(bb) = Km; UPmax(bb) = pmax; tau = tf
    pu = UPmax(malg_i) * tf * (frp/(frp + KP(malg_i))) * (Kq /(Kq + pu - IPmin(malg_i)))

    print *,'P_UP',pu

    !-- Update the internal phosphorous store
    _FLUX_VAR_B_(data%id_malg_ip(malg_i)) =  _FLUX_VAR_B_(data%id_malg_ip(malg_i)) &
                                          +  pu * malg

    ! CALL EnforceIntNutLimits(IPmin(bb),IPmax(bb),malg_ip,malg,"P")


    !--------------------------------------------------------------------------!
    !-- Get the INTERNAL NITROGEN stores for the macroalgae groups.            !
    !   Recall that internal nutrient is in g and must be converted to g N/g C !
    !   by division by macroalgae biomass for the nitrogen limitation          !

    !-- Get the internal nitrogen ratio.
    IF(malg>zero_) THEN
      nu = malg_in/malg
    ELSEWHERE
      nu = zero_
    END WHERE

    print *,'IN/M',nu

    !-- Macroalgae nitrogen uptake
    nu = UNmax(bb) * WQ%sed(SedCells,SED(M(macIndex))) * tf * (INmax(bb) - nu) &
       / (INmax(bb)-INmin(bb)) * (WQ%a3d(BotCells,NO3) + WQ%a3d(BotCells,NH4)) &
       / (WQ%a3d(BotCells,NO3) + WQ%a3d(BotCells,NH4) + KN(bb))

    !-- Time series quantities
    IF(macGroup == ms) THEN
      map = DDT * nu(vs2)
      mar = DDT * rf(vs2) * sf(vs2) * WQ%sed(SedCells(vs2),SED(INm(macIndex)))
      mag = DDT * mf(vs2) * WQ%sed(SedCells(vs2),SED(INm(macIndex)))
    END IF

    !-- Update the internal nitrogen store
    !_FLUX_VAR_B_(data%id_malg_in(malg_i)) =  _FLUX_VAR_B_(data%id_malg_in(malg_i)) &
    !                                      +  nu * malg
    _STATE_VAR_S_(data%id_malg_in(malg_i)) =  0.9 * INmax(bb) * malg

    !-- Limit to minimum and maximum ratios.
    ! CALL EnforceIntNutLimits(INmin(bb),INmax(bb),malg_in,malg,"N")


    !--------------------------------------------------------------------------!
    !-- If base is dead, then slough live cladophora into the water column     !
    pf_MB = zero_

    !-- first need to examine light & temp at the bottom of the bed
    temp = AvgTemp/35.0

   !lght = MIN(600.0,macroPAR_Bot)/1235.0
    lght = MIN(600.0, AvgLight * EXP(-macroExt * macroHgt) )/1235.0

    pplt = a01                             &
         + a02 * temp                      &
         + a03 * lght                      &
         + a04 * temp * temp               &
         + a05 * temp * lght               &
         + a06 * lght * lght               &
         + a07 * temp * temp * temp        &
         + a08 * temp * temp * lght        &
         + a09 * temp * lght * lght        &
         + a10 * lght * lght * lght        &
         + a11 * temp * temp * temp * temp &
         + a12 * temp * temp * temp * lght &
         + a13 * temp * temp * lght * lght &
         + a14 * temp * lght * lght * lght &
         + a15 * lght * lght * lght * lght

    pf_MB = (data%Vmax(malg_i) * pplt * 5.43) * sf  -  rf

    print *,'pfMB',pf_MB


    IF(malg > 0.1) THEN
      !SloughTrigger(malg_i) = SloughTrigger(malg_i) + pf_MB * DDT
      _FLUX_VAR_B_(data%id_slough_trig) = _FLUX_VAR_B_(data%id_slough_trig) + pf_MB
    ENDIF
    print *,'SLTG',SloughTrigger(malg_i)

    sf = one_


! Uptake / Excretion feedbacks  C/N/P/DO

! Slough variable update feedback


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   DO malg_i=1,data%num_malg

      fI   = photosynthesis_irradiance(data%malgdata(malg_i)%lightModel, &
                    data%malgdata(malg_i)%I_K, data%malgdata(malg_i)%I_S, par, extc, Io, dz)
      fT   = 1.

      primprod(malg_i) = data%malgdata(malg_i)%R_growth * fI * fT

      ! Respiration and general metabolic loss
      respiration(malg_i) = bio_respiration(data%malgdata(malg_i)%R_resp, data%malgdata(malg_i)%theta_resp,temp)

      ! Salinity stress effect on respiration
      fSal = 1.0 ! phyto_salinity(data%malgdata,malg_i,salinity)
      respiration(malg_i) = respiration(malg_i) * fSal

      malg_flux = (primprod(malg_i) - respiration(malg_i)) *  malg

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     _FLUX_VAR_B_(data%id_malg(malg_i)) = _FLUX_VAR_B_(data%id_malg(malg_i)) + malg_flux

     _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + malg
   ENDDO

   ! Export diagnostic variables
   _DIAG_VAR_S_(data%id_gpp  ) = SUM(primprod)*secs_per_day
   _DIAG_VAR_S_(data%id_diag_par )  = par
   _DIAG_VAR_S_(data%id_p2r )  = SUM(primprod(:)) / SUM(respiration(:))

!   ! Update critical shear stress in this location (if feedback is on)
!   if (simFeedbacks .AND. data%id_tau >0) then
!     _DIAG_VAR_S_(data%id_tau)  = (data%taumax(matz))*(1. -  MIN(malg/max_malg,1.))
!   endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END SUBROUTINE aed2_calculate_benthic_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_slough_cladophora(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update cladophora biomass due to sloughing if stress is enough
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_cladophora_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: tau
   ! State
   AED_REAL :: malg
!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT. data%simSloughing) RETURN


   ! Retrieve current environmental conditions for the cell.
   bottom_stress = _STATE_VAR_S_(data%id_taub)   ! local bottom shear stress
   bottom_stress = MIN(bottom_stress, 100.)


    !--------------------------------------------------------------------------!
    !-- Update moving average for stress
    !   (averaged over the past 2 hrs)

    !_DIAG_VAR_S_(data%id_tau_avg) = _DIAG_VAR_S_(data%id_tau_avg) * (1-(DDT/StrAvgTime)) + bottom_stress *(DDT/StrAvgTime)
    AvgStress = bottom_stress ! _DIAG_VAR_S_(data%id_tau_avg)

    StrAvgTime = MIN(StrAvgTime+DDT,1.0/12.0)

   ! Retrieve current (local) state variable values.
   DO malg_i=1,data%num_malg
     malg = _STATE_VAR_(data%id_malg(malg_i))
     slough_trigger = _STATE_VAR_S_(data%id_slough_trig(malg_i))
     ! Check if growth phase; if so clip to 0
     IF(slough_trigger > zero_) _STATE_VAR_S_(data%id_slough_trig(malg_i)) = zero_
     ! Slough off weakened filaments (those with cumulative respiration excess)
     IF(slough_trigger < -0.4) THEN
       hf = one_
       _STATE_VAR_S_(data%id_slough_trig(malg_i)) = zero_
     ELSE
       hf = zero_
     ENDIF

     ! Sloughing of even healthy filaments due to high shear stresses
     IF(malg>0.1 .AND. slough_trigger<-0.001 .AND. AvgStress>data%PCm(malg_i))THEN

      ! Lss = Lmax (T/Tcrit) * (X/Xmax)

      !hf = 0.242 * exp(-0.3187 * (DWQ%LZ(BotCells)+DWQ%DZ(BotCells)))          &
      !   * AvgStress / PCm(macGroup)                                 &
      !   * WQ%sed(SedCells,SED(M(macIndex)))/ X_maxp

      hf = 0.0010 * (AvgStress - data%PCm(malg_i)) / one_  * malg/ X_maxp
      hf = hf * DDT    !??

     ENDIF

    print *,'hf :',hf,AvgStress

    IF (hf>0.99) hf = 0.99

    ! Update the malg and slough variables with following the slough event
    _STATE_VAR_S_(data%id_slough_trig(malg_i))    = slough_trigger
  ENDDO
END SUBROUTINE aed2_slough_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_cladophora(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to cladophora biomass
!
!  WARNING - THIS IS ADDDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_cladophora_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: dz,malg
!
!-------------------------------------------------------------------------------
!BEGIN

   DO malg_i=1,data%num_malg
      ! Retrieve current (local) state variable values
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5
      malg = _STATE_VAR_S_(data%id_malg(malg_i))! cladophora group i

      ! Self-shading depending on amount of carbon in water volume
      extinction = extinction + (data%malgdata(malg_i)%KePHY * (malg/dz) )
   ENDDO

END SUBROUTINE aed2_light_extinction_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_bio_drag_cladophora(data,column,layer_idx,drag)
!-------------------------------------------------------------------------------
! Get the effect of cladophora biomass on benthic drag
!
!  WARNING - THIS IS ADDDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_cladophora_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   AED_REAL :: dz,malg
!
!-------------------------------------------------------------------------------
!BEGIN

  DO malg_i=1,data%num_malg
      ! Retrieve current (local) state variable values
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5
      malg = _STATE_VAR_S_(data%id_malg(malg_i))! cladophora group i

      ! Self-shading depending on amount of carbon in water volume
      drag = drag + (data%malgdata(malg_i)%KCD * (malg/dz) )
   ENDDO
END SUBROUTINE aed2_bio_drag_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed2_cladophora
