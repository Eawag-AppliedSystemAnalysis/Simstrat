!###############################################################################
!#                                                                             #
!# aed2_phytoplankton.F90                                                      #
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
!# Created August 2011                                                         #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#define _MOB_CONST_ 0
#define _MOB_TEMP_  1
#define _MOB_CALC_  2


MODULE aed2_phytoplankton
!-------------------------------------------------------------------------------
!  aed2_phytoplankton --- phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,ONLY : find_free_lun, &
                       exp_integral,  &
                       aed2_bio_temp_function, &
                       fTemp_function
   USE aed2_phyto_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_phytoplankton_data_t
!

   TYPE,extends(aed2_model_data_t) :: aed2_phytoplankton_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_p(:)
      INTEGER,ALLOCATABLE :: id_in(:)
      INTEGER,ALLOCATABLE :: id_ip(:)
      INTEGER,ALLOCATABLE :: id_rho(:)
      INTEGER :: id_Pexctarget,id_Pmorttarget,id_Pupttarget(1:2)
      INTEGER :: id_Nexctarget,id_Nmorttarget,id_Nupttarget(1:4)
      INTEGER :: id_Cexctarget,id_Cmorttarget,id_Cupttarget
      INTEGER :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      INTEGER :: id_DOupttarget
      INTEGER :: id_par, id_tem, id_sal, id_dz, id_extc
      INTEGER :: id_I_0
      INTEGER :: id_GPP, id_NCP, id_PPR, id_NPR, id_dPAR
      INTEGER :: id_TPHY, id_TCHLA, id_TIN, id_TIP
      INTEGER :: id_NUP, id_PUP, id_CUP
      INTEGER,ALLOCATABLE :: id_NtoP(:)
      INTEGER,ALLOCATABLE :: id_fT(:), id_fI(:), id_fNit(:), id_fPho(:), id_fSil(:), id_fSal(:)

      !# Model parameters
      INTEGER  :: num_phytos
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: phytos
      ! LOGICAL  :: do_exc,do_mort,do_upt, do_N2uptake
      LOGICAL  :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL  :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL  :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL  :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER  :: nnup, npup
      AED_REAL :: dic_per_n

     CONTAINS
         PROCEDURE :: define            => aed2_define_phytoplankton
         PROCEDURE :: calculate         => aed2_calculate_phytoplankton
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_phytoplankton
         PROCEDURE :: mobility          => aed2_mobility_phytoplankton
         PROCEDURE :: light_extinction  => aed2_light_extinction_phytoplankton
!        PROCEDURE :: delete            => aed2_delete_phytoplankton

   END TYPE

   AED_REAL :: dtlim = 0.9 * 3600
   LOGICAL  :: extra_debug = .false.

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_phytoplankton_load_params(data, dbase, count, list, w_model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_phytoplankton_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: w_model(*)
!
!LOCALS
   INTEGER  :: status
!  AED_REAL :: p_initial=0.
!  AED_REAL :: p0=0.0225
!  AED_REAL :: w_p=-1.157407e-05
!  AED_REAL :: i_min=25.
!  AED_REAL :: rmax=1.157407e-05
!  AED_REAL :: alpha=0.3
!  AED_REAL :: rpn=1.157407e-07
!  AED_REAL :: rpdu=2.314814e-07
!  AED_REAL :: rpdl=1.157407e-06

   INTEGER  :: i,tfil
   AED_REAL :: minNut

   AED_REAL, parameter :: secs_pr_day = 86400.
   TYPE(phyto_nml_data) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /phyto_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD', iostat=status)
    IF (status /= 0) STOP 'Cannot open phyto_data namelist file'
    read(tfil,nml=phyto_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist phyto_data'

    data%num_phytos = count
    ALLOCATE(data%phytos(count))
    ALLOCATE(data%id_p(count))
    ALLOCATE(data%id_in(count))
    ALLOCATE(data%id_ip(count))
    ALLOCATE(data%id_rho(count))
    ALLOCATE(data%id_NtoP(count))
    IF (extra_debug) THEN
       ALLOCATE(data%id_fT(count))
       ALLOCATE(data%id_fI(count))
       ALLOCATE(data%id_fNit(count))
       ALLOCATE(data%id_fPho(count))
       ALLOCATE(data%id_fSil(count))
       ALLOCATE(data%id_fSal(count))
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%phytos(i)%p_name       = pd(list(i))%p_name
       data%phytos(i)%p0           = pd(list(i))%p0
       data%phytos(i)%w_p          = pd(list(i))%w_p/secs_pr_day
       data%phytos(i)%w_model      = w_model(i)
       data%phytos(i)%Xcc          = pd(list(i))%Xcc
       data%phytos(i)%R_growth     = pd(list(i))%R_growth/secs_pr_day
       data%phytos(i)%fT_Method    = pd(list(i))%fT_Method
       data%phytos(i)%theta_growth = pd(list(i))%theta_growth
       data%phytos(i)%T_std        = pd(list(i))%T_std
       data%phytos(i)%T_opt        = pd(list(i))%T_opt
       data%phytos(i)%T_max        = pd(list(i))%T_max
       data%phytos(i)%lightModel   = pd(list(i))%lightModel
       data%phytos(i)%I_K          = pd(list(i))%I_K
       data%phytos(i)%I_S          = pd(list(i))%I_S
       data%phytos(i)%KePHY        = pd(list(i))%KePHY
       data%phytos(i)%f_pr         = pd(list(i))%f_pr
       data%phytos(i)%R_resp       = pd(list(i))%R_resp/secs_pr_day
       data%phytos(i)%theta_resp   = pd(list(i))%theta_resp
       data%phytos(i)%k_fres       = pd(list(i))%k_fres
       data%phytos(i)%k_fdom       = pd(list(i))%k_fdom
       data%phytos(i)%salTol       = pd(list(i))%salTol
       data%phytos(i)%S_bep        = pd(list(i))%S_bep
       data%phytos(i)%S_maxsp      = pd(list(i))%S_maxsp
       data%phytos(i)%S_opt        = pd(list(i))%S_opt
       data%phytos(i)%simDINUptake = pd(list(i))%simDINUptake
       data%phytos(i)%simDONUptake = pd(list(i))%simDONUptake
       data%phytos(i)%simNFixation = pd(list(i))%simNFixation
       data%phytos(i)%simINDynamics= pd(list(i))%simINDynamics
       data%phytos(i)%N_o          = pd(list(i))%N_o
       data%phytos(i)%K_N          = pd(list(i))%K_N
       data%phytos(i)%X_ncon       = pd(list(i))%X_ncon
       data%phytos(i)%X_nmin       = pd(list(i))%X_nmin
       data%phytos(i)%X_nmax       = pd(list(i))%X_nmax
       data%phytos(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_pr_day
       data%phytos(i)%k_nfix       = pd(list(i))%k_nfix
       data%phytos(i)%R_nfix       = pd(list(i))%R_nfix/secs_pr_day
       data%phytos(i)%simDIPUptake = pd(list(i))%simDIPUptake
       data%phytos(i)%simIPDynamics= pd(list(i))%simIPDynamics
       data%phytos(i)%P_0          = pd(list(i))%P_0
       data%phytos(i)%K_P          = pd(list(i))%K_P
       data%phytos(i)%X_pcon       = pd(list(i))%X_pcon
       data%phytos(i)%X_pmin       = pd(list(i))%X_pmin
       data%phytos(i)%X_pmax       = pd(list(i))%X_pmax
       data%phytos(i)%R_puptake    = pd(list(i))%R_puptake/secs_pr_day
       data%phytos(i)%simSiUptake  = pd(list(i))%simSiUptake
       data%phytos(i)%Si_0         = pd(list(i))%Si_0
       data%phytos(i)%K_Si         = pd(list(i))%K_Si
       data%phytos(i)%X_sicon      = pd(list(i))%X_sicon

       ! Register group as a state variable
       data%id_p(i) = aed2_define_variable(                         &
                              TRIM(data%phytos(i)%p_name),                     &
                              'mmol/m**3',                                     &
                              'phytoplankton '//TRIM(data%phytos(i)%p_name),   &
                              pd(list(i))%p_initial,                           &
                              minimum=pd(list(i))%p0,                          &
                              mobility = data%phytos(i)%w_p)

       IF (data%phytos(i)%w_model == _MOB_CALC_) THEN
          ! Register rho group as a state variable
          data%id_rho(i) = aed2_define_variable(                    &
                              TRIM(data%phytos(i)%p_name)//'_rho',             &
                              'mmol/m**3',                                     &
                        'phytoplankton '//TRIM(data%phytos(i)%p_name)//'_rho', &
                              pd(list(i))%w_p,                                 &
                              minimum=zero_,                                   &
                              mobility = data%phytos(i)%w_p)

       ENDIF
       IF (data%phytos(i)%simINDynamics /= 0) THEN
          IF(data%phytos(i)%simINDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_ncon
          ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_nmin
          ENDIF
          ! Register IN group as a state variable
          data%id_in(i) = aed2_define_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IN',              &
                              'mmol/m**3',                                     &
                         'phytoplankton '//TRIM(data%phytos(i)%p_name)//'_IN', &
                              pd(list(i))%p_initial*data%phytos(i)%X_ncon,     &
                              minimum=minNut,                                  &
                              mobility = data%phytos(i)%w_p)

       ENDIF
       IF (data%phytos(i)%simIPDynamics /= 0) THEN
          IF(data%phytos(i)%simIPDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_pcon
          ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_pmin
          ENDIF
          ! Register IP group as a state variable
          data%id_ip(i) = aed2_define_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IP',              &
                              'mmol/m**3',                                     &
                         'phytoplankton '//TRIM(data%phytos(i)%p_name)//'_IP', &
                              pd(list(i))%p_initial*data%phytos(i)%X_pcon,     &
                              minimum=minNut,                                  &
                              mobility = data%phytos(i)%w_p)

       ENDIF

       data%id_NtoP(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_NtoP','mmol/m**3', 'INi/IPi')

       IF (extra_debug) THEN
          data%id_fT(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fT', 'mmol/m**3', 'fT')
          data%id_fI(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fI', 'mmol/m**3', 'fI')
          data%id_fNit(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fNit', 'mmol/m**3', 'fNit')
          data%id_fPho(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fPho', 'mmol/m**3', 'fPho')
          data%id_fSil(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSil', 'mmol/m**3', 'fSil')
          data%id_fSal(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSal', 'mmol/m**3', 'fSal')
       ENDIF
    ENDDO
END SUBROUTINE aed2_phytoplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_phytoplankton(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the phytoplankton biogeochemical model
!
!  Here, the aed2_p_m namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_phytoplankton_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_phytos
   INTEGER  :: the_phytos(MAX_PHYTO_TYPES)
   INTEGER  :: w_model(MAX_PHYTO_TYPES)
   CHARACTER(len=64)  :: p_excretion_target_variable=''
   CHARACTER(len=64)  :: p_mortality_target_variable=''
   CHARACTER(len=64)  :: p1_uptake_target_variable=''
   CHARACTER(len=64)  :: p2_uptake_target_variable=''
   CHARACTER(len=64)  :: n_excretion_target_variable=''
   CHARACTER(len=64)  :: n_mortality_target_variable=''
   CHARACTER(len=64)  :: n1_uptake_target_variable=''
   CHARACTER(len=64)  :: n2_uptake_target_variable=''
   CHARACTER(len=64)  :: n3_uptake_target_variable=''
   CHARACTER(len=64)  :: n4_uptake_target_variable=''
   CHARACTER(len=64)  :: c_excretion_target_variable=''
   CHARACTER(len=64)  :: c_mortality_target_variable=''
   CHARACTER(len=64)  :: c_uptake_target_variable=''
   CHARACTER(len=64)  :: do_uptake_target_variable=''
   CHARACTER(len=64)  :: si_excretion_target_variable=''
   CHARACTER(len=64)  :: si_mortality_target_variable=''
   CHARACTER(len=64)  :: si_uptake_target_variable=''
   CHARACTER(len=128) :: dbase='aed2_phyto_pars.nml'


   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   AED_REAL           :: zerolimitfudgefactor = 0.9 * 3600
   NAMELIST /aed2_phytoplankton/ num_phytos, the_phytos, w_model,              &
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    si_excretion_target_variable,si_mortality_target_variable, &
                      si_uptake_target_variable,                               &
                    dbase, zerolimitfudgefactor, extra_debug
!-----------------------------------------------------------------------
!BEGIN
   w_model = _MOB_CONST_
   ! Read the namelist
   read(namlst,nml=aed2_phytoplankton,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_phytoplankton'
   dtlim = zerolimitfudgefactor

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   CALL aed2_phytoplankton_load_params(data, dbase, num_phytos, the_phytos, w_model)

   CALL aed2_bio_temp_function(data%num_phytos,             &
                              data%phytos%theta_growth,     &
                              data%phytos%T_std,            &
                              data%phytos%T_opt,            &
                              data%phytos%T_max,            &
                              data%phytos%aTn,              &
                              data%phytos%bTn,              &
                              data%phytos%kTn,              &
                              data%phytos%p_name)


   ! Register link to nutrient pools, if variable names are provided in namelist.
   data%do_Pexc = p_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Pexctarget = aed2_locate_variable( p_excretion_target_variable)
   ENDIF
   data%do_Nexc = n_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Nexctarget = aed2_locate_variable( n_excretion_target_variable)
   ENDIF
   data%do_Cexc = c_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Cexctarget = aed2_locate_variable( c_excretion_target_variable)
   ENDIF
   data%do_Siexc = si_excretion_target_variable .NE. ''
   IF (data%do_Siexc) THEN
     data%id_Siexctarget = aed2_locate_variable( si_excretion_target_variable)
   ENDIF

   data%do_Pmort = p_mortality_target_variable .NE. ''
   IF (data%do_Pmort) THEN
     data%id_Pmorttarget = aed2_locate_variable( p_mortality_target_variable)
   ENDIF
   data%do_Nmort = n_mortality_target_variable .NE. ''
   IF (data%do_Nmort) THEN
     data%id_Nmorttarget = aed2_locate_variable( n_mortality_target_variable)
   ENDIF
   data%do_Cmort = c_mortality_target_variable .NE. ''
   IF (data%do_Cmort) THEN
     data%id_Cmorttarget = aed2_locate_variable( c_mortality_target_variable)
   ENDIF
   data%do_Simort = si_mortality_target_variable .NE. ''
   IF (data%do_Simort) THEN
     data%id_Simorttarget = aed2_locate_variable( si_mortality_target_variable)
   ENDIF

   data%npup = 0
   IF (p1_uptake_target_variable .NE. '') data%npup = 1
   IF (p2_uptake_target_variable .NE. '') data%npup = 2
   data%do_Puptake = .FALSE.
   IF (data%npup>0) data%do_Puptake=.TRUE.
   IF (data%do_Puptake) THEN
     IF (data%npup>0) data%id_Pupttarget(1) = aed2_locate_variable( p1_uptake_target_variable); ifrp=1
     IF (data%npup>1) data%id_Pupttarget(2) = aed2_locate_variable( p2_uptake_target_variable); idop=2
   ENDIF
   data%nnup = 0
   IF (n1_uptake_target_variable .NE. '') data%nnup = 1
   IF (n2_uptake_target_variable .NE. '') data%nnup = 2
   IF (n3_uptake_target_variable .NE. '') data%nnup = 3
   IF (n4_uptake_target_variable .NE. '') data%nnup = 4
   data%do_Nuptake = .false.
   IF (data%nnup>0) data%do_Nuptake=.true.
   IF (data%do_Nuptake) THEN
     IF (data%nnup>0) data%id_Nupttarget(1) = aed2_locate_variable( n1_uptake_target_variable); ino3=1
     IF (data%nnup>1) data%id_Nupttarget(2) = aed2_locate_variable( n2_uptake_target_variable); inh4=2
     IF (data%nnup>2) data%id_Nupttarget(3) = aed2_locate_variable( n3_uptake_target_variable); idon=3
     IF (data%nnup>3) data%id_Nupttarget(4) = aed2_locate_variable( n4_uptake_target_variable); in2 =4
   ENDIF
   data%do_Cuptake = c_uptake_target_variable .NE. ''
   IF (data%do_Cuptake) THEN
     data%id_Cupttarget = aed2_locate_variable( c_uptake_target_variable)
   ENDIF
   data%do_DOuptake = do_uptake_target_variable .NE. ''
   IF (data%do_DOuptake) THEN
     data%id_DOupttarget = aed2_locate_variable( do_uptake_target_variable)
   ENDIF
   data%do_Siuptake = si_uptake_target_variable .NE. ''
   IF (data%do_Siuptake) THEN
     data%id_Siupttarget = aed2_locate_variable( si_uptake_target_variable)
   ENDIF

   ! Register diagnostic variables
   data%id_GPP = aed2_define_diag_variable('GPP','mmol/m**3',  'gross primary production')
   data%id_NCP = aed2_define_diag_variable('NCP','mmol/m**3',  'net community production')
   data%id_PPR = aed2_define_diag_variable('PPR','mmol/m**3/d','gross primary production rate')
   data%id_NPR = aed2_define_diag_variable('NPR','mmol/m**3/d','net community production rate')

   data%id_NUP = aed2_define_diag_variable('NUP','mmol/m**3/d','nitrogen uptake')
   data%id_PUP = aed2_define_diag_variable('PUP','mmol/m**3/d','phosphorous uptake')
   data%id_CUP = aed2_define_diag_variable('CUP','mmol/m**3/d','carbon uptake')

   data%id_dPAR = aed2_define_diag_variable('PAR','W/m**2',  'photosynthetically active radiation')
   data%id_TCHLA = aed2_define_diag_variable('TCHLA','ug/L', 'Total Chlorophyll-a')
   data%id_TPHY = aed2_define_diag_variable('TPHYS','ug/L',  'Total Phytoplankton')
   data%id_TIN = aed2_define_diag_variable('IN','ug/L',      'Total Chlorophyll-a')
   data%id_TIP = aed2_define_diag_variable('IP','ug/L',      'Total Chlorophyll-a')

   ! Register environmental dependencies
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_par = aed2_locate_global('par')
   data%id_I_0 = aed2_locate_global_sheet('par_sf')
   data%id_dz = aed2_locate_global('layer_ht')
   data%id_extc = aed2_locate_global('extc_coef')
END SUBROUTINE aed2_define_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_phytoplankton(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: phy, tphy, tin, tip, tchla
   AED_REAL :: INi, IPi
   AED_REAL :: pup
   AED_REAL :: no3up,nh4up
   AED_REAL :: cup, rsiup
   AED_REAL :: temp, par, Io, salinity, extc, dz
   AED_REAL :: primprod(data%num_phytos), exudation(data%num_phytos), &
               a_nfix(data%num_phytos), respiration(data%num_phytos)
   AED_REAL :: cuptake(data%num_phytos), cexcretion(data%num_phytos), cmortality(data%num_phytos)
   AED_REAL :: nuptake(data%num_phytos,1:4), nexcretion(data%num_phytos), nmortality(data%num_phytos)
   AED_REAL :: puptake(data%num_phytos,1:2), pexcretion(data%num_phytos), pmortality(data%num_phytos)
   AED_REAL :: siuptake(data%num_phytos), siexcretion(data%num_phytos), simortality(data%num_phytos)
   AED_REAL :: fT, fNit, fPho, fSil, fI, fXl, fSal, PNf
   AED_REAL :: upTot

   INTEGER  :: phy_i,c
   AED_REAL :: flux, available

! MH to fix
!  AED_REAL :: dt = 3600. ! just for now, hard code it
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)    ! local temperature
   salinity = _STATE_VAR_(data%id_sal)! local salinity
   par = _STATE_VAR_(data%id_par)     ! local photosynthetically active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation

   pup = 0.
   ! Retrieve current (local) state variable values.
   IF (data%do_Puptake)  pup = _STATE_VAR_(data%id_Pupttarget(1))

   no3up = 0.
   nh4up = 0.
   IF (data%do_Nuptake) THEN
       no3up = _STATE_VAR_(data%id_Nupttarget(1))
       nh4up = _STATE_VAR_(data%id_Nupttarget(2))
   ENDIF
   cup = 0.
   IF (data%do_Cuptake)  cup = _STATE_VAR_(data%id_Cupttarget)
   rsiup = 0.
   IF (data%do_Siuptake)  rsiup = _STATE_VAR_(data%id_Siupttarget)

   tphy = 0.0
   tchla = 0.0
   tin  = 0.0
   tip  = 0.0

   INi = 0.
   IPi = 0.

   DO phy_i=1,data%num_phytos

!     ! Available nutrients - this is set to the concentration of nutrients from
!     ! other aed2_modules unless they are limiting
!     no3up = sum(nuptake(phy_i, 1) / dtlim
!     nh4up = nuptake(phy_i, 2) / dtlim
!     pup   = puptake(phy_i, 1) / dtlim
!     cup   = cuptake(phy_i)    / dtlim
!     rsiup = siuptake(phy_i)   / dtlim

      primprod(phy_i)    = zero_
      exudation(phy_i)   = zero_
      a_nfix(phy_i)      = zero_
      respiration(phy_i) = zero_

      cuptake(phy_i)     = zero_
      cexcretion(phy_i)  = zero_
      cmortality(phy_i)  = zero_
      nuptake(phy_i,:)   = zero_
      nexcretion(phy_i)  = zero_
      nmortality(phy_i)  = zero_
      puptake(phy_i,:)   = zero_
      pexcretion(phy_i)  = zero_
      pmortality(phy_i)  = zero_

      ! Retrieve this phytoplankton group
      phy = _STATE_VAR_(data%id_p(phy_i))

      ! Get the temperature limitation function
      fT = fTemp_function(data%phytos(phy_i)%fT_Method,    &
                          data%phytos(phy_i)%T_max,        &
                          data%phytos(phy_i)%T_std,        &
                          data%phytos(phy_i)%theta_growth, &
                          data%phytos(phy_i)%aTn,          &
                          data%phytos(phy_i)%bTn,          &
                          data%phytos(phy_i)%kTn,temp)

      ! Get the light and nutrient limitation.
      ! NITROGEN.
      fNit = 0.0
      IF(data%phytos(phy_i)%simINDynamics /= 0) THEN
         ! IN variable available
         INi = _STATE_VAR_(data%id_in(phy_i))
      ELSE
         ! Assumed constant IN:
         INi = phy*data%phytos(phy_i)%X_ncon
      END IF

      ! Estimate fN limitation from IN or ext N value
      IF(data%phytos(phy_i)%simINDynamics > 1) THEN
         IF (phy > data%phytos(phy_i)%p0) THEN
            fNit = INi / phy
            fNit = phyto_fN(data%phytos,phy_i,IN=fNit)
         ENDIF
         IF (phy > zero_ .AND. phy <= data%phytos(phy_i)%p0) THEN
            fNit = phyto_fN(data%phytos,phy_i,din=no3up+nh4up)
         ENDIF
      ELSE
         fNit = phyto_fN(data%phytos,phy_i,din=no3up+nh4up)
      ENDIF
      IF (data%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixer: apply no N limitation. N Fixation ability
         ! depends on DIN concentration
         a_nfix = (one_ - fNit)
         fNit = one_
      ENDIF


      ! PHOSPHOROUS.
      fPho = zero_
      IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! IP variable available
         IPi = _STATE_VAR_(data%id_ip(phy_i))
      ELSE
         ! Assumed constant IP:
         IPi = phy*data%phytos(phy_i)%X_pcon
      END IF

      ! Estimate fP limitation from IP or ext P value
      IF (data%phytos(phy_i)%simIPDynamics > 1) THEN
         IF (phy > data%phytos(phy_i)%p0) THEN
            fPho = IPi / phy
            fPho = phyto_fP(data%phytos,phy_i,IP=fPho)
         ENDIF
         IF (phy > zero_ .AND. phy <= data%phytos(phy_i)%p0) THEN
            fPho = phyto_fP(data%phytos,phy_i,frp=pup)
         ENDIF
      ELSE
         fPho = phyto_fP(data%phytos,phy_i,frp=pup)
      ENDIF

      ! SILICA.
      fSil = phyto_fSi(data%phytos,phy_i,rsiup)


      ! LIGHT
      extc = _STATE_VAR_(data%id_extc)
      ! dz = 0.5     !MH: to fix
      dz = _STATE_VAR_(data%id_dz)
      fI = phyto_light(data%phytos, phy_i, par, extc, Io, dz)
      ! fI = 0.1


      ! METAL AND TOXIC EFFECTS
      fXl = 1.0

      ! Primary production rate
      primprod(phy_i) = data%phytos(phy_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl

      ! Adjust primary production rate for nitrogen fixers
      IF (data%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixing species, and the growth rate to  must be reduced
         ! to compensate for the increased metabolic cost of this process
         primprod(phy_i) = primprod(phy_i) * (data%phytos(phy_i)%k_nfix + &
                           (1.0-a_nfix(phy_i))*(1.0-data%phytos(phy_i)%k_nfix))
      ENDIF


      ! Respiration and general metabolic loss

      respiration(phy_i) = phyto_respiration(data%phytos,phy_i,temp)

      ! Salinity stress effect on respiration
      fSal =  phyto_salinity(data%phytos,phy_i,salinity)
      respiration(phy_i) = respiration(phy_i) * fSal

      ! photo-exudation
      exudation(phy_i) = primprod(phy_i)*data%phytos(phy_i)%f_pr

      ! Limit respiration if at the min biomass to prevent
      ! leak in the C mass balance
      IF (phy <= data%phytos(phy_i)%p0) THEN
         respiration(phy_i) = zero_
         exudation(phy_i) = zero_
      ENDIF

      ! write(*,"(4X,'limitations (fT,fI,fN,fP,fSi,Io, par, mu): ',9F9.2)")fT,fI,fNit,fPho,fSil,Io,par,primprod*secs_pr_day


      ! Carbon uptake and excretion

      cuptake(phy_i)    = -primprod(phy_i) * phy
      cexcretion(phy_i) = (data%phytos(phy_i)%k_fdom*(1.0-data%phytos(phy_i)%k_fres)*respiration(phy_i)+exudation(phy_i)) * phy
      cmortality(phy_i) = ((1.0-data%phytos(phy_i)%k_fdom)*(1.0-data%phytos(phy_i)%k_fres)*respiration(phy_i)) * phy

      ! Nitrogen uptake and excretion

      CALL phyto_internal_nitrogen(data%phytos,phy_i,data%do_N2uptake,phy,INi,primprod(phy_i),&
                             fT,no3up,nh4up,a_nfix(phy_i),respiration(phy_i),exudation(phy_i),PNf,&
                                   nuptake(phy_i,:),nexcretion(phy_i),nmortality(phy_i))

      ! Phosphorus uptake and excretion

      CALL phyto_internal_phosphorus(data%phytos,phy_i,data%npup,phy,IPi,primprod(phy_i),&
                                 fT,pup,respiration(phy_i),exudation(phy_i),&
                                         puptake(phy_i,:),pexcretion(phy_i),pmortality(phy_i))

      ! Silica uptake and excretion

      IF (data%phytos(phy_i)%simSiUptake > 0) THEN
         siuptake(phy_i)    =-data%phytos(phy_i)%X_sicon * primprod(phy_i) * phy
         siexcretion(phy_i) = data%phytos(phy_i)%X_sicon * (data%phytos(phy_i)%k_fdom*respiration(phy_i)+exudation(phy_i)) * phy
         simortality(phy_i) = data%phytos(phy_i)%X_sicon * ((1.0-data%phytos(phy_i)%k_fdom)*respiration(phy_i)) * phy
      ELSE
         siuptake(phy_i)    = zero_
         siexcretion(phy_i) = zero_
         simortality(phy_i) = zero_
      ENDIF

      ! Diagnostic info

      _DIAG_VAR_(data%id_NtoP(phy_i)) =  INi/IPi

      IF (extra_debug) THEN
         _DIAG_VAR_(data%id_fT(phy_i)) =  fT
         _DIAG_VAR_(data%id_fI(phy_i)) =  fI
         _DIAG_VAR_(data%id_fNit(phy_i)) =  fNit
         _DIAG_VAR_(data%id_fPho(phy_i)) =  fPho
         _DIAG_VAR_(data%id_fSil(phy_i)) =  fSil
         _DIAG_VAR_(data%id_fSal(phy_i)) =  fSal
      ENDIF
   END DO


   !-----------------------------------------------------------------
   ! Check uptake values for availability to prevent -ve numbers

   ! pup   - p available
   ! no3up - no3 available
   ! nh4up - nh4 available
   ! cup   - c available
   ! rsiup - Si available

   IF (data%do_Puptake) THEN
      upTot = sum(puptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= pup ) THEN
         DO phy_i=1,data%num_phytos
            puptake(phy_i,1) = (pup*0.99/dtlim) * (puptake(phy_i,1)/upTot)
         ENDDO
      ENDIF
   ENDIF

   IF (data%do_Nuptake) THEN
      upTot = sum(nuptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= no3up ) THEN
         DO phy_i=1,data%num_phytos
            nuptake(phy_i,1) = (no3up*0.99/dtlim) * (nuptake(phy_i,1)/upTot)
         ENDDO
      ENDIF

      upTot = sum(nuptake(:,2))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= nh4up ) THEN
         DO phy_i=1,data%num_phytos
            nuptake(phy_i,2) = (nh4up*0.99/dtlim) * (nuptake(phy_i,2)/upTot)
         ENDDO
      ENDIF
   ENDIF
   IF (data%do_Cuptake) THEN
      upTot = sum(cuptake)*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= cup ) THEN
         DO phy_i=1,data%num_phytos
            cuptake(phy_i) = (cup*0.99/dtlim) * (cuptake(phy_i)/upTot)
         ENDDO
      ENDIF
   ENDIF
!  IF (data%do_DOuptake) THEN
!     !
!  ENDIF
   IF (data%do_Siuptake) THEN
      upTot = sum(siuptake)*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= rsiup ) THEN
         DO phy_i=1,data%num_phytos
            siuptake(phy_i) = (rsiup*0.99/dtlim) * (siuptake(phy_i)/upTot)
         ENDDO
      ENDIF
   ENDIF

   DO phy_i=1,data%num_phytos

      ! Retrieve this phytoplankton group
      !-----------------------------------------------------------------
      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Phytoplankton production / losses
      phy = _STATE_VAR_(data%id_p(phy_i))
      flux = (primprod(phy_i) - respiration(phy_i) - exudation(phy_i)) * phy
      available = MAX(zero_, phy - data%phytos(phy_i)%p0)
      IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
      _FLUX_VAR_(data%id_p(phy_i)) = _FLUX_VAR_(data%id_p(phy_i)) + ( flux)

      IF (data%phytos(phy_i)%simINDynamics /= 0) THEN
         ! _FLUX_VAR_(data%id_in(phy_i)) = _FLUX_VAR_(data%id_in(phy_i)) + ( (-sum(nuptake) - nexcretion(phy_i) - nmortality(phy_i) )*INi )
         INi = _STATE_VAR_(data%id_in(phy_i))
         flux = (-sum(nuptake(phy_i,:)) - nexcretion(phy_i) - nmortality(phy_i) )
         available = MAX(zero_, INi - data%phytos(phy_i)%X_nmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_in(phy_i)) = _FLUX_VAR_(data%id_in(phy_i)) + ( flux)
      ENDIF
      IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! _FLUX_VAR_(data%id_ip(phy_i)) = _FLUX_VAR_(data%id_ip(phy_i)) + ( (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) ) )
         IPi = _STATE_VAR_(data%id_ip(phy_i))
         flux = (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) )
         available = MAX(zero_, IPi - data%phytos(phy_i)%X_pmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_ip(phy_i)) = _FLUX_VAR_(data%id_ip(phy_i)) + ( flux)
      ENDIF

      ! Now manage uptake of nutrients, CO2 and DO - these cumulative fluxes already limited above loop
      IF (data%do_Puptake) THEN
         DO c = 1,data%npup
            _FLUX_VAR_(data%id_Pupttarget(c)) = _FLUX_VAR_(data%id_Pupttarget(c)) + ( puptake(phy_i,c))
         ENDDO
      ENDIF
      IF (data%do_Nuptake) THEN
         DO c = 1,data%nnup
            _FLUX_VAR_(data%id_Nupttarget(c)) = _FLUX_VAR_(data%id_Nupttarget(c)) + ( nuptake(phy_i,c))
         ENDDO
      ENDIF
      IF (data%do_Cuptake) THEN
         _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) + (  cuptake(phy_i) - respiration(phy_i)*data%phytos(phy_i)%k_fres*phy )
      ENDIF
      IF (data%do_DOuptake) THEN
         _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) + ( -cuptake(phy_i) + respiration(phy_i)*data%phytos(phy_i)%k_fres*phy )
      ENDIF
      IF (data%do_Siuptake) THEN
         _FLUX_VAR_(data%id_Siupttarget) = _FLUX_VAR_(data%id_Siupttarget) + ( siuptake(phy_i))
      ENDIF
      ! Now manage mortality contributions to POM
      IF (data%do_Pmort) THEN
         _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + (pmortality(phy_i))
      ENDIF
      IF (data%do_Nmort) THEN
         _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + (nmortality(phy_i))
      ENDIF
      IF (data%do_Cmort) THEN
         _FLUX_VAR_(data%id_Cmorttarget) = _FLUX_VAR_(data%id_Cmorttarget) + (cmortality(phy_i))
      ENDIF
      IF (data%do_Simort) THEN
         _FLUX_VAR_(data%id_Simorttarget) = _FLUX_VAR_(data%id_Simorttarget) + (simortality(phy_i))
      ENDIF
      ! Now manage excretion/exudation contributions to DOM
      IF (data%do_Pexc) THEN
         _FLUX_VAR_(data%id_Pexctarget) = _FLUX_VAR_(data%id_Pexctarget) + (pexcretion(phy_i))
      ENDIF
      IF (data%do_Nexc) THEN
         _FLUX_VAR_(data%id_Nexctarget) = _FLUX_VAR_(data%id_Nexctarget) + (nexcretion(phy_i))
      ENDIF
      IF (data%do_Cexc) THEN
         _FLUX_VAR_(data%id_Cexctarget) = _FLUX_VAR_(data%id_Cexctarget) + (cexcretion(phy_i))
      ENDIF
      IF (data%do_Siexc) THEN
         _FLUX_VAR_(data%id_Siexctarget) = _FLUX_VAR_(data%id_Siexctarget) + (siexcretion(phy_i))
      ENDIF

      !-----------------------------------------------------------------
      ! export diagnostic variables

      ! Total phytoplankton carbon
      tphy = tphy + phy

      ! Total chlorophyll-a
      IF (data%phytos(phy_i)%Xcc > 0.1) THEN
        ! Assume Xcc (mol C/ mol chla) is a constant
        tchla = tchla + ( phy / data%phytos(phy_i)%Xcc ) * 12.0
      ELSE
        ! Use dynamic equation (Eq 13: of Baklouti, Cloern et al. 1995)
        ! theta = 1/Xcc [mg Chl (mg C)1] = 0.003 + 0.0154  e^0.050T  e^0.059E mu
        tchla = tchla + ( phy * (0.003 + 0.0154 * exp(0.050*temp) * exp(0.059*par) &
                        * primprod(phy_i)))
      ENDIF

      ! Total internal nutrients
      tin = tin + INi
      tip = tip + IPi

   ENDDO

   _DIAG_VAR_(data%id_GPP) =  sum(primprod)
   _DIAG_VAR_(data%id_NCP) =  sum(primprod - respiration)
   _DIAG_VAR_(data%id_PPR) =  sum(primprod) * secs_pr_day
   _DIAG_VAR_(data%id_NPR) =  sum(primprod - respiration) * secs_pr_day

   _DIAG_VAR_(data%id_NUP) =  sum(nuptake)
   _DIAG_VAR_(data%id_PUP) =  sum(puptake)
   _DIAG_VAR_(data%id_CUP) =  sum(cuptake)

   _DIAG_VAR_(data%id_dPAR) =  par
   _DIAG_VAR_(data%id_TCHLA) =  tchla
   _DIAG_VAR_(data%id_TPHY) =  tphy
   _DIAG_VAR_(data%id_TIN) =  tin
   _DIAG_VAR_(data%id_TIP) =  tip


END SUBROUTINE aed2_calculate_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_phytoplankton(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of phytoplankton.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: phy        ! State
   INTEGER  :: phy_i
   AED_REAL :: phy_flux

   ! Parameters
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN


   DO phy_i=1,data%num_phytos
      ! Retrieve current (local) state variable values.
      phy = _STATE_VAR_(data%id_p(phy_i))! phytoplankton

      phy_flux = zero_  !data%phytos(phy_i)%w_p*MAX(phy,zero_)

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     ! Transfer sediment flux value to AED2.
     _FLUX_VAR_(data%id_p(phy_i)) = _FLUX_VAR_(data%id_p(phy_i)) + (phy_flux)

   ENDDO

END SUBROUTINE aed2_calculate_benthic_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_mobility_phytoplankton(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: temp
   AED_REAL :: w_rho
   AED_REAL :: phy
   INTEGER  :: phy_i
!
!-------------------------------------------------------------------------------
!BEGIN

    temp = _STATE_VAR_(data%id_tem)

   DO phy_i=1,data%num_phytos
      SELECT CASE (data%phytos(phy_i)%w_model)
         CASE (_MOB_CONST_)
            w_rho = data%phytos(phy_i)%w_p
         CASE (_MOB_TEMP_)
            w_rho = data%phytos(phy_i)%w_p * temp !# MH to fix
         CASE (_MOB_CALC_)
            !# MH to complete
             phy = _STATE_VAR_(data%id_p(phy_i))! phytoplankton
             w_rho = data%phytos(phy_i)%w_p
             w_rho = _STATE_VAR_(data%id_rho(phy_i))
         CASE DEFAULT
             STOP
      END SELECT

      mobility(data%id_p(phy_i)) =  w_rho
   ENDDO

END SUBROUTINE aed2_mobility_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_phytoplankton(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: phy
   INTEGER  :: phy_i
!
!-----------------------------------------------------------------------
!BEGIN

   DO phy_i=1,data%num_phytos
      ! Retrieve current (local) state variable values.
      phy = _STATE_VAR_(data%id_p(phy_i))! phytoplankton

      ! Self-shading with explicit contribution from background phytoplankton concentration.
      extinction = extinction + (data%phytos(phy_i)%KePHY*phy)

   ENDDO


END SUBROUTINE aed2_light_extinction_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_phytoplankton
