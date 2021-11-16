!###############################################################################
!#                                                                             #
!# aed_pathogens.F90                                                           #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2017 - 2021 -  The University of Western Australia               #
!#                                                                             #
!#   AED2+ is free software: you can redistribute it and/or modify             #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED2+ is distributed in the hope that it will be useful,                  #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#   For use in a commercial seeting please contact the authors.               #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Originally created July 2012 by Matthew Hipsey, UWA                         #
!# Follow updates @ https://github.com/AquaticEcoDynamics/libaed2-plus         #
!#                                                                             #
!# For background to the module, please refer to the original implmentation    #
!# published here:                                                             #
!#                                                                             #
!# Hipsey, M.R., Antenucci, J.P. and Brookes, J.D., 2008. A generic,           #
!#       process based model of microbial pollution in aquatic systems.        #
!#       Water Resources Research, 44(7).                                      #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_pathogens
!-------------------------------------------------------------------------------
!  aed_pathogens --- pathogen biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_pathogens_data_t
!

   !----------------------------------------------------------------------------
!  %% NAMELIST    %%  pathogen_param_t
   TYPE pathogen_param_t
      CHARACTER(64) :: p_name
      AED_REAL      :: coef_grwth_uMAX                     !-- Max growth rate at 20C
      AED_REAL      :: coef_grwth_Tmin, coef_grwth_Tmax    !-- Tmin and Tmax, f(T)
      AED_REAL      :: coef_grwth_T1, coef_grwth_T2        !-- coef_grwth_T1  and  coef_grwth_T2
      AED_REAL      :: coef_grwth_Kdoc                     !-- Half-saturation for growth, coef_grwth_Kdoc
      AED_REAL      :: coef_grwth_ic                       !-- coef_grwth_ic
      AED_REAL      :: coef_mort_kd20                      !-- Mortality rate (Dark death rate) @ 20C and 0 psu
      AED_REAL      :: coef_mort_theta                     !-- Temperature multiplier for mortality: coef_mort_theta
      AED_REAL      :: coef_mort_c_SM, coef_mort_alpha, coef_mort_beta  !-- Salinity effect on mortality
      AED_REAL      :: coef_mort_c_PHM, coef_mort_K_PHM, coef_mort_delta_M  !-- pH effect on mortality
      AED_REAL      :: coef_mort_fdoc                      !-- Fraction of mortality back to doc
      AED_REAL      :: coef_light_kb_vis, coef_light_kb_uva, coef_light_kb_uvb !-- Light inactivation
      AED_REAL      :: coef_light_cSb_vis, coef_light_cSb_uva, coef_light_cSb_uvb !-- Salinity effect on light inactivation
      AED_REAL      :: coef_light_kDOb_vis, coef_light_kDOb_uva, coef_light_kDOb_uvb !-- DO effect on light
      AED_REAL      :: coef_light_cpHb_vis, coef_light_cpHb_uva, coef_light_cpHb_uvb !-- pH effect on light inactivation
      AED_REAL      :: coef_light_KpHb_vis, coef_light_KpHb_uva, coef_light_KpHb_uvb !-- pH effect on light inactivation
      AED_REAL      :: coef_light_delb_vis, coef_light_delb_uva, coef_light_delb_uvb !-- exponent for pH effect on light inactivation
      AED_REAL      :: coef_pred_kp20, coef_pred_theta_P   !-- Loss rate due to predation and temp multiplier
      AED_REAL      :: coef_sett_fa                        !-- Attached fraction in water column
      AED_REAL      :: coef_sett_w_path      !-- Sedimentation velocity (m/d) at 20C (-ve means down) for NON-ATTACHED orgs
      AED_REAL      :: coef_resus_epsilonP   !-- Pathogen resuspension rate
      AED_REAL      :: coef_resus_tauP_0     !-- Critical shear stress for organism resuspension
   END TYPE
!  %% END NAMELIST    %%  pathogen_param_t

!  TYPE pathogen_data
!     ! General Attributes
!     TYPE(pathogen_param_t) :: par
!  END TYPE

   TYPE,extends(aed_model_data_t) :: aed_pathogens_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_pf(:), id_pd(:)            ! Column ID of pathogens (alive and dead)
      INTEGER,ALLOCATABLE :: id_pa(:)                      ! Column ID of pathogens attached to ss
      INTEGER,ALLOCATABLE :: id_ps(:)                      ! Column ID of pathogens in sediment
      INTEGER,ALLOCATABLE :: id_ss(:)                      ! Column ID of ss if chosen

      ! Diagnostic IDs for processes
      INTEGER,ALLOCATABLE :: id_growth(:), id_mortality(:), id_sunlight(:), id_grazing(:), id_total(:)

      INTEGER  :: id_oxy, id_pH,  id_doc, id_tss           ! Dependency ID
      INTEGER  :: id_tem, id_sal                           ! Environmental IDs (3D)
      INTEGER  :: id_par, id_nir, id_uva, id_uvb           ! Environmental IDs (3D)
      INTEGER  :: id_I_0                                   ! Environmental ID (2D)
      INTEGER  :: resuspension
      INTEGER  :: id_epsilon, id_taub

      !# Model parameters
      INTEGER  :: num_pathogens
      TYPE(pathogen_param_t),DIMENSION(:),ALLOCATABLE :: pathogens
      INTEGER  :: num_ss
      AED_REAL :: tau_0_min, kTau_0
      AED_REAL,ALLOCATABLE :: epsilon(:), tau_0(:), tau_r(:), Ke_ss(:)
      AED_REAL,ALLOCATABLE :: epsilonP(:), tauP_0(:)
      AED_REAL,DIMENSION(:),ALLOCATABLE :: ss_set, ss_tau, ss_ke
      LOGICAL :: sim_sedorgs
      AED_REAL :: att_ts

     CONTAINS
         PROCEDURE :: define            => aed_define_pathogens
         PROCEDURE :: calculate         => aed_calculate_pathogens
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_pathogens
         PROCEDURE :: mobility          => aed_mobility_pathogens
         PROCEDURE :: light_extinction  => aed_light_extinction_pathogens
!        PROCEDURE :: delete            => aed_delete_pathogens
   END TYPE

!MODULE GLOBALS
   LOGICAL  :: extra_diag = .false.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_pathogens(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the pathogen biogeochemical model
!
!  Here, the aed_p_m namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_pathogens_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i

!  %% NAMELIST   %%  /aed_pathogens/
!  %% Last Checked 20/08/2021
   INTEGER  :: num_pathogens
   INTEGER  :: the_pathogens(MAX_PATHO_TYPES)
   INTEGER  :: num_ss = 0
   AED_REAL :: ss_set(MAX_PATHO_TYPES)=zero_
   AED_REAL :: ss_tau(MAX_PATHO_TYPES)=one_
   AED_REAL :: ss_ke(MAX_PATHO_TYPES) =zero_
   AED_REAL :: ss_initial = zero_
   AED_REAL :: epsilon(100)
   AED_REAL :: tau_0(100)
   AED_REAL :: tau_r(100)
   AED_REAL :: tau_0_min
   AED_REAL :: Ktau_0
   AED_REAL :: att_ts = (1./(secs_per_day*7.)) ! attachment fraction re-equilibrates over a week
   INTEGER  :: resuspension
   LOGICAL  :: sim_sedorgs = .FALSE.
   CHARACTER(len=64)  :: oxy_variable = ''
   CHARACTER(4) :: trac_name
   CHARACTER(len=128) :: dbase='aed_pathogen_pars.nml'

! From Module Globals
!  LOGICAL  :: extra_diag = .false.      !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_pathogens/

   NAMELIST /aed_pathogens/ num_pathogens, the_pathogens, resuspension, &
            num_ss, ss_set, ss_tau, ss_ke, sim_sedorgs, oxy_variable,    &
            epsilon, tau_0, tau_0_min, Ktau_0, dbase, extra_diag, att_ts,&
            diag_level
!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_pathogens initialization"

   ! Read the namelist
   read(namlst,nml=aed_pathogens,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_pathogens'

   IF ( extra_diag ) diag_level = 10

   data%tau_0_min = tau_0_min
   data%Ktau_0 = Ktau_0
   data%att_ts = att_ts

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   CALL aed_pathogens_load_params(data, dbase, num_pathogens, the_pathogens)

   IF ( num_ss > 0 ) THEN
      ALLOCATE(data%id_ss(num_ss))
      ALLOCATE(data%ss_set(num_ss)) ; data%ss_set(1:num_ss) = ss_set(1:num_ss)/secs_per_day
      ALLOCATE(data%ss_ke(num_ss))  ; data%ss_ke(1:num_ss)  = ss_ke(1:num_ss)
      ALLOCATE(data%ss_tau(num_ss)) ; data%ss_tau(1:num_ss) = ss_tau(1:num_ss)

      ALLOCATE(data%epsilon(num_ss)); data%epsilon(1:num_ss) = epsilon(1:num_ss)
      ALLOCATE(data%tau_0(num_ss))  ; data%tau_0(1:num_ss)   = tau_0(1:num_ss)

      trac_name = 'ss0'
      ! Register state variables
      DO i=1,num_ss
         trac_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_ss(i) = aed_define_variable(TRIM(trac_name),'g/m**3','path ss', &
                                                  ss_initial,minimum=zero_)
      ENDDO
   ENDIF
   IF ( resuspension == 2 ) THEN
      data%id_epsilon =  aed_define_sheet_diag_variable('epsilon','g/m**2/s', 'Resuspension rate')
   ENDIF

   ! Register state dependancies
   data%id_tss=-1 ;  data%id_doc=-1 ;  data%id_pH=-1;  data%id_oxy=-1
   IF (oxy_variable .NE. '') data%id_oxy = aed_locate_variable(oxy_variable)

   ! Register environmental dependencies
   data%id_tem = aed_locate_global('temperature')
   data%id_sal = aed_locate_global('salinity')
   data%id_par = aed_locate_global('par')
   data%id_nir = aed_locate_global('nir')
   data%id_uva = aed_locate_global('uva')
   data%id_uvb = aed_locate_global('uvb')
   data%id_I_0 = aed_locate_sheet_global('par_sf')
   IF ( resuspension > 0 ) &
      data%id_taub = aed_locate_sheet_global('taub')

   data%resuspension = resuspension
END SUBROUTINE aed_define_pathogens
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION load_csv(dbase, pd)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(pathogen_param_t) :: pd(MAX_PHYTO_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0
   LOGICAL :: meh
   INTEGER :: ret = 0
!
!BEGIN
!-------------------------------------------------------------------------------
   unit = aed_csv_read_header(dbase, csvnames, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   ALLOCATE(values(nccols))

   DO WHILE ( aed_csv_read_row(unit, values) )
      DO ccol=2,nccols
         pd(ccol)%p_name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('coef_grwth_uMAX')     ; pd(ccol)%coef_grwth_uMAX     = extract_double(values(ccol))
            CASE ('coef_grwth_Tmin')     ; pd(ccol)%coef_grwth_Tmin     = extract_double(values(ccol))
            CASE ('coef_grwth_Tmax')     ; pd(ccol)%coef_grwth_Tmax     = extract_double(values(ccol))
            CASE ('coef_grwth_T1')       ; pd(ccol)%coef_grwth_T1       = extract_double(values(ccol))
            CASE ('coef_grwth_T2')       ; pd(ccol)%coef_grwth_T2       = extract_double(values(ccol))
            CASE ('coef_grwth_Kdoc')     ; pd(ccol)%coef_grwth_Kdoc     = extract_double(values(ccol))
            CASE ('coef_grwth_ic')       ; pd(ccol)%coef_grwth_ic       = extract_double(values(ccol))
            CASE ('coef_mort_kd20')      ; pd(ccol)%coef_mort_kd20      = extract_double(values(ccol))
            CASE ('coef_mort_theta')     ; pd(ccol)%coef_mort_theta     = extract_double(values(ccol))
            CASE ('coef_mort_c_SM')      ; pd(ccol)%coef_mort_c_SM      = extract_double(values(ccol))
            CASE ('coef_mort_alpha')     ; pd(ccol)%coef_mort_alpha     = extract_double(values(ccol))
            CASE ('coef_mort_beta')      ; pd(ccol)%coef_mort_beta      = extract_double(values(ccol))
            CASE ('coef_mort_c_PHM')     ; pd(ccol)%coef_mort_c_PHM     = extract_double(values(ccol))
            CASE ('coef_mort_K_PHM')     ; pd(ccol)%coef_mort_K_PHM     = extract_double(values(ccol))
            CASE ('coef_mort_delta_M')   ; pd(ccol)%coef_mort_delta_M   = extract_double(values(ccol))
            CASE ('coef_mort_fdoc')      ; pd(ccol)%coef_mort_fdoc      = extract_double(values(ccol))
            CASE ('coef_light_kb_vis')   ; pd(ccol)%coef_light_kb_vis   = extract_double(values(ccol))
            CASE ('coef_light_kb_uva')   ; pd(ccol)%coef_light_kb_uva   = extract_double(values(ccol))
            CASE ('coef_light_kb_uvb')   ; pd(ccol)%coef_light_kb_uvb   = extract_double(values(ccol))
            CASE ('coef_light_cSb_vis')  ; pd(ccol)%coef_light_cSb_vis  = extract_double(values(ccol))
            CASE ('coef_light_cSb_uva')  ; pd(ccol)%coef_light_cSb_uva  = extract_double(values(ccol))
            CASE ('coef_light_cSb_uvb')  ; pd(ccol)%coef_light_cSb_uvb  = extract_double(values(ccol))
            CASE ('coef_light_kDOb_vis') ; pd(ccol)%coef_light_kDOb_vis = extract_double(values(ccol))
            CASE ('coef_light_kDOb_uva') ; pd(ccol)%coef_light_kDOb_uva = extract_double(values(ccol))
            CASE ('coef_light_kDOb_uvb') ; pd(ccol)%coef_light_kDOb_uvb = extract_double(values(ccol))
            CASE ('coef_light_cpHb_vis') ; pd(ccol)%coef_light_cpHb_vis = extract_double(values(ccol))
            CASE ('coef_light_cpHb_uva') ; pd(ccol)%coef_light_cpHb_uva = extract_double(values(ccol))
            CASE ('coef_light_cpHb_uvb') ; pd(ccol)%coef_light_cpHb_uvb = extract_double(values(ccol))
            CASE ('coef_light_KpHb_vis') ; pd(ccol)%coef_light_KpHb_vis = extract_double(values(ccol))
            CASE ('coef_light_KpHb_uva') ; pd(ccol)%coef_light_KpHb_uva = extract_double(values(ccol))
            CASE ('coef_light_KpHb_uvb') ; pd(ccol)%coef_light_KpHb_uvb = extract_double(values(ccol))
            CASE ('coef_light_delb_vis') ; pd(ccol)%coef_light_delb_vis = extract_double(values(ccol))
            CASE ('coef_light_delb_uva') ; pd(ccol)%coef_light_delb_uva = extract_double(values(ccol))
            CASE ('coef_light_delb_uvb') ; pd(ccol)%coef_light_delb_uvb = extract_double(values(ccol))
            CASE ('coef_pred_kp20')      ; pd(ccol)%coef_pred_kp20      = extract_double(values(ccol))
            CASE ('coef_pred_theta_P')   ; pd(ccol)%coef_pred_theta_P   = extract_double(values(ccol))
            CASE ('coef_sett_fa')        ; pd(ccol)%coef_sett_fa        = extract_double(values(ccol))
            CASE ('coef_sett_w_path')    ; pd(ccol)%coef_sett_w_path    = extract_double(values(ccol))
            CASE ('coef_resus_epsilonP') ; pd(ccol)%coef_resus_epsilonP = extract_double(values(ccol))
            CASE ('coef_resus_tauP_0')   ; pd(ccol)%coef_resus_tauP_0   = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit)
   !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pathogens_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pathogens_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in) :: count
   INTEGER,INTENT(in) :: list(*)
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,tfil
   AED_REAL :: minPath

   TYPE(pathogen_param_t) :: pd(MAX_PATHO_TYPES)
   NAMELIST /pathogen_data/ pd    ! %% pathogen_param_t - see above
!-------------------------------------------------------------------------------
!BEGIN
    minPath = 1e-5
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, pd)
       CASE (NML_TYPE)
           tfil = find_free_lun()
           open(tfil,file=dbase, status='OLD',iostat=status)
           IF (status /= 0) STOP 'Error opening namelist pathogen_data'
           read(tfil,nml=pathogen_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist pathogen_data'


    data%num_pathogens = count
    ALLOCATE(data%pathogens(count))
    ALLOCATE(data%id_pf(count))
    ALLOCATE(data%id_pd(count))
    ALLOCATE(data%id_pa(count))
    IF (data%sim_sedorgs) THEN
       ALLOCATE(data%id_ps(count))
    ENDIF
    ALLOCATE(data%id_total(count))
    IF ( diag_level >= 10 ) THEN
       ALLOCATE(data%id_growth(count))
       ALLOCATE(data%id_sunlight(count))
       ALLOCATE(data%id_mortality(count))
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%pathogens(i) = pd(list(i))
       data%pathogens(i)%coef_sett_w_path = data%pathogens(i)%coef_sett_w_path/secs_per_day
       data%pathogens(i)%coef_mort_kd20 = data%pathogens(i)%coef_mort_kd20/secs_per_day

       ! Register group as a state variable
       data%id_pf(i) = aed_define_variable(                                  &
                             TRIM(data%pathogens(i)%p_name)//'_f',            &
                             'orgs/m**3', 'pathogen alive',                   &
                             minPath,                                         &
                            ! pd(list(i))%p_initial,                          &
                             minimum=minPath,                                 &
                             !minimum=pd(list(i))%p0,                         &
                             mobility = data%pathogens(i)%coef_sett_w_path)


       ! Check if we need to registrer a variable for the attached fraction
       IF (data%pathogens(i)%coef_sett_fa > zero_) THEN
         data%id_pa(i) = aed_define_variable(                                &
                             TRIM(data%pathogens(i)%p_name)//'_a',            &
                             'orgs/m**3', 'pathogen attached',                &
                             minPath,                                         &
                            ! pd(list(i))%p_initial,                          &
                             minimum=minPath,                                 &
                             !minimum=pd(list(i))%p0,                         &
                             mobility = data%pathogens(i)%coef_sett_w_path)
       ENDIF

      !IF (data%pathogens(i)%p_name == 'crypto') THEN
         ! Register a state variable for dead fraction
         data%id_pd(i) = aed_define_variable(                               &
                             TRIM(data%pathogens(i)%p_name)//'_d',           &
                             'orgs/m**3', 'pathogen dead',                   &
                             zero_,                                          &
                             zero_,                                          &
                             mobility = data%pathogens(i)%coef_sett_w_path)

      !ENDIF


       IF (data%sim_sedorgs) THEN
          data%id_ps(i) = aed_define_sheet_variable( TRIM(data%pathogens(i)%p_name)//'_s', 'orgs/m2', 'pathogens in sediment')
          PRINT *,'WARNING: sim_sedorgs is not complete, and sediment population is not growing or resuspending'
       ENDIF

       data%id_total(i) = aed_define_diag_variable( TRIM(data%pathogens(i)%p_name)//'_t', 'orgs/m3', 'total')
       IF ( diag_level >= 10 ) THEN
          data%id_growth(i) = aed_define_diag_variable( TRIM(data%pathogens(i)%p_name)//'_g', 'orgs/m3/day', 'growth')
          data%id_sunlight(i) = aed_define_diag_variable( TRIM(data%pathogens(i)%p_name)//'_l', 'orgs/m3/day', 'sunlight')
          data%id_mortality(i) = aed_define_diag_variable( TRIM(data%pathogens(i)%p_name)//'_m', 'orgs/m3/day', 'mortality')
       ENDIF
   ENDDO
END SUBROUTINE aed_pathogens_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_pathogens(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of pathogen biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pathogens_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: pth_f, pth_a, pth_d
   AED_REAL :: temp,salinity,oxy,pH,doc
   AED_REAL :: Io,par,uva,uvb
   AED_REAL :: growth,light,mortality, predation, attachment
   AED_REAL :: f_AOC,f_pH,f_DO,phi,lightBW,phstar,att_frac

   INTEGER  :: pth_i

!-------------------------------------------------------------------------------
!BEGIN
   pth_d = 0.0 ; pth_a = 0.0

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)     ! local temperature
   salinity = _STATE_VAR_(data%id_sal) ! local salinity
   IF (data%id_oxy>0) THEN  ! & use_oxy
      oxy = _STATE_VAR_(data%id_oxy)   ! local oxygen
   ELSE
      oxy = 10.0 !mg/L
   ENDIF

   !doc = _STATE_VAR_(data%id_doc)    ! local DOC
   !ph = _STATE_VAR_(data%id_ph)      ! local pH
   phstar = 0.0                       ! abs(ph-7.)

   ! Get light bandwidth intensities
   Io = _STATE_VAR_S_(data%id_I_0)    ! surface short wave radiation
   par = _STATE_VAR_(data%id_par)     ! local photosynthetically active radiation (45% of sw)
   IF ( data%id_uva > 0 ) THEN
      uva = _STATE_VAR_(data%id_uva)
   ELSE
      uva = (par/0.45)*0.03           ! uva is 3% of sw (Kirk 1994)
   ENDIF
   IF ( data%id_uvb > 0 ) THEN
      uvb = _STATE_VAR_(data%id_uvb)
   ELSE
      uvb = (par/0.45)*0.003          ! uvb is 0.3% of sw
   ENDIF

   DO pth_i=1,data%num_pathogens

      ! Retrieve this pathogen group
      pth_f = _STATE_VAR_(data%id_pf(pth_i))
      IF (data%pathogens(pth_i)%coef_sett_fa > zero_) THEN
        pth_a = _STATE_VAR_(data%id_pa(pth_i))
      END IF

      growth    = zero_
      predation = zero_

      ! Natural mortality (as impacted by T, S, pH; see Hipsey et al 2008)
      f_AOC = 1.0 ! aoc / (K_AOC + aoc)
      f_pH  = 1.0 ! + c_PH * ( pH_star**delta / (pH_star**delta+K_PH**delta) )
      mortality = data%pathogens(pth_i)%coef_mort_kd20  &
                + (data%pathogens(pth_i)%coef_mort_c_SM*salinity**data%pathogens(pth_i)%coef_mort_alpha) &
                * ((1.0-f_AOC)**data%pathogens(pth_i)%coef_mort_beta) * f_pH
      mortality = mortality * (data%pathogens(pth_i)%coef_mort_theta**(temp-20.0))


      ! Sunlight inactivation (as impacted by S, DO and pH; see Hipsey et al 2008)
      light     = zero_
      lightBW   = zero_
      phi  = 1e-6  ! Convert J to MJ as kb is in m2/MJ)
      ! Visible
      f_DO = oxy / (data%pathogens(pth_i)%coef_light_kDOb_vis + oxy)
      f_pH = 1.0
           !(1.0+coef_light_cpHb_vis*(pH_star**coef_light_delb_vis/(coef_light_KpHb_vis**coef_light_delb_vis+pH_star**coef_light_delb_vis)))
      lightBW = phi * (data%pathogens(pth_i)%coef_light_kb_vis + data%pathogens(pth_i)%coef_light_cSb_vis*salinity)
      lightBW = lightBW * par * f_pH * f_DO
      light     = light + lightBW
      ! UV-A
      f_DO = oxy / (data%pathogens(pth_i)%coef_light_kDOb_uva + oxy)
      f_pH = 1.0
           !(1.0 + coef_light_cpHb_uva*(pH_star**coef_light_delb_uva/(coef_light_KpHb_uva**coef_light_delb_uva+pH_star**coef_light_delb_uva)))
      lightBW = phi * (data%pathogens(pth_i)%coef_light_kb_uva + data%pathogens(pth_i)%coef_light_cSb_uva*salinity)
      lightBW = lightBW * uva * f_pH * f_DO
      light     = light + lightBW
      ! UV-B
      f_DO = oxy / (data%pathogens(pth_i)%coef_light_kDOb_uvb + oxy)
      f_pH = 1.0
           !(1.0 + coef_light_cpHb_uvb*(pH_star**coef_light_delb_uvb/(coef_light_KpHb_uvb**coef_light_delb_uvb+pH_star**coef_light_delb_uvb)))
      lightBW = phi * (data%pathogens(pth_i)%coef_light_kb_uvb + data%pathogens(pth_i)%coef_light_cSb_uvb*salinity)
      lightBW = lightBW * uvb * f_pH * f_DO
      light   = light + lightBW

      ! Attachment of free orgs to particles (as impacted by SS and desired attachment ratio)
      attachment = zero_
      IF (data%pathogens(pth_i)%coef_sett_fa > zero_ .AND. (pth_f+pth_a) > 1e-2) THEN
         ! First check if ratio at last time step is less than desired (ie coef_sett_fa)
         att_frac = pth_a/(pth_a+pth_f)
         IF (att_frac<data%pathogens(pth_i)%coef_sett_fa) THEN
            ! Assume rate of attachment is slow based on att_ts (orgs/m3/s)
            ! Small fraction, so the "attachment deficit" is redistributed each time step
            attachment = data%att_ts * (data%pathogens(pth_i)%coef_sett_fa*pth_a - pth_f)
            IF(attachment>zero_ .AND. attachment > _FLUX_VAR_(data%id_pf(pth_i)))THEN
              ! proposed attachment is more than is available.
              attachment = 0.5 * _FLUX_VAR_(data%id_pf(pth_i))
            ELSEIF(attachment<zero_ .AND. attachment > _FLUX_VAR_(data%id_pa(pth_i)))THEN
              ! proposed de-attachment is more than is available.
              attachment = 0.5 * _FLUX_VAR_(data%id_pa(pth_i))
            ENDIF
         ENDIF
      ENDIF

      !-----------------------------------------------------------------
      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Pathogen production / losses
      _FLUX_VAR_(data%id_pf(pth_i)) = _FLUX_VAR_(data%id_pf(pth_i)) +          &
                                      ( (growth - light - mortality - predation)*pth_f - attachment)
      _FLUX_VAR_(data%id_pd(pth_i)) = _FLUX_VAR_(data%id_pd(pth_i)) +          &
                                      ( ( light + mortality + predation)*pth_f)
      ! Now the separate attached pathogen fraction
      IF (data%pathogens(pth_i)%coef_sett_fa > zero_) THEN
         _FLUX_VAR_(data%id_pa(pth_i)) = _FLUX_VAR_(data%id_pa(pth_i)) +      &
                                      ( (growth - light/2. - mortality )*pth_a + attachment )
         _FLUX_VAR_(data%id_pd(pth_i)) = _FLUX_VAR_(data%id_pd(pth_i)) +      &
                                      ((light/2. + mortality)*pth_a)
      ENDIF


      !-----------------------------------------------------------------
      ! SET DIAGNOSTICS
      _DIAG_VAR_(data%id_total(pth_i)) =  pth_f + pth_a + pth_d  ! orgs/m3/s
      IF ( diag_level >= 10 ) THEN
         _DIAG_VAR_(data%id_growth(pth_i)) =  growth*(pth_f + pth_a)  ! orgs/m3/s
         _DIAG_VAR_(data%id_sunlight(pth_i)) =  light*pth_f + (light/2.)*pth_a  ! orgs/m3/s
         _DIAG_VAR_(data%id_mortality(pth_i)) =  mortality*(pth_f + pth_a)  ! orgs/m3/s
      ENDIF
   ENDDO
END SUBROUTINE aed_calculate_pathogens
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_pathogens(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of pathogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pathogens_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: ss       ! State
   INTEGER  :: pth_i, ss_i
   AED_REAL :: resus_flux, sett_flux, kin_flux
   AED_REAL :: pth_sed, pth_wat_free, pth_wat_att

   AED_REAL :: bottom_stress, dummy_tau

   ! Parameters
!
!-------------------------------------------------------------------------------
!BEGIN
   resus_flux = 0.0   !## CAB [-Wmaybe-uninitialized]

!  vel = _STATE_VAR_(data%id_vel)

   IF ( data%resuspension  > 0) THEN
      bottom_stress = _STATE_VAR_S_(data%id_taub)
      bottom_stress = MIN(bottom_stress, 100.)
   ENDIF

   ! SEDIMENT PARTICULATE RESUSPENSION - WHAT PATHOGENS ARE ATTACHED TO
   IF (data%num_ss>0) THEN
      DO ss_i=1,data%num_ss
         ss = _STATE_VAR_(data%id_ss(ss_i))     ! ss conc of ith group

         IF ( data%resuspension > 0 ) THEN
            IF (bottom_stress > data%tau_0(ss_i)) THEN
               resus_flux = data%epsilon(ss_i) * ( bottom_stress - data%tau_0(ss_i)) / data%tau_r(ss_i)
            ELSE
               resus_flux = 0.
            ENDIF
         ENDIF

         _FLUX_VAR_(data%id_ss(ss_i)) = _FLUX_VAR_(data%id_ss(ss_i)) + resus_flux
      ENDDO
   ENDIF

   ! SEDIMENT ORGANISMS - KINETICS AND RESUSPENSION
   IF (data%sim_sedorgs) THEN

      ! Dynamic sediment pool of organisms that increase and decrease
      sett_flux = zero_

      DO pth_i=1,data%num_pathogens
         ! Retrieve current (local) state variable values.
         pth_wat_free = _STATE_VAR_(data%id_pf(pth_i))     ! pathogen (benthic water - free)
         IF (data%pathogens(pth_i)%coef_sett_fa>zero_) THEN
            pth_wat_att  = _STATE_VAR_(data%id_pa(pth_i))  ! pathogen (benthic water - attached)
         END IF

         pth_sed = _STATE_VAR_S_(data%id_ps(pth_i)) ! pathogen (sediment pool)


         ! Compute the resuspension flux from the sediment to water
         resus_flux = zero_  !data%pathogens(pth_i)%w_p*MAX(pth,zero_)

         ! Compute the settling flux from the  water to sediment
         sett_flux = data%pathogens(pth_i)%coef_sett_w_path * pth_wat_free
         IF (data%pathogens(pth_i)%coef_sett_fa>zero_)  THEN
            sett_flux = sett_flux + data%ss_set(1)  * pth_wat_att
         ENDIF

         ! Grwoth/death
         kin_flux = 0.0;

         ! Org flux to / from the sediment (orgs/m2/s)
         _FLUX_VAR_B_(data%id_ps(pth_i)) = _FLUX_VAR_B_(data%id_ps(pth_i)) + sett_flux - resus_flux + kin_flux

         ! Add to respective pools in water (free/attached)
         _FLUX_VAR_(data%id_pf(pth_i)) = _FLUX_VAR_(data%id_pf(pth_i)) +      &
                                                        (resus_flux)*(1.-data%pathogens(pth_i)%coef_sett_fa)
         IF (data%pathogens(pth_i)%coef_sett_fa>zero_) THEN
            _FLUX_VAR_(data%id_pa(pth_i)) = _FLUX_VAR_(data%id_pa(pth_i)) +      &
                                                        (resus_flux)*data%pathogens(pth_i)%coef_sett_fa
         ENDIF
      ENDDO
   ELSE

      ! Unresolved sediment pool of organisms, but we will still predict generic resuspension flux
      DO pth_i=1,data%num_pathogens

         ! Compute the resuspension flux from the sediment to water
         IF ( data%resuspension > 0 ) THEN
            resus_flux = 0.
            IF (bottom_stress > data%tauP_0(pth_i)) THEN
               ! Pathogens need to resuspend
               resus_flux = data%epsilonP(pth_i) * ( bottom_stress - data%tauP_0(pth_i)) / data%tau_r(1)

              ! Add to respective pools in water (free/attached)
              IF (bottom_stress > data%tau_0(1) .AND. data%pathogens(pth_i)%coef_sett_fa>zero_) THEN
                ! Attached and free organisms lifting in association with uplifted ss
                _FLUX_VAR_(data%id_pf(pth_i)) = _FLUX_VAR_(data%id_pf(pth_i)) +      &
                                                        (resus_flux)*(1.-data%pathogens(pth_i)%coef_sett_fa)
                _FLUX_VAR_(data%id_pa(pth_i)) = _FLUX_VAR_(data%id_pa(pth_i)) +      &
                                                        (resus_flux)*data%pathogens(pth_i)%coef_sett_fa
              ELSE
                ! Free organisms only (surfical organisms lifting before sediment uplift)
                _FLUX_VAR_(data%id_pf(pth_i)) = _FLUX_VAR_(data%id_pf(pth_i)) + (resus_flux)
              ENDIF
            ENDIF
         ENDIF

      ENDDO
   ENDIF

END SUBROUTINE aed_calculate_benthic_pathogens
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_pathogens(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pathogens_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: temp, vel
   INTEGER  :: ss_i,pth_i
!
!-------------------------------------------------------------------------------
!BEGIN
   temp = _STATE_VAR_(data%id_tem)

   ! First set velocity for free pathogen groups
   DO pth_i=1,data%num_pathogens
      mobility(data%id_pf(pth_i)) =  data%pathogens(pth_i)%coef_sett_w_path
   ENDDO

   ! Compute settling rate of particles
   DO ss_i=1,data%num_ss
      ! Update the settling rate and assign to mobility array
      mobility(data%id_ss(ss_i)) = data%ss_set(ss_i)  !Stokes()
   ENDDO

   ! Set velocity of attached pathogens, if simulated.
   !##NB Currently attached all assume on SS1
   DO pth_i=1,data%num_pathogens
      IF (data%pathogens(pth_i)%coef_sett_fa>zero_)  THEN
         mobility(data%id_pa(pth_i)) =  data%ss_set(1)
      ENDIF
   ENDDO
END SUBROUTINE aed_mobility_pathogens
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_pathogens(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to ss variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pathogens_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: ss
   INTEGER  :: ss_i
!
!-----------------------------------------------------------------------
!BEGIN
   DO ss_i=1,data%num_ss
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(ss_i)) ! ss conc form last timestep

      ! Self-shading with explicit contribution from background phytoplankton concentration.
      ! m^-1 = (m^-1)/(g/m3) * (g/m3)
      extinction = extinction + (data%ss_ke(ss_i)*ss)
   ENDDO
END SUBROUTINE aed_light_extinction_pathogens
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_pathogens
