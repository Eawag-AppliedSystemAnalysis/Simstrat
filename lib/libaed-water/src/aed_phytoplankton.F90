!###############################################################################
!#                                                                             #
!# aed_phytoplankton.F90                                                       #
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
!# Created August 2011                                                         #
!#                                                                             #
!#  Track changes on GitHub @ https://github.com/AquaticEcoDynamics/libaed-water
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |   ______     | || |  ____  ____  | || |  ____  ____  | |         !
!         | |  |_   __ \   | || | |_   ||   _| | || | |_  _||_  _| | |         !
!         | |    | |__) |  | || |   | |__| |   | || |   \ \  / /   | |         !
!         | |    |  ___/   | || |   |  __  |   | || |    \ \/ /    | |         !
!         | |   _| |_      | || |  _| |  | |_  | || |    _|  |_    | |         !
!         | |  |_____|     | || | |____||____| | || |   |______|   | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

MODULE aed_phytoplankton
!-------------------------------------------------------------------------------
!  aed_phytoplankton --- phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util
   USE aed_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private

   PUBLIC aed_phytoplankton_data_t


   TYPE,extends(aed_model_data_t) :: aed_phytoplankton_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_p(:)
      INTEGER,ALLOCATABLE :: id_in(:)
      INTEGER,ALLOCATABLE :: id_ip(:)
      INTEGER,ALLOCATABLE :: id_rho(:)
      INTEGER,ALLOCATABLE :: id_NtoP(:)
      INTEGER,ALLOCATABLE :: id_vvel(:)
      INTEGER,ALLOCATABLE :: id_fT(:), id_fI(:), id_fNit(:), &
                             id_fPho(:), id_fSil(:), id_fSal(:)
      INTEGER :: id_Pexctarget,id_Pmorttarget,id_Pupttarget(1:2)
      INTEGER :: id_Nexctarget,id_Nmorttarget,id_Nupttarget(1:4)
      INTEGER :: id_Cexctarget,id_Cmorttarget,id_Cupttarget
      INTEGER :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      INTEGER :: id_DOupttarget, id_l_resus, id_Psed_phy
      INTEGER :: id_par, id_I_0, id_extc, id_sedzone
      INTEGER :: id_tem, id_sal, id_dz, id_dens
      INTEGER :: id_GPP, id_NCP, id_PPR, id_NPR, id_dPAR
      INTEGER :: id_TPHY, id_TCHLA, id_TIN, id_TIP
      INTEGER :: id_MPB, id_d_MPB, id_d_BPP, id_d_BCP, id_d_mpbv
      INTEGER :: id_NUP, id_NUP2, id_PUP, id_CUP

      !# Model parameters
      INTEGER  :: num_phytos
      TYPE(phyto_data_t),DIMENSION(:),ALLOCATABLE :: phytos
      ! LOGICAL  :: do_exc,do_mort,do_upt, do_N2uptake
      LOGICAL  :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL  :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL  :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL  :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER  :: do_mpb, n_zones
      AED_REAL :: R_mpbg, R_mpbr, I_Kmpb, mpb_max, theta_mpb_growth, theta_mpb_resp
      INTEGER  :: nnup, npup
      AED_REAL :: dic_per_n
      AED_REAL :: min_rho,max_rho
      AED_REAL,ALLOCATABLE :: resuspension(:), active_zones(:)

     CONTAINS
         PROCEDURE :: define            => aed_define_phytoplankton
         PROCEDURE :: calculate         => aed_calculate_phytoplankton
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_phytoplankton
         PROCEDURE :: mobility          => aed_mobility_phytoplankton
         PROCEDURE :: light_extinction  => aed_light_extinction_phytoplankton
!        PROCEDURE :: delete            => aed_delete_phytoplankton

   END TYPE

! MODULE GLOBALS
   AED_REAL :: dtlim = 0.9 * 3600
   LOGICAL  :: extra_diag = .false.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS


!###############################################################################
INTEGER FUNCTION load_csv(dbase,pd)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(phyto_param_t) :: pd(MAX_PHYTO_TYPES)
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
            CASE ('p0')            ; pd(ccol)%p0            = extract_double(values(ccol))
            CASE ('w_p')           ; pd(ccol)%w_p           = extract_double(values(ccol))
            CASE ('Xcc')           ; pd(ccol)%Xcc           = extract_double(values(ccol))
            CASE ('R_growth')      ; pd(ccol)%R_growth      = extract_double(values(ccol))
            CASE ('fT_Method')     ; pd(ccol)%fT_Method     = extract_integer(values(ccol))
            CASE ('theta_growth')  ; pd(ccol)%theta_growth  = extract_double(values(ccol))
            CASE ('T_std')         ; pd(ccol)%T_std         = extract_double(values(ccol))
            CASE ('T_opt')         ; pd(ccol)%T_opt         = extract_double(values(ccol))
            CASE ('T_max')         ; pd(ccol)%T_max         = extract_double(values(ccol))
            CASE ('lightModel')    ; pd(ccol)%lightModel    = extract_integer(values(ccol))
            CASE ('I_K')           ; pd(ccol)%I_K           = extract_double(values(ccol))
            CASE ('I_S')           ; pd(ccol)%I_S           = extract_double(values(ccol))
            CASE ('KePHY')         ; pd(ccol)%KePHY         = extract_double(values(ccol))
            CASE ('f_pr')          ; pd(ccol)%f_pr          = extract_double(values(ccol))
            CASE ('R_resp')        ; pd(ccol)%R_resp        = extract_double(values(ccol))
            CASE ('theta_resp')    ; pd(ccol)%theta_resp    = extract_double(values(ccol))
            CASE ('k_fres')        ; pd(ccol)%k_fres        = extract_double(values(ccol))
            CASE ('k_fdom')        ; pd(ccol)%k_fdom        = extract_double(values(ccol))
            CASE ('salTol')        ; pd(ccol)%salTol        = extract_integer(values(ccol))
            CASE ('S_bep')         ; pd(ccol)%S_bep         = extract_double(values(ccol))
            CASE ('S_maxsp')       ; pd(ccol)%S_maxsp       = extract_double(values(ccol))
            CASE ('S_opt')         ; pd(ccol)%S_opt         = extract_double(values(ccol))
            CASE ('simDINUptake')  ; pd(ccol)%simDINUptake  = extract_integer(values(ccol))
            CASE ('simDONUptake')  ; pd(ccol)%simDONUptake  = extract_integer(values(ccol))
            CASE ('simNFixation')  ; pd(ccol)%simNFixation  = extract_integer(values(ccol))
            CASE ('simINDynamics') ; pd(ccol)%simINDynamics = extract_integer(values(ccol))
            CASE ('N_o')           ; pd(ccol)%N_o           = extract_double(values(ccol))
            CASE ('K_N')           ; pd(ccol)%K_N           = extract_double(values(ccol))
            CASE ('X_ncon')        ; pd(ccol)%X_ncon        = extract_double(values(ccol))
            CASE ('X_nmin')        ; pd(ccol)%X_nmin        = extract_double(values(ccol))
            CASE ('X_nmax')        ; pd(ccol)%X_nmax        = extract_double(values(ccol))
            CASE ('R_nuptake')     ; pd(ccol)%R_nuptake     = extract_double(values(ccol))
            CASE ('k_nfix')        ; pd(ccol)%k_nfix        = extract_double(values(ccol))
            CASE ('R_nfix')        ; pd(ccol)%R_nfix        = extract_double(values(ccol))
            CASE ('simDIPUptake')  ; pd(ccol)%simDIPUptake  = extract_integer(values(ccol))
            CASE ('simIPDynamics') ; pd(ccol)%simIPDynamics = extract_integer(values(ccol))
            CASE ('P_0')           ; pd(ccol)%P_0           = extract_double(values(ccol))
            CASE ('K_P')           ; pd(ccol)%K_P           = extract_double(values(ccol))
            CASE ('X_pcon')        ; pd(ccol)%X_pcon        = extract_double(values(ccol))
            CASE ('X_pmin')        ; pd(ccol)%X_pmin        = extract_double(values(ccol))
            CASE ('X_pmax')        ; pd(ccol)%X_pmax        = extract_double(values(ccol))
            CASE ('R_puptake')     ; pd(ccol)%R_puptake     = extract_double(values(ccol))
            CASE ('simSiUptake')   ; pd(ccol)%simSiUptake   = extract_integer(values(ccol))
            CASE ('Si_0')          ; pd(ccol)%Si_0          = extract_double(values(ccol))
            CASE ('K_Si')          ; pd(ccol)%K_Si          = extract_double(values(ccol))
            CASE ('X_sicon')       ; pd(ccol)%X_sicon       = extract_double(values(ccol))

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
SUBROUTINE aed_phytoplankton_load_params(data, dbase, count, list, settling, resuspension)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_phytoplankton_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: settling(*)
   AED_REAL,INTENT(in)         :: resuspension(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(phyto_param_t) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /phyto_data/ pd     ! %% type phyto_param_t - see aed_bio_utils
!-------------------------------------------------------------------------------
!BEGIN
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, pd)
       CASE (NML_TYPE)
           tfil = find_free_lun()
           open(tfil,file=dbase, status='OLD',iostat=status)
           IF (status /= 0) STOP 'Cannot open phyto_data namelist file'
           read(tfil,nml=phyto_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist phyto_data'

    data%num_phytos = count
    ALLOCATE(data%phytos(count))
    ALLOCATE(data%id_p(count)) ; data%id_p(:) = 0
    ALLOCATE(data%id_in(count)) ; data%id_in(:) = 0
    ALLOCATE(data%id_ip(count)) ; data%id_ip(:) = 0
    ALLOCATE(data%id_rho(count)) ; data%id_rho(:) = 0
    ALLOCATE(data%id_NtoP(count)) ; data%id_NtoP(:) = 0
    IF ( diag_level >= 10 ) THEN
       ALLOCATE(data%id_fT(count)) ; data%id_fT(:) = 0
       ALLOCATE(data%id_fI(count)) ; data%id_fI(:) = 0
       ALLOCATE(data%id_fNit(count)) ; data%id_fNit(:) = 0
       ALLOCATE(data%id_fPho(count)) ; data%id_fPho(:) = 0
       ALLOCATE(data%id_fSil(count)) ; data%id_fSil(:) = 0
       ALLOCATE(data%id_fSal(count)) ; data%id_fSal(:) = 0
       ALLOCATE(data%id_vvel(count)) ; data%id_vvel(:) = 0
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%phytos(i)%p_name       = pd(list(i))%p_name
       data%phytos(i)%p0           = pd(list(i))%p0
       data%phytos(i)%w_p          = pd(list(i))%w_p/secs_per_day
       data%phytos(i)%settling     = settling(i)
       data%phytos(i)%resuspension = resuspension(i)
       data%phytos(i)%Xcc          = pd(list(i))%Xcc
       data%phytos(i)%R_growth     = pd(list(i))%R_growth/secs_per_day
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
       data%phytos(i)%R_resp       = pd(list(i))%R_resp/secs_per_day
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
       data%phytos(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_per_day
       data%phytos(i)%k_nfix       = pd(list(i))%k_nfix
       data%phytos(i)%R_nfix       = pd(list(i))%R_nfix/secs_per_day
       data%phytos(i)%simDIPUptake = pd(list(i))%simDIPUptake
       data%phytos(i)%simIPDynamics= pd(list(i))%simIPDynamics
       data%phytos(i)%P_0          = pd(list(i))%P_0
       data%phytos(i)%K_P          = pd(list(i))%K_P
       data%phytos(i)%X_pcon       = pd(list(i))%X_pcon
       data%phytos(i)%X_pmin       = pd(list(i))%X_pmin
       data%phytos(i)%X_pmax       = pd(list(i))%X_pmax
       data%phytos(i)%R_puptake    = pd(list(i))%R_puptake/secs_per_day
       data%phytos(i)%simSiUptake  = pd(list(i))%simSiUptake
       data%phytos(i)%Si_0         = pd(list(i))%Si_0
       data%phytos(i)%K_Si         = pd(list(i))%K_Si
       data%phytos(i)%X_sicon      = pd(list(i))%X_sicon

       data%phytos(i)%c1           = 0.0124/60.   ! From Chung et al (2014)
       data%phytos(i)%c3           = 0.0230/60.   !  "
       data%phytos(i)%f1           = 0.675        ! Ross and Sharples (2007)
       data%phytos(i)%f2           = 0.750        !  "
       data%phytos(i)%d_phy        = 1e-5

       ! Register group as a state variable
       data%id_p(i) = aed_define_variable(                                    &
                              TRIM(data%phytos(i)%p_name),                     &
                              'mmol/m**3',                                     &
                              'phytoplankton '//TRIM(data%phytos(i)%p_name),   &
                              pd(list(i))%p_initial,                           &
                              minimum=pd(list(i))%p0,                          &
                              maximum=1e4,                          &
                              mobility = data%phytos(i)%w_p)

       ! Register rho (density) group as a state variable, if required
       IF (data%phytos(i)%settling == _MOB_STOKES_) THEN
          data%id_rho(i) = aed_define_variable(                               &
                              TRIM(data%phytos(i)%p_name)//'_rho',             &
                              'kg/m**3',                                       &
                        'phytoplankton '//TRIM(data%phytos(i)%p_name)//'_rho', &
                              (data%min_rho+data%max_rho)/2.,                  &
                              minimum=data%min_rho,                            &
                              maximum=data%max_rho,                            &
                              mobility = data%phytos(i)%w_p)
       ENDIF

       ! Register internal nitrogen group as a state variable, if required
       IF (data%phytos(i)%simINDynamics /= 0) THEN
          IF(data%phytos(i)%simINDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_ncon
          ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_nmin
          ENDIF
          ! Register IN group as a state variable
          data%id_in(i) = aed_define_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IN',              &
                              'mmol/m**3',                                     &
                         'phytoplankton '//TRIM(data%phytos(i)%p_name)//'_IN', &
                              pd(list(i))%p_initial*data%phytos(i)%X_ncon,     &
                              minimum=minNut,                                  &
                              mobility = data%phytos(i)%w_p)

       ELSE
         IF (data%phytos(i)%settling == _MOB_MOTILE_) THEN
             data%phytos(i)%settling = _MOB_CONST_
             print *,'Motility can not be simulated as no IN'
         ENDIF
       ENDIF

       ! Register internal phosphorus group as a state variable, if required
       IF (data%phytos(i)%simIPDynamics /= 0) THEN
          IF(data%phytos(i)%simIPDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_pcon
          ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_pmin
          ENDIF
          ! Register IP group as a state variable
          data%id_ip(i) = aed_define_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IP',              &
                              'mmol/m**3',                                     &
                         'phytoplankton '//TRIM(data%phytos(i)%p_name)//'_IP', &
                              pd(list(i))%p_initial*data%phytos(i)%X_pcon,     &
                              minimum=minNut,                                  &
                              mobility = data%phytos(i)%w_p)
       ENDIF

       ! Group specific diagnostic variables
       data%id_NtoP(i) = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_NtoP','-', 'internal n:p ratio')
       IF ( diag_level >= 10 ) THEN
          data%id_fI(i)   = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fI', '-', 'fI (0-1)')
          data%id_fNit(i) = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fNit', '-', 'fNit (0-1)')
          data%id_fPho(i) = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fPho', '-', 'fPho (0-1)')
          data%id_fSil(i) = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSil', '-', 'fSil (0-1)')
          data%id_fT(i)   = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fT', '-', 'fT (>0)')
          data%id_fSal(i) = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSal', '-', 'fSal (>1)')
          ! Register vertical velocity diagnostic, where relevant
          IF (data%phytos(i)%settling == _MOB_STOKES_ .OR. &
                                   data%phytos(i)%settling == _MOB_MOTILE_) THEN
            data%id_vvel(i) = aed_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_vvel', 'm/day', 'vertical velocity')
          ENDIF
       ENDIF
    ENDDO
END SUBROUTINE aed_phytoplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_phytoplankton(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the phytoplankton biogeochemical model
!
!  Here, the aed_p_m namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_phytoplankton_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status,i

!  %% NAMELIST
!  %% Last Checked 20/08/2021
   ! Default settings
   INTEGER            :: num_phytos = 0
   INTEGER            :: the_phytos(MAX_PHYTO_TYPES) = 0
   INTEGER            :: settling(MAX_PHYTO_TYPES)     = _MOB_CONST_
   AED_REAL           :: resuspension(MAX_PHYTO_TYPES) = 0.
   CHARACTER(len=64)  :: resus_link='NCS_resus'
   CHARACTER(len=64)  :: p_excretion_target_variable='OGM_dop'
   CHARACTER(len=64)  :: p_mortality_target_variable='OGM_pop'
   CHARACTER(len=64)  :: p1_uptake_target_variable='PHS_frp'
   CHARACTER(len=64)  :: p2_uptake_target_variable=''
   CHARACTER(len=64)  :: n_excretion_target_variable='OGM_don'
   CHARACTER(len=64)  :: n_mortality_target_variable='OGM_pon'
   CHARACTER(len=64)  :: n1_uptake_target_variable='NIT_nit'
   CHARACTER(len=64)  :: n2_uptake_target_variable='NIT_amm'
   CHARACTER(len=64)  :: n3_uptake_target_variable=''
   CHARACTER(len=64)  :: n4_uptake_target_variable=''
   CHARACTER(len=64)  :: c_excretion_target_variable='OGM_doc'
   CHARACTER(len=64)  :: c_mortality_target_variable='OGM_poc'
   CHARACTER(len=64)  :: c_uptake_target_variable=''
   CHARACTER(len=64)  :: do_uptake_target_variable='OXY_oxy'
   CHARACTER(len=64)  :: si_excretion_target_variable=''
   CHARACTER(len=64)  :: si_mortality_target_variable=''
   CHARACTER(len=64)  :: si_uptake_target_variable=''
   CHARACTER(len=128) :: dbase='aed_phyto_pars.nml'
   AED_REAL           :: zerolimitfudgefactor = 0.9 * 3600
   AED_REAL           :: R_mpbg = 0.
   AED_REAL           :: R_mpbr = 0.
   AED_REAL           :: I_Kmpb = 100.
   AED_REAL           :: mpb_max = 1000.
   AED_REAL           :: theta_mpb_growth = 1.05
   AED_REAL           :: theta_mpb_resp   = 1.05
   AED_REAL           :: min_rho = 900.
   AED_REAL           :: max_rho = 1200.
   INTEGER            :: do_mpb = 0
   INTEGER            :: n_zones = 0
   AED_REAL           :: active_zones(1000) = 0

!  From module globals
   LOGICAL  :: extra_debug = .false.      !## Obsolete Use diag_level = 10
!  LOGICAL  :: extra_diag = .false.      !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST

   NAMELIST /aed_phytoplankton/ num_phytos, the_phytos, settling, resuspension,&
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    si_excretion_target_variable,si_mortality_target_variable, &
                      si_uptake_target_variable,                               &
                    dbase, zerolimitfudgefactor, extra_debug, extra_diag,      &
                    do_mpb, R_mpbg, R_mpbr, I_Kmpb, mpb_max, min_rho, max_rho, &
                    resus_link, n_zones, active_zones, diag_level,             &
                    theta_mpb_growth,theta_mpb_resp
!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_phytoplankton initialization"

   ! Read the namelist, and set module parameters
   read(namlst,nml=aed_phytoplankton,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_phytoplankton'
   dtlim = zerolimitfudgefactor
   IF( extra_debug ) extra_diag = .true.       ! legacy use of extra_debug
   IF ( extra_diag ) diag_level = 10

   ! Set module parameters
   data%min_rho = min_rho ; data%max_rho = max_rho
   data%do_mpb = do_mpb
   data%R_mpbg = R_mpbg/secs_per_day   ; data%R_mpbr = R_mpbr/secs_per_day
   data%I_Kmpb = I_Kmpb   ; data%mpb_max = mpb_max
   data%theta_mpb_growth = theta_mpb_growth   ; data%theta_mpb_resp = theta_mpb_resp
   ALLOCATE(data%resuspension(num_phytos)); data%resuspension = resuspension(1:num_phytos)
   data%n_zones = n_zones
   IF( n_zones>0 ) THEN
     ALLOCATE(data%active_zones(n_zones))
     DO i=1,n_zones
       data%active_zones(i) = active_zones(i)
     ENDDO
   ENDIF

   ! Store species parameter values in our own derived type
   ! Note: all rates are provided in values per day,
   !     and are converted in here to values per second.
   CALL aed_phytoplankton_load_params(data,dbase,num_phytos,the_phytos,settling,resuspension)

   CALL aed_bio_temp_function( data%num_phytos,                                &
                               data%phytos%theta_growth,                       &
                               data%phytos%T_std,                              &
                               data%phytos%T_opt,                              &
                               data%phytos%T_max,                              &
                               data%phytos%aTn,                                &
                               data%phytos%bTn,                                &
                               data%phytos%kTn,                                &
                               data%phytos%p_name)

   ! Register microphytbenthos as a state variable
   IF (data%do_mpb>0) THEN
     data%id_mpb = aed_define_sheet_variable(  'mpb',                          &
                                               'mmol/m**2',                    &
                                               'microphytobenthos biomass',    &
                                                0.001,                         &
                                                minimum=0.001)
   ENDIF

   ! Register link to nutrient pools, if variable names are provided in namelist.
   data%do_Pexc = p_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Pexctarget = aed_locate_variable(p_excretion_target_variable)
   ENDIF
   data%do_Nexc = n_excretion_target_variable .NE. ''
   IF (data%do_Nexc) THEN
     data%id_Nexctarget = aed_locate_variable(n_excretion_target_variable)
   ENDIF
   data%do_Cexc = c_excretion_target_variable .NE. ''
   IF (data%do_Cexc) THEN
     data%id_Cexctarget = aed_locate_variable(c_excretion_target_variable)
   ENDIF
   data%do_Siexc = si_excretion_target_variable .NE. ''
   IF (data%do_Siexc) THEN
     data%id_Siexctarget = aed_locate_variable(si_excretion_target_variable)
   ENDIF

   data%do_Pmort = p_mortality_target_variable .NE. ''
   IF (data%do_Pmort) THEN
     data%id_Pmorttarget = aed_locate_variable(p_mortality_target_variable)
   ENDIF
   data%do_Nmort = n_mortality_target_variable .NE. ''
   IF (data%do_Nmort) THEN
     data%id_Nmorttarget = aed_locate_variable(n_mortality_target_variable)
   ENDIF
   data%do_Cmort = c_mortality_target_variable .NE. ''
   IF (data%do_Cmort) THEN
     data%id_Cmorttarget = aed_locate_variable(c_mortality_target_variable)
   ENDIF
   data%do_Simort = si_mortality_target_variable .NE. ''
   IF (data%do_Simort) THEN
     data%id_Simorttarget = aed_locate_variable(si_mortality_target_variable)
   ENDIF

   data%npup = 0
   IF (p1_uptake_target_variable .NE. '') data%npup = 1
   IF (p2_uptake_target_variable .NE. '') data%npup = 2
   data%do_Puptake = .FALSE.
   IF (data%npup>0) data%do_Puptake=.TRUE.
   IF (data%do_Puptake) THEN
   ! IF (data%npup>0) THEN ; data%id_Pupttarget(1) = aed_locate_variable(p1_uptake_target_variable); ifrp=1 ; ENDIF
   ! IF (data%npup>1) THEN ; data%id_Pupttarget(2) = aed_locate_variable(p2_uptake_target_variable); idop=2 ; ENDIF
     !# CAB ifrp and idop are now constants in bio_utils
     IF (data%npup>0) data%id_Pupttarget(ifrp) = aed_locate_variable(p1_uptake_target_variable)
     IF (data%npup>1) data%id_Pupttarget(idop) = aed_locate_variable(p2_uptake_target_variable)
   ENDIF
   data%nnup = 0
   IF (n1_uptake_target_variable .NE. '') data%nnup = 1
   IF (n2_uptake_target_variable .NE. '') data%nnup = 2
   IF (n3_uptake_target_variable .NE. '') data%nnup = 3
   IF (n4_uptake_target_variable .NE. '') data%nnup = 4
   data%do_Nuptake = .false.
   IF (data%nnup>0) data%do_Nuptake=.true.
   IF (data%do_Nuptake) THEN
   ! IF (data%nnup>0) THEN ; data%id_Nupttarget(1) = aed_locate_variable( n1_uptake_target_variable); ino3=1 ; ENDIF
   ! IF (data%nnup>1) THEN ; data%id_Nupttarget(2) = aed_locate_variable( n2_uptake_target_variable); inh4=2 ; ENDIF
   ! IF (data%nnup>2) THEN ; data%id_Nupttarget(3) = aed_locate_variable( n3_uptake_target_variable); idon=3 ; ENDIF
   ! IF (data%nnup>3) THEN ; data%id_Nupttarget(4) = aed_locate_variable( n4_uptake_target_variable); in2 =4 ; ENDIF
     !# CAB ino3, inh4, idon and in2 are now constants in bio_utils
     IF (data%nnup>0) data%id_Nupttarget(ino3) = aed_locate_variable( n1_uptake_target_variable)
     IF (data%nnup>1) data%id_Nupttarget(inh4) = aed_locate_variable( n2_uptake_target_variable)
     IF (data%nnup>2) data%id_Nupttarget(idon) = aed_locate_variable( n3_uptake_target_variable)
     IF (data%nnup>3) data%id_Nupttarget(in2)  = aed_locate_variable( n4_uptake_target_variable)
   ENDIF
   data%do_Cuptake = c_uptake_target_variable .NE. ''
   IF (data%do_Cuptake) THEN
     data%id_Cupttarget = aed_locate_variable( c_uptake_target_variable)
   ENDIF
   data%do_DOuptake = do_uptake_target_variable .NE. ''
   IF (data%do_DOuptake) THEN
     data%id_DOupttarget = aed_locate_variable( do_uptake_target_variable)
   ENDIF
   data%do_Siuptake = si_uptake_target_variable .NE. ''
   IF (data%do_Siuptake) THEN
     data%id_Siupttarget = aed_locate_variable( si_uptake_target_variable)
   ENDIF

   !-- sedimentation link variable
   data%id_Psed_phy = aed_define_diag_variable('Psed_phy','mmol/m**2/s','PHY sedimentation')

   !-- resuspension link variable
   IF ( .NOT.resus_link .EQ. '' ) THEN
      data%id_l_resus  = aed_locate_sheet_variable(TRIM(resus_link))
   ELSE
      data%id_l_resus = 0
      data%resuspension = 0.
   ENDIF

   ! Register diagnostic variables
   data%id_GPP = aed_define_diag_variable('GPP','mmol/m**3/d',  'gross primary production')
   data%id_NCP = aed_define_diag_variable('NCP','mmol/m**3/d',  'net community production')
   IF ( diag_level >= 10 ) data%id_PPR = aed_define_diag_variable('PPR','-','phytoplankton p/r ratio (gross)')
   IF ( diag_level >= 10 ) data%id_NPR = aed_define_diag_variable('NPR','-','phytoplankton p/r ratio (net)')

   data%id_NUP = aed_define_diag_variable('NUP_no3','mmol/m**3/d','nitrogen (NO3) uptake')
   data%id_NUP2= aed_define_diag_variable('NUP_nh4','mmol/m**3/d','nitrogen (NH4) uptake')
   data%id_PUP = aed_define_diag_variable('PUP','mmol/m**3/d','phosphorous uptake')
   data%id_CUP = aed_define_diag_variable('CUP','mmol/m**3/d','carbon uptake')

   IF ( diag_level >= 10 ) data%id_dPAR = aed_define_diag_variable('PAR','W/m**2', 'photosynthetically active radiation')
   data%id_TCHLA = aed_define_diag_variable('TCHLA','ug/L', 'total chlorophyll-a')
   IF ( diag_level >= 10 ) data%id_TPHY = aed_define_diag_variable('TPHYS','mmol/m**3', 'total phytoplankton')
   data%id_TIN = aed_define_diag_variable('IN','mmol/m**3', 'total phy nitrogen')
   data%id_TIP = aed_define_diag_variable('IP','mmol/m**3', 'total phy phosphorus')
   IF(do_mpb>0) THEN
     data%id_d_MPB = aed_define_sheet_diag_variable('MPB','mmol/m**2', 'microphytobenthos density')
     data%id_d_BPP = aed_define_sheet_diag_variable('BPP','mmol/m**2/d', 'benthic gross productivity')
     data%id_d_BCP = aed_define_sheet_diag_variable('BCP','mmol/m**2/d', 'benthic net productivity')
     data%id_d_mpbv= aed_define_sheet_diag_variable('MPBV','mmol/m**2/d', 'mpb vertical exchange')
   ENDIF

   ! Register environmental dependencies
   data%id_tem = aed_locate_global('temperature')
   data%id_sal = aed_locate_global('salinity')
   data%id_par = aed_locate_global('par')
   data%id_I_0 = aed_locate_sheet_global('par_sf')
   data%id_dz = aed_locate_global('layer_ht')
   data%id_extc = aed_locate_global('extc_coef')
   data%id_dens = aed_locate_global('density')
   data%id_sedzone = aed_locate_sheet_global('sed_zone')

END SUBROUTINE aed_define_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_phytoplankton(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Right hand sides of phytoplankton biogeochemical model
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
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
   AED_REAL :: upTot,net_cuptake

   INTEGER  :: phy_i,c
   AED_REAL :: flux, available

!------------------------------------------------------------------------------+
!BEGIN
   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)      ! local temperature
   salinity = _STATE_VAR_(data%id_sal)  ! local salinity
   par = _STATE_VAR_(data%id_par)       ! local photosynth. active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation

   ! Retrieve current (local) state variable values.
   pup = 0.
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

   !---------------------------------------------------------------------------+
   DO phy_i=1,data%num_phytos

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

      !------------------------------------------------------------------------+
      ! Get the temperature limitation function
      fT = fTemp_function(data%phytos(phy_i)%fT_Method,    &
                          data%phytos(phy_i)%T_max,        &
                          data%phytos(phy_i)%T_std,        &
                          data%phytos(phy_i)%theta_growth, &
                          data%phytos(phy_i)%aTn,          &
                          data%phytos(phy_i)%bTn,          &
                          data%phytos(phy_i)%kTn,temp)

      !------------------------------------------------------------------------+
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
      fI = photosynthesis_irradiance(data%phytos(phy_i)%lightModel, &
               data%phytos(phy_i)%I_K, data%phytos(phy_i)%I_S, par, extc, Io, dz)
      ! fI = 0.1


      ! METAL AND TOXIC EFFECTS
      fXl = 1.0

      !------------------------------------------------------------------------+
      ! Primary production rate
      primprod(phy_i) = data%phytos(phy_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl

      ! Adjust primary production rate for nitrogen fixers
      IF (data%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixing species, and the growth rate to  must be reduced
         ! to compensate for the increased metabolic cost of this process
         primprod(phy_i) = primprod(phy_i) * (data%phytos(phy_i)%k_nfix + &
                           (1.0-a_nfix(phy_i))*(1.0-data%phytos(phy_i)%k_nfix))
      ENDIF

      !------------------------------------------------------------------------+
      ! Respiration and general metabolic loss
      respiration(phy_i) = bio_respiration(data%phytos(phy_i)%R_resp,data%phytos(phy_i)%theta_resp,temp)

      ! Salinity stress effect on respiration (or growth)
      fSal =  phyto_salinity(data%phytos,phy_i,salinity)
      IF( data%phytos(phy_i)%salTol >= 4) THEN
        primprod(phy_i) = primprod(phy_i) * fSal   ! growth limtation rather than mortality enhancement
      ELSE
        respiration(phy_i) = respiration(phy_i) * fSal
      ENDIF

      ! photo-exudation
      exudation(phy_i) = primprod(phy_i)*data%phytos(phy_i)%f_pr

      ! Limit respiration if at the min biomass to prevent
      ! leak in the C mass balance
      IF (phy <= data%phytos(phy_i)%p0) THEN
         respiration(phy_i) = zero_
         exudation(phy_i) = zero_
      ENDIF

      ! write(*,"(4X,'limitations (fT,fI,fN,fP,fSi,Io, par, mu): ',9F9.2)")fT,fI,fNit,fPho,fSil,Io,par,primprod*secs_per_day

      !------------------------------------------------------------------------+
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

      !------------------------------------------------------------------------+
      ! Diagnostic info
      _DIAG_VAR_(data%id_NtoP(phy_i)) =  INi/IPi

      IF ( diag_level >= 10 ) THEN
         _DIAG_VAR_(data%id_fT(phy_i))   =  fT
         _DIAG_VAR_(data%id_fI(phy_i))   =  fI
         _DIAG_VAR_(data%id_fNit(phy_i)) =  fNit
         _DIAG_VAR_(data%id_fPho(phy_i)) =  fPho
         _DIAG_VAR_(data%id_fSil(phy_i)) =  fSil
         _DIAG_VAR_(data%id_fSal(phy_i)) =  fSal
      ENDIF
   END DO


   !---------------------------------------------------------------------------+
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

   !---------------------------------------------------------------------------+
   ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER
   net_cuptake = zero_
   DO phy_i=1,data%num_phytos

      !------------------------------------------------------------------------+
      !# PHYTOPLANKTON PRODUCTION & RESPIRATION
      phy = _STATE_VAR_(data%id_p(phy_i))
      flux = (primprod(phy_i) - respiration(phy_i) - exudation(phy_i)) * phy
      available = MAX(zero_, phy - data%phytos(phy_i)%p0)
      IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
      _FLUX_VAR_(data%id_p(phy_i)) = _FLUX_VAR_(data%id_p(phy_i)) + (flux)

      !------------------------------------------------------------------------+
      !# PHYTOPLANKTON INTERNAL NITROGEN
      IF (data%phytos(phy_i)%simINDynamics /= 0) THEN
         ! _FLUX_VAR_(data%id_in(phy_i)) = _FLUX_VAR_(data%id_in(phy_i)) +     &
         !      ( (-sum(nuptake) - nexcretion(phy_i) - nmortality(phy_i) )*INi )
         INi = _STATE_VAR_(data%id_in(phy_i))
         flux = (-sum(nuptake(phy_i,:)) - nexcretion(phy_i) - nmortality(phy_i))
         available = MAX(zero_, INi - data%phytos(phy_i)%X_nmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_in(phy_i)) = _FLUX_VAR_(data%id_in(phy_i)) + (flux)
      ENDIF

      !------------------------------------------------------------------------+
      !# PHYTOPLANKTON INTERNAL PHOSPHORUS
      IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! _FLUX_VAR_(data%id_ip(phy_i)) = _FLUX_VAR_(data%id_ip(phy_i)) +     &
         !  ( (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i)) )
         IPi = _STATE_VAR_(data%id_ip(phy_i))
         flux = (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i))
         available = MAX(zero_, IPi - data%phytos(phy_i)%X_pmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_ip(phy_i)) = _FLUX_VAR_(data%id_ip(phy_i)) + (flux)
      ENDIF

      !------------------------------------------------------------------------+
      !# PHYTOPLANKTON CELL DENSITY
      IF ( data%id_rho(phy_i)>0 ) THEN
         ! density increases during carbohydrate creation (daytime)
         flux = zero_
         IF( par>zero_ ) THEN
           flux = data%phytos(phy_i)%c1 * &
              (one_ - EXP(-par/data%phytos(phy_i)%I_K) ) - data%phytos(phy_i)%c3
         ELSE
           ! darkness
           flux = -data%phytos(phy_i)%c3
         ENDIF
        _FLUX_VAR_(data%id_rho(phy_i)) = _FLUX_VAR_(data%id_rho(phy_i)) + flux
         ! check maximum/minimum density are not exceeded
         IF( _STATE_VAR_(data%id_rho(phy_i))>data%max_rho ) THEN
            _FLUX_VAR_(data%id_rho(phy_i)) =zero_
            _STATE_VAR_(data%id_rho(phy_i))=data%max_rho
         ENDIF
         IF( _STATE_VAR_(data%id_rho(phy_i))<data%min_rho ) THEN
             _FLUX_VAR_(data%id_rho(phy_i)) =zero_
             _STATE_VAR_(data%id_rho(phy_i))=data%min_rho
         ENDIF
      ENDIF

      !------------------------------------------------------------------------+
      ! BIOGEOCHEMICAL FEEDBACKS
      ! Now manage uptake of nutrients, CO2 and DO - these cumulative fluxes
      ! are already limited above loop
      IF (data%do_Puptake) THEN
         DO c = 1,data%npup
            _FLUX_VAR_(data%id_Pupttarget(c)) =     &
                              _FLUX_VAR_(data%id_Pupttarget(c)) + ( puptake(phy_i,c))
         ENDDO
      ENDIF
      IF (data%do_Nuptake) THEN
         DO c = 1,data%nnup
            _FLUX_VAR_(data%id_Nupttarget(c)) =     &
                               _FLUX_VAR_(data%id_Nupttarget(c)) + ( nuptake(phy_i,c))
         ENDDO
      ENDIF
      IF (data%do_Cuptake) THEN
         _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) +           &
                       (  cuptake(phy_i) - respiration(phy_i)*data%phytos(phy_i)%k_fres*phy )
      ENDIF
      net_cuptake = net_cuptake + (cuptake(phy_i) - respiration(phy_i)*data%phytos(phy_i)%k_fres*phy)
      IF (data%do_DOuptake) THEN
         _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) +         &
                   ( -cuptake(phy_i) + respiration(phy_i)*data%phytos(phy_i)%k_fres*phy )
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

      !------------------------------------------------------------------------+
      ! UPDATE DIAGNOSTIC VARIABLES

      ! total phytoplankton carbon
      tphy = tphy + phy

      ! total chlorophyll-a
      IF (data%phytos(phy_i)%Xcc > 0.1) THEN
        !-- assume Xcc (mol C/ mol chla) is a constant
        tchla = tchla + ( phy / data%phytos(phy_i)%Xcc ) * 12.0
      ELSE
        !-- use dynamic equation (Eq 13: of Baklouti, Cloern et al. 1995)
        !-- theta = 1/Xcc [mg Chl (mg C)1] = 0.003 + 0.0154  e^0.050T  e^0.059E mu
        tchla = tchla + ( phy * (0.003 + 0.0154 * exp(0.050*temp) * exp(0.059*par) &
                        * primprod(phy_i)))
      ENDIF

      ! total internal nutrients
      tin = tin + INi
      tip = tip + IPi
   ENDDO

   !---------------------------------------------------------------------------+
   ! Set diagnostic arrays for combined assemblage properties
   _DIAG_VAR_(data%id_GPP) =  sum(cuptake)*secs_per_day
   _DIAG_VAR_(data%id_NCP) =  net_cuptake*secs_per_day
   IF ( diag_level >= 10 ) _DIAG_VAR_(data%id_PPR) =  -999. !sum(cuptake) / ( sum(cuptake) - net_cuptake)
   IF ( diag_level >= 10 ) _DIAG_VAR_(data%id_NPR) =  -999. !net_cuptake / ( sum(cuptake) - net_cuptake)
   _DIAG_VAR_(data%id_NUP) =  sum(nuptake(:,1))*secs_per_day
   _DIAG_VAR_(data%id_NUP2)=  sum(nuptake(:,2))*secs_per_day
   _DIAG_VAR_(data%id_PUP) =  sum(puptake)*secs_per_day
   _DIAG_VAR_(data%id_CUP) =  sum(cuptake)*secs_per_day

   IF ( diag_level >= 10 ) _DIAG_VAR_(data%id_dPAR) =  par
   _DIAG_VAR_(data%id_TCHLA)=  tchla
   IF ( diag_level >= 10 ) _DIAG_VAR_(data%id_TPHY) =  tphy
   _DIAG_VAR_(data%id_TIN)  =  tin
   _DIAG_VAR_(data%id_TIP)  =  tip

END SUBROUTINE aed_calculate_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_phytoplankton(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Calculate sedimentation of phytoplankton.
! Everything in units per surface area (not volume!) per time.
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: phy_i
   AED_REAL :: phy,mpb,temp,extc,par,dz,Io,fI,matz        ! State
   AED_REAL :: Fsed_phy,Psed_phy,mpb_flux,mpb_prod,mpb_resp
!
!------------------------------------------------------------------------------+
!BEGIN
   ! Loop through pelagic plankton groups (CURRENTLY DISABLED)
   !DO phy_i=1,data%num_phytos
   !    ! Retrieve current (local) state variable values.
   !    phy = _STATE_VAR_(data%id_p(phy_i))! phytoplankton
   !    phy_flux = zero_  ! no groups currently accumulate biomass in the sediment
   !    ! Set bottom fluxes for the pelagic (change per surface area per second)
   !    !_FLUX_VAR_(data%id_p(phy_i)) = _FLUX_VAR_(data%id_p(phy_i)) + phy_flux
   !    !_FLUX_VAR_B_(data%id_p(phy_i)) = _FLUX_VAR_B_(data%id_p(phy_i)) - phy_flux
   !ENDDO

   ! Process microphytobenthos (MPB)
   IF ( data%do_mpb>0 ) THEN
     ! Get local conditions
     matz = _STATE_VAR_S_(data%id_sedzone) ! local benthic type
     mpb  = _STATE_VAR_S_(data%id_mpb)     ! local mpb density
     temp = _STATE_VAR_(data%id_tem)       ! local temperature
     extc = _STATE_VAR_(data%id_extc)      ! cell extinction
     dz   = _STATE_VAR_(data%id_dz)        ! cell depth
     par  = _STATE_VAR_(data%id_par)       ! local photosynt. active radiation
     Io   = _STATE_VAR_S_(data%id_I_0)     ! surface short wave radiation

     ! Get sedimentation flux (mmmol/m2/s) loss into the benthos (diagnostic was set in mobility)
     Psed_phy = _DIAG_VAR_(data%id_Psed_phy)

     ! Compute photosynthesis and respiration
     fI = photosynthesis_irradiance(3,data%I_Kmpb,data%I_Kmpb,par,extc,Io,dz)
     mpb_prod = data%R_mpbg*fI*(data%theta_mpb_growth**(temp-20.))*(1.-(MIN(mpb,data%mpb_max)/data%mpb_max))
     mpb_resp = (data%R_mpbr*(data%theta_mpb_resp**(temp-20.)))           !*(( (mpb-mpb_min)/(data%mpb_max-mpb_min) )
     mpb_flux = (mpb_prod-mpb_resp)*mpb

     ! Update the MPB biomass, and O2/CO2 fluxes (mmol/m2/day)
     _FLUX_VAR_B_(data%id_mpb) = _FLUX_VAR_B_(data%id_mpb) + mpb_flux + Psed_phy
     IF (data%do_DOuptake) THEN
        _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) + mpb_flux
     ENDIF
     IF (data%do_Cuptake) THEN
        _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) - mpb_flux
     ENDIF
     ! A quick and dirty nutrient uptake by MPB; needs cleaning to account for limitation and excretion of DOM
     IF (data%do_Nuptake) THEN
        _FLUX_VAR_(data%id_Nupttarget(1)) = &
                           _FLUX_VAR_(data%id_Nupttarget(1)) - mpb_flux * (16./106.) *0.5
        _FLUX_VAR_(data%id_Nupttarget(2)) = &
                           _FLUX_VAR_(data%id_Nupttarget(2)) - mpb_flux * (16./106.) *0.5
     ENDIF
     IF (data%do_Puptake) THEN
        _FLUX_VAR_(data%id_Pupttarget(1)) = _FLUX_VAR_(data%id_Pupttarget(1)) - mpb_flux * (1./106.)
     ENDIF

     ! Resuspension (a simple assumption here)
     Fsed_phy = zero_
     IF ( data%n_zones > 0 ) THEN
        IF( in_zone_set(matz,data%active_zones) .AND. data%id_l_resus > 0 ) THEN
           Fsed_phy = _DIAG_VAR_S_(data%id_l_resus) * data%resuspension(1)
        ENDIF
     ENDIF
     _FLUX_VAR_B_(data%id_mpb) = _FLUX_VAR_B_(data%id_mpb) - Fsed_phy
     _FLUX_VAR_(data%id_p(1)) = _FLUX_VAR_(data%id_p(1)) + Fsed_phy

     ! Update the diagnostic variables
     _DIAG_VAR_S_(data%id_d_mpb) = mpb
     _DIAG_VAR_S_(data%id_d_bpp) =(mpb_prod) * mpb * secs_per_day
     _DIAG_VAR_S_(data%id_d_bcp) = mpb_flux * secs_per_day
     _DIAG_VAR_S_(data%id_d_mpbv)=(Psed_phy - Fsed_phy) * secs_per_day
   ENDIF
END SUBROUTINE aed_calculate_benthic_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_phytoplankton(data,column,layer_idx,mobility)
!------------------------------------------------------------------------------+
! Get the vertical movement values for phytoplankton cells
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: temp, par, rho_p, Io
   AED_REAL :: vvel
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: IN, IC, Q, Qmax
   INTEGER  :: phy_i
!
!------------------------------------------------------------------------------+
!BEGIN
   _DIAG_VAR_(data%id_Psed_phy) = zero_

   DO phy_i=1,data%num_phytos
      SELECT CASE (data%phytos(phy_i)%settling)

         CASE ( _MOB_OFF_ )
            ! disable settling by setting vertical velocity to 0
            vvel = zero_

         CASE ( _MOB_CONST_ )
            ! constant settling velocity using user provided value
            vvel = data%phytos(phy_i)%w_p

         CASE ( _MOB_TEMP_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_dens)
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%phytos(phy_i)%w_p*mu20*pw / ( mu*pw20 )

         CASE ( _MOB_STOKES_ )
            ! settling velocity based on Stokes Law calculation and cell density
            pw = _STATE_VAR_(data%id_dens)       ! water density
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)                 ! water dynamic viscosity
            IF( data%id_rho(phy_i)>0 ) THEN
              rho_p = _STATE_VAR_(data%id_rho(phy_i))  ! cell density
            ELSE
              rho_p = data%phytos(phy_i)%rho_phy
            ENDIF
            vvel = -9.807*(data%phytos(phy_i)%d_phy**2.)*( rho_p-pw ) / ( 18.*mu )

          CASE ( _MOB_MOTILE_ )
             ! vertical velocity based on motility and behaviour of phyto group
             ! modelled as in Ross and Sharples (2007)
             par = _STATE_VAR_(data%id_par)     ! local photosynthetically active radiation
             Io = _STATE_VAR_S_(data%id_I_0)    ! surface short wave radiation

             vvel = zero_
             IC = _STATE_VAR_(data%id_p(phy_i))
             IN = _STATE_VAR_(data%id_in(phy_i))
             Q = IN/IC
             Qmax = data%phytos(phy_i)%X_nmax
             vvel = zero_
             IF( Q<(data%phytos(phy_i)%f1*Qmax) ) THEN
               !-- nutrient starvation, swim down
               vvel = -data%phytos(phy_i)%w_p
             ELSEIF(Q>(data%phytos(phy_i)%f2*Qmax)) THEN
               !-- nutrient replete
               IF( par > data%phytos(phy_i)%I_K ) THEN
                 !-- high light, swim up
                 vvel = data%phytos(phy_i)%w_p
               ELSE
                 !-- low light
                 vvel = zero_
               ENDIF
            ELSE
               vvel = zero_
            ENDIF
            IF (par>0.95*Io) vvel = zero_   ! already at the surface

         CASE DEFAULT
            ! unknown settling/migration option selection
            vvel =  zero_

      END SELECT
      ! set global mobility array
      mobility(data%id_p(phy_i)) = vvel
      IF( diag_level >= 10  .AND. data%id_vvel(phy_i)>0) &
              _DIAG_VAR_(data%id_vvel(phy_i)) = vvel * secs_per_day
      ! set sedimentation flux (mmmol/m2) for later use/reporting
      _DIAG_VAR_(data%id_Psed_phy) =   &
              _DIAG_VAR_(data%id_Psed_phy) + vvel*_STATE_VAR_(data%id_p(phy_i))
    ENDDO
END SUBROUTINE aed_mobility_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_phytoplankton(data,column,layer_idx,extinction)
!------------------------------------------------------------------------------+
! Get the light extinction coefficient due to biogeochemical variables
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phytoplankton_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: phy
   INTEGER  :: phy_i
!
!------------------------------------------------------------------------------+
!BEGIN

   DO phy_i=1,data%num_phytos
      ! Retrieve current (local) state variable values.
      phy = _STATE_VAR_(data%id_p(phy_i))! phytoplankton

      ! Self-shading with contribution from this phytoplankton concentration.
      extinction = extinction + (data%phytos(phy_i)%KePHY*phy)
   ENDDO
END SUBROUTINE aed_light_extinction_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phytoplankton
