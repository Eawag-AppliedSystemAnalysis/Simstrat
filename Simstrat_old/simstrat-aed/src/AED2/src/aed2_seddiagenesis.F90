!###############################################################################
!#                                                                             #
!# aed2_seddiagenesis.F90                                                      #
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
!# Created May 2012                                                            #
!#                                                                             #
!###############################################################################


#include "aed2.h"


#define _MAX_LINKS_    700

#ifdef SINGLE
#  define SED_REAL REAL(SEDP)
#  define PO4AdsorptionFractionL PO4AdsorptionFractionLocal
#else
#  define PO4AdsorptionFractionL PO4AdsorptionFraction
#endif

!
MODULE aed2_seddiagenesis

   USE aed2_core
   USE aed2_util,        ONLY : PO4AdsorptionFraction
   USE aed2_gctypes
   USE aed2_gclib,       ONLY : AED_GC_Input,printTables
   USE aed2_sedeqsolver, ONLY : ConfigEquilibriumSolver,                  &
                                InitialiseGCProperties,                   &
                                UpdateEquilibration
   USE aed2_sedcandi,    ONLY : ConfigureCANDI, InitialiseCANDI, doCANDI, &
                                dfluxes, aed2_sed_candi_t

   IMPLICIT NONE

      PRIVATE
!
      PUBLIC aed2_seddiagenesis_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_seddiagenesis_data_t
      !# Variable identifiers
      INTEGER  :: id_zone, id_temp, id_salt, id_par
      INTEGER  :: id_Fsed_oxy, id_Fsed_rsi, id_Fsed_amm, id_Fsed_nit
      INTEGER  :: id_Fsed_frp, id_Fsed_pon, id_Fsed_don
      INTEGER  :: id_Fsed_pop, id_Fsed_dop, id_Fsed_poc, id_Fsed_doc
      INTEGER  :: id_Fsed_dic, id_Fsed_ch4, id_Fsed_feii

      INTEGER, ALLOCATABLE :: id_bw(:)
      INTEGER, ALLOCATABLE :: id_df(:),id_pf(:)

      INTEGER  :: id_zones
!     INTEGER  :: id_tot_sed

      !# Model parameters
      INTEGER  :: sed_modl, n_zones, nSedLayers, startSedCol, endSedCol

      TYPE(aed2_sed_candi_t),POINTER :: candi
      AED_REAL, DIMENSION(:,:,:), ALLOCATABLE :: SedLayerData

     CONTAINS
         PROCEDURE :: define            => aed2_define_seddiagenesis
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_seddiagenesis
!        PROCEDURE :: mobility          => aed2_mobility_seddiagenesis
!        PROCEDURE :: light_extinction  => aed2_light_extinction_seddiagenesis
!        PROCEDURE :: delete            => aed2_delete_seddiagenesis

   END TYPE

   CHARACTER(len=40) :: variables(_MAX_LINKS_),      &
                        water_link(_MAX_LINKS_),     &
                        diss_flux_link(_MAX_LINKS_), &
                        part_sed_link(_MAX_LINKS_)

   AED_REAL          :: default_vals(_MAX_LINKS_)
   AED_REAL          :: initial_vals(_MAX_LINKS_)

   AED_REAL,DIMENSION(:,:), ALLOCATABLE :: dynamic_swibc
   AED_REAL,DIMENSION(:), ALLOCATABLE :: bctime


   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: listSedDissEqVars,listSedPartEqVars
   INTEGER, ALLOCATABLE :: dColMap(:)
   INTEGER, ALLOCATABLE :: pColMap(:)
   INTEGER  :: nSedDissEqVars, nSedPartEqVars


   TYPE(AEDConstDiagenesisType), DIMENSION(:), ALLOCATABLE :: diagenesis

   INTEGER  :: fgcount, updateStep, nvars

   INTEGER, PARAMETER :: sedfid=101010

   AED_REAL :: time
  !INTEGER  :: substep, thisStep
   INTEGER  :: thisStep

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_seddiagenesis(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_seddiagenesis_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS

   AED_REAL          :: FsedA_initial=0.01
   AED_REAL          :: FsedN_initial=0.01
   CHARACTER(len=64) :: sediment_model=''
   INTEGER  :: i, status



!  AED_REAL :: Fsed_oxy, Fsed_rsi, Fsed_amm, Fsed_nit, Fsed_frp, &
!              Fsed_pon, Fsed_don, Fsed_pop, Fsed_dop, &
!              Fsed_poc, Fsed_doc, Fsed_dic

   NAMELIST /aed2_sediment/ sediment_model
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_sediment,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_seddiagenesis'

   data%sed_modl = -1

   IF ( sediment_model .EQ. "DYNAMIC" ) THEN
      data%sed_modl = 3

   ELSEIF ( sediment_model .EQ. "DYNAMIC2D" ) THEN
      data%id_zones = aed2_locate_global( 'sed_zone')
      data%sed_modl = 4

   ELSE
      print *,'Not supported'
      STOP

   ENDIF

   CALL load_sed_candi_data(data,namlst)

   ! Register state dependent variables
   ! NOTE the benthic=.true. argument, which specifies the variable is benthic.

   ! Bottom water conc variables set from all modules - we need these for bottom water BC
   ALLOCATE(data%id_bw(nvars))
   DO i =1,nvars
     IF(TRIM(water_link(i))/='') THEN
       data%id_bw(i) = aed2_locate_variable( TRIM(water_link(i)))
     ENDIF
   ENDDO
   ! Dissolved sed flux rates back to all modules - we need these to flux back to water
   ALLOCATE(data%id_df(nvars))
   DO i =1,nvars
     IF(TRIM(diss_flux_link(i))/='') THEN
       data%id_df(i) = aed2_locate_sheet_variable( TRIM(diss_flux_link(i)))
     ENDIF
   ENDDO
   ! Particulate flux rates from other modules - these are deposition to sediment
   ALLOCATE(data%id_pf(nvars))
   DO i =1,nvars
     IF(TRIM(part_sed_link(i))/='') THEN
       data%id_pf(i) = aed2_locate_sheet_variable( TRIM(part_sed_link(i)))
     ENDIF
   ENDDO

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global( 'temperature')
   data%id_salt = aed2_locate_global( 'salinity')
   data%id_par = aed2_locate_global(  'par')

   ! Now initialise it
   CALL initialise_dynamic_sediment(data,namlst,default_vals)


CONTAINS
!===============================================================================


   !############################################################################
   SUBROUTINE load_sed_candi_data(data,namlst)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      CLASS (aed2_seddiagenesis_data_t),INTENT(inout) :: data
      INTEGER,INTENT(in)                              :: namlst
   !
   !LOCALS
      LOGICAL,DIMENSION(:),ALLOCATABLE           :: dummySolid
      CHARACTER(LEN=64),DIMENSION(:),ALLOCATABLE :: dummyName
     !CHARACTER(len=64) :: sediment_model=''

      INTEGER  :: n_zones,i, nSedCols, parCounter, disCounter, var
      INTEGER  :: num_components, num_minerals
      CHARACTER(len=64) :: dis_components(MAX_GC_COMPONENTS)
      CHARACTER(len=64) :: the_minerals(MAX_GC_MINERALS)
      INTEGER, ALLOCATABLE :: dColMapTmp(:)
      INTEGER, ALLOCATABLE :: pColMapTmp(:)

      !-- Timey things --!
        INTEGER  :: timeswitch
        REAL     :: num_days
        REAL     :: fluxon
        REAL     :: fluxoff
        REAL     :: substep
        REAL     :: driverDT

      ! Stoichiometric coefficients
      ! OM model 1
        REAL     :: fracCPL                  ! Dan added
        REAL     :: fracCPR                  ! Dan added
        REAL     :: fracCPspecial            ! Dan added
        REAL     :: fracNPL                  ! Dan added
        REAL     :: fracNPR                  ! Dan added
        REAL     :: fracNPspecial            ! Dan added
        REAL     :: fracPPL                  ! Dan added
        REAL     :: fracPPR                  ! Dan added
        REAL     :: fracPPspecial            ! Dan added

        REAL     :: fracOAc                  ! Dan added
        REAL     :: fracH2                   ! Dan added

        REAL     :: fracCDHyd
        REAL     :: fracNDHyd
        REAL     :: fracPDHyd
        REAL     :: fracCOAc
        REAL     :: fracNOAc
        REAL     :: fracPOAc
        REAL     :: fracCH2
        REAL     :: fracNH2
        REAL     :: fracPH2

      !-- Initial & Boundary Conditions  --!
        INTEGER  :: ibc2
        INTEGER  :: ibbc
        INTEGER  :: startSteady
        REAL     :: flux_scale
        REAL     :: POMVR
        REAL     :: InitSedDepth
        CHARACTER(LEN=5) :: OMInitMethodL
        REAL     :: OM_topL
        REAL     :: OM_minL
        REAL     :: OM_cfL
        REAL     :: InitMinDepthL
        CHARACTER(LEN=5) :: OMInitMethodR
        REAL     :: OM_topR
        REAL     :: OM_minR
        REAL     :: OM_cfR
        REAL     :: InitMinDepthR
        REAL     :: InMinDep
        INTEGER  :: OutputUnits

        !-- other
        INTEGER  :: nZones, the_zones(10)
        CHARACTER(len=128) :: zone_pars

        INTEGER  :: Bsolidswitch

 NAMELIST /aed2_sed_candi/ timeswitch, num_days, fluxon, fluxoff, substep, driverDT, &

                           fracCPL, fracCPR, fracCPspecial,                    &
                           fracNPL, fracNPR, fracNPspecial,                    &
                           fracPPL, fracPPR, fracPPspecial,                    &
                           fracOAc, fracH2,                                    &
                           fracCDHyd, fracNDHyd, fracPDHyd,                    &
                           fracCOAc, fracNOAc, fracPOAc,                       &
                           fracCH2, fracNH2, fracPH2,                          &

                           ibc2, ibbc, startSteady, flux_scale, POMVR, InitSedDepth, &
                           OMInitMethodL, OM_topL, OM_minL, OM_cfL, InitMinDepthL,   &
                           OMInitMethodR, OM_topR, OM_minR, OM_cfR, InitMinDepthR,   &
                           InMinDep, OutputUnits,                                    &

                           ! time stuff
                           timeswitch, num_days, fluxon, fluxoff, substep,     &
                           driverDT,                                           &

                           ! var links etc
                           variables, default_vals, initial_vals,              &
                           water_link, diss_flux_link, part_sed_link,          &

                           ! zone stuff
                           nZones, the_zones, zone_pars, Bsolidswitch


   !----------------------------------------------------------------------------
   !BEGIN
      variables = ''

      ! Read the namelist
      read(namlst,nml=aed2_sed_candi,iostat=status)
      IF (status /= 0) STOP 'Error reading namelist aed2_sed_candi'

      IF (substep < 0) substep = 24   ! Assuming AED2 is running hourly

      data%n_zones = 1 ! n_zones

      ALLOCATE(diagenesis(data%n_zones))

      DO i=1,data%n_zones
      ! OM model 1
        diagenesis(i)%stcoef%fracCPL       = fracCPL !(i)                 ! Dan added
        diagenesis(i)%stcoef%fracCPR       = fracCPR !(i)                 ! Dan added
        diagenesis(i)%stcoef%fracCPspecial = fracCPspecial !(i)           ! Dan added
        diagenesis(i)%stcoef%fracNPL       = fracNPL !(i)                 ! Dan added
        diagenesis(i)%stcoef%fracNPR       = fracNPR !(i)                 ! Dan added
        diagenesis(i)%stcoef%fracNPspecial = fracNPspecial !(i)           ! Dan added
        diagenesis(i)%stcoef%fracPPL       = fracPPL !(i)                 ! Dan added
        diagenesis(i)%stcoef%fracPPR       = fracPPR !(i)                 ! Dan added
        diagenesis(i)%stcoef%fracPPspecial = fracPPspecial !(i)           ! Dan added
      ! OM model 2
      ! Nothing.
      ! OM model 3
        diagenesis(i)%stcoef%fracCDHyd    = fracCDHyd
        diagenesis(i)%stcoef%fracNDHyd    = fracNDHyd
        diagenesis(i)%stcoef%fracPDHyd    = fracPDHyd
        diagenesis(i)%stcoef%fracCOAc     = fracCOAc
        diagenesis(i)%stcoef%fracNOAc     = fracNOAc
        diagenesis(i)%stcoef%fracPOAc     = fracPOAc
        diagenesis(i)%stcoef%fracCH2      = fracCH2
        diagenesis(i)%stcoef%fracNH2      = fracNH2
        diagenesis(i)%stcoef%fracPH2      = fracPH2

        diagenesis(i)%BC%ibc2          = ibc2 !(i)
        diagenesis(i)%BC%ibbc          = ibbc !(i)
        diagenesis(i)%BC%startSteady   = startSteady !(i)
        diagenesis(i)%BC%flux_scale    = flux_scale !(i)
        diagenesis(i)%BC%POMVR         = POMVR !(i)
        diagenesis(i)%BC%InitSedDepth  = InitSedDepth !(i)
        diagenesis(i)%BC%OMInitMethodL = OMInitMethodL !(i)
        diagenesis(i)%BC%OM_topL       = OM_topL !(i)
        diagenesis(i)%BC%OM_minL       = OM_minL !(i)
        diagenesis(i)%BC%OM_cfL        = OM_cfL !(i)
        diagenesis(i)%BC%InitMinDepthL = InitMinDepthL !(i)
        diagenesis(i)%BC%OMInitMethodR = OMInitMethodR !(i)
        diagenesis(i)%BC%OM_topR       = OM_topR !(i)
        diagenesis(i)%BC%OM_minR       = OM_minR !(i)
        diagenesis(i)%BC%OM_cfR        = OM_cfR !(i)
        diagenesis(i)%BC%InitMinDepthR = InitMinDepthR !(i)
        diagenesis(i)%BC%InMinDep      = InMinDep !(i)
        diagenesis(i)%BC%OutputUnits   = OutputUnits !(i)

        diagenesis(i)%time%timeswitch    = timeswitch !(i)
        diagenesis(i)%time%num_days      = num_days !(i)
        diagenesis(i)%time%fluxon        = fluxon !(i)
        diagenesis(i)%time%fluxoff       = fluxoff  !(i)
        diagenesis(i)%time%substep       = substep  !(i)
        diagenesis(i)%time%driverDT      = driverDT  !(i)
      ENDDO

      data%nSedLayers = diagenesis(1)%param%maxnpts

      WRITE(*,"(6X,'Configuring CANDI...')")
      DO i =1,_MAX_LINKS_
       IF(TRIM(variables(i))/='') THEN
        nvars = i
       END IF
      END DO
      !nvars = 38

      ! Variables not requested by user that are compulsory, and added to the end
      variables(nvars+1) = 'pH'
      variables(nvars+2) = 'pe'
      variables(nvars+3) = 'ubalchg'

      WRITE(*,"(8X,'Num sediment biogeochemical variables =',I5)")nvars


      ALLOCATE( dColMapTmp(nvars) )  ! Column numbers of dissolved vars that are participating in eq solver
      ALLOCATE( pColMapTmp(nvars) )  ! Column numbers of pure phase vars that are participating in eq solver
      ALLOCATE( dummySolid(nvars+3)) ! T/F Flag for all vars for CANDI to know if its particulate
      dummySolid = .false.

      !-- Loop through and assign solid status, and geochemical status, if relevant
      i = 0;
      disCounter = 1
      parCounter = 1
      DO var = 1, nvars

        ! manganese iv 1
        IF(TRIM(variables(var)) == 'MnIV' .OR. TRIM(variables(var)) == 'mniv') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'MnIV'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! oxygen 2
        ELSEIF(TRIM(variables(var)) == 'oxy' .OR. TRIM(variables(var)) == 'O2') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'O2'

        ! sulfate 3
        ELSEIF(TRIM(variables(var)) == 'SO4' .OR. TRIM(variables(var)) == 'so4') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'SO4'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! phosphate 4
        ELSEIF(TRIM(variables(var)) == 'PO4' .OR. TRIM(variables(var)) == 'frp') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'PO4'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! ammonium 5
        ELSEIF(TRIM(variables(var)) == 'NH4' .OR. TRIM(variables(var)) == 'amm') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'NH4'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! methane 6
        ELSEIF(TRIM(variables(var)) == 'CH4' .OR. TRIM(variables(var)) == 'ch4') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'CH4'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! dic 7
        ELSEIF(TRIM(variables(var)) == 'DIC' .OR. TRIM(variables(var)) == 'dic') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'DIC'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! sulfide 8
        ELSEIF(TRIM(variables(var)) == 'H2S' .OR. TRIM(variables(var)) == 'h2s') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'H2S'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! adsorped PO4 9
        ELSEIF(TRIM(variables(var)) == 'PIP' .OR. TRIM(variables(var)) == 'frpads' .OR. TRIM(variables(var)) == 'pip') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! adsorbed N? (Dan) 10
        ELSEIF(TRIM(variables(var)) == 'PIN' .OR. TRIM(variables(var)) == 'nh4ads' .OR. TRIM(variables(var)) == 'pin') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! dissolved organic carbon 11
        ELSEIF(TRIM(variables(var)) == 'DOCl' .OR. TRIM(variables(var)) == 'docl') THEN
          dummySolid(var) = .false.;  i=i+1;

        ! dissolved organic carbon 12
        ELSEIF(TRIM(variables(var)) == 'DOCr' .OR. TRIM(variables(var)) == 'docr') THEN
          dummySolid(var) = .false.;  i=i+1;

        ! dissolved organic nitrogen 13
        ELSEIF(TRIM(variables(var)) == 'DONl' .OR. TRIM(variables(var)) == 'donl') THEN
          dummySolid(var) = .false.;  i=i+1;

        ! dissolved organic nitrogen 14
        ELSEIF(TRIM(variables(var)) == 'DONr' .OR. TRIM(variables(var)) == 'donr') THEN
          dummySolid(var) = .false.;  i=i+1;

        ! dissolved organic phosphorus 15
        ELSEIF(TRIM(variables(var)) == 'DOPl' .OR. TRIM(variables(var)) == 'dopl') THEN
          dummySolid(var) = .false.;  i=i+1;

        ! dissolved organic phosphorus 16
        ELSEIF(TRIM(variables(var)) == 'DOPr' .OR. TRIM(variables(var)) == 'dopr') THEN
          dummySolid(var) = .false.;  i=i+1;

        ! particulate organic carbon 17
        ELSEIF(TRIM(variables(var)) == 'POCl' .OR. TRIM(variables(var)) == 'pocl') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! particulate organic carbon 18
        ELSEIF(TRIM(variables(var)) == 'POCr' .OR. TRIM(variables(var)) == 'pocr') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! particulate organic nitrogen 19
        ELSEIF(TRIM(variables(var)) == 'PONl' .OR. TRIM(variables(var)) == 'ponl') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! particulate organic nitroegn 20
        ELSEIF(TRIM(variables(var)) == 'PONr' .OR. TRIM(variables(var)) == 'ponr') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! particulate organic phosphorus 21
        ELSEIF(TRIM(variables(var)) == 'POPl' .OR. TRIM(variables(var)) == 'popl') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! particulate organic phosphorus 22
        ELSEIF(TRIM(variables(var)) == 'POPr' .OR. TRIM(variables(var)) == 'popr') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! microphytobenthos 23
        ELSEIF(TRIM(variables(var)) == 'chla' .OR. TRIM(variables(var)) == 'Chla') THEN
          dummySolid(var) = .true.;  i=i+1;

        ! benthic grazers 24
        ELSEIF(TRIM(variables(var)) == 'zoo' .OR. TRIM(variables(var)) == 'zoop') THEN
          dummySolid(var) = .true.;  i=i+1;

          !-----------------------------------------------------------!
                  !Dan adds new variables!
          !POML 25
        ELSEIF(TRIM(variables(var)) == 'POML' .OR. TRIM(variables(var)) == 'poml') THEN
          dummySolid(var) = .true.;  i=i+1;
            ! POMR 26
        ELSEIF(TRIM(variables(var)) == 'POMR' .OR. TRIM(variables(var)) == 'pomr') THEN
          dummySolid(var) = .true.;  i=i+1;
          !DOMR 27
        ELSEIF(TRIM(variables(var)) == 'DOMR' .OR. TRIM(variables(var)) == 'domr') THEN
        !       dis_components(disCounter) = 'DOMR'
          dummySolid(var) = .false.;  i=i+1;
            ! dHyd 28
        ELSEIF(TRIM(variables(var)) == 'dhyd' .OR. TRIM(variables(var)) == 'DHyd') THEN
        !       dis_components(disCounter) = 'DHYD'
          dummySolid(var) = .false.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'POM1' .OR. TRIM(variables(var)) == 'pom1') THEN

          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'POM2' .OR. TRIM(variables(var)) == 'pom2') THEN

          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'POM3' .OR. TRIM(variables(var)) == 'pom3') THEN

          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'POM4' .OR. TRIM(variables(var)) == 'pom4') THEN

          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'OAc' .OR. TRIM(variables(var)) == 'oac') THEN

          dummySolid(var) = .false.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'H2' .OR. TRIM(variables(var)) == 'h2') THEN

          dummySolid(var) = .false.;  i=i+1;
        ELSEIF(TRIM(variables(var)) == 'Necromass' .OR. TRIM(variables(var)) == 'Necromass') THEN

          dummySolid(var) = .true.;  i=i+1;

        ELSEIF(TRIM(variables(var)) == 'BAer' .OR. TRIM(variables(var)) == 'baer') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BAer'
        END IF

        ELSEIF(TRIM(variables(var)) == 'BDen' .OR. TRIM(variables(var)) == 'bden') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BDen'
        END IF

        ELSEIF(TRIM(variables(var)) == 'BMan' .OR. TRIM(variables(var)) == 'bman') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BMan'
        END IF

        ELSEIF(TRIM(variables(var)) == 'BIro' .OR. TRIM(variables(var)) == 'biro') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BIro'
        END IF

        ELSEIF(TRIM(variables(var)) == 'BSul' .OR. TRIM(variables(var)) == 'bsul') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BSul'
        END IF

        ELSEIF(TRIM(variables(var)) == 'BMet' .OR. TRIM(variables(var)) == 'bmet') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BMet'
        END IF

        ELSEIF(TRIM(variables(var)) == 'BFer' .OR. TRIM(variables(var)) == 'bfer') THEN
        IF(Bsolidswitch==1)THEN
          dummySolid(var) = .true.;  i=i+1;
        ELSEIF(Bsolidswitch==2)THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'BFer'
        END IF
          !-----------------------------------------------------------!

        ! manganese
        ELSEIF(TRIM(variables(var)) == 'MnII' .OR. TRIM(variables(var)) == 'mnii') THEN

          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'MnII'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! nitrate
        ELSEIF(TRIM(variables(var)) == 'NO3' .OR. TRIM(variables(var)) == 'nit') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'NO3'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! manganese
        ELSEIF(TRIM(variables(var)) == 'MnO2A' .OR. TRIM(variables(var)) == 'mno2' .OR. TRIM(variables(var)) == 'mno2a') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'MnO2A'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! manganese
        ELSEIF(TRIM(variables(var)) == 'MnO2B' .OR. TRIM(variables(var)) == 'mno2b') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'MnO2B'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! mnco3
        ELSEIF(TRIM(variables(var)) == 'MnCO3' .OR. TRIM(variables(var)) == 'mnco3') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'MnCO3'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! iron (reduced)
        ELSEIF(TRIM(variables(var)) == 'FeII' .OR. TRIM(variables(var)) == 'feii') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'FeII'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! iron (oxidised)
        ELSEIF(TRIM(variables(var)) == 'FeIII' .OR. TRIM(variables(var)) == 'feiii') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'FeIII'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! particulate iron (amorphous)
        ELSEIF(TRIM(variables(var)) == 'FeOH3A' .OR. TRIM(variables(var)) == 'feoh3' .OR. TRIM(variables(var)) == 'feoh3a') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'FeOH3A'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

          ! particulate iron (crystalline)
        ELSEIF(TRIM(variables(var)) == 'FeOH3B' .OR. TRIM(variables(var)) == 'feoh3b') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'FeOH3B'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! iron monosulfide (MBO)
        ELSEIF(TRIM(variables(var)) == 'FeS' .OR. TRIM(variables(var)) == 'fes') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'FeS'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! pyrite
        ELSEIF(TRIM(variables(var)) == 'FeS2' .OR. TRIM(variables(var)) == 'pyrite' .OR. TRIM(variables(var)) == 'fes2') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'FeS2'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! feco3
        ELSEIF(TRIM(variables(var)) == 'FeCO3' .OR. TRIM(variables(var)) == 'feco3') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'FeCO3'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! ca
        ELSEIF(TRIM(variables(var)) == 'Ca' .OR. TRIM(variables(var)) == 'ca') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'Ca'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! caco3
        ELSEIF(TRIM(variables(var)) == 'CaCO3' .OR. TRIM(variables(var)) == 'calcite' .OR. TRIM(variables(var)) == 'caco3') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'CaCO3'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

        ! zn
        ELSEIF(TRIM(variables(var)) == 'Zn' .OR. TRIM(variables(var)) == 'zn') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'Zn'
          dColMapTmp(disCounter) = var ; disCounter=disCounter+1

        ! ZnS
        ELSEIF(TRIM(variables(var)) == 'ZnS' .OR. TRIM(variables(var)) == 'zns') THEN
          dummySolid(var) = .true.;  i=i+1;
          the_minerals(parCounter) = 'ZnS'
          pColMapTmp(parCounter) = var ; parCounter=parCounter+1

          !POML
        ELSEIF(TRIM(variables(var)) == 'POMspecial' .OR. TRIM(variables(var)) == 'pomspecial') THEN
          dummySolid(var) = .true.;  i=i+1;
                 ! N2
        ELSEIF(TRIM(variables(var)) == 'N2' .OR. TRIM(variables(var)) == 'n2') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'N2'

        ELSEIF(TRIM(variables(var)) == 'proton' .OR. TRIM(variables(var)) == 'proton') THEN
          dummySolid(var) = .false.;  i=i+1;
          dis_components(disCounter) = 'proton'


!        ! pH
!        ELSEIF(TRIM(variables(var)) == 'Ca' .OR. TRIM(variables(var)) == 'ca') THEN
!
!         dummySolid(var) = .false.;  i=i+1;

        ELSE

          print *,'CANDI-AED Setup - NOTE unknown variable: ',TRIM(variables(var))
          print *,' -> assuming it is dissolved and non-reactive'
          dummySolid(var) = .false.

        END IF

      END DO

      !----------------------------------------------------------------------------
      !-- Configure the main CANDI-AED model configuration
      data%candi => ConfigureCANDI(diagenesis(1),&
                          data%nSedLayers,&
                          nvars+3,&
                          variables(1:nvars+3),&
                          dummySolid(1:nvars+3),&
                          nSedCols)
      data%startSedCol = 0
      data%endSedCol   = nSedCols-1

      !-- Allocate space component objects
      ALLOCATE(data%SedLayerData(data%n_zones,data%nSedLayers,data%startSedCol:data%endSedCol))
      data%SedLayerData = 0.0

      time = 0.0
      thisStep = 0


     !----------------------------------------------------------------------------
     ! Configure the geochemical solver for sediment vars
     CALL AED_GC_Input('aed2_geochem_pars.dat')

     num_components = disCounter -1
     num_minerals   = parCounter -1
     IF(num_components<0)num_components=0
     IF(num_minerals<0)num_minerals=0

     !num_minerals = 6  ! Currently minerals in the Equil Solver are set to be ignored...

     CALL ConfigEquilibriumSolver( num_components,  num_minerals,                &
                                   dis_components(1:num_components),             &
                                   the_minerals(1:num_minerals),                 &
                                   nSedDissEqVars, nSedPartEqVars,               &
                                   listSedDissEqVars,listSedPartEqVars)

     !data%num_comp = nSedDissEqVars
     !data%num_mins = nSedPartEqVars

     WRITE(*,"(8X,'Num vars for Eq solver =',I5,I5)")nSedDissEqVars,nSedPartEqVars

     ALLOCATE( dColMap(nSedDissEqVars) )
     ALLOCATE( pColMap(nSedPartEqVars) )

     dColMap(1:num_components) = dColMapTmp(1:num_components) -1
     pColMap(1:num_minerals)   = pColMapTmp(1:num_minerals)   -1

     dColMap(num_components+1) = nvars+1 -1
     dColMap(num_components+2) = nvars+2 -1
     dColMap(num_components+3) = nvars+3 -1

     print *,'dcolmap',dColMap
     print *,'pcolmap',pColMap

     !----------------------------------------------------------------------------
     ! Check BCs

     ! If we are running with a dynamically varying external sed-water interface BC
     IF(IBC2 == 10) THEN
       CALL AED_SWI_Input('aed2_sediment_swibc.dat')
     ENDIF

      WRITE(*,"(6X,'CANDI configured OK')")

   END SUBROUTINE load_sed_candi_data
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE aed2_define_seddiagenesis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE initialise_dynamic_sediment(data,namlst,default_vals)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_seddiagenesis_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in)                             :: namlst
   AED_REAL,DIMENSION(:),INTENT(inout)            :: default_vals
!
!LOCALS
   AED_REAL,DIMENSION(:),   ALLOCATABLE   :: tt
   REAL(SEDP),DIMENSION(:), ALLOCATABLE   :: bottomConcs
   REAL(SEDP),DIMENSION(:), ALLOCATABLE   :: sedimentConcs
!
   INTEGER  :: allocStatus, openStatus, i, zone, kk, var
   AED_REAL,DIMENSION(nSedDissEqVars)     :: dissConcs
   AED_REAL,DIMENSION(nSedPartEqVars)     :: partConcs

   CHARACTER(LEN=(16+12+4)) :: FileName
!
!-------------------------------------------------------------------------------
!BEGIN
    WRITE(*,"(4X,'Initialising sediment diagenesis model (CANDI-AED)...')")

    !-- Initalisations and temporary allocations
    fgcount = 0

    ALLOCATE(bottomConcs(data%startSedCol:data%endSedCol))
    ALLOCATE(sedimentConcs(data%startSedCol:data%endSedCol))
    ALLOCATE(tt(data%nSedLayers))
    tt=15.0

    DO zone = 1,data%n_zones         ! DO Horse

       WRITE(*,"(6X,'Call initialiseCANDI @ zone: ',I5,'/',I5)")zone,data%n_zones

       !-- Prepare bottom conc information for CANDI initialisation
       DO var = data%startSedCol,data%endSedCol         ! DO Donkey

           bottomConcs(var)   = default_vals(var+1) /1e3  ! mmol/m3 -> mmol/L
           IF(data%candi%isSolid(var)) THEN
              ! Assumign solids are fluxing, init val all the way to the top
              bottomConcs(var)   = initial_vals(var+1) /1e3  ! mmol/m3 -> mmol/L
           END IF
           sedimentConcs(var) = initial_vals(var+1) /1e3  ! mmol/m3 -> mmol/L

       ENDDO         ! ENDDO Donkey
       DO var = data%endSedCol-2,data%endSedCol     ! DO Rabbit
           bottomConcs(var)   = default_vals(var+1) ! pH,pe,uchg
           sedimentConcs(var) = initial_vals(var+1)
print *, 'sedimentConcs', sedimentConcs
       bottomConcs((nvars+1)-1)     =   8.0             !pH
       sedimentConcs((nvars+1)-1)  =    8.0
       bottomConcs((nvars+2)-1)     =   8.0             !pe
       sedimentConcs((nvars+2)-1)  =    8.0
       bottomConcs((nvars+3)-1)     =   0 ! -0.056      !ubalchg
       sedimentConcs((nvars+3)-1)  =    0 ! -0.056
       ENDDO            ! ENDDO Rabbit
       ! For the fish project, the field data shows a consistent water pH of 8.0
       !-- Initialise diagenesis model
        !print *, 'bottomConcs', bottomConcs
       CALL initialiseCANDI(data%candi,bottomConcs,sedimentConcs)
       !This function lives in aed2_sedcandi
       !-- IF INITIALIZING WE MUST PERFORM AN INITIAL EQUILIBRATION
       DO kk = 1,data%nSedLayers      !DO Unicorns
         dissConcs = data%candi%Y(dColMap,kk)
         partConcs = data%candi%Y(pColMap,kk)
        !print *, 'dissConcs', dissConcs
        !print *,'bottomConcs,',bottomConcs
        !print *,'partConcs,',partConcs
         CALL InitialiseGCProperties( dissConcs, partConcs, 0 )
         data%candi%Y(dColMap,kk) = dissConcs
         data%candi%Y(pColMap,kk) = partConcs
         ENDDO                          !ENDDO Unicorns
!       !-- Store initialised "Y" array into SedimentData for this zone
!       DO kk = 1,data%nSedLayers
!          data%SedLayerData(zone,kk,:) =  Y(:,kk)
       ENDDO            ! END DO Donkey
       !-- Open some output text files to put all this info
       DO var = data%startSedCol,data%endSedCol
         FileName = "results/candi_aed/"//TRIM(variables(var+1))//".sed"   ! Needs zone dir
         OPEN(UNIT = sedfid+var, FILE = TRIM(FileName), &
             STATUS = "REPLACE", ACTION = "WRITE", IOSTAT = openStatus)
         IF (openStatus /= 0) THEN
           PRINT *,' initialise_dynamic_sediment(): Error on attempt to open sed file:',FileName
           PRINT *,'   Make sure the directory ./results/candi_aed/ exists'
           STOP    " PROGRAM STOPPED"
         END IF
         !-- Now print file header containing info and layer details
         WRITE(sedfid+var, &
           "('%! AED SedDiagenesis Output: ',A40 )")FileName
         WRITE(sedfid+var, &
            "('%!   1st column is time (years), remaining ',I3,' columns correspond to sed layers')")data%nSedLayers
         WRITE(sedfid+var, &
           "('%!   1st row is layer depths (cm)')")
         WRITE(sedfid+var,"(E13.5, 100E11.4)")0.0,data%candi%rpar(:)
       ENDDO
       DO var = 1,SIZE(partConcs)
         FileName = "results/candi_aed/IAP_"//TRIM(listSedPartEqVars(var))//".sed"   ! Needs zone dir
         OPEN(UNIT = sedfid+(nvars+3)+var+1, FILE = TRIM(FileName), &
             STATUS = "REPLACE", ACTION = "WRITE", IOSTAT = openStatus)
         IF (openStatus /= 0) THEN
           PRINT *,' initialise_dynamic_sediment(): Error on attempt to open sed file:',FileName
           PRINT *,'   Make sure the directory ./results/candi_aed/ exists'
           STOP    " PROGRAM STOPPED"
         END IF
         !-- Now print file header containing info and layer details
         WRITE(sedfid+(nvars+3)+var+1, &
           "('%! AED SedDiagenesis IAP Output: ',A40 )")FileName
         WRITE(sedfid+(nvars+3)+var+1, &
            "('%!   1st column is time (years), remaining ',I3,' columns correspond to sed layers')")data%nSedLayers
         WRITE(sedfid+(nvars+3)+var+1, &
           "('%!   1st row is layer depths (cm)')")
         WRITE(sedfid+(nvars+3)+var+1,"(E13.5, 100E11.4)")0.0,data%candi%rpar(:)
       ENDDO
       FileName = "results/candi_aed/swi_fluxes.sed"   ! need a file for dissolved fluxes
       OPEN(UNIT = sedfid+(nvars+3)+SIZE(partConcs)+10, FILE = TRIM(FileName), &
           STATUS = "REPLACE", ACTION = "WRITE", IOSTAT = openStatus)
       IF (openStatus /= 0) THEN
         PRINT *,' initialise_dynamic_sediment(): Error on attempt to open sed file:',FileName
         PRINT *,'   Make sure the directory ./results/candi_aed/ exists'
         STOP    " PROGRAM STOPPED"
       END IF
       !-- Now print file header containing info and layer details
       WRITE(sedfid+(nvars+3)+SIZE(partConcs)+10, &
         "('%! AED SedDiagenesis Sed-Water Interface FLUX Rates: ',A40 )")FileName
       WRITE(sedfid+(nvars+3)+SIZE(partConcs)+10, &
          "('%!   1st column is time (years), remaining ',I3,' columns correspond to variables')")nvars
!    ENDDO
    DEALLOCATE(sedimentConcs)
    DEALLOCATE(bottomConcs)
    DEALLOCATE(tt)
    WRITE(*,"(6X,'CANDI initialised OK')")
!    RETURN
!99 CALL fatal_error('aed2_sediment_init','Error reading namelist initialise_dynamic_sediment')
END SUBROUTINE initialise_dynamic_sediment
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed2_calculate_benthic_seddiagenesis(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED sediment.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_seddiagenesis_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
 ! AED_REAL :: steady, conc, flux, factor, default_bottom_ubalchg, default_bottom_pH
   DOUBLETYPE :: steady
   AED_REAL :: conc, flux, factor, default_bottom_ubalchg, default_bottom_pH
   INTEGER  :: zone,layer, var, rxn_mode
   AED_REAL,DIMENSION(nSedDissEqVars)     :: dissConcs
   AED_REAL,DIMENSION(nSedPartEqVars)     :: partConcs, IAPtmp, KIAPtmp, QIAPtmp
   AED_REAL,DIMENSION(SIZE(default_vals)) :: dynamic_vals
   AED_REAL,DIMENSION(SIZE(default_vals)) :: swi_flux
   AED_REAL :: PO4dis,PO4par,pH

   AED_REAL :: Kpo4p, Kadsratio, Qmax   ! CAB Added

   TYPE(aed2_sed_candi_t),POINTER :: candi
!
!-------------------------------------------------------------------------------
!BEGIN

   candi => data%candi

   !
   thisStep = thisStep+1
   IF(thisStep<candi%time%substep) THEN
     RETURN
   END IF

   time = time + candi%deltaT
   !print *, 'time', time
   DO zone = 1,data%n_zones
       WRITE(*,"(8X,'Updating sediment profiles @ zone: ',I5,'/',I5,' :',F8.5,' years')")zone,data%n_zones,time


       !-- Retrieve current environmental conditions for the cell.
       candi%TEMP = _STATE_VAR_(data%id_temp) ! local temperature
       candi%SALT = _STATE_VAR_(data%id_salt) ! local salinity
       candi%PRES = 1.025E-3 * 1.0




      IF(candi%BC%IBC2==10) THEN
         dynamic_vals = GetDynamicBCVals(time)
      END IF


       !-- Prepare bottom water concentration for boundary specification
       DO var = data%startSedCol,data%endSedCol-3

         ! Only set dissoved species concs here
         IF(.NOT. candi%isSolid(var)) THEN

           IF(TRIM(water_link(var+1))/='') THEN
              conc = _STATE_VAR_(data%id_bw(var+1))
             candi%Y(var,1) =  conc /1e3

           ELSE
             IF(candi%BC%IBC2==10) THEN
               ! Dynamic BW conc assigned as read in
               candi%Y(var,1) =  dynamic_vals(var+1) /1e3
             ELSE
               ! Constant BW conc at default value assumed
                           candi%Y(var,1) =  default_vals(var+1) /1e3

             END IF

           ENDIF
         END IF

       ENDDO


       default_bottom_pH = 8.00
       dissConcs = candi%Y(dColMap,1)
       partConcs = candi%Y(pColMap,1)
       dissConcs(SIZE(dissConcs)-2)=default_bottom_pH

       CALL InitialiseGCProperties( dissConcs, partConcs, 0 )

       default_bottom_ubalchg = dissConcs(SIZE(dissConcs))
!print *, 'dissConcs  ===========================', dissConcs
       candi%Y((nvars+1)-1,1) = default_bottom_pH
       candi%Y((nvars+3)-1,1) = default_bottom_ubalchg

       !-- Prepare particulate flux amounts for boundary specification

       IF(candi%BC%IBC2 == 0) THEN
         !-- Get sediment surface particulate flux from deposition amounts
         DO var = data%startSedCol,data%endSedCol-3
            ! umol/cm2/yr (from g/ts)
            candi%PartFluxes(var) = 0.0

            ! Only flux dissoved species concs here
            IF(candi%isSolid(var)) THEN
              ! Search for linked var (most likely from aed2_sedflux)
              IF(TRIM(part_sed_link(var+1))/='') THEN

                !-- Convert from mmol/m2/day -> umol/cm2/yr ????????

                factor = 1e-3 * 365.25 / (100.*100.)  ! may need seconds here if AED2 AED is in "s"

                 flux = _STATE_VAR_S_(data%id_pf(var+1))
                candi%PartFluxes(var) = flux/factor

              ELSE
                ! No linked var so just use default_var
                candi%PartFluxes(var) =  default_vals(var+1) !/1e3

              ENDIF
            END IF
         END DO

       ELSE IF(candi%BC%IBC2 == 1) THEN
         !-- Fix surface organic particulates at initial value
         DO var = data%startSedCol,data%endSedCol-3
           ! Only set particulate species CONCS here
           IF(candi%isSolid(var)) THEN
             candi%Y(var,1) =  default_vals(var+1) /1e3
                                         ENDIF
         ENDDO

       ELSE IF(candi%BC%IBC2 == 2) THEN
         !-- Set prescribed flux rate here or set from default_vals
         DO var = data%startSedCol,data%endSedCol  !DO King Kong
           ! Only set particulate species FLUX here
           IF(candi%isSolid(var)) THEN
                        ! Optional place to hard-code fluxes here
                        ! The point of this section is to have a simple code for making the flux
                        ! off, then on, then off again.
                        ! I did it to simulate Abrolhos Islands, add a fish farm, remove the fish farm
                        ! If you don't like this bit just change timeswitch to any number other than 1 in fabm.nml
                        IF (candi%time%timeswitch == 1) THEN   !
                                        IF (variables(var+1)=='pomspecial') THEN ! If it is 'pomspecial'
                                        !IF(var==45) THEN ! If it is 'pomspecial'
                                                        IF(time<candi%time%fluxon .OR. time>candi%time%fluxoff) THEN
                                                          candi%PartFluxes(var)=0.00E+00
                                                        ELSE
                                                          candi%PartFluxes(var)= default_vals(var+1) * 1e3 / 1e4 ! (CANDI expects umol/cm2/yr)
                                                        !!ELSEIF(time>time%fluxon .AND. time>time%fluxoff) THEN
                                                        !  !PartFluxes(var)=0.00E+00
                                                        !ELSE
                                                        !PartFluxes(var)=0.00E+00
                                                        ENDIF
!                                                       print*, 'Flux', PartFluxes(var)
                                        ELSE !Else if not pomspecial
                                          candi%PartFluxes(var) =  default_vals(var+1) * 1e3 / 1e4  ! (CANDI expects umol/cm2/yr)

                                          !PRINT *,'PartFlux_doCANDI',var,default_vals(var+1)

                                        ENDIF  !End whether it is 'pomspecial'
                        ENDIF !End if timeswitch = 1
           ENDIF !ENDIF is Solid
           ENDDO             !ENDDO King Kong

                ! PartFluxes(var) = X
       ELSE IF(candi%BC%IBC2 == 10) THEN
         !-- Set prescribed flux rate here or set from default_vals
         DO var = data%startSedCol,data%endSedCol
           ! Only set particulate species FLUX here
           IF(candi%isSolid(var)) THEN
             candi%PartFluxes(var) =  dynamic_vals(var+1) * 1e3 / 1e4  ! (CANDI expects umol/cm2/yr)
           ENDIF
         ENDDO
         ! Optional place to hard-code fluxes here
         ! PartFluxes(var) = X
       END IF


       !-- Calculate advection velocities for solids (wvel) and liquids (uvel)
       DO layer = 1,data%nSedLayers
         candi%wvel(layer)   = candi%param%w00*(1.0 - candi%param%p00)/(candi%ps(layer))
         candi%uvel(layer)   = candi%param%w00*candi%param%p00/candi%poros(layer)
         candi%kg0var(layer) = 0.01*exp(-candi%rpar(layer)**2.0/(2.0*candi%param%x2))  &   !!&&???
                       + 0.16*(candi%param%x2/candi%param%w00+candi%rpar(layer)/candi%wvel(layer))**(-0.95)

       END DO


       !-- Set uniform starting value if running in steady mode
       IF(candi%iSTEADY == 1) THEN
         DO layer = 2,data%nSedLayers
             candi%Y(candi%hco3y,layer) = candi%Y(candi%hco3y,1)
             candi%Y(candi%hsy,layer)   = 0.0
         END DO
       END IF

!
      !-- Perform thermodynamic equilibration prior to diagenesis equations
       DO layer = 2,data%nSedLayers

         dissConcs = candi%Y(dColMap,layer)
         WHERE(dissConcs(1:SIZE(dissConcs)-1)<0.0)dissConcs(1:SIZE(dissConcs)-1)=0.0
         dissConcs(SIZE(dissConcs)-1) = 8.0
         partConcs = candi%Y(pColMap,layer)
         WHERE(partConcs<0.0)partConcs=0.0

         IAPtmp = 0.0
         KIAPtmp = 0.0
         QIAPtmp = 0.0

         ! Do geochemical equilibration
         IF (rxn_mode==0)THEN
                 !'UpdateEquilibration' not called
         ELSEIF (rxn_mode==1)THEN
           CALL UpdateEquilibration(dissConcs, partConcs, concMode=0, inTemp=REAL(candi%TEMP), stoEq=.true., IAP=IAPtmp)
         ELSEIF (rxn_mode==2)THEN
           CALL UpdateEquilibration(dissConcs, partConcs, concMode=0, inTemp=REAL(candi%TEMP), stoEq=.false., IAP=IAPtmp)
         ELSEIF (rxn_mode==3)THEN
           CALL UpdateEquilibration(dissConcs, partConcs, concMode=0, inTemp=REAL(candi%TEMP), stoEq=.false.,IAP=IAPtmp)
         ENDIF ! End if rxn_mode


         !CALL UpdateEquilibration(dissConcs, partConcs, concMode=0, inTemp=REAL(candi%TEMP), stoEq=.true.)
         !CALL UpdateEquilibration(dissConcs, partConcs, concMode=0, inTemp=REAL(candi%TEMP), stoEq=.false.,&
         !IAP=IAPtmp, KIAP=KIAPtmp, QIAP=QIAPtmp)

         candi%Y(dColMap,layer) = dissConcs
         candi%Y(pColMap,layer) = partConcs

         candi%IAP(pColMap,layer) = IAPtmp
         candi%KIAP(pColMap,layer) = KIAPtmp
         candi%QIAP(pColMap,layer) = QIAPtmp


         IF(diagenesis(zone)%param%PO4AdsorptionModel==1 .OR. diagenesis(zone)%param%PO4AdsorptionModel==2) THEN

           Kpo4p = diagenesis(zone)%param%Kpo4p ; Kadsratio = diagenesis(zone)%param%Kadsratio ; Qmax = diagenesis(zone)%param%Qmax
          ! PO4 adsorption based on pH and Fe particulate concentration
           IF(diagenesis(zone)%param%ads_use_pH) THEN
             pH = candi%Y((nvars+1)-1,layer)     ! Watch out! pH is first column after nvars
print *, 'pH', pH
             CALL PO4AdsorptionFractionL(diagenesis(zone)%param%PO4AdsorptionModel, &
                                      candi%Y(candi%po4ly,layer)+candi%Y(candi%po4sy,layer),       &  ! Check my units please :)
                                      candi%Y(candi%feohy,layer)+candi%Y(candi%feohby,layer),      &  ! Check my units please :)
                                      Kpo4p,Kadsratio,Qmax, &
                                      PO4dis,PO4par,                       &
                                      thepH=pH)

           ELSE
             CALL PO4AdsorptionFractionL(diagenesis(zone)%param%PO4AdsorptionModel, &
                                        candi%Y(candi%po4ly,layer)+candi%Y(candi%po4sy,layer),       &  ! Check my units please :)
                                      candi%Y(candi%feohy,layer)+candi%Y(candi%feohby,layer),      &  ! Check my units please :)
                                      Kpo4p,Kadsratio,Qmax, &
                                      PO4dis,PO4par)
           ENDIF ! End if 'ads_use_pH' is on

           candi%Y(candi%po4ly,layer) = PO4dis
           candi%Y(candi%po4sy,layer) = PO4par

         END IF ! End if adsorption model is on

    END DO

    !   DO layer = 1,data%nSedLayers
    !  !   print *,'Y',Y(:,layer)
    !   END DO
       !-- Run main diagenesis model (kinetic components and transport)
       CALL doCANDI(candi,steady)

       ! Check columns (except last 3) for -ve
       WHERE(candi%Y(1:nvars,:)<1e-20)candi%Y(1:nvars,:)=0.00

    !   !-- Update main sediment data store with new CANDI data (Y)
    !   DO layer = 1,data%nSedLayers
    !    ! data%SedLayerData(zone,layer,:) = Y(:,layer)
    !   END DO

       !-- Write species concs to results/candi/*.sed files */
       DO var = data%startSedCol,data%endSedCol
          WRITE(sedfid+var,"(E13.5, 100E14.5)")time, candi%Y(var,:)
          if(var==data%endSedCol) THEN
            !print *,'sss',Y(var,:)
            !print *,'sss',Y(var-1,:)
            !pause
            END IF
       END DO
       !-- Write particualte IAP to results/candi/*.sed files */
       DO var = 1,SIZE(IAPtmp)
          WRITE(sedfid+(nvars+3)+var+1,"(E13.5, 100E14.5)")time, candi%IAP(pColMap(var),:)
       END DO

       !-----------------------------------------------------------------------
       !-- Set the flux of the dissolved species between water & sediment
       DO var = data%startSedCol,data%endSedCol-3

         ! Only flux dissoved species concs here
         IF(.NOT. candi%isSolid(var)) THEN
           ! SEarch for linked var (most likely from aed2_sedflux)
           IF(TRIM(diss_flux_link(var+1))/='') THEN

             !-- Convert from 1000 mmol/cm2/year -> g/m2/day ????????
             !-- (the extra factor of 1000 comes in when calculating flux from
             !-- mol/L concs; 1L=1000cm3)
             ! Dan:
             !       dfluxes is in 1e6 umol/cm2/y
             !       I want mmol/m2/y .: factor should = 1e-3 (1e6 umol->mmol) * 100*100 (cm2->m2)
             factor = 1e-6 *100.*100./365.25
             !factor = 1e-3 *100.*100.     ! Now factor = 1e1; 1e3/1e4 = 1e-1; you could divide PartFluxes by factor
             flux = -dfluxes(candi,2,var) * factor

             _STATE_VAR_S_(data%id_df(var+1)) =  flux
           ENDIF
           swi_flux(var) = dfluxes(candi,2,var)*(1e-3*100.*100.)  ! This unit conversion by Dan, December 2014
         ELSE
           swi_flux(var) = candi%PartFluxes(var)/(1e3/1e4)
         END IF
       END DO
       ! Output file designed to output in and out fluxes at the interface
       WRITE(sedfid+(nvars+3)+SIZE(IAPtmp)+10,"(E13.5, 100E14.5)")time, swi_flux(1:nvars)

   END DO
   thisStep = 0

END SUBROUTINE aed2_calculate_benthic_seddiagenesis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_seddiagenesis_update_state_benthic(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_seddiagenesis_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, tss

   ! State
   DOUBLETYPE,   DIMENSION(8)  :: dissConcs
   DOUBLETYPE,   DIMENSION(2)  :: partConcs

   ! Temporary variables
   INTEGER  :: i
!
!-------------------------------------------------------------------------------
!BEGIN

   !Check if the interval is right (sub-timestepping)
!   IF(updateStep<data%speciation_dt) THEN
!     updateStep = updateStep+1
!     RETURN
!   END IF
   updateStep = 0


   ! Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature


END SUBROUTINE aed2_seddiagenesis_update_state_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION GetDynamicBCVals(theTime) RESULT(theVals)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL, INTENT(IN) :: theTime
   AED_REAL,DIMENSION(SIZE(default_vals)):: theVals
!
!LOCALS
   ! Temporary variables
   INTEGER  :: i

!
!-------------------------------------------------------------------------------
!BEGIN


   theVals(1:nvars) = dynamic_swibc(1,1:nvars)



END FUNCTION GetDynamicBCVals
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE AED_SWI_Input(fileName)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(LEN=*), INTENT(IN) :: fileName
!
!LOCALS
   ! Temporary variables
   INTEGER  :: i,j,ncols,ntimes,openStatus


!
!-------------------------------------------------------------------------------
!BEGIN

      PRINT *,'       Prescribed particulate fluxes set from aed2_sediment_swibc.dat ... '
      OPEN(UNIT = 1000051, FILE = TRIM(fileName), STATUS = "OLD", &
                                 ACTION = "READ", IOSTAT = openStatus)

      READ (1000051,*)
      READ (1000051,*)ncols
      READ (1000051,*)ntimes

      IF(ncols/=nvars) THEN
        print *,'STOOPPOP'
        STOP
      ENDIF

      print *,' ->ncols,ntimes',ncols,ntimes
      ALLOCATE(dynamic_swibc(ntimes,ncols))
      ALLOCATE(bctime(ntimes))

      DO i = 1,ntimes
        READ (1000051,*) bctime(i),(dynamic_swibc(i,j),j=1,ncols)
      END DO



END SUBROUTINE AED_SWI_Input
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef SINGLE
!###############################################################################
SUBROUTINE PO4AdsorptionFractionLocal(PO4AdsorptionModel, &
                                 PO4tot_,            &
                                 ParticleConc_,      &
                                 Kpo4p,K,Qm,         &
                                 PO4dis,PO4par,      &
                                 thepH)
!-------------------------------------------------------------------------------
! Routine to compute fraction of PO4 adsorped to sediment/particulate concentration
!-------------------------------------------------------------------------------
    INTEGER,  INTENT(IN)  :: PO4AdsorptionModel
    SED_REAL, INTENT(IN)  :: PO4tot_, ParticleConc_
    AED_REAL, INTENT(IN)  :: Kpo4p,K,Qm
    AED_REAL, INTENT(OUT) :: PO4dis, PO4par
    AED_REAL, INTENT(IN), OPTIONAL :: thepH

!-------------------------------------------------------------------------------
!LOCALS
    AED_REAL :: buffer, f_pH, pH
    AED_REAL :: PO4tot, ParticleConc
    AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10

!
!-------------------------------------------------------------------------------
!BEGIN
   PO4dis   = zero_
   PO4par   = zero_
   buffer   = zero_
   f_pH     = one_

   ! calculate the total possible PO4 for sorption, and solids
   PO4tot        = MAX(one_e_neg_ten, PO4tot_ )       ! Co in Chao (mg)
   ParticleConc  = MAX(one_e_neg_ten, ParticleConc_ ) ! s in Chao  (mg = mol/L * g/mol * mg/g)


   IF(PO4AdsorptionModel == 1) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Ji 2008:
     !
     ! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.

     PO4par = (Kpo4p*ParticleConc) / (one_+Kpo4p*ParticleConc) * PO4tot
     PO4dis = one_ / (one_+Kpo4p*ParticleConc) * PO4tot


   ELSEIF(PO4AdsorptionModel == 2) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Chao et al. 2010:
     !
     ! Chao, X. et al. 2010. Three-dimensional numerical simulation of
     !   water quality and sediment associated processes with application
     !   to a Mississippi delta lake. J. Environ. Manage. 91 p1456-1466.
     IF(PRESENT(thepH)) THEN
       pH = thepH
       IF(thepH > 11.) pH = 11.0
       IF(thepH < 3.)  pH = 3.0

       ! -0.0094x2 + 0.0428x + 0.9574
       ! (ursula.salmon@uwa.edu.au: fPH for PO4 sorption to Fe in Mine Lakes)
       f_pH = -0.0094*pH*pH + 0.0428*pH + 0.9574
     ELSE
       f_pH = one_
     END IF

     ! calculate particulate fraction based on quadratic solution

     ! Chao Eq 16
     buffer = SQRT(((PO4tot+(1./K)-(ParticleConc*Qm*f_pH)))**2. + (4.*f_pH*ParticleConc*Qm/K))
     PO4par  = 0.5 * ((PO4tot+(1./K)+(ParticleConc*Qm*f_pH))  - buffer  )

     ! Check for stupid solutions
     IF(PO4par > PO4tot) PO4par = PO4tot
     IF(PO4par < zero_) PO4par = zero_

     ! Now set dissolved portion
     PO4dis = PO4tot - PO4par

   ELSE
     !-----------------------------------------------------
     ! No model is selected

     PO4dis = PO4tot
     PO4par = zero_

   END IF

END SUBROUTINE PO4AdsorptionFractionLocal
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

END MODULE aed2_seddiagenesis
