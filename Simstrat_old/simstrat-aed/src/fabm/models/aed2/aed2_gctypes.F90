!###############################################################################
!#                                                                             #
!# aed2_gctypes.F90                                                            #
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
!# Created March 2012                                                          #
!# Amended March 2014 (Dan)                                                    #
!###############################################################################

#include "aed2.h"


MODULE aed2_gctypes
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  PUBLIC  ! all

  integer, parameter :: dble_prec = SELECTED_REAL_KIND(15,307)
  integer, parameter :: GCHP = dble_prec
  integer, parameter :: GCHM_STR_LEN = 32

  INTEGER,   PARAMETER :: IVOID    = -9999
  REAL,      PARAMETER :: VOID     = -9999.0
  CHARACTER(*), PARAMETER :: VOID_STR = "VOID_STRING"



  INTEGER,  PARAMETER  ::  SEDP = GCHP


!  ===============================
  TYPE AEDConstDiagenesisType
     REAL                      :: db0
     INTEGER  :: imix
     REAL                      :: xs
     REAL                      :: x1
     REAL                      :: x2
     INTEGER  :: irrg
     REAL                      :: alpha0
     REAL                      :: xirrig
     REAL                      :: ventflow
     REAL                      :: w00
     REAL                      :: p0
     REAL                      :: p00
     REAL                      :: bp
     INTEGER  :: torteq
     REAL                      :: an, aa, ab
     REAL                      :: xl
     INTEGER  :: maxnpts
     INTEGER  ::numOM
     LOGICAL :: simDOM
     INTEGER  :: OMapproach
     INTEGER  :: OMModel                                    ! Dan added
     INTEGER  :: FTemswitch                                ! Dan added
     INTEGER  :: FTswitch                                     ! Dan added
     INTEGER  :: FBIOswitch                                 ! Dan added
     INTEGER  :: FINswitch                                   ! Dan added
     INTEGER  :: FInO2OnlySwitch
     INTEGER  :: FOMswitch                                   ! Dan added
     INTEGER  :: Bsolidswitch                                   ! Dan added
     INTEGER    :: timeswitch                                  ! Dan added
     REAL       :: num_days                                  ! Dan added
     REAL       :: fluxon                                                   ! Dan added
     REAL       :: fluxoff                                          ! Dan added
     REAL       :: substep                                          ! Dan added
     REAL       :: driverDT                                         ! Dan added

     LOGICAL :: simMnFe, simFeS, simX, simCaCO3, simFeCO3, simMnCO3  ,simNPAds
 !-- OMModel 1
     REAL                      :: poml2dic                 ! Dan added
     REAL                      :: pomr2dic                 ! Dan added
     REAL                      :: pomspecial2dic           ! Dan added
 !-- OMModel 2
     REAL                      :: docl2dic
     REAL                      :: donl2din
     REAL                      :: dopl2dip
     REAL                      :: pocl2docl
     REAL                      :: ponl2donl
     REAL                      :: popl2dopl
     REAL                      :: docr2docl
     REAL                      :: donr2donl
     REAL                      :: dopr2dopl
     REAL                      :: pocr2docr
     REAL                      :: ponr2donr
     REAL                      :: popr2dopr
     REAL                      :: pocvr2docr
     REAL                      :: ponvr2donr
     REAL                      :: popvr2dopr
 !-- OMModel 3
     REAL                      :: dHyd2dic                 ! Dan added
     !REAL                      :: dFer2dic                 ! Dan added
     !REAL                      :: poml2dic                ! Dan added
     REAL                      :: domr2dic                 ! Dan added
     !REAL                      :: dHyd2dFer                ! Dan added
     REAL                      :: dHyd2domr                ! Dan added
     REAL                      :: domr2pomr                ! Dan added
     REAL                      :: poml2doml                ! Dan added
     !REAL                     :: dFer2domr                ! Dan added
! Stoichiometric coefficients
! OM model 1
     REAL                      :: fracCPL                  ! Dan added
     REAL                      :: fracCPR                  ! Dan added
     REAL                      :: fracCPspecial            ! Dan added
     REAL                      :: fracNPL                  ! Dan added
     REAL                      :: fracNPR                  ! Dan added
     REAL                      :: fracNPspecial            ! Dan added
     REAL                      :: fracPPL                  ! Dan added
     REAL                      :: fracPPR                  ! Dan added
     REAL                      :: fracPPspecial            ! Dan added
! OM model 2
! Nothing.
! OM model 3
     REAL                      :: fracOAc                  ! Dan added
     REAL                      :: fracH2                   ! Dan added

     REAL                      :: fracCDHyd
     REAL                      :: fracNDHyd
     REAL                      :: fracPDHyd
     REAL                      :: fracCOAc
     REAL                      :: fracNOAc
     REAL                      :: fracPOAc
     REAL                      :: fracCH2
     REAL                      :: fracNH2
     REAL                      :: fracPH2

     REAL                      :: FTR                      !Dan added
     REAL                      :: FTT                      !Dan added
     REAL                      :: deltaGATP                !Dan added
     REAL                      :: dG0FerDHyd
     REAL                      :: dG0AerDHy
     REAL                      :: dG0AerOAc
     REAL                      :: dG0DenDHy
     REAL                      :: dG0DenOAc
     REAL                      :: dG0DenH2
     REAL                      :: dG0ManOAc
     REAL                      :: dG0IroOAc
     REAL                      :: dG0IroH2
     REAL                      :: dG0SulOAc
     REAL                      :: dG0SulH2
     REAL                      :: dG0MetOAc
     REAL                      :: dG0MetH2

     REAL                      :: kgrowthFer
     REAL                      :: kgrowthAer
     REAL                      :: kgrowthDen
     REAL                      :: kgrowthMan
     REAL                      :: kgrowthIro
     REAL                      :: kgrowthSul
     REAL                      :: kgrowthMet
     REAL                      :: kdeathFer
     REAL                      :: kdeathAer
     REAL                      :: kdeathDen
     REAL                      :: kdeathMan
     REAL                      :: kdeathIro
     REAL                      :: kdeathSul
     REAL                      :: kdeathMet

     REAL                      :: kHyd1
     REAL                      :: kHyd2
     REAL                      :: kHyd3
     REAL                      :: kHyd4
     REAL                      :: kHydN
     REAL                      :: CellWeight
     REAL                      :: fuse
     REAL                      :: Tiny
     REAL                      :: Temporary_proton
     REAL                      :: BMax
     REAL                      :: e
     REAL                      :: F
     REAL                      :: n
     REAL                      :: dPsi
     REAL                      :: KDHyd
     REAL                      :: KOAc
     REAL                      :: KH2
     REAL                      :: YDHyAer
     REAL                      :: YDHyFer
     REAL                      :: YDenDHy
     REAL                      :: YAerOAc
     REAL                      :: YDenOAc
     REAL                      :: YDenH2
     REAL                      :: YManOAc
     REAL                      :: YIroOAc
     REAL                      :: YIroH2
     REAL                      :: YSulOAc
     REAL                      :: YSulH2
     REAL                      :: YMetOAc
     REAL                      :: YMetH2

     REAL                      :: knh4ox
     REAL                      :: ktsno3
     REAL                      :: ktsox
     REAL                      :: kmnox
     REAL                      :: kmnno3
     REAL                      :: kfeox
     REAL                      :: kfesox
     REAL                      :: kfes2ox
     !REAL                      :: kfe1ox
     !REAL                      :: kfe2ox
     REAL                      :: kfeAge
     REAL                      :: kmnAge
     REAL                      :: kfeno3
     !REAL                      :: kfe2no3
     REAL                      :: kfesmn
     REAL                      :: kfesfe
     REAL                      :: kfesppt
     REAL                      :: ktsfe
     REAL                      :: kpyrite
     REAL                      :: kmnfe
     REAL                      :: ktsmn
     REAL                      :: kch4ox
     REAL                      :: kch4so4
     REAL                      :: ko2, kpo2
     REAL                      :: kno3, kpno3
     REAL                      :: kmn, kpmn
     REAL                      :: kfe, kpfe
     REAL                      :: kso4, kpso4
     REAL                      :: lo2, lpo2                 ! Dan added
     REAL                      :: lno3, lpno3               ! Dan added
     REAL                      :: lmn, lpmn                 ! Dan added
     REAL                      :: lfe, lpfe                 ! Dan added
     REAL                      :: lso4, lpso4               ! Dan added
     REAL                      :: kanh4, kapo4
     INTEGER  :: ibc2
     INTEGER  :: startSteady
     REAL                      :: flux_scale
     REAL                      :: POMVR
     CHARACTER(LEN=5)          :: OMInitMethodL
     CHARACTER(LEN=5)          :: OMInitMethodR
     REAL                      :: OM_topL, OM_minL, OM_cfL, InitMinDepthL
     REAL                      :: OM_topR, OM_minR, OM_cfR, InitMinDepthR
     REAL                      :: kSidppt
     REAL                      :: kCalppt
     REAL                      :: kRodppt
     REAL                      :: kMnO2Appt
     REAL                      :: kMnO2Bppt
     REAL                      :: kFeOHAppt
     REAL                      :: kFeOHBppt
     CHARACTER(LEN=5)          :: Xname
     REAL                      :: Xmk, Xfl, Xfm
     REAL                      :: kXSppt
     INTEGER  :: ibbc
     INTEGER  :: rxn_mode, OutputUnits
     REAL                      :: InitSedDepth
     INTEGER  :: PO4AdsorptionModel
     REAL                      :: KPO4p
     REAL                      :: Kadsratio
     REAL                      :: Qmax
     LOGICAL                   :: ads_use_pH
  END TYPE AEDConstDiagenesisType
!-------------------------------------------------------------------------------

!  ===============================
  TYPE AEDChmConstType
     LOGICAL                     :: TRA_PARS
     LOGICAL                     :: DOX_PARS
     LOGICAL                     :: CH4_PARS
     LOGICAL                     :: BAC_PARS
     LOGICAL                     :: CNP_PARS
     LOGICAL                     :: IRN_PARS
     LOGICAL                     :: MAN_PARS
     LOGICAL                     :: ASN_PARS
     LOGICAL                     :: SUL_PARS
     LOGICAL                     :: ALM_PARS
     LOGICAL                     :: ZNC_PARS
     LOGICAL                     :: GCH_DATA
     LOGICAL                     :: ORG_PARS
#if 0
     TYPE (AEDConstCNPType)      :: CNP
     TYPE (AEDConstFeType)       :: Iron
     TYPE (AEDConstMnType)       :: Manganese
     TYPE (AEDConstAsType)       :: Arsenic
     TYPE (AEDConstSuType)       :: Sulfur
     TYPE (AEDConstCH4Type)      :: Methane
     TYPE (AEDConstAlType)       :: Aluminium
     TYPE (AEDConstZnType)       :: Zinc
     TYPE (AEDConstDOType)       :: DissolvedOxygen
     TYPE (AEDConstPOCType)      :: Particulate_Organic_Carbon
     TYPE (AEDConstColourType)   :: Tracer
     TYPE (AEDConstBacType)      :: Bacteria
     TYPE (AEDConstOrgChmType)   :: OrgComp
#endif
  END TYPE AEDChmConstType


  ! -------------------------------------------------------------
  ! ---------- 1. Geochemical component structure, Xj  ----------

  TYPE gcUnknowns
    !-- Required at data prep time:
    CHARACTER(GCHM_STR_LEN)      :: EltName   =  VOID_STR !Element Name
    CHARACTER(GCHM_STR_LEN)      :: CompName  =  VOID_STR !Component Name
    INTEGER  :: CompType  =  IVOID    !ComponentType
    INTEGER  :: CompIndex =  IVOID    !Component Number, j
    DOUBLETYPE                   :: MolWeight =  VOID     !Molecular Weight (u)
    INTEGER  :: Charge    =  IVOID    !Component  charge
    ! -- Required at run-time:
    INTEGER  :: WQindex   =  IVOID    !Index in AED WQ arrays
    INTEGER  :: EQindex   =  IVOID    !Index in the inequality array
    DOUBLETYPE                   :: Value     =  VOID     !Number of Moles, Xj
    DOUBLETYPE                   :: Total     =  VOID     !Total Moles, Tj
    DOUBLETYPE                   :: Delta     =  VOID     !dX
    TYPE(gcSpecies),     POINTER :: Master    => NULL()
    TYPE(PurePhaseInfo), POINTER :: PPData    => NULL()   !Only required if Type = PUREPHASE
  END TYPE gcUnknowns

  ! ------------------------------------------------------
  ! ---------- 2. Aqueous species structure, Ci ----------

  TYPE gcSpecies
    CHARACTER(GCHM_STR_LEN)      :: Name         =  VOID_STR !Species Name
    INTEGER  :: SpeciesIndex =  IVOID    !Species Number, i
    DOUBLETYPE                   :: Moles        =  VOID     !Ci
    INTEGER  :: Charge       =  IVOID    !qi
    DOUBLETYPE                   :: logKat25     =  VOID     !log K(T=25C)
    DOUBLETYPE                   :: logK         =  VOID     !log K; K = equilib. const.
    DOUBLETYPE                   :: deltaH       =  VOID     !dH
    INTEGER  :: Gflag        =  IVOID    !Gamma method
    DOUBLETYPE                   :: logMoles     =  VOID     !log(ni)
    DOUBLETYPE                   :: logActivity  =  VOID     !ln(ai)
    DOUBLETYPE                   :: logGamma     =  VOID     !logG
    DOUBLETYPE                   :: delGamma     =  VOID     !d/dG
    DOUBLETYPE                   :: H2OStoich    =  VOID     !
    DOUBLETYPE                   :: HStoich      =  VOID     !
    ! -- Stoichiometric coeffs:
    REAL,  DIMENSION(:), POINTER :: Stoich => NULL()     !stoichiometric coeffs
  END TYPE gcSpecies
  ! ------------------------------------------------------
  ! ---------- 3. Pure phase -----------------------------

  ! -- Subtype for gcUnknows:
  TYPE PurePhaseInfo
    DOUBLETYPE                   :: logKat25    = VOID  !log K(T=25C)
    DOUBLETYPE                   :: logK        = VOID  !log K; K = equilib. const.
    DOUBLETYPE                   :: deltaH      = VOID  !
    DOUBLETYPE                   :: Diameter    = VOID  !diameter of particle  [m]
    DOUBLETYPE                   :: Moles       = VOID  !quantity of material
    DOUBLETYPE                   :: SettlingVel = VOID  !settling vel     [m s^-1]
    ! -- Stoichiometric coeffs:
    REAL,  DIMENSION(:), POINTER :: Stoich => NULL()    !stoichiometric coeffs
  END TYPE PurePhaseInfo
!!##/GEOCH_MOD


!#----------------- sed_candi types ----------------------------------
!
  ! maximal value for prescribed corgflux values
  INTEGER,    PARAMETER :: MAXCORG=2000

  ! readspecies=maximal number of included in CCANDI
  !        (for read in STEUER.f/STEUER.DAT)
  INTEGER,    PARAMETER :: readspecies = 26

TYPE aed2_sed_candi_t
  !-- Module switches
  LOGICAL  :: simCaXCO3    = .FALSE.
  LOGICAL  :: simCaCO3     = .FALSE.
  LOGICAL  :: simFeCO3     = .FALSE.
  LOGICAL  :: simMnCO3     = .FALSE.
  LOGICAL  :: simFeII      = .FALSE.
  LOGICAL  :: simMnFe      = .FALSE.
  LOGICAL  :: simFeS       = .FALSE.
  LOGICAL  :: simCaCO3C12  = .FALSE.
  LOGICAL  :: simC12       = .FALSE.
  LOGICAL  :: simTracer    = .FALSE.
  LOGICAL  :: simRefOM     = .FALSE.
  LOGICAL  :: simX         = .FALSE.
  LOGICAL  :: simXS        = .FALSE.
  LOGICAL  :: simXO        = .FALSE.
  LOGICAL  :: simAdsp      = .FALSE.

  INTEGER  :: eqmethod, mod_kgvar, mod_mnfe,                                   &
              mod_feii, mod_tracer, mod_po4solid,                              &
              mod_fes, mod_caxco3,mod_caco3,                                   &
              mod_C12, mod_CaCO3C12
  INTEGER  :: OMModel                                                 !Dan added
  INTEGER  :: OMapproach                                              !Dan added
  INTEGER  :: FTemswitch                                              !Dan added
  INTEGER  :: FBIOswitch                                              !Dan added
  INTEGER  :: FTswitch                                                !Dan added
  INTEGER  :: FINswitch                                               !Dan added
  INTEGER  :: timeswitch                                              !Dan added
  INTEGER  :: FOMswitch                                               !Dan added
  INTEGER  :: FInO2OnlySwitch
  INTEGER  :: BsolidSwitch
  REAL     :: num_days                                                !Dan added
  REAL     :: fluxon                                                  !Dan added
  REAL     :: fluxoff                                                 !Dan added
  REAL     :: substep                                                 !Dan added
  REAL     :: driverDT                                                !Dan added
  INTEGER  :: kgmod       = 0           ! OM variables
  INTEGER  :: MnFemod     = 0           ! Mn & Fe asscoiated variables
  INTEGER  :: alkmod      = 0           ! Alkalinity assoc. variables
  INTEGER  :: rmmod       = 0           !
  INTEGER  :: FeIImod     = 0           ! FeII associated variables
  INTEGER  :: Tracermod   = 0           ! Num of "non-reactive" vars
  INTEGER  :: FeSmod      = 0           ! FeS(Pyrite) variable
  INTEGER  :: PO4smod     = 0           ! Solid PO4 variable
  INTEGER  :: CaXCO3mod   = 0           ! Mn & Fe carbonate variables
  INTEGER  :: CaCO3mod    = 0           ! Calcite variables
  INTEGER  :: FeCO3mod    = 0           ! Siderite variables
  INTEGER  :: MnCO3mod    = 0           ! Rhodocrosite variables
  INTEGER  :: C12mod      = 0           ! Variables assoc. with C12
  INTEGER  :: CaCO3C12mod = 0           ! Variables assoc. with C12
  INTEGER  :: xmod        = 0           ! Variables assoc. with X

  !-- Column inidicies
  INTEGER  :: O2y
  INTEGER  :: NO3y
  INTEGER  :: SO4y
  INTEGER  :: PO4ly
  INTEGER  :: NH4y
  INTEGER  :: CH4y
  INTEGER  :: DICy
  INTEGER  :: protony
  INTEGER  :: tbohy
  INTEGER  :: HSy
  INTEGER  :: CaY
  INTEGER  :: Xy
  INTEGER  :: N2y
  INTEGER  :: OAcy
  INTEGER  :: H2y
  INTEGER  :: HCO3y
  INTEGER  :: DOCLy
  INTEGER  :: POCLy
  INTEGER  :: DOCRy
  INTEGER  :: POCRy
  INTEGER  :: DONLy
  INTEGER  :: PONLy
  INTEGER  :: DONRy
  INTEGER  :: PONRy
  INTEGER  :: DOPLy
  INTEGER  :: POPLy
  INTEGER  :: DOPRy
  INTEGER  :: POPRy
  INTEGER  :: POXy
  INTEGER  :: DOXy

  INTEGER  :: DOMRy           !Dan added
  INTEGER  :: POMLy           !Dan added
  INTEGER  :: POMRy           !Dan added
  INTEGER  :: POMspecialy     !Dan added
  INTEGER  :: dhydy           !Dan added
  INTEGER  :: POM1y           !Dan added
  INTEGER  :: POM2y           !Dan added
  INTEGER  :: POM3y           !Dan added
  INTEGER  :: POM4y           !Dan added
  INTEGER  :: BAery
  INTEGER  :: BDeny
  INTEGER  :: BMany
  INTEGER  :: BIroy
  INTEGER  :: BSuly
  INTEGER  :: BMety
  INTEGER  :: BFery
  INTEGER  :: Necromassy

  INTEGER  :: MnO2y
  INTEGER  :: FeOHy
  INTEGER  :: MnO2By
  INTEGER  :: FeOHBy
  INTEGER  :: MnIVy
  INTEGER  :: MnIIy
  INTEGER  :: FeIIy
  INTEGER  :: FeIIIy
  INTEGER  :: FeSY
  INTEGER  :: FeS2Y
  INTEGER  :: PO4sy
  INTEGER  :: NH4sy
  INTEGER  :: XSy
  INTEGER  :: araY
  INTEGER  :: calY
  INTEGER  :: sidY
  INTEGER  :: rodY
  INTEGER,           DIMENSION(:), ALLOCATABLE :: TracerY
  CHARACTER(LEN=12), DIMENSION(:), ALLOCATABLE :: CANDIVarNames
  LOGICAL,           DIMENSION(:), ALLOCATABLE :: isSolid


  !-- Configuration constants
  INTEGER  :: MAXNPTS                              !Deepest sediment layer
  INTEGER  :: nSPECIES
  INTEGER  :: maxneq
  INTEGER  :: MU
  INTEGER  :: ML
  INTEGER  :: LRW
  INTEGER  :: LIW

  ! Model parameters
  REAL(SEDP) :: corg(MAXCORG)
  REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: DIFFC
  ! Reaction rates
  REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: reac
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGC
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGN
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGP
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROXDHYD
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROXDOMR

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FNO3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FMnO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FFeOH
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FSO4
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FMet

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNH4OX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnOX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeOX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSOX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RCH4OX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSOX
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeS2OX

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNH4NO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnNO3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeNO3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSNO3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeMnA
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeMnB
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSMnA
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSMnB
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSMnA
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSMnB
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSFeA
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSFeB
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSFeA
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSFeB
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RCH4SO4

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnAge
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeOHAppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeOHBppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeAge
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPyrite
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RXSppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RSidppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RRodppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RCalppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnO2Appt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnO2Bppt
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPO4ads
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNH4ads

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPOM1
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPOM2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPOM3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPOM4
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RPOMspecial
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNecro
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RAerDHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RDenO2DHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RDenNO3DHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFerDHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RDHyd

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathFer
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathAer
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathDen
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathMan
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathIro
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathSul
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathMet
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathTot

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RAerOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RDenOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RManOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RIroOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RSulOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMetOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RAerH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RDenH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RManH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RIroH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RSulH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMetH2

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGFerDHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTFerDHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGAerOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTAerOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGDenOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTDenOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGDenH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTDenH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGManOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTManOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGIroOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTIroOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGIroH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTIroH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGSulOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTSulOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGSulH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTSulH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGMetOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTMetOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGMetH2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTMetH2

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: Btot
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FBHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FDHyd
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FOAc
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FH2
    !Dan added
  ! Reaction constants
  ! OM model 1
  REAL(SEDP) :: poml2dic, pomr2dic, pomspecial2dic              !Dan added
  ! OM model 2
  REAL(SEDP) :: docl2dic,   donl2din,   dopl2dip
  REAL(SEDP) :: pocl2docl,  ponl2donl,  popl2dopl
  REAL(SEDP) :: docr2docl,  donr2donl,  dopr2dopl
  REAL(SEDP) :: pocr2docr,  ponr2donr,  popr2dopr
  REAL(SEDP) :: pocvr2docr, ponvr2donr, popvr2dopr
  ! OM model 3
  REAL(SEDP) :: domr2pomr                                       !Dan added
  REAL(SEDP) :: poml2doml                                       !Dan added
  REAL(SEDP) :: domr2dic                                        !Dan added
  !REAL(SEDP) :: poml2dic                                       !Dan added
  !REAL(SEDP) :: pomr2dic                                       !Dan added
  ! FT numbers
  REAL(SEDP) :: FTR                                             !Dan added
  REAL(SEDP) :: FTT                                             !Dan added
  REAL(SEDP) :: deltaGATP                                       !Dan added
  REAL(SEDP) :: dGmp                                            !Dan added
  REAL(SEDP) :: BMax                                            !Dan added

  REAL(SEDP) :: O2perdhyd                                       !Dan added
  REAL(SEDP) :: no3perdhyd                                      !Dan added
  REAL(SEDP) :: MnO2perdhyd                                     !Dan added
  REAL(SEDP) :: FeOHperdhyd                                     !Dan added
  REAL(SEDP) :: SO4perdhyd                                      !Dan added
!
  REAL(SEDP) :: fracOAc
  REAL(SEDP) :: fracH2

  REAL(SEDP) :: YDHyAer
  REAL(SEDP) :: YDHyFer
  REAL(SEDP) :: YDenDHy
  REAL(SEDP) :: YAerOAc
  REAL(SEDP) :: YDenOAc
  REAL(SEDP) :: YDenH2
  REAL(SEDP) :: YManOAc
  REAL(SEDP) :: YIroOAc
  REAL(SEDP) :: YIroH2
  REAL(SEDP) :: YSulOAc
  REAL(SEDP) :: YSulH2
  REAL(SEDP) :: YMetOAc
  REAL(SEDP) :: YMetH2

  REAL(SEDP) ::dG0FerDHyd
  REAL(SEDP) ::dG0AerDHy
  REAL(SEDP) ::dG0AerOAc
  REAL(SEDP) ::dG0DenDHy
  REAL(SEDP) ::dG0DenOAc
  REAL(SEDP) ::dG0DenH2
  REAL(SEDP) ::dG0ManOAc
  REAL(SEDP) ::dG0IroOAc
  REAL(SEDP) ::dG0IroH2
  REAL(SEDP) ::dG0SulOAc
  REAL(SEDP) ::dG0SulH2
  REAL(SEDP) ::dG0MetOAc
  REAL(SEDP) ::dG0MetH2

  REAL(SEDP) ::e
  REAL(SEDP) ::F
  REAL(SEDP) ::n
  REAL(SEDP) ::dPsi

  REAL(SEDP) ::KDHyd
  REAL(SEDP) ::KOAc
  REAL(SEDP) ::KH2

  REAL(SEDP) :: kgrowthFer
  REAL(SEDP) :: kgrowthAer
  REAL(SEDP) :: kgrowthDen
  REAL(SEDP) :: kgrowthMan
  REAL(SEDP) :: kgrowthIro
  REAL(SEDP) :: kgrowthSul
  REAL(SEDP) :: kgrowthMet

  REAL(SEDP) ::kdeathFer
  REAL(SEDP) ::kdeathAer
  REAL(SEDP) ::kdeathDen
  REAL(SEDP) ::kdeathMan
  REAL(SEDP) ::kdeathIro
  REAL(SEDP) ::kdeathSul
  REAL(SEDP) ::kdeathMet

  REAL(SEDP) ::kHyd1
  REAL(SEDP) ::kHyd2
  REAL(SEDP) ::kHyd3
  REAL(SEDP) ::kHyd4
  REAL(SEDP) ::kHydN


!  REAL(SEDP) ::domr2dic
!  REAL(SEDP) ::domr2pomr
!  REAL(SEDP) ::poml2doml
  REAL(SEDP) ::fuse
  REAL(SEDP) ::CellWeight
  REAL(SEDP) ::Tiny
  REAL(SEDP) ::Temporary_proton
    ! Stoichiometric coefficients
  ! OM model 1
  REAL(SEDP) :: fracCPL, fracCPR, fracCPspecial                 !Dan added
  REAL(SEDP) :: fracNPL, fracNPR, fracNPspecial                 !Dan added
  REAL(SEDP) :: fracPPL, fracPPR, fracPPspecial                 !Dan added

  ! OM model 2
  ! Nothing.

  ! OM model 3
  REAL(SEDP) :: fracCDHyd, fracNDHyd, fracPDHyd                 !Dan added
  REAL(SEDP) :: fracCOAc , fracNOAc , fracPOAc                  !Dan added
  REAL(SEDP) :: fracCH2  , fracNH2  , fracPH2                   !Dan added
  REAL(SEDP) :: kO2,   kpO2,   lO2,   lpO2
  REAL(SEDP) :: kNO3,  kpNO3,  lNO3,  lpNO3
  REAL(SEDP) :: kMnO2, kpMnO2, lMnO2, lpMnO2
  REAL(SEDP) :: kFeOH, kpFeOH, lFeOH, lpFeOH
  REAL(SEDP) :: kSO4,  kpSO4,  lSO4,  lpSO4
  REAL(SEDP) :: kCH4,  kpCH4 ! I think this parameter is not used - Dan

  REAL(SEDP) :: kNH4OX
  REAL(SEDP) :: kMnOX
  REAL(SEDP) :: kFeOX
  REAL(SEDP) :: kTSOX
  REAL(SEDP) :: kCH4OX
  REAL(SEDP) :: kFeSOX
  REAL(SEDP) :: kFeS2OX

  REAL(SEDP) :: kNH4NO2
  REAL(SEDP) :: kMnNO3
  REAL(SEDP) :: kFeNO3
  REAL(SEDP) :: kTSNO3
  REAL(SEDP) :: kMnFe
  REAL(SEDP) :: kTSMn
  REAL(SEDP) :: kFeSMn
  REAL(SEDP) :: kTSFe
  REAL(SEDP) :: kFeSFe
  REAL(SEDP) :: kCH4SO4

  REAL(SEDP) :: kMnAge
  REAL(SEDP) :: kFeOHAppt
  REAL(SEDP) :: kFeOHBppt
  REAL(SEDP) :: kFeAge
  REAL(SEDP) :: kFeSppt
  REAL(SEDP) :: kPyrite
  REAL(SEDP) :: kXSppt
  REAL(SEDP) :: kSidppt
  REAL(SEDP) :: kRodppt
  REAL(SEDP) :: kCalppt
  REAL(SEDP) :: kMnO2Appt
  REAL(SEDP) :: kMnO2Bppt
  REAL(SEDP) :: kPO4ads
  REAL(SEDP) :: kNH4ads

  REAL(SEDP) :: KSPFES,KSPPO4
  REAL(SEDP) :: K(12)

  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: kg0var
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: kgpo4up

  REAL(SEDP) :: sRNHOX,  sRSOX, sRMNOX, sRFEOX, sRMNFE, sRMNO2TS
  REAL(SEDP) :: sRCH4OX, sRCH4SO4
  REAL(SEDP) :: sRNO3TS, sRTSFE3
  REAL(SEDP) :: sRDCO3,  sRDHCO3, sRDCO2
  REAL(SEDP) :: sRFESFE3,sRFESMN4,sRFESOX
  REAL(SEDP) :: sRMnNO3, sRFENO3
  REAL(SEDP) :: spo4sppt,spo4sdis
  REAL(SEDP) :: sRCALDIS,sRCALPPT,sRARADIS, sRARAPPT
  REAL(SEDP) :: sRMnPPT, sRMnDIS, sRFePPT,  sRFeDIS
  REAL(SEDP) :: sRFe1OX, sRFe2OX, sRFe1NO3, sRFe2NO3

  REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: Y
  REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: Ydot
  REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: Ytemp
  REAL, DIMENSION(:,:), ALLOCATABLE :: IAP
  REAL, DIMENSION(:,:), ALLOCATABLE :: KIAP
  REAL, DIMENSION(:,:), ALLOCATABLE :: QIAP

  REAL(SEDP) :: co2suc,hco3suc,hssuc,count

  REAL(SEDP) :: W00,ventflow
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: pps
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: psp
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: poros
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: t2
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dpdx
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dt2dx
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: uvel
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: wvel
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: bioturb
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dbdx
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: ps
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: cirrig
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dff
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dudx
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: poros_bg
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: poros_dot

  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: bottom, top_bound, FixedBottomConc
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dh
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dh2
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dhr
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dhf
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: RPAR
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: rtol,atol

  REAL(SEDP) :: MNO20,FEOH0, PB0, ASH0, Fe10,Fe20,CaCO30
  REAL(SEDP) :: O2L,NO3L,SO4L,TPO4L,TNH4L,TH2SL,MNIIL,FEIIL
  REAL(SEDP) :: TCO2L,CH4L,MNO2L,FEOHL,G0L,G1L,G2L,FESL,CAL,CLL
  REAL(SEDP) :: PBL,ASHL,Fe1L,Fe2L
  REAL(SEDP) :: G00,G10,G20
  REAL(SEDP) :: O2W,NO3W,SO4W,TPO4W,TNH4W,TH2SW
  REAL(SEDP) :: MNIIW,FEIIW,TCO2W,ALKW,CH4W,CAW,TBOH4W,CLW
  REAL(SEDP) :: CO2W, CO3W, HCO3W, HSW, H2SW
  REAL(SEDP) :: fluxirr, fluxdiff, fluxadv

  INTEGER  :: imix
  REAL(SEDP) :: XL,XLdouble,TDH
  REAL(SEDP) :: ALPHA0,XIRRIG,porosmin
  REAL(SEDP) :: DB0,XS,X1,X2,nld,lld
  REAL(SEDP) :: P0,P00,BP
  REAL(SEDP) :: AN, aa, ab
  REAL(SEDP) :: PO4salpha, po4sat,kpo4up
  INTEGER  :: torteq
  INTEGER  :: IBBC
  REAL(SEDP), DIMENSION(:), ALLOCATABLE :: PartFluxes

  REAL(SEDP) :: FG0,FG1,FG2,fash,fcaco3
  REAL(SEDP) :: fluxo2,fluxno3,fluxso4,fluxtpo4,fluxtnh4
  REAL(SEDP) :: fluxmnii,fluxfeii,fluxco2,fluxhco3,fluxco3
  REAL(SEDP) :: fluxh2s,fluxhs,fluxboh3,fluxboh4,fluxth2s
  REAL(SEDP) :: fluxca
  REAL(SEDP) :: fluxtco2,fluxch4,fluxhco3Y,fluxco2Y,fluxco3Y
  REAL(SEDP) :: fluxh2sY,fluxhsY
  REAL(SEDP) :: flux_scale

  INTEGER :: npt,oldday,randini
  CHARACTER (LEN=40) :: FileCon, FileRat, FilePar, FileSpe, FileDanPar, FileRat2
  CHARACTER (LEN=40) :: fname(0:readspecies-1), corgfile,FileAuf, FileratOAc
  REAL(SEDP) :: temp, salt, PRES, tout, deltat,xMn, xFe
  INTEGER :: iso4, irrg, iSTEADY, IBC2, startSteady, iPRINT
  INTEGER :: aufflagout, useVODE, hofmu, rxn_mode
  LOGICAL :: restartCANDI
  INTEGER :: hofw(0:readspecies-1)
  INTEGER :: writestep,corginp,ashflux,caco3flux,job
  REAL(SEDP) :: D13CCH4SW,D13CCH4DF,D13CTCO2SW,D13CTCO2DF,D13CGI
  REAL(SEDP) :: alphaM, alphaOM, alphaAOM, PDB
  REAL(SEDP) :: Eco2hco3, Eco3hco3
  REAL(SEDP) :: D13CCAL, D13CARA, alphaCAL, alphaARA
  REAL(SEDP) :: SC,SN,SP   ! OM Stoichiometry
  REAL(SEDP) :: MK,FL,FM   ! Mineral Metal Stoichiometry
  REAL(SEDP) :: eqtol,lost,zeitschritt,day
  REAL(SEDP) :: DF(28)

  REAL(SEDP) :: POMVR, InitSedDepth
  REAL(SEDP) :: InitMinDepthL, OM_topL, OM_minL, OM_cfL
  REAL(SEDP) :: InitMinDepthR, OM_topR, OM_minR, OM_cfR
  CHARACTER(LEN=4) :: OMInitMethodL, OMInitMethodR
  CHARACTER(LEN=4) :: Xname
  INTEGER :: OutputUnits
END TYPE


END MODULE aed2_gctypes
!------------------------------------------------------------------------------!
