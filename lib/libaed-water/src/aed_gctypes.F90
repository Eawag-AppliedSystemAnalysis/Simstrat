!###############################################################################
!#                                                                             #
!# aed_gctypes.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2012 - 2021 -  The University of Western Australia               #
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
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created March 2012                                                          #
!# Amended March 2014 (Dan)                                                    #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_gctypes
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  PUBLIC  ! all

  INTEGER, PARAMETER      :: dble_prec    = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER      :: GCHP         = dble_prec
  INTEGER, PARAMETER      :: GCHM_STR_LEN = 32

  INTEGER,   PARAMETER    :: IVOID        = -9999
  REAL,      PARAMETER    :: VOID         = -9999.0
  CHARACTER(*), PARAMETER :: VOID_STR     = "VOID_STRING"

  INTEGER,  PARAMETER     :: SEDP         = GCHP

  ! -------------------------------------------------------------
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


  ! -------------------------------------------------------------
  ! ---------- 1. Geochemical component structure, Xj  ----------

  TYPE gcUnknowns
    !-- Required at data prep time:
    CHARACTER(GCHM_STR_LEN)      :: EltName   =  VOID_STR !Element Name
    CHARACTER(GCHM_STR_LEN)      :: CompName  =  VOID_STR !Component Name
    INTEGER                      :: CompType  =  IVOID    !ComponentType
    INTEGER                      :: CompIndex =  IVOID    !Component Number, j
    DOUBLETYPE                   :: MolWeight =  VOID     !Molecular Weight (u)
    INTEGER                      :: Charge    =  IVOID    !Component  charge
    ! -- Required at run-time:
    INTEGER                      :: WQindex   =  IVOID    !Index in AED WQ arrays
    INTEGER                      :: EQindex   =  IVOID    !Index in the inequality array
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
    INTEGER                      :: SpeciesIndex =  IVOID    !Species Number, i
    DOUBLETYPE                   :: Moles        =  VOID     !Ci
    INTEGER                      :: Charge       =  IVOID    !qi
    DOUBLETYPE                   :: logKat25     =  VOID     !log K(T=25C)
    DOUBLETYPE                   :: logK         =  VOID     !log K; K = equilib. const.
    DOUBLETYPE                   :: deltaH       =  VOID     !dH
    INTEGER                      :: Gflag        =  IVOID    !Gamma method
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

#if 0
!#----------------- sed_candi types ----------------------------------
!#----------------- sed_candi types ----------------------------------
!#----------------- sed_candi types ----------------------------------
!
#ifdef SINGLE
# define SED_REAL REAL
#else
# define SED_REAL REAL(SEDP)
#endif

  ! ------------------------------------------------------
   TYPE aed_sed_candi_time_t
      !-- Timey wimey things --!
      INTEGER  :: timeswitch       !# 0       ! If timeswitch = 1 then hard coded pom flux will be switched
      REAL     :: num_days         !# Dan added
      REAL     :: firststeps       !# 5       ! If time<fluxon then flux will be zero, if time>fluxon then flux will be default_vals
      REAL     :: fluxon           !# 5       ! If time<fluxon then flux will be zero, if time>fluxon then flux will be default_vals
      REAL     :: fluxoff          !# 15      ! If time>fluxoff then flux will be zero
      REAL     :: justwaitabit     !# 15      ! If time>fluxoff then flux will be zero
      REAL     :: substep          !# 240     ! Hourly if glm is hourly, i.e. 3600 s; CANDI all units are y
      REAL     :: substep_0
      REAL     :: substep_1
      REAL     :: substep_2
      REAL     :: substep_3
      REAL     :: driverDT         !# 3600.   ! s Should be 3600 if hourly
         END TYPE aed_sed_candi_time_t

  ! ------------------------------------------------------
   TYPE aed_sed_candi_bc_t
      !-- Initial & Boundary Conditions  --!
      INTEGER  :: ibc2             !#   2
      INTEGER  :: ibbc             !#   1
      INTEGER  :: startSteady      !#  1
      SED_REAL :: flux_scale       !#  1.0
      SED_REAL :: POMVR            !#  0.3 ! 0.0210
   !  INTEGER  :: Irrigswitch      !# 1
      SED_REAL :: InitSedDepth     !# 2.0e-10
      CHARACTER(LEN=4) :: OMInitMethodL    !# CO_I
      SED_REAL :: OM_topL          !#  1.0
      SED_REAL :: OM_minL          !#  0.9
      SED_REAL :: OM_cfL           !#  0.6
      SED_REAL :: InitMinDepthL    !# 5
      CHARACTER(LEN=4) :: OMInitMethodR    !# EX_I
      SED_REAL :: OM_topR          !#  1.0
      SED_REAL :: OM_minR          !#  0.9
      SED_REAL :: OM_cfR           !#  0.6
      SED_REAL :: InitMinDepthR    !# 5
      SED_REAL :: InMinDep         !# 99.
      SED_REAL :: OutputUnits      !#  0
   END TYPE aed_sed_candi_bc_t

  ! ------------------------------------------------------
   TYPE aed_sed_candi_param_t
      !-- Sediment Physical Properties --!
      SED_REAL :: db0              !#  20.00
      INTEGER  :: imix             !#  1
      SED_REAL :: xs               !#  20.0E+0
      SED_REAL :: x1               !#  15.0e-0
      SED_REAL :: x2               !#  35.0e-0
      INTEGER  :: irrg             !#  1
      SED_REAL :: alpha0           !#  30.00  ! y^-1
      INTEGER  :: xirrig           !#  7   ! cm
      SED_REAL :: poreflux         !#  -15.0E-05 !
      SED_REAL :: w00              !#  0.15 ! cm y^-1 WVC 1996
      LOGICAL  :: w00h00_0n
      LOGICAL  :: pf_on
      SED_REAL :: p0               !#  0.9 !
      SED_REAL :: p00              !#  0.6 !
      SED_REAL :: density          !#  2.5 ! kg L^-1
      SED_REAL :: bp               !#  0.6
      INTEGER  :: torteq           !#  2
      SED_REAL :: an               !#  2.14
      SED_REAL :: aa               !#  3.79
      SED_REAL :: ab               !#  2.02
      SED_REAL :: xl               !#  40.0
      INTEGER  :: maxnpts          !#  100
      INTEGER  :: job              !#  1 ! Job 0 = even grid; job 1 = exponential

      !-- Biogeochemical Configuration --!
      INTEGER  :: numOM            !# 3
      LOGICAL  :: simDOM           !# .true.
      INTEGER  :: OMapproach       !# 1
      INTEGER  :: OMModel          !# 1
      INTEGER  :: FTemswitch       !# 2
      INTEGER  :: fancyNswitch     !# 1
      INTEGER  :: FTswitch
      INTEGER  :: FBIOswitch       !# 1
      INTEGER  :: FINswitch        !# 2
      INTEGER  :: FInO2Onlyswitch  !# 2
      INTEGER  :: FOMswitch        !# 1
      INTEGER  :: Bsolidswitch     !# 1

      !# mmolLsolid !umolg  !mM ! percentsolids ! mM ! percentsolidsmM ! percentsolids ! mM ! percentsolids
      CHARACTER(len=20) :: SolidInitialUnit

!     SED_REAL :: SolidFluxUnit    !#  mMm2y ! percentsolids ! mMm2y ! percentsolids
!     SED_REAL :: SolidOutputUnit  !#  mmolLsolids !mM

      LOGICAL  :: VCW              !# .false.
      LOGICAL  :: simMnFe          !# .true.
      LOGICAL  :: simFeS           !# .true.
      LOGICAL  :: simX             !# .false.
      LOGICAL  :: simCaCO3         !# .false.
      LOGICAL  :: simFeCO3         !# .false.
      LOGICAL  :: simMnCO3         !# .false.
      LOGICAL  :: simNPAds         !# .false.
      INTEGER  :: rxn_mode         !#  2  ! Flag for whether IAP/K is used for mineral precip
      LOGICAL  :: ads_use_pH       !# .true.

      !-- Reaction rate constants --!
        !-- Primary reactions --!
        !OMModel 1
      SED_REAL :: poml2dic         !# 0.5
      SED_REAL :: pomr2dic         !# 0.0003        ! From the Moon
      SED_REAL :: pomspecial2dic   !# 2.00        ! From Brigolin
      SED_REAL :: R0               !# 135.0       ! mmol L^-1 y^-1Van Cappellen and Wang 96
      SED_REAL :: VCWBeta          !# 0.208       ! cm^-1 Van Cappellen and Wang 96
         !OMModel 2
      SED_REAL :: docl2dic         !#  0.0 ! 85.001  !y^-1 From Norlem
      SED_REAL :: donl2din         !#  0.0 ! 100.001 !y^-1 From Norlem
      SED_REAL :: dopl2dip         !#  0.0 ! 100.001 !y^-1 From Norlem
      SED_REAL :: pocl2docl        !#  0.0 ! 20.0 !y^-1 From Norlem
      SED_REAL :: ponl2donl        !#  0.0 ! 25.0 !y^-1 From Norlem
      SED_REAL :: popl2dopl        !#  0.0 ! 25.0 !y^-1 From Norlem
      SED_REAL :: docr2docl        !#  0.0 ! 0.20 !y^-1 From Norlem
      SED_REAL :: donr2donl        !#  0.0 ! 0.20 !y^-1 From Norlem
      SED_REAL :: dopr2dopl        !#  0.0 ! 0.20 !y^-1  From Norlem
      SED_REAL :: pocr2docr        !#  0.0 ! 0.03 !y-1 From Norlem
      SED_REAL :: ponr2donr        !#  0.0 ! 0.03 !y^-1 From Norlem
      SED_REAL :: popr2dopr        !#  0.0 ! 0.03 !y^-1 From Norlem
      SED_REAL :: pocvr2docr       !#  0.0 ! 1.0e-7 !y^-1 From Norlem
      SED_REAL :: ponvr2donr       !#  0.0 ! 1.0e-7 !y^-1 From Norlem
      SED_REAL :: popvr2dopr       !#  0.0 ! 1.0e-7 !y^-1 From Norlem
         !OMModel 3
      SED_REAL :: kgrowthFer       !#  5.0E1 ! 055.000 ! y^-1
      SED_REAL :: kgrowthAer       !#  5.0E1 !  055.000 ! y^-1
      SED_REAL :: kgrowthDen       !#  5.0E1 !  055.000 ! y^-1
      SED_REAL :: kgrowthMan       !#  5.0E1 !  055.000 ! y^-1
      SED_REAL :: kgrowthIro       !#  5.0E1 !  055.000 ! y^-1
      SED_REAL :: kgrowthSul       !#  5.0E1 !  055.000 ! y^-1
      SED_REAL :: kgrowthMet       !#  5.0E1 !  055.000 ! y^-1

      SED_REAL :: kdeathFer        !#  5.0E-1 ! y^-1
      SED_REAL :: kdeathAer        !#  5.0E-1  ! y^-1
      SED_REAL :: kdeathDen        !#  5.0E-1  ! y^-1
      SED_REAL :: kdeathMan        !#  5.0E-1  ! y^-1
      SED_REAL :: kdeathIro        !#  5.0E-1  ! y^-1
      SED_REAL :: kdeathSul        !#  5.0E-1  ! y^-1
      SED_REAL :: kdeathMet        !#  5.0E-1  ! y^-1

      SED_REAL :: kHyd1            !# 0.95  ! y^-1
      SED_REAL :: kHyd2            !# 1. ! y^-1
      SED_REAL :: kHyd3            !# 0.1  ! y^-1
      SED_REAL :: kHyd4            !# 0.01  ! y^-1
      SED_REAL :: kHydN            !# 0.  ! y^-1

      SED_REAL :: domr2dic         !# 0.01              ! From the moon
      SED_REAL :: domr2pomr        !# 0.001             ! From the moon
      SED_REAL :: poml2doml        !# 1.0               ! From the moon

      SED_REAL :: BMax             !# 5.50E-00 ! 1e12 cells / CellWeight 113 g mol^-1
                                               ! * CellBiomass 5e-13 g cell^-1
                                               ! * 1e3 mmol mol^-1 = = 4
      SED_REAL :: BMin             !# 0.050E-00
      SED_REAL :: fuse             !# 0.1
      SED_REAL :: Cellweight       !# 113. ! g mol^-1
      SED_REAL :: Tiny             !# 1.0000000E-08
      SED_REAL :: Temporary_proton !# 123456789

      SED_REAL :: dnra               !# Fraction of NO2 that converts to NH4, the remaining fraction goes to N2

      !-- Stoichiometric constants for OMModel 1 --!
    !--C--!  For the VCW or Canavan stoichiometry, this is like "x"
           SED_REAL :: xlab        !# 140
           SED_REAL :: xref        !# 140
           SED_REAL :: xspecial    !# 140!9.09
    !--N--! For the VCW or Canavan stoichiometry, this is like "y"
           SED_REAL :: ylab        !# 10.04545
           SED_REAL :: yref        !# 10.04545
           SED_REAL :: yspecial    !# 10 !0.76
    !--P--!   For the VCW or Canavan stoichiometry, this is like "z"
           SED_REAL :: zlab        !# 0.5
           SED_REAL :: zref        !# 0.5
           SED_REAL :: zspecial    !# 0.5
    !--XMetal--!   For the VCW or Canavan stoichiometry, this is like "z"
           SED_REAL :: XMetal_lab  !# 0.01
           SED_REAL :: XMetal_ref  !# 0.01
           SED_REAL :: XMetal_special   !# 1

      !-- Stoichiometric constants for OMModel 3 --!
      SED_REAL :: xPOM1            !# 108.
      SED_REAL :: yPOM1            !# 16.
      SED_REAL :: zPOM1            !# 1.
      SED_REAL :: XMetal_POM1      !# 1
      SED_REAL :: xPOM2            !# 108.
      SED_REAL :: yPOM2            !# 16.
      SED_REAL :: zPOM2            !# 1.
      SED_REAL :: XMetal_POM2      !# 1
      SED_REAL :: xPOM3            !# 108.
      SED_REAL :: yPOM3            !# 16.
      SED_REAL :: zPOM3            !# 1.
      SED_REAL :: XMetal_POM3      !# 1
      SED_REAL :: xPOM4            !# 108.
      SED_REAL :: yPOM4            !# 16.
      SED_REAL :: zPOM4            !# 1.
      SED_REAL :: XMetal_POM4      !# 1
      SED_REAL :: xDHyd            !# 6.
      SED_REAL :: yDHyd            !# 0.
      SED_REAL :: zDHyd            !# 0.
      SED_REAL :: XMetal_DHyd      !# 1
      SED_REAL :: xOAc             !# 2.
      SED_REAL :: yOAc             !# 0.
      SED_REAL :: zOAc             !# 0.
      SED_REAL :: xH2              !# 1.
      SED_REAL :: yH2              !# 0.
      SED_REAL :: zH2              !# 0.

      !-- Secondary reactions --!
      SED_REAL :: kNH4OX           !#  5.0e3   ! mM^-1 y^-1 From  WVC 96 k11
      SED_REAL :: knh4ox_to_N2O    !#
      SED_REAL :: kTSNO3           !#  0.0    ! 8.0e1 ! setting it to close to kTSOX
      SED_REAL :: kTSOX            !#  1.6e2  ! mM^-1 y^-1 From  WVC 96 k12
      SED_REAL :: kMnOX            !#  0.0    ! mM^-1 y^-1 From  WVC 96
      SED_REAL :: kMnadsOx         !#  5.0e3  ! mM^-1 y^-1 From  WVC 96 k7
      SED_REAL :: kMnNO3           !#  0.0    ! Turn it off to match VCW 1.0e+2 !mM^-1 y^-1 From Norlem !1e+5! uM^-1 y^-1
      SED_REAL :: kFeOX            !#  1.4e5  ! mM^-1 y^-1 From  WVC 96 k8
      SED_REAL :: kFeadsOX         !#  5.0e4  ! mM^-1 y^-1 From  WVC 96 k9
      SED_REAL :: kFeSOX           !#  2.2E+2 ! mM^-1 y^-1  WVC 96 k15
      SED_REAL :: kFeS2OX          !#  1.0    ! 1.0e-4                   ! Moon? !changed by GK to 1
      SED_REAL :: kMnAge           !#  0.0    ! 3.0e-2                   ! mM^-1 y^-1 From Norlem
      SED_REAL :: kFeAge           !#  0.00   ! From the Moon 0.6        ! y^-1 from Reed ! 3.0e-2 ! mM^-1 y^-1 From Norlem
      SED_REAL :: kFeNO3           !#  0.0 ! Off to match VCW1.0e+1
            !Tuning it to reduce NO3 in the surface layers ! 1.0e+2   ! mM^-1 y^-1 From Norlem !1e+5! uM^-1 y^-1 From Norlem

      SED_REAL :: kFeSMn       !#  0.0 !5.0e-9    ! From the Moon  ! mM^-1 y^-1 From Norlem    !7e+5! uM^-1 y^-1 From Norlem
      SED_REAL :: kFeSFe       !#  0.0 !1.0E-9    ! From the Moon
      SED_REAL :: kFeSppt      !#  1.5e-1   ! mM^-1 y^-1 *density in code. From  WVC 96
      SED_REAL :: kFeSdis      !#  1.0e-3 ! mol g^-1 y^-1 VCW k23
      SED_REAL :: kTSFe        !#  8.0e0  ! mM^-1 y^-1 From  WVC 96 k14
      SED_REAL :: kPyrite      !#  2.3e-4 ! Turn it off to match VCW  0e-4 ! Moon 1.0e2 ! g mmol^-1 y^-1 ~From Katsev et al. 2013
      SED_REAL :: kMnFe        !#  3.0e3   ! mM^-1 y^-1 From  WVC 96  k10
      SED_REAL :: kTSMn        !#  2.0e1   ! mM^-1 y^-1 From  WVC 96 k13
      SED_REAL :: kCH4OX       !#  1.0e7   ! mM^-1 y^-1 From  WVC 96 k16
      SED_REAL :: kCH4SO4      !#  1.0e+1
                        ! mM^-1 y^-1 VCW k17 !From Brigolin 2009 ! 1.0 !mM^-1 y^-1 From Norlem  !1.0E+03! uM^-1 y^-1 From Norlem
      SED_REAL :: kSidppt      !#  4.5e+2 ! mol g^-1 y^-1 *density in code.  5.4e-1        ! mM^-1 y^-1 From  Thullner 2005
      SED_REAL :: kSiddis      !#  2.5e-1 ! y^-1
      SED_REAL :: kRodppt      !#  1e+2 ! mol g^-1 y^-1 VCW k21 *density in code. 6.0e-1 !Tuning 6.0e-2     ! From Thullner 2005
      SED_REAL :: kRoddis          !#  2.5e-1 ! y^-1 VCW k-21
      SED_REAL :: kCalppt          !#  0.0 ! Turn it off to match VCW 1e-4         ! From the moon
      SED_REAL :: kMnO2Appt        !#  0.0 ! Turn it off to match VCW 1e2          ! From the moon
      SED_REAL :: kMnO2Bppt        !#  0.0e-2 ! From the moon
      SED_REAL :: kFeOHAppt        !#  0.0 ! Turn it off to match VCW 0.1          ! From the moon
      SED_REAL :: kFeOHBppt        !#  0.0e-2 ! From the moon
      SED_REAL :: knh4no2
      SED_REAL :: kno2o2

   !------------------------------------------------------------------------
   !##########!           FT setup  v        ################################
   !------------------------------------------------------------------------
     SED_REAL :: FTR               !# 8.31e-3   ! kJ k^-1 mol^-1
     SED_REAL :: FTT               !# 298.0         ! K
     SED_REAL :: deltaGATP         !# 45.0      ! kJ mol^-1
     SED_REAL :: e                 !# 2.71828182845904523536028747135266249775724709369995
     SED_REAL :: F                 !# 96.48534      ! kJ (volt gram equivalent)^-1
     SED_REAL :: n                 !# 1.            ! electron
     SED_REAL :: dPsi              !# 0.120         ! volt
          ! FT energies kJ mol substrate^-1
     SED_REAL :: dG0FerDHyd        !# -30.      ! Dale et al. 2008a. 24*1.25 kJ e-^-1
     SED_REAL :: dG0AerDHy         !# -2883.    ! #Roden and Jin 2011
     SED_REAL :: dG0AerOAc         !# -847.0    ! #Roden and Jin 2011

     SED_REAL :: dG0DenDHy         !# -2774.    ! #Roden and Jin 2011
     SED_REAL :: dG0DenOAc         !# -813.0    ! #Roden and Jin 2011
     SED_REAL :: dG0DenH2          !# -226.0    ! #Roden and Jin 2011
     SED_REAL :: dG0ManOAc         !# -625.0    ! #Roden and Jin 2011
     SED_REAL :: dG0ManH2

     SED_REAL :: dG0IroOAc         !# -736.6   ! #Roden and Jin 2011 |# -61     #kJ mol OAc^-1 Bethke 2011
     SED_REAL :: dG0IroH2          !# -230.7   ! #Roden and Jin 2011 |# -30      #kJ mol H2^-1  Bethke 2011
     SED_REAL :: dG0SulOAc         !# -64.7    ! #Roden and Jin 2011 |# -6*8    #kJ mol OAc^-1 Bethke 2011
     SED_REAL :: dG0SulH2          !# -38.8    ! #Roden and Jin 2011 |# -19*8    #kJ mol H2^-1  Bethke 2011
     SED_REAL :: dG0MetOAc         !# -31.7    ! #Roden and Jin 2011 |# -3.88*8 #kJ mol OAc^-1 Bethke 2011
     SED_REAL :: dG0MetH2          !# -31.3    ! #Roden and Jin 2011 |# -16.92*8 #kJ mol H2^-1  Bethke 2011

     SED_REAL :: YDHyAer   !# 68.41 ! /CellWeight !# Roden and Jin 2011
                                    ! #0.515375 / CellWeight*GlucoseWeight # mol cells (mol glucose)^-1
                                    !# (Not sure where this came from)
      SED_REAL :: YDHyFer     !# 11.95 ! /CellWeight # Roden and Jin 2011 | 0.2 #% Watson 2003
      SED_REAL :: YDenDHy     !# 67.10 ! /CellWeight # Roden and Jin 2011 | #0.333 / CellWeight*GlucoseWeight
                                                                                  !# mol cells (mol glucose)^-1
                                                                                  !# Not sure where this came from
      SED_REAL :: YAerOAc     !# 17.25 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011
      SED_REAL :: YDenOAc     !# 16.82 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011
      SED_REAL :: YDenH2      !# 04.49 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011
      SED_REAL :: YManH2      !#    CAB - Not Used?
      SED_REAL :: YManOAc     !# 14.14 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011
      SED_REAL :: YIroOAc     !# 18.70 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011 |  #9.0/CellWeight #Watson 2003
      SED_REAL :: YIroH2      !# 04.54 ! /CellWeight #gB/mol H2  / g/mol cell Roden and Jin 2011 |  #9.0/CellWeight #Watson 2003
      SED_REAL :: YSulOAc     !# 02.03 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011 |  #4.3/CellWeight #Watson 2003
      SED_REAL :: YSulH2      !# 01.15 ! /CellWeight #gB/mol H2  / g/mol cell Roden and Jin 2011 |  #4.3/CellWeight #Watson 2003
      SED_REAL :: YMetOAc     !# 01.02 ! /CellWeight #gB/mol OAc / g/mol cell Roden and Jin 2011 |  #2.0/CellWeight #Watson 2003
      SED_REAL :: YMetH2      !# 00.94 ! /CellWeight #gB/mol H2  / g/mol cell Roden and Jin 2011 |  #2.0/CellWeight #Watson 2003

   !------------------------------------------------------------------------
   !##########!              FT setup ^           ###########################
   !------------------------------------------------------------------------
      ! FTEA and FIN constants
          ! Approach 1
      SED_REAL :: kO2         !# 20.0E-3  ! WVC96 10.0E-03  ! mM from Megonigal  !  0.04 !mM From Norlem   !40! uM from Norlem
      SED_REAL :: kpO2        !# 20.0E-3  ! WVC96 10.0E-03  ! mM from Megonigal  ! 0.04  !mM From Norlem   !40! uM from Norlem
      SED_REAL :: kin_denitrat ! This was k1o2
      SED_REAL :: kin_denitrit ! This was k1o2
      SED_REAL :: kpart_ammox
      SED_REAL :: kin_deamm
      SED_REAL :: kin_denitrous
      SED_REAL :: klim_denitrat
      !# 3.0E-3   ! WVC96 27.0E-03  ! mM from Megonigal  ! 2.5e-2 !mM From Norlem  !25! uM from Norlem
      SED_REAL :: klim_denitrit        !# From the moon
      SED_REAL :: kpart_denitrit       !# From the moon
      SED_REAL :: klim_denitrous        !# From the moon
      SED_REAL :: kpNO3       !# 3.0E-3   ! WVC96 27.0E-03  ! mM from Megonigal  ! 2.5e-2 !mM From Norlem  !25! uM from Norlem
      SED_REAL :: kpNO2       !#
      SED_REAL :: kpN2O       !#
      SED_REAL :: kMnO2    !# 16.0E-0  ! WVC96 13.3E-03  ! mM from Brigolin et al. 2009
                                                         ! 1.4 !mmol/g From Norlem   !1400! umol/g from Norlem
      SED_REAL :: kpMnO2   !# 16.0E-0  ! WVC96 13.3E-03  ! mM from Brigolin et al. 2009
                                                         ! 1.4 !mmol/g From Norlem   !1400! umol/g from Norlem
      SED_REAL :: kFeOH    !# 100.0E-0 ! WVC96 87.0E-03  ! mM from Brigolin et al. 2009
                                                         ! 18.0 !mmol/g From Norlem  !18000! umol/g from Norlem
      SED_REAL :: kpFeOH   !# 100.0E-0 ! WVC96 87.0E-03  ! mM from Brigolin et al. 2009
                                                         ! 18.0 !mmol/g From Norlem  !18000! umol/g from Norlem
      SED_REAL :: kSO4     !# 1600.0E-3 ! WVC96 1180.0E-03 ! mM from Brigolin et al. 2009
                                                         ! 0.8 !mM From Norlem   !800! uM from Norlem
      SED_REAL :: kpSO4    !# 1600.0E-3 ! WVC96 1180.0E-03 ! mM from Brigolin et al. 2009
                                                         ! 0.8  !mM From Norlem  !800! uM from Norlem

      SED_REAL :: kpO2NO3     !#   3E-3
      SED_REAL :: kpO2NO2
      SED_REAL :: kpO2N2O
      SED_REAL :: kpO2MnO2    !#       200E-0
      SED_REAL :: kpO2FeOH    !#       200E-0
      SED_REAL :: kpO2SO4     !#       2E-3
      SED_REAL :: kpO2CH4     !#       2E-3

          ! Approach 2
      SED_REAL :: lO2         !#  20e-3     ! mM From Thuller 2005
      SED_REAL :: lpO2        !#  20e-3     ! mM From Thuller 2005
      SED_REAL :: lNO3        !#  5e-3      ! mM From Thuller 2005
      SED_REAL :: lpNO3       !#  5e-3      ! mM From Thuller 2005
      SED_REAL :: lNO2        !#  From the moon
      SED_REAL :: lpNO2       !#
      SED_REAL :: lN2O        !#
      SED_REAL :: lpN2O       !#

      SED_REAL :: lMnO2            !#  16e-0     ! umol g sediment^-1
                                   !# = umol cm^-3 sediment - need to normalize by density, porosity

      SED_REAL :: lpMnO2           !#  16e-0     ! umol g sediment^-1
      SED_REAL :: lFeOH            !#  100e-0    ! umol g sediment^-1
      SED_REAL :: lpFeOH           !#  100e-0    ! umol g sediment^-1

      SED_REAL :: lSO4             !#  1600e-3   !mM From Thuller 2005
      SED_REAL :: lpSO4            !#  1600e-3       !mM From Thuller 2005

      SED_REAL :: lpo2no3          !#  20.0e-3   !mmol L^-1 from Thullner 2005
      SED_REAL :: lpo2no2          !#
      SED_REAL :: lpo2n2o          !#
      SED_REAL :: lpo2mno2         !#  20.0e-3   !200e-3 !20.0e-3 ! 200.0e-3     !mmol L^-1 from Thullner 2005
      SED_REAL :: lpo2feoh         !#  20.0e-3   !200.0e-3 !20.0e-3 ! 200.0e-3   !mmol L^-1 from Thullner 2005
      SED_REAL :: lpo2so4          !#  20.0e-3   !2.0e-3 !20.0e-3 ! 2.0e-3   !mmol L^-1 from Thullner 2005
      SED_REAL :: lpo2ch4          !#  20.0e-3   !2.0e-3 !20.0e-3 ! 2.0e-3   !mmol L^-1 from Thullner 2005

      ! FOM
      SED_REAL :: KDHyd            !#  1.0E-00   !From the moon ! 1.0E+00 mM Dale 2008
      SED_REAL :: KOAc             !#  1.0E-2    ! mM thesis compilation table
      SED_REAL :: KH2              !#  1.0E-2    ! mM thesis compilation table
      SED_REAL :: fracOAc          !#  0.666687
      SED_REAL :: fracH2           !#  0.333333

        ! Geochemcial setup
      SED_REAL :: kNH4Ads          !# 1.4       ! 1.4
      SED_REAL :: VCWSb            !# 30        ! micromol g^-1
      SED_REAL :: GammaS           !# 0.11      ! 0.11 to 0.18
      SED_REAL :: KSMnads          !# 3.5
      SED_REAL :: KSFeads          !# 3.7
      SED_REAL :: kPO4Ads          !#  3.0 ! L mg^-1 From Norlem
      INTEGER  :: PO4AdsorptionModel !# 1
      SED_REAL :: KPO4p            !# 1.05 ! From aed_phosphorus
      SED_REAL :: Kadsratio        !# 1.05 ! From aed_phosphorus
      SED_REAL :: Qmax             !# 1.05 ! From aed_phosphorus
      SED_REAL :: KPSid            !# -8.4
      SED_REAL :: KPRod            !# -8.5
      SED_REAL :: KPFeS            !#  2.2
      CHARACTER(len=5) :: Xname    !#  "Zn"
      SED_REAL :: Xmk              !# 0.06
      SED_REAL :: Xfl              !# 0.065
      SED_REAL :: Xfm              !# 0.0001
      SED_REAL :: kXSppt           !# 0.0

      SED_REAL :: kmn
      SED_REAL :: kpmn
      SED_REAL :: kfe
      SED_REAL :: kpfe
      SED_REAL :: lmn              ! Dan added
      SED_REAL :: lpmn             ! Dan added
      SED_REAL :: lfe              ! Dan added
      SED_REAL :: lpfe             ! Dan added
      SED_REAL :: kanh4
      SED_REAL :: kapo4
   END TYPE aed_sed_candi_param_t

  ! ------------------------------------------------------
   TYPE aed_candi_stioc_coefs_t
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
   END TYPE aed_candi_stioc_coefs_t

!#----------------- sed_candi types ----------------------------------
!
  ! maximal value for prescribed corgflux values
  INTEGER,    PARAMETER :: MAXCORG=2000

  ! readspecies=maximal number of included in CCANDI
  !        (for read in STEUER.f/STEUER.DAT)
  INTEGER,    PARAMETER :: readspecies = 26

  TYPE aed_sed_candi_t
     TYPE(aed_sed_candi_time_t)  :: time
     TYPE(aed_sed_candi_bc_t)    :: BC
     TYPE(aed_sed_candi_param_t) :: param
     TYPE(aed_candi_stioc_coefs_t) :: stcoef

!    !-- Extra Module switches
     LOGICAL  :: simCaXCO3    = .FALSE.
     LOGICAL  :: simFeII      = .FALSE.
     LOGICAL  :: simCaCO3C12  = .FALSE.
     LOGICAL  :: simC12       = .FALSE.
     LOGICAL  :: simTracer    = .FALSE.
     LOGICAL  :: simRefOM     = .FALSE.
     LOGICAL  :: simXS        = .FALSE.
     LOGICAL  :: simXO        = .FALSE.
     LOGICAL  :: simAdsp      = .FALSE.

     INTEGER  :: kgmod       = 0          ! OM variables
     INTEGER  :: MnFemod     = 0          ! Mn & Fe asscoiated variables
     INTEGER  :: alkmod      = 0          ! Alkalinity assoc. variables
     INTEGER  :: rmmod       = 0          !
     INTEGER  :: FeIImod     = 0          ! FeII associated variables
     INTEGER  :: Tracermod   = 0          ! Num of "non-reactive" vars
     INTEGER  :: FeSmod      = 0          ! FeS(Pyrite) variable
     INTEGER  :: PO4smod     = 0          ! Solid PO4 variable
     INTEGER  :: CaXCO3mod   = 0          ! Mn & Fe carbonate variables
     INTEGER  :: CaCO3mod    = 0          ! Calcite variables
     INTEGER  :: FeCO3mod    = 0          ! Siderite variables
     INTEGER  :: MnCO3mod    = 0          ! Rhodocrosite variables
     INTEGER  :: C12mod      = 0          ! Variables assoc. with C12
     INTEGER  :: CaCO3C12mod = 0          ! Variables assoc. with C12
     INTEGER  :: xmod        = 0          ! Variables assoc. with X
!
!    !-- Column inidicies
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
     INTEGER  :: N2Oy
     INTEGER  :: NO2y
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
     INTEGER  :: MPBy

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
     INTEGER  :: Feadsy
     INTEGER  :: Mnadsy
     INTEGER  :: FeSY
     INTEGER  :: FeS2Y
     INTEGER  :: PO4sy
     INTEGER  :: NH4sy
     INTEGER  :: XSy
     INTEGER  :: araY
     INTEGER  :: calY
     INTEGER  :: sidY
     INTEGER  :: rodY

     INTEGER  :: FIN_O2Y

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

!    ! Model parameters
!    REAL(SEDP) :: corg(MAXCORG)
     REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: DIFFC
!    ! Reaction rates
     REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: reac
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGC
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGN
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGP
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RGX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROXDHYD
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROXDOMR
!
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FNO3
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FNO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FN2O
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FMnO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FFeOH
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FSO4
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FMet
!
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNH4OX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: rnh4ox_to_N2O
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnOX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnadsOX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeOX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeadsox
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RTSOX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RCH4OX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeSOX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeS2OX
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RXSOx
!
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: rnh4no2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNO2O2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RN2OO2
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
!
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

!# CAB new
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: QSid
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: deltaSid
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: deltadisSid
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: QRod
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: deltaRod
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: deltadisRod
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: QFeS
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: deltaFeS
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: deltadisFeS
!# CAB end

!
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

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNO3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RNO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RN2O
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMnO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RSO4
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RFeOH
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RMet
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RCO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RCH4
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMNO3
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMNO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMN2O
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMMnO2
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMSO4
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMFeOH
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMMet
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMMnii
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: ROMFeii
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: TP
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: TN
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: TOC

  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: k_iron
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: kp_iron
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: k_manganese
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: kp_manganese
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: l_manganese
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: lp_manganese
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: l_iron
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: lp_iron


  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: TerminalOxidation
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: Nrelease
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: NH4release
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: NO3release
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: Prelease
  REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: XMetalrelease

     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathFer
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathAer
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathDen
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathMan
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathIro
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathSul
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathMet
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: RdeathTot
!
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
!
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGFerDHyd
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTFerDHyd
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGAerOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTAerOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGDenOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTDenOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGDenH2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTDenH2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGManOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: dGManH2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTManOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTManH2
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
!
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: Btot
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: Bsubtot
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FBHyd
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FDHyd
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FOAc
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FH2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FBMax

!# CAB new
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_NO3
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_NO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_N2O
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_MnO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_FeOH
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_SO4
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_CH4

     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2NO3
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2NO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2N2O
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2MnO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2FeOH
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2SO4
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FIN_O2CH4

     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_O2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_NO3
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_NO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_N2O
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_MnO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_FeOH
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_SO4
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTEA_CH4

     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_O2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_NO3
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_NO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_N2O
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_MnO2
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_FeOH
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_SO4
     REAL(SEDP), DIMENSION(:),   ALLOCATABLE :: FTem_Met
!# CAB end

     REAL(SEDP) :: dGmp                                            !Dan added


!
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: kg0var
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: kgpo4up
!
     REAL(SEDP) :: sRNHOX,  sRSOX, sRMNOX, sRFEOX, sRMNFE, sRMNO2TS
     REAL(SEDP) :: sRCH4OX, sRCH4SO4
     REAL(SEDP) :: sRNO3TS, sRTSFE3
     REAL(SEDP) :: sRDCO3,  sRDHCO3, sRDCO2
     REAL(SEDP) :: sRFESFE3,sRFESMN4,sRFESOX
     REAL(SEDP) :: sRMnNO3, sRFENO3
     REAL(SEDP) :: sRFe1OX, sRFe2OX, sRFe1NO3, sRFe2NO3
!
     REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: Y
     REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: Ydot
     REAL(SEDP), DIMENSION(:,:), ALLOCATABLE :: Ytemp
     REAL, DIMENSION(:,:), ALLOCATABLE :: IAP
     REAL, DIMENSION(:,:), ALLOCATABLE :: KIAP
     REAL, DIMENSION(:,:), ALLOCATABLE :: QIAP
     REAL, DIMENSION(:,:), ALLOCATABLE :: morevarlist

     !
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
!
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: bottom, top_bound, FixedBottomConc
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dh
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dh2
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dhr
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: dhf
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: RPAR
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: rtol,atol
!
     REAL(SEDP) :: fluxirr, fluxdiff, fluxadv

     SED_REAL :: xldouble, tdh
     SED_REAL :: porosmin
     SED_REAL :: nld, lld
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: PartFluxes
     REAL(SEDP), DIMENSION(:), ALLOCATABLE :: deepconcs

     REAL(SEDP) :: FG0,FG1,FG2,fash,fcaco3

     INTEGER :: npt,oldday,randini
     CHARACTER (LEN=40) :: FileCon, FileRat, FilePar, FileSpe, FileDanPar, FileRat2
     CHARACTER (LEN=40) :: fname(0:readspecies-1), corgfile,FileAuf, FileratOAc

     SED_REAL :: temp, salt, PRES, tout, deltat
     INTEGER  :: iso4, iSTEADY, iPRINT
     INTEGER  :: aufflagout, useVODE, hofmu
     LOGICAL  :: restartCANDI
     INTEGER  :: writestep, corginp, job
     REAL(SEDP) :: SC,SN,SP   ! OM Stoichiometry
     REAL(SEDP) :: MK,FL,FM   ! Mineral Metal Stoichiometry
     REAL(SEDP) :: eqtol,lost,zeitschritt,day
     REAL(SEDP) :: DF(28)
     REAL(SEDP) :: daysFromStart
     INTEGER    :: thisStep
  END TYPE aed_sed_candi_t


!  ===============================
  TYPE AEDConstDiagenesisType
     TYPE(aed_sed_candi_time_t)  :: time
     TYPE(aed_sed_candi_bc_t)    :: bc
     TYPE(aed_sed_candi_param_t) :: param
     TYPE(aed_candi_stioc_coefs_t) :: stcoef
     !???
     REAL                      :: Xmk, Xfl, Xfm
  END TYPE AEDConstDiagenesisType
!-------------------------------------------------------------------------------
#endif

END MODULE aed_gctypes
!------------------------------------------------------------------------------!
