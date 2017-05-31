!###############################################################################
!#                                                                             #
!# aed2_sedeqsolver.F90                                                        #
!#                                                                             #
!# Calculate the updated values for geochemical species, including             #
!#  metals, cations and anions based on geochemical equilibrium                #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
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
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"
#include "aed2_debug.h"

#define DP 8

MODULE aed2_sedeqsolver
!##############################################################################!

  USE aed2_core  ! temporarily
! USE DataTypes
  USE aed2_gctypes

  USE aed2_gclib,only : antiLOG,NumOfComps,NumOfSpecies,GetSpeciesName,         &
                       GetCompCompName, GetCompInfo, GetPhaseInfo, GetSpeciesInfo


 IMPLICIT NONE


 PRIVATE
 PUBLIC :: ConfigEquilibriumSolver,                                            &
           InitialiseGCProperties,                                             &
           UpdateEquilibration,                                                &
           reportGeochemConfig,                                                &
           GetListOfGeochemDiagnostics,                                        &
      !    returnGCDerivedVector,                                              &
           allComponents, nComponents,                                         &
           PUREPHASE,                                                          &
           nFixedDICHMVars,                                                    &
           simManganRedox,                                                     &
           simArsenicRedox,                                                    &
           simSulfurRedox,                                                     &
           simCarbonRedox,                                                     &
           simIronRedox,                                                       &
           FEII, FEIII, MNII, MNIV, ALIII, NiII,                               &
           ASIII, ASV, CDII, ZNII, PBII, SO4, H2S,                             &
           chargeBalCol, CO3, CH4


!------------------------------------------------------------------------------!
! Module level declarations                                                    !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! CAB - my quick haque additions
!
 INTEGER, PARAMETER :: MAX_LABEL_LENGTH    = 32
 INTEGER, PARAMETER :: idl=10    ! Character length of keywords
 REAL, PARAMETER  ::  WaterFreezingInKelvins = 273.15  ! !##get ref
 REAL, PARAMETER  ::  GasConstant = 8.3143             ! Drever, 1997, pg19

 LOGICAL :: simC_DIC      ! Simulate Inorganic Carbon in C Cycle

!------------------------------------------------------------------------------!

 !-- Component Subset required for this simulation
 TYPE (gcUnknowns), DIMENSION(:), ALLOCATABLE, TARGET  :: allComponents
 TYPE (gcSpecies),  DIMENSION(:), ALLOCATABLE, TARGET  :: allSpecies

 !-- Some Shortcut Pointers
 TYPE (gcUnknowns), POINTER                            :: ionStrength
 TYPE (gcUnknowns), POINTER                            :: activityWater
 TYPE (gcUnknowns), POINTER                            :: hydrogenIon
 TYPE (gcUnknowns), POINTER                            :: massWater
 TYPE (gcUnknowns), POINTER                            :: massHydrogen

 !-- Module Constants
 INTEGER, PARAMETER :: nCompulsoryUnknowns = 5      ! MU,AH2O,CB,MH,MH2O
 !-- Component Types
 INTEGER, PARAMETER :: PUREPHASE   = 1
 INTEGER, PARAMETER :: MOLEBLNCE   = 2
 INTEGER, PARAMETER :: IONSTNGTH   = 3
 INTEGER, PARAMETER :: CHARGEBAL   = 4
 INTEGER, PARAMETER :: ACTIVYH2O   = 5
 INTEGER, PARAMETER :: MASSOXYGN   = 6
 INTEGER, PARAMETER :: MASSHYDGN   = 7
 !-- Species Gamma Calculation Method
 INTEGER, PARAMETER :: UNCHGD      = 0
 INTEGER, PARAMETER :: DAVIES      = 1
 INTEGER, PARAMETER :: WATEQDH     = 2
 !-- Input concentration mode
 INTEGER, PARAMETER :: MMOLPERL    = 0
 INTEGER, PARAMETER :: MMOLPERM3   = 2
 INTEGER, PARAMETER :: MGPERL      = 1
 !-- Output Verbosity
 INTEGER, PARAMETER :: baseVerb    = 0
 INTEGER  :: verbosity   = 1
 !-- General
 DOUBLETYPE, PARAMETER :: gc_zero =  0.0_GCHP
 DOUBLETYPE, PARAMETER :: gc_one  =  1.0_GCHP
 DOUBLETYPE, PARAMETER :: minLogActivity = -30.0

 !-- Module Variables
 DOUBLETYPE, DIMENSION(:,:), ALLOCATABLE :: inequalityArray
 DOUBLETYPE :: Waq         = gc_one
 DOUBLETYPE :: molWgtH2O   = 0.0180160000
 DOUBLETYPE :: cellTemp    = 20.00
 DOUBLETYPE :: cellSal     = 0.00
 INTEGER  :: residColumn     = 0
 INTEGER  :: chargeBalCol    = 0
 INTEGER  :: nComponents     = 0
 INTEGER  :: nSpecies        = 0
 INTEGER  :: nOPTEqs         = 0
 INTEGER  :: nEQLEqs         = 0
 INTEGER  :: nINQEqs         = 0
 INTEGER  :: nMBEqs          = 0
 INTEGER  :: nPPEqs          = 0
 INTEGER  :: nFixedDICHMVars = 0
 LOGICAL :: pHisFixed       = .TRUE.
 LOGICAL :: peisFixed       = .FALSE.
 LOGICAL :: SolveWithPP     = .FALSE.

 LOGICAL :: simManganRedox  = .FALSE.
 LOGICAL :: simArsenicRedox  = .FALSE.
 LOGICAL :: simIronRedox    = .FALSE.
 LOGICAL :: simSulfurRedox  = .FALSE.
 LOGICAL :: simCarbonRedox  = .FALSE.
 INTEGER  :: FEII, FEIII, MNII, MNVII, ALIII, MNIV, NiII
 INTEGER  :: ZNII, ASIII, ASV, PBII, CDII, H2S, SO4, CH4, CO3

 LOGICAL :: solvePPtoEquilibrium

 !-- Derived Variables
 INTEGER  :: nGCDerivedVars = 0
 INTEGER, DIMENSION(:),ALLOCATABLE        :: derivedGCList
 DOUBLETYPE, DIMENSION(:), ALLOCATABLE    :: derivedGCVals
!##############################################################################!


CONTAINS


!##############################################################################!
! ConfigureGeochem:                                                            !
!                                                                              !
! Main setup routine for geochemistry solver.                                  !
!------------------------------------------------------------------------------!
 SUBROUTINE ConfigEquilibriumSolver( nUserReqComponents,  nUserReqPurePhases,  &
                              listofCompNames,     listofPPNames,              &
                              nDissTransportables, nPartTransportables,        &
                              listDissTransVars,   listPartTransVars)          !
   !-- Incoming                                                                !
   INTEGER,                         INTENT(IN)    :: nUserReqComponents        !
   INTEGER,                         INTENT(IN)    :: nUserReqPurePhases        !
   CHARACTER(LEN=64), DIMENSION(:), INTENT(IN)    :: listofCompNames           !
   CHARACTER(LEN=64), DIMENSION(:), INTENT(IN)    :: listofPPNames             !
   !-- Outgoing                                                                !
   INTEGER,                         INTENT(INOUT) :: nDissTransportables       !
   INTEGER,                         INTENT(INOUT) :: nPartTransportables       !
   CHARACTER(LEN=64), ALLOCATABLE, DIMENSION(:), INTENT(OUT)   :: listDissTransVars       !
   CHARACTER(LEN=64), ALLOCATABLE, DIMENSION(:), INTENT(OUT)   :: listPartTransVars       !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "ConfigEquilibriumSolver"            !
   TYPE(gcUnknowns)        :: dummyComp                                        !
   INTEGER  :: status                                           !
   INTEGER  :: unknownsThatAreAlreadyVars                       !
   INTEGER  :: i, cIndex, wqColCounter                          !
                                                                               !
   !---------------------------------------------------------------------------!

   !-- Lets begin the setup
   nFixedDICHMVars = 0 !nDissTransportables

   !-------
   !-- Calculate the total number of geochemical components
   !-- This is not the size of DICHM - it is the number of
   !-- Components to be solved for in this module,
   !-- i.e. nComponents
   !-------
   nComponents = nUserReqComponents &
               + nUserReqPurePhases &
               + nCompulsoryUnknowns


   !-- Allocate space component objects
   ALLOCATE(allComponents(nComponents),STAT=status)
   CheckAllocStatus(status,THIS_PROC,"allComponents")


   !-------
   !-- In the WQ array, DICHM MUST have the *essential* inorganic
   !-- nutrients: DO, PO4, NH4, NO3 and SiO2
   !-- If these are also to be included in the speciation calcs,
   !-- then we must make sure we dont count them twice
   !-------
   unknownsThatAreAlreadyVars = 0
   !DO i = 1,nUserReqComponents
   !  IF(TRIM(listofCompNames(i)) == "DO"  .OR. &
   !     TRIM(listofCompNames(i)) == "PO4" .OR. &
   !     TRIM(listofCompNames(i)) == "NH4" .OR. &
   !     TRIM(listofCompNames(i)) == "NO3" .OR. &
   !     TRIM(listofCompNames(i)) == "SiO2") THEN
   !
   !     unknownsThatAreAlreadyVars = unknownsThatAreAlreadyVars + 1
   !  END IF
   !END DO

   !-------
   !-- Check with the database to see the components requested
   !-- by the user (in listofCompNames read in from *.con) are
   !-- valid and present in the chem database (*.chm). If they
   !-- are present, then populate the local Component space.
   !-------
   wqColCounter = 0! nDissTransportables
   DO i = 1,nUserReqComponents
     allComponents(i)%EltName    = TRIM(listofCompNames(i))

     CALL GetCompInfo(allComponents(i))   !!##GEOCH_MOD

     allComponents(i)%wqIndex = GetDissChemIndex(allComponents(i)%EltName,       &
                                wqColCounter)
     allComponents(i)%CompIndex = i
     allComponents(i)%CompType  = MOLEBLNCE

     ! DIC is a special case: If this has been selected
     ! then we must perform some other
     IF(TRIM(allComponents(i)%EltName) == "DIC") THEN
       !DICColNum = allComponents(i)%wqIndex
       simC_DIC  = .TRUE.
     END IF
   END DO

   !-------
   !-- The geochem routine requires some compulsory unknowns
   !-- to be included - these are added after the user requested
   !-- unknowns.
   !-------

   !-- Compulsory Unknown 1: Ionic Strength
   allComponents(nUserReqComponents+1)%CompType   =  IONSTNGTH
   allComponents(nUserReqComponents+1)%CompName   = "IONSTNGTH"
   allComponents(nUserReqComponents+1)%EltName    = "IONSTNGTH"
   allComponents(nUserReqComponents+1)%CompIndex  = nUserReqComponents+1
   allComponents(nUserReqComponents+1)%wqIndex    = 0 ! NOT TRANSPORTABLE
   ionStrength => allComponents(nUserReqComponents+1)

   !-- Compulsory Unknown 2: Activity of Water
   allComponents(nUserReqComponents+2)%CompType   =  ACTIVYH2O
   allComponents(nUserReqComponents+2)%CompName   = "H2O"
   allComponents(nUserReqComponents+2)%EltName    = "ACTIVYH2O"
   allComponents(nUserReqComponents+2)%CompIndex  = nUserReqComponents+2
   allComponents(nUserReqComponents+2)%wqIndex    = 0 ! NOT TRANSPORTABLE
   activityWater => allComponents(nUserReqComponents+2)

   !-- Compulsory Unknown 3: Charge Balance (lnaH+)
   allComponents(nUserReqComponents+3)%CompType   =  CHARGEBAL
   allComponents(nUserReqComponents+3)%CompName   = "H+"
   allComponents(nUserReqComponents+3)%EltName    = "pH"
   allComponents(nUserReqComponents+3)%CompIndex  = nUserReqComponents+3
   allComponents(nUserReqComponents+3)%wqIndex    = &
            GetDissChemIndex(allComponents(nUserReqComponents+3)%EltName, &
                           wqColCounter)
   hydrogenIon   => allComponents(nUserReqComponents+3)


   !-- Compulsory Unknown 4: Mass Hydrogen (lnae-)
   allComponents(nUserReqComponents+4)%CompType   =  MASSHYDGN
   allComponents(nUserReqComponents+4)%CompName   = "e-"
   allComponents(nUserReqComponents+4)%EltName    = "pe"
   allComponents(nUserReqComponents+4)%CompIndex = nUserReqComponents+4
   allComponents(nUserReqComponents+4)%wqIndex = &
            GetDissChemIndex(allComponents(nUserReqComponents+4)%EltName, &
                           wqColCounter)
   massHydrogen  => allComponents(nUserReqComponents+4)

   !-- Compulsory Unknown 5: Mass Water
   allComponents(nUserReqComponents+5)%CompType   = MASSOXYGN
   allComponents(nUserReqComponents+5)%CompName   = "MASSOXYGN"
   allComponents(nUserReqComponents+5)%EltName    = "MASSOXYGN"
   allComponents(nUserReqComponents+5)%CompIndex  = nUserReqComponents+5
   allComponents(nUserReqComponents+5)%wqIndex    = 0 ! NOT TRANSPORTABLE
   massWater     => allComponents(nUserReqComponents+5)

   !-- Special transportable so we can remember the conc of unbalaned charge
   chargeBalCol = GetDissChemIndex("ubalchg   ",wqColCounter)


   !-------
   !-- OK, now we can calculate the size of DICHM - the dissolved
   !-- transported components. This is not the same as
   !-- nComponents! Now that we know what is being simulated
   !-- we can return to WQ_config() the list of names to be
   !-- included in WQ3D%a3d
   !-------

   nDissTransportables = wqColCounter
   !-- Allocate and populate list of DICHM vars
   ALLOCATE(listDissTransVars(nDissTransportables),STAT=status)
   listDissTransVars = ""
   listDissTransVars = GetDissChemNames(nDissTransportables, &
              allComponents(1:nUserReqComponents+nCompulsoryUnknowns))


   DO i = 1,nDissTransportables
     PRINT *,'  Diss Names:',i,listDissTransVars(i)
   END DO


   !-------
   !-- Now we know what the components to be included in the simulation
   !-- are, we can sort through the species database and find the
   !-- ones that should be simulated. The simulated subset is
   !-- placed in allSpecies(:), and all relevant pointers are setup
   !-------


   CALL createSpeciesSubSetforSim( &
           allComponents(1:nUserReqComponents+nCompulsoryUnknowns))


   !-------
   !-- The pure phase unknowns must now be taken care of
   !-- and PICHM setup (similar to above). There are currently
   !-- no compulsory PICHM vars, so this is easy.
   !-------

   ALLOCATE(dummyComp%PPData,STAT=status)
   CheckAllocStatus(status,THIS_PROC,"dummyComp%PPData")

   ALLOCATE(dummyComp%PPData%stoich(NumOfComps()),STAT=status)
   CheckAllocStatus(status,THIS_PROC,"dummyComp%PPData%Stoich")
   dummyComp%PPData%stoich = gc_zero

   wqColCounter = 0            !nPartTransportables
   DO i = 1,nUserReqPurePhases
     cIndex = nUserReqComponents + nCompulsoryUnknowns + i

     dummyComp%EltName = TRIM(listOfPPNames(i))

     CALL GetPhaseInfo(dummyComp)

     dummyComp%CompIndex = cIndex
     allComponents(cIndex) = setupPhaseStoich(dummyComp,     &
                 allComponents(1:nUserReqComponents+nCompulsoryUnknowns))

     allComponents(cIndex)%wqIndex = GetPartChemIndex(wqColCounter)
     allComponents(cIndex)%CompType  = PUREPHASE

   END DO

   !-------
   !-- Allocate and populate list of PICHM vars
   !-------
   nPartTransportables = wqColCounter

   ALLOCATE(listPartTransVars(nPartTransportables),STAT=status)
   listPartTransVars = ""
   listPartTransVars = GetPPChemNames(nPartTransportables, &
         allComponents( (nUserReqComponents+nCompulsoryUnknowns+1) :    &
                         nComponents) )

   DO i = 1,nPartTransportables
     PRINT *,'  Part Names:',i,listPartTransVars(i)
   END DO

   !solvePPtoEquilibrium = .TRUE.
   solvePPtoEquilibrium = .FALSE.

   !-------
   !-- Update species stoich array with relevant pure phase
   !-- data.
   !-------


   !-------
   !-- Reset special CompName's for internal use
   !-------
   allComponents(nUserReqComponents+2)%CompName   = "ACTIVYH2O"
   allComponents(nUserReqComponents+3)%CompName   = "CHARGEBAL"
   allComponents(nUserReqComponents+4)%CompName   = "MASSHYDGN"


   !-------
   !-- Now define the dimensions of the reaction set
   !-------
   nMBEqs  = nUserReqComponents
   nPPEqs  = nUserReqPurePhases


   Waq       = gc_one
   molWgtH2O = 0.0180160000

   IF(verbosity > 0) THEN
     print *,'  Aqueous Species List: ----------'
!     DO i = 1,nSpecies
!       PRINT *,'   Species:',i, TRIM(allSpecies(i)%Name),': ',                 &
!                                INT(allSpecies(i)%Stoich)
!     END DO
     print *,'--------------------------------'
     IF(verbosity > 2) THEN
       CALL outputSpeciesData(allSpecies)
     END IF
   END IF

   !-------
   !-- Now find out if we have any redox pairs
   !-------
   CALL setupWCRedoxConfiguration( allComponents(1:nUserReqComponents) )


 Call ReportGeochemConfig(6)

 END SUBROUTINE ConfigEquilibriumSolver                                               !
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
! InitialiseGCProperties:                                                      !
!                                                                              !
! Based on user input initial conditions for geochemical variables, this       !
! rotuine will solve for equilibrium conditions. It assumes pH is fixed        !
! (ie. given based on measurement), and solves for the charge imbalance, which !
! is stored in the transportable variable CHGBAL.                              !
! No mineral phases are included.                                              !
!------------------------------------------------------------------------------!
 SUBROUTINE InitialiseGCProperties(dissConcs, partConcs, concMode)             !
   !-- Incoming                                                                !
 !  AED_REAL,    INTENT(IN)    :: inTemp           ! Temperature   !
   INTEGER, INTENT(IN)                     :: concMode         ! Conc units    !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER                 :: THIS_PROC = "InitialiseGCProperties"
   INTEGER,      DIMENSION(:), ALLOCATABLE :: componentList                    !
!  DOUBLETYPE,   DIMENSION(:)  :: partConcs, dissConcs             !
   AED_REAL,   DIMENSION(:)  :: partConcs, dissConcs             !
   DOUBLETYPE                              :: fff, weight                      !
   LOGICAL                                 :: REPEAT, FAIL,failed              !
   INTEGER  :: cellIndex,                       &
                                              componentIndex,                  &
                                              nComps2Process,                  &
                                              iter                             !
                                                                               !
   !---------------------------------------------------------------------------!

   !WRITE(*,"(8X,'Initialise geochemical properties...')")

    ! verbosity = 0

     cellTemp  = 20.


     !---------------
     !-- This initial round requires that we hold pH const and
     !-- ignore pure phase so that we can deterimine the charge
     !-- balance and remember it for future calculations
     pHisFixed = .TRUE.
     peisFixed = .TRUE.
     SolveWithPP = .FALSE.

     nComps2Process = nMBEqs + 2     ! Mole-Balance Eqs + ionicStrength + aH2O
     ALLOCATE(componentList(nComps2Process))


     ! Need to set component Eq Indicies
     CALL seteqIndiciesforSpeciation(allComponents, nComps2Process,componentList)
     nOPTEqs = 1                            ! Number of Rows in A matrix
     nEQLEqs = nComps2Process -  nOPTEqs    ! Number of Rows in C matrix
     nINQEqs = 0                            ! Number of Rows in E matrix


     !-- 1. Set start/total moles of components to be the WQ3F value and
     !--    change units as necessary
     CALL setComponentTotalConc(dissConcs,partConcs,allComponents, concMode)

     ionStrength%Total = initialIonicStrength(allComponents)

     CALL updateLogK(allSpecies,allComponents,cellTemp)

     CALL updateGammas(allSpecies,ionStrength%Total)

     REPEAT = .TRUE.
     FAIL   = .FALSE.

     iter = 0
     DO WHILE (REPEAT)
       iter = iter + 1

       IF (verbosity > 10) THEN
         PRINT *,'Beginning set iteration: ',iter
       END IF
       IF (iter == 10) THEN
         PRINT *, 'Did not converge in set iteration',iter
         fail = .TRUE.
       END IF


       CALL updateSpeciesMoles(allSpecies, allComponents)


       CALL moleBalanceSums(allSpecies, allComponents)

       REPEAT = .FALSE.
       DO componentIndex = 1, nComps2Process
         IF(allComponents(componentList(componentIndex))%CompType == MOLEBLNCE) THEN

                IF (ABS(allComponents(componentList(componentIndex))%Total) < 1e-30) THEN
                  allComponents(componentList(componentIndex))%Total = gc_zero
                END IF

                fff = ABS(allComponents(componentList(componentIndex))%Value)

          IF (fff == gc_zero .AND. allComponents(componentList(componentIndex))%Total == gc_zero) THEN
            allComponents(componentList(componentIndex))%Master%logActivity = minLogActivity
            CYCLE

          ELSE IF (fff == gc_zero) THEN
            REPEAT = .TRUE.

            allComponents(componentList(componentIndex))%Master%logActivity =  &
            allComponents(componentList(componentIndex))%Master%logActivity +5.0

          ELSE IF (fail .AND. fff < 1.5 * ABS(allComponents(componentList(componentIndex))%Total)) THEN
            CYCLE
          ELSE IF (fff > 1.5 * ABS(allComponents(componentList(componentIndex))%Total) .OR. &
                    fff < 1e-5 * ABS(allComponents(componentList(componentIndex))%Total) ) THEN

            IF (fff < (1e-5 * ABS(allComponents(componentList(componentIndex))%Total))) THEN
              weight = 0.3
            ELSE
              weight = 1.0
            END IF

            IF (allComponents(componentList(componentIndex))%Total <= gc_zero) THEN
              allComponents(componentList(componentIndex))%Master%logActivity = minLogActivity
            ELSE
              REPEAT = .TRUE.
              allComponents(componentList(componentIndex))%Master%logActivity = &
              allComponents(componentList(componentIndex))%Master%logActivity + &
              weight * LOG10(ABS(allComponents(componentList(componentIndex))%Total / &
              allComponents(componentList(componentIndex))%Value))
            END IF

          END IF
         END IF
       END DO
     END DO

     ionStrength%Total = updateIonicStrength(allSpecies)

     CALL updateGammas(allSpecies,ionStrength%Total)


     !-- 2. Equilibrate water column cells
     CALL solveEquilibriumReactions( allSpecies,   allComponents,              &
                                     componentList,nComps2Process, failed )

     IF(verbosity > 2) THEN
       DO componentIndex = 1,nComps2Process
         print *,'After Aq Stage',allcomponents(componentIndex)%CompName,      &
         allcomponents(componentIndex)%Total,allcomponents(componentIndex)%Value
       END DO
     END IF

     DEALLOCATE(componentList)




     !---------------
     !-- 3. Update state variable arrays with new component estimates
     !      and revert to original units
     CALL updateWQArray(dissConcs,partConcs,allComponents, concMode)

     !WQ(cellIndex,DICHM(:)) = REAL(dissConcs, r_wq)
     !WQ(cellIndex,PICHM(:)) = REAL(partConcs, r_wq)

     !-- Charge Balance Condition
     dissConcs(chargeBalCol) = hydrogenIon%Value


     !-- 4. Now update derived variable arrays
     IF(concMode /= MMOLPERL) THEN
       !-- Not CANDI calling, so update away
       DO iter = 1, nGCDerivedVars
         IF(derivedGCList(iter) > 0) THEN
           derivedGCVals(iter) = allSpecies(derivedGCList(iter))%Moles
         ELSE IF (derivedGCList(iter) == -999) THEN
           derivedGCVals(iter) = calcpCO2(allSpecies, cellTemp, cellSal)
         END IF
       END DO
     END IF


   verbosity = baseVerb


 END SUBROUTINE InitialiseGCProperties
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
! UpdateEquilibration:                                                        !
!                                                                              !
! Main run-time routine to conduct equilibration on supplied cells             !
! Incoming cells may be from water column or sediment, this rotuine is generic !
!                                                                              !
!------------------------------------------------------------------------------!
 SUBROUTINE UpdateEquilibration(dissConcs, partConcs, concMode, inTemp,       &
 stoEq, IAP, KIAP, QIAP, upDerv)
!  DOUBLETYPE,   DIMENSION(:)  :: partConcs, dissConcs             !
   AED_REAL,   DIMENSION(:)  :: partConcs, dissConcs             !
   INTEGER,                  INTENT(IN)              :: concMode  ! g/m3 or mol/L      !
   REAL,    INTENT(IN),    OPTIONAL :: inTemp    ! Temperature        !
   LOGICAL,                  INTENT(IN),    OPTIONAL :: stoEq     ! solvePPtoEquil     !
   AED_REAL, DIMENSION(:), INTENT(INOUT), OPTIONAL :: IAP
   AED_REAL, DIMENSION(:), INTENT(INOUT), OPTIONAL :: KIAP !
   AED_REAL, DIMENSION(:), INTENT(INOUT), OPTIONAL :: QIAP !
   LOGICAL,                  INTENT(IN),    OPTIONAL :: upDerv    ! Update derived     !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER                 :: THIS_PROC = "Equilibrate"        !
   INTEGER,      DIMENSION(:), ALLOCATABLE :: componentList                    !
   DOUBLETYPE                              :: fff, weight, TotalH, TotalO      !
   LOGICAL                                 :: REPEAT, FAIL, failed ,updateDerivedVars            !
   INTEGER  :: cellIndex,                       &
                                              componentIndex,                  &
                                              nComps2Process,                  &
                                              sIndex,                          &
                                              iter,                            &
                                              nFailedCells                     !
                                                                               !
   !---------------------------------------------------------------------------!
   !  WRITE(*,"(6X,'Chemistry equilibration')")

   IF(PRESENT(stoEq)) THEN
     solvePPtoEquilibrium = stoEq
   END IF
   IF(PRESENT(upDerv)) THEN
     updateDerivedVars    = upDerv
   ELSE
     updateDerivedVars    = .FALSE.
   END IF


     nFailedCells = 0

     !verbosity = baseVerb

     cellSal   = 0.0
     IF(PRESENT(inTemp)) THEN
       cellTemp  = REAL(inTemp,GCHP)
     ELSE
       cellTemp  = 20.0
     END IF

     !---------------
     !-- We have already held pH const and calc'd the charge bal
     !-- in initialiseGeochem. So now we must estimate the pH whilst
     !-- keeping the charge balance const. Lets do this without
     !-- pure phases first.

     pHisFixed = .FALSE.
     peisFixed = .TRUE.
     SolveWithPP = .FALSE.

     !-- Before pure phases come in, we just include the mole balances
     !-- ionic strength, activity of H2O and pH
     nComps2Process = nMBEqs + 3
     ALLOCATE(componentList(nComps2Process))


     ! Need to set component Eq Indicies
     CALL seteqIndiciesforSpeciation(allComponents,nComps2Process,componentList)
     nOPTEqs = 1                                    ! Number of Rows in A matrix
     nEQLEqs = nComps2Process -  nOPTEqs            ! Number of Rows in C matrix
     nINQEqs = 0                                    ! Number of Rows in E matrix


     !-- 1. Set start/total moles of components to be the WQ3F value and
     !--    change units as necessary
     CALL setComponentTotalConc(dissConcs,partConcs,allComponents,concMode)

     ionStrength%Total = initialIonicStrength(allComponents)
     hydrogenIon%Total = dissConcs(chargeBalCol)

     CALL updateLogK(allSpecies,allComponents,cellTemp)

     CALL updateGammas(allSpecies,ionStrength%Total)

     REPEAT = .TRUE.
     FAIL   = .FALSE.

     iter = 0
     DO WHILE (REPEAT)
       iter = iter + 1

       IF (verbosity > 10) THEN
         print *,'Beginning set iteration: ',iter
       END IF
       IF (iter == 10) THEN
         print *, 'Did not converge in set iteration',iter
         fail = .TRUE.
       END IF

       CALL updateSpeciesMoles(allSpecies, allComponents)

       CALL moleBalanceSums(allSpecies, allComponents)

       REPEAT = .FALSE.
       DO componentIndex = 1, nComps2Process
         IF(allComponents(componentList(componentIndex))%CompType==MOLEBLNCE) THEN

           IF(ABS(allComponents(componentList(componentIndex))%Total) <1e-30) THEN
             allComponents(componentList(componentIndex))%Total = gc_zero
           END IF

           fff = ABS(allComponents(componentList(componentIndex))%Value)
           IF (fff == gc_zero .AND.                                            &
               allComponents(componentList(componentIndex))%Total==gc_zero) THEN

         allComponents(componentList(componentIndex))%Master%logActivity   &
                  = minLogActivity
             CYCLE

           ELSE IF (fff == gc_zero) THEN

         REPEAT = .TRUE.
             allComponents(componentList(componentIndex))%Master%logActivity = &
             allComponents(componentList(componentIndex))%Master%logActivity+5.0

           ELSE IF (fail .AND. fff < 1.5 * ABS(allComponents(componentList(componentIndex))%Total)) THEN

             CYCLE

           ELSE IF (fff > 1.5 * ABS(allComponents(componentList(componentIndex))%Total) .OR. &

         fff < 1e-5 * ABS(allComponents(componentList(componentIndex))%Total) ) THEN
             IF (fff < (1e-5 * ABS(allComponents(componentList(componentIndex))%Total))) THEN
               weight = 0.3
             ELSE
               weight = 1.0
             END IF

             IF (allComponents(componentList(componentIndex))%Total <= gc_zero) THEN
               allComponents(componentList(componentIndex))%Master%logActivity = minLogActivity
             ELSE
               REPEAT = .TRUE.
               allComponents(componentList(componentIndex))%Master%logActivity = &
               allComponents(componentList(componentIndex))%Master%logActivity + &
               weight * LOG10(ABS(allComponents(componentList(componentIndex))%Total / &
               allComponents(componentList(componentIndex))%Value))
             END IF
           END IF
         END IF
       END DO

     END DO

     ionStrength%Total = updateIonicStrength(allSpecies)

     CALL updateGammas(allSpecies,ionStrength%Total)

     !-- 2. Equilibrate cells/layers
     CALL solveEquilibriumReactions(allSpecies,allComponents,componentList,    &
                                   nComps2Process,failed)

     IF(verbosity > 2) THEN
       DO componentIndex = 1,nComps2Process
        print *,'   After Aq Stage',allcomponents(componentIndex)%CompName,    &
           allcomponents(componentIndex)%Total,allcomponents(componentIndex)%Value
       END DO
     END IF

     DEALLOCATE(componentList)




     IF(nPPEqs > 0 .AND. solvePPtoEquilibrium) THEN

       IF(verbosity > 1) THEN
         PRINT *,'Starting Batch Rxn:'
       END IF

       !---------------
       ! -- Add pure phases and perform full solution
       ! --

       !-- Update switches
       pHisFixed   = .FALSE.
       peisFixed   = .FALSE.
       SolveWithPP   = .TRUE.

       !-- Count unknown equations for inequality solver
       nOPTEqs = nPPEqs                             ! Number of Rows in A matrix
       nEQLEqs = nMBEqs + nCompulsoryUnknowns       ! Number of Rows in C matrix
       nINQEqs = nPPEqs                             ! Number of Rows in E matrix

       ! n = PP + MB's + MU + aH2O + CB + MH + MH2O
       nComps2Process = nPPEqs + nMBEqs + nCompulsoryUnknowns

        !-- Update component Eqn Indicies and populate componentList
       ALLOCATE(componentList(nComps2Process))
       CALL seteqIndiciesforBatchRxn(allComponents,nComps2Process,componentList)

       activityWater%Total = (antiLOG(activityWater%master%logActivity,10.0_GCHP) - 1.0)

       massWater%Total     = massWater%Value
       hydrogenIon%Total   = hydrogenIon%Value
       molWgtH2O           = 1.0/massWater%Total
       ionStrength%Total   = updateIonicStrength(allSpecies)

       TotalH = gc_zero
       TotalO = gc_zero
       DO sIndex = 1,nSpecies
          TotalH = TotalH+ allSpecies(sIndex)%Hstoich  *allSpecies(sIndex)%Moles
          TotalO = TotalO+ allSpecies(sIndex)%H2Ostoich*allSpecies(sIndex)%Moles
       END DO
       TotalH = TotalH + 2*massWater%Total

       massHydrogen%Total  = TotalH-2*TotalO

       !-- 2. Equilibrate water column cells
       CALL solveEquilibriumReactions(allSpecies,allComponents,componentList,  &
                                      nComps2Process,failed)

       DEALLOCATE(componentList)
     END IF


     !---------------
     !-- 3. Update state variable arrays with new component estimates
     ! and revert to original units
     IF(.NOT. failed) THEN

       CALL updateWQArray(dissConcs,partConcs,allComponents,concMode)

       WHERE(partConcs<gc_zero)partConcs=gc_zero

      ! WQ3F(cellIndex,DICHM(:)) = REAL(dissConcs, r_wq)
      ! WQ3F(cellIndex,PICHM(:)) = REAL(partConcs, r_wq)

       !-- 4. Now update derived variable arrays
       IF(updateDerivedVars) THEN
        DO iter = 1, nGCDerivedVars
         IF(derivedGCList(iter) > 0) THEN
           derivedGCVals(iter) = allSpecies(derivedGCList(iter))%Moles
         ELSE IF (derivedGCList(iter) <= -1 .AND. derivedGCList(iter) >= -5  ) THEN
           derivedGCVals(iter) = calcSIforPP(-derivedGCList(iter))
         ELSE IF (derivedGCList(iter) == -998) THEN
           derivedGCVals(iter) = gc_zero
         ELSE IF (derivedGCList(iter) == -999) THEN
           derivedGCVals(iter) = calcpCO2(allSpecies, cellTemp, cellSal)
         END IF
        END DO
       END IF
     ELSE
       ! Non convergances
       PRINT *,'Convergence failure in cell: '

       nFailedCells = nFailedCells + 1
       IF(updateDerivedVars) THEN
        DO iter = 1, nGCDerivedVars
         IF (derivedGCList(iter) == -998) THEN
           derivedGCVals(iter) = 1.0
         END IF
        END DO
       END IF
     END IF

     !-- 5. If pure phases are simulated an optional array
     !   to return the IAP to mediate the kinetic reactions is possible
     IF(nPPEqs > 0 .AND. PRESENT(IAP)) THEN
         DO iter = 1,nPPEqs
           IAP(iter) = CalcIAPforPP( nComponents - nPPEqs + iter )
           ! Note: this assumes PP's are the last in the component array
         END DO
     END IF

     IF(nPPEqs > 0 .AND. PRESENT(KIAP)) THEN
        DO iter = 1,nPPEqs
          KIAP(iter) = CalcIAPforPP2( nComponents - nPPEqs + iter )
          ! Note: this assumes PPs are the last in the component array
        END DO
     END IF

     IF(nPPEqs > 0 .AND. PRESENT(QIAP)) THEN
        DO iter = 1,nPPEqs
          QIAP(iter) = CalcIAPforPP3( nComponents - nPPEqs + iter )
          ! Note: this assumes PPs are the last in the component array
        END DO
     END IF



!   IF(nFailedCells > 0) THEN
!     PRINT *,' >---'
!     PRINT *,' Number of cells that failed converging: ',nFailedCells,' (of ',SIZE(WQ3F,1),')'
!     PRINT *,'  Cell molalities are being set to the previous time step.'
!     PRINT *,' ---<'
!     nFailedCells = 0
!   END IF


 END SUBROUTINE UpdateEquilibration                                           !
!------------------------------------------------------------------------------!







!------------------------------------------------------------------------------!
! SolveEquilibriumReactions                                                    !
!                                                                              !
! Numerical solution to reaction set defined in the above routines. Solution   !
! is based on a iterative Newton-Rhapson scheme and the Simplex matrix solver  !
!                                                                              !
!------------------------------------------------------------------------------!
 SUBROUTINE solveEquilibriumReactions(species, components,                     &
                                      cList,   nComps2Process, failed)         !
   !-- Incoming                                                                !
   TYPE (gcSpecies),   DIMENSION(:)  :: species                                !
   TYPE (gcUnknowns),  DIMENSION(:)  :: components                             !
   INTEGER, INTENT(IN)               :: nComps2Process                         !
   INTEGER, DIMENSION(:), INTENT(IN) :: cList                                  !
   LOGICAL, INTENT(OUT)              :: failed                                 !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER           :: THIS_PROC = "solveEquilibriumReactions"!
   INTEGER  :: Status     ! memory allocation status  !
   INTEGER  :: componentIndex, miscIndex              !
   INTEGER  :: nItern, ineqCntr                       !
   INTEGER  :: inKODE, outKODE, KODE                  !
   INTEGER  :: inITER, outITER, ITER                  !
   INTEGER  :: ITMAX = 150                            !
   LOGICAL                           :: converged                              !
   DOUBLETYPE                        :: ERROR                                  !
   DOUBLETYPE                        :: max, min                               !
   DOUBLETYPE,   PARAMETER           :: convTolerance = 1e-12   ! Convg toler  !
   DOUBLETYPE,   PARAMETER           :: min_value     = 1e-10   ! For Scaling  !
   DOUBLETYPE,   PARAMETER           :: TOLER         = 1e-14   ! For Simplex  !
   DOUBLETYPE,   DIMENSION(nOPTEqs+nEQLEqs+nINQEqs)  :: RES                    !
   DOUBLETYPE,   DIMENSION(nComps2Process+2)         :: deltaConc              !
   DOUBLETYPE,   DIMENSION(nComps2Process)           :: normal                 !
   !---------------------------------------------------------------------------!
   !-- Initial items
   IF(nPPEqs > 1) THEN
     ITMAX = 500
   ELSE
     ITMAX = 150
   END IF

   !-- Allocate space for the inequality solver array
   !-- (This includes jacobian sums + residuals + inequality constraints)
   ALLOCATE(inequalityArray(nOPTEqs+nEQLeqs+nINQEqs+2,nComps2Process+2),       &
                                                      STAT=Status)
   CheckAllocStatus(Status,THIS_PROC,"inequalityArray")
   inequalityArray = gc_zero

   residColumn = nComps2Process+1
   nItern      = 0
   converged   = .FALSE.

   inKODE      = 1
   outKODE     = 0
   KODE        = 0

   inITER      = 80
   outITER     = 80
   ITER        = 80


   !---------------------------------------------------------------------------!
   !-- Enter the Newton-Rhapson iterative sequence
   DO WHILE(.NOT.converged .AND. nItern < ITMAX)

     IF (verbosity > 1) THEN
       print *,'---',nItern
     END IF

     !-- Initialise
     !--
     inequalityArray = gc_zero
     deltaConc(:)    = gc_zero
     RES(:)          = gc_zero

     !-- Calculate residuals
     !-- Y = AX - T
     converged = .TRUE.
     DO componentIndex = 1,nComps2Process
       inequalityArray(components(componentIndex)%eqIndex,residColumn)       = &
                                     components(cList(componentIndex))%Total - &
                                     components(cList(componentIndex))%Value
       IF(ABS(inequalityArray(components(componentIndex)%eqIndex,residColumn))>&
                                                             convTolerance) THEN
         !-- Pure phase
         IF(components(cList(componentIndex))%CompType == PUREPHASE) THEN
           !-- May be no moles and undersaturated:
           IF(components(cList(componentIndex))%ppData%moles <= gc_zero .AND.  &
              components(cList(componentIndex))%Total -                        &
              components(cList(componentIndex))%Value > 1e-8 ) THEN

           ELSE
             converged = .FALSE.
           END IF
         ELSE
           !-- Non pure phase component has residual > tolerance
           converged = .FALSE.
         END IF
       END IF
     END DO

     IF (verbosity > 5) THEN
       print *,'RESIDUALS---',nItern
       DO componentIndex = 1,nComps2Process
         print *,' ',components(cList(componentIndex))%CompName,               &
                 components(cList(componentIndex))%Total,                      &
                 components(cList(componentIndex))%Value,                      &
                 inequalityArray(components(componentIndex)%eqIndex,residColumn)
       END DO
     END IF


     !-- Make Newton-Rhapson Jacobian matrix for new estimates
     !-- Z = dY/dX = aaC/X
     !-- Z.dX = Y

     inequalityArray(1:nOPTEqs+nEQLeqs,1:nComps2Process) = &
                  buildJacobian(components(cList), species, nComps2Process)


     !-- Add and inequality constraints: 1 for each pure phase
     IF(nINQEqs > 0) THEN

       ineqCntr = 0
       DO componentIndex = 1,nComps2Process
         IF(components(cList(componentIndex))%CompType==PUREPHASE) THEN

          ineqCntr = ineqCntr + 1

          IF(components(cList(componentIndex))%ppData%moles <= gc_zero) THEN
            !-- del = -1 will ensure the solution can only precipitate
            deltaConc(nEQLeqs+ineqCntr) = -gc_one
          ELSE
            !-- moles of PP left, set and inequality to stop dissolution>moles
            inequalityArray(nOPTEqs+nEQLeqs+ineqCntr,nEQLeqs+ineqCntr) = gc_one
            inequalityArray(nOPTEqs+nEQLeqs+ineqCntr,residColumn) =            &
                                components(cList(componentIndex))%ppData%moles
          END IF

          !-- No moles and undersaturated
          IF (components(cList(componentIndex))%Total -                        &
              components(cList(componentIndex))%Value > 1e-8 .AND.             &
              components(cList(componentIndex))%ppData%moles <= gc_zero) THEN

            inequalityArray(nOPTEqs+nEQLeqs+ineqCntr,residColumn) = gc_zero
            res((components(cList(componentIndex))%eqIndex))      = gc_zero

          ELSE
            res((components(cList(componentIndex))%eqIndex)) = gc_one
          END IF

         END IF
       END DO
     END IF


     !-- Normalize columns if necessary
     normal = gc_one

     DO componentIndex = 1, nComps2Process

       max = gc_zero
       IF (components(cList(componentIndex))%CompType == MOLEBLNCE) THEN

         DO miscIndex = 1, nComps2Process
           IF (ABS(inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex)) > max) THEN
             max = ABS(inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex))
             IF (max > min_value) THEN
               CYCLE
             END IF
           END IF
         END DO

         IF (max == gc_zero) THEN
           miscIndex  = components(cList(componentIndex))%eqIndex
           inequalityArray(miscIndex,componentIndex) = 1e-5*components(cList(componentIndex))%Total
           max = 1e-5*components(cList(componentIndex))%Total
         END IF

       END IF

       IF (components(cList(componentIndex))%CompType == MASSHYDGN) THEN
         !/* make absolute value of diagonal at least 1e-12 */

         min = 1e-25
         miscIndex  = components(cList(componentIndex))%eqIndex
         inequalityArray(miscIndex,componentIndex) =                           &
                                 inequalityArray(miscIndex,componentIndex) + min
         IF (ABS(inequalityArray(miscIndex,componentIndex)) < min) THEN
           inequalityArray(miscIndex,componentIndex) = min
         END IF

         max = gc_zero
         DO miscIndex = 1, nComps2Process
           IF(components(cList(miscIndex))%CompType == MOLEBLNCE   .OR.        &
              components(cList(miscIndex))%CompType == MASSHYDGN   .OR.        &
              components(cList(miscIndex))%CompType == MASSOXYGN) THEN

             IF(ABS(inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex)) > max) THEN
               max = ABS(inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex))
               IF (max > min_value) THEN
                 CYCLE
               END IF
             END IF
           END IF
         END DO
       END IF

       IF(max > gc_zero .AND. max < min_value) THEN

         !-- Scaling column
         DO miscIndex = 1, nComps2Process
           inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex) = &
           inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex) * &
           min_value / max
         END DO
         normal(componentIndex) =  min_value / max

       END IF
     END DO

     !-- No moles of pure phase and undersaturated, so force mass transfer to zero
     IF(nPPEqs > 0) THEN
         DO componentIndex = 1,nComps2Process
           IF(components(cList(componentIndex))%CompType==PUREPHASE) THEN
             IF(components(cList(componentIndex))%ppData%moles <= gc_zero .AND.&
                components(cList(componentIndex))%Total -  &
                  components(cList(componentIndex))%Value > 1e-8 ) THEN

               DO miscIndex = 1, nComps2Process
                 inequalityArray(components(cList(miscIndex))%eqIndex,         &
                                               componentIndex) = gc_zero
               END DO
             END IF

           END IF
         END DO
     END IF


     IF(verbosity > 10) THEN
       DO miscIndex = 1,nOPTEqs+nEQLEqs+nINQEqs
         print *,'Q: ',inequalityArray(miscIndex,:)
       END DO
       print *,'deltaconc1',deltaConc(1:nComps2Process)
       print *,'RES:',res
     END IF


     !-- Solve for dX
     KODE = inKODE
     ITER = inITER
     CALL simplexAlgorithm(nOPTEqs,nEQLEqs,nINQEqs,nComps2Process,& ! K,L,M,N
                           nOPTEqs+nEQLEqs+nINQEqs,               & ! KLMD
                           nOPTEqs+nEQLEqs+nINQEqs+nComponents,   & ! NKLMD
                           inequalityArray,                       & ! Q
                           KODE, TOLER,                           & ! KODE,TOLER
                           ITER, deltaConc, RES, ERROR)             ! ITER,X,RES,ERR
     outKODE = KODE
     outITER = ITER


     !-- For any columns that were normalised, we must re-scale the
     !-- respective delta value
     deltaConc(1:nComps2Process) = deltaConc(1:nComps2Process) * normal


     IF(verbosity > 10) THEN
       print *,'outKODE',outKODE
       print *,'deltaconc2',deltaConc(1:nComps2Process)
       print *,'RES:',res
     END IF


     !-- Rescale columns of array  ! Not sure why - Who cares?
     DO componentIndex = 1, nComps2Process
       IF (ABS(normal(componentIndex) - gc_one) > 0.0001) THEN
         DO miscIndex = 1, nComps2Process
          inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex)=&
          inequalityArray(components(cList(miscIndex))%eqIndex,componentIndex)/&
          normal(componentIndex)
         END DO
       END IF
     END DO


     !-- Now we must check the mineral transfer rates and update
     !--     the unknowns with the scaled delta values
     CALL updateUnknownsWithdX(deltaConc, components, cList)


     !-- Update Gammas
     CALL updateGammas(species, ionStrength%Total)


     !-- Update species molalities
     CALL updateSpeciesMoles(species, components)


     !-- Recalculate mole balance sums based on the new species concs
     CALL moleBalanceSums(species, components)


     !-- Lastly, perform a final check for convergence in case
     !--     something went wrong
     DO componentIndex = 1,nComps2Process
       ! Calculate residual
       inequalityArray(components(componentIndex)%eqIndex,residColumn) =       &
                        components(cList(componentIndex))%Total -              &
                        components(cList(componentIndex))%Value

       IF(components(componentIndex)%CompType /= PUREPHASE) THEN
        IF(ABS(inequalityArray(components(componentIndex)%eqIndex,residColumn))&
                                                            >convTolerance) THEN
          converged = .FALSE.
        END IF
       END IF
     END DO
     IF (verbosity > 5) THEN
       print *,'RESIDUALS2---',nItern
       DO componentIndex = 1,nComps2Process
         print *,' ',components(cList(componentIndex))%CompName,               &
                     components(cList(componentIndex))%Total,                  &
                     components(cList(componentIndex))%Value,                  &
                 inequalityArray(components(componentIndex)%eqIndex,residColumn)
       END DO
     END IF


     nItern = nItern + 1


   END DO


   !---------------------------------------------------------------------------!
   !-- Now the Newton-Rhapson sequence is finished, interogate and             !
   !   convergence issues                                                      !

   failed = .FALSE.

   IF(nItern == ITMAX) THEN
     PRINT *,'Non-convergence by ITMAX: ',nItern

     CALL outputComponentData(components)

     DO componentIndex = 1,nComps2Process
       IF(ABS(inequalityArray(components(componentIndex)%eqIndex,residColumn))>&
                                                             convTolerance) THEN

         print *,'Component Type & Name: ',                                    &
                                components(cList(componentIndex))%CompType,    &
                components(cList(componentIndex))%CompName

         IF(components(cList(componentIndex))%CompType == PUREPHASE)  THEN
!           print *,' PurePhase non-convergence after ',nItern,' steps'
         print *,' ->PP: ',  componentIndex,   cList(componentIndex),                  &
                components(cList(componentIndex))%ppdata%moles,&
                inequalityArray(components(componentIndex)%eqIndex,residColumn)

           IF(components(cList(componentIndex))%ppdata%moles < 1e-20) THEN
             print *,'  # moles <1e-20, so proceeding'
           ELSE
             print *,'  # moles >1e-20, covnvergence = FAILED'
             failed = .TRUE.
           END IF

         ELSE IF(components(cList(componentIndex))%CompType == MASSHYDGN) THEN
!           print *,'just MASSHYDGN, dont worry so much ',nItern
         ELSE IF(components(cList(componentIndex))%CompType == ACTIVYH2O) THEN
!           print *,'just ACTIVYH2O, dont worry so much ',nItern
         ELSE IF(components(cList(componentIndex))%CompType == MASSOXYGN) THEN
!           print *,'just MASSOXYGN, dont worry so much ',nItern
         ELSE
           ! Else we must have not converged yet
!            print *,'Convergence failure ',nItern
            failed = .TRUE.
         END IF
       END IF
     END DO
   END IF

   DEALLOCATE(inequalityArray)


 END SUBROUTINE solveEquilibriumReactions                                      !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! UpdateLogK                                                                   !
!                                                                              !
! Routine that updates the logK values for each species according to Tempture  !
! and their user supplied deltaH values                                        !
!------------------------------------------------------------------------------!
 SUBROUTINE updateLogK(specs, comps, cellTemp)                                 !
   !-- Incoming                                                                !
   TYPE (gcSpecies),  DIMENSION(:), INTENT(INOUT) :: specs                     !
   TYPE (gcUnknowns), DIMENSION(:), INTENT(INOUT) :: comps                     !
   DOUBLETYPE,                     INTENT(IN)    :: cellTemp                  !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER                        :: THIS_PROC = "updateLogK"  !
   DOUBLETYPE                                    :: tempK                     !
   INTEGER  :: nSpec, nComp              !
   INTEGER  :: sIndex, cIndex            !


   nSpec = SIZE(specs)
   nComp = SIZE(comps)

   tempK = cellTemp + WaterFreezingInKelvins

   DO sIndex = 1, nSpec
     specs(sIndex)%logK = ( specs(sIndex)%deltaH/(LOG(10.0)*GasConstant/1e3) ) &
                        * ( (1.0/298.15) - (1.0/tempK) )
     specs(sIndex)%logK = specs(sIndex)%logKat25 + specs(sIndex)%logK
   END DO

   DO cIndex = 1, nComp
     IF(comps(cIndex)%CompType == PUREPHASE) THEN
       comps(cIndex)%PPData%logK = (                                           &
                    comps(cIndex)%PPData%deltaH/(LOG(10.0)*GasConstant/1e3) )  &
                    * ( (1.0/298.15) - (1.0/tempK) )
       comps(cIndex)%PPData%logK  = comps(cIndex)%PPData%logKat25              &
                                  + comps(cIndex)%PPData%logK
     END IF

   END DO


 END SUBROUTINE updateLogK
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! UpdateGammas                                                                 !
!                                                                              !
! Routine that recalculates activity coefficients (gamma values) according to  !
! temperature and ionic strength                                               !
!------------------------------------------------------------------------------!
 SUBROUTINE updateGammas(specs, ionicStrength)                                 !
   !-- Incoming                                                                !
   TYPE (gcSpecies), DIMENSION(:), INTENT(INOUT) :: specs                      !
   DOUBLETYPE,                    INTENT(IN)    :: ionicStrength              !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER                       :: THIS_PROC = "updateGammas" !
   DOUBLETYPE                                   :: tempK                      !
   DOUBLETYPE                                   :: Aconst, Bconst             !
   DOUBLETYPE                                   :: d1, s1, s2, s3             !
   INTEGER  :: nSpec                      !
   INTEGER  :: sIndex                     !
                                                                               !

   nSpec = SIZE(specs)

   tempK = cellTemp + WaterFreezingInKelvins

   s1 = 374.11 - cellTemp

   s2 = s1**(1.0/3.0)

   s3 = gc_one + 0.1342489*s2 - 3.946263e-03*s1

   s3 = s3 / &
  (3.1975 - 0.3151548*s2 - 1.203374e-03*s1 + 7.48908e-13*(s1*s1*s1*s1))

   s3 = SQRT(s3)

   IF (tempK >= 373.15) THEN
     d1 = 5321.0/tempK+ 233.76- tempK*(tempK*(8.292e-07*tempK-1.417e-03)+0.9297)
   ELSE
     d1 = 2727.586 + 0.6224107*tempK - 466.9151*LOG(tempK) - 52000.87/tempK
   END IF

   d1 = SQRT(d1*tempK)
   Aconst = 1824827.7 * s3 / (d1 * d1 * d1)
   Bconst = 50.2905   * s3 / d1

   DO sIndex = 1, nSpec

     IF(specs(sIndex)%gFlag == UNCHGD) THEN
       !-- Uncharged species use the Setchenow equation, with b = 0.1
       specs(sIndex)%logGamma = 0.1000 * ionicStrength
       specs(sIndex)%delGamma = 0.1000 * LOG(REAL(10.0,DP)) * specs(sIndex)%Moles

     ELSE IF(specs(sIndex)%gFlag == DAVIES) THEN
       !-- Gamma calculated using the Davies equation
       specs(sIndex)%logGamma = gammaDavies(specs(sIndex)%charge,ionicStrength,Aconst)
       specs(sIndex)%delGamma = specs(sIndex)%Moles * &
               delgammaDavies(specs(sIndex)%charge,ionicStrength,Aconst)

     ELSE IF(specs(sIndex)%gFlag == WATEQDH) THEN
       !-- Gamma calculated using the Extended/WATEQ Debye-Huckel equation

       specs(sIndex)%logGamma = gammaWATEQDH(specs(sIndex)%charge,ionicStrength,Aconst,Bconst)
       specs(sIndex)%delGamma = specs(sIndex)%Moles * &
               delgammaWATEQDH(specs(sIndex)%charge,ionicStrength,Aconst,Bconst)
     ELSE
       specs(sIndex)%logGamma  = gc_zero
       specs(sIndex)%delGamma  = gc_zero

     END IF

   END DO


 END SUBROUTINE updateGammas
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
 FUNCTION gammaDavies(speciesCharge, ionicStrength, AA) RESULT(logGamma)       !
   !-- Incoming                                                                !
   INTEGER, INTENT(IN)      :: speciesCharge
   DOUBLETYPE, INTENT(IN)  :: ionicStrength
   DOUBLETYPE, INTENT(IN)  :: AA
   !-- Outgoing
   DOUBLETYPE  :: logGamma

   logGamma = -AA*(REAL(speciesCharge,GCHP)**2.0)                              &
           *((SQRT(ionicStrength)/(1.0+SQRT(ionicStrength)))-0.30*ionicStrength)

 END FUNCTION gammaDavies
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
 FUNCTION delgammaDavies(speciesCharge, ionicStrength, AA) RESULT(delGamma)    !
   !-- Incoming                                                                !
   INTEGER, INTENT(IN)      :: speciesCharge
   DOUBLETYPE, INTENT(IN)  :: ionicStrength
   DOUBLETYPE, INTENT(IN)  :: AA
   !-- Outgoing
   DOUBLETYPE  :: delGamma

   delGamma =  - (LOG(REAL(10.0,DP)) * AA  * (REAL(speciesCharge,GCHP))**2.0)  &
            * (1.0/(2.0*SQRT(ionicStrength)*(SQRT(ionicStrength)+1.0)**2.0)-0.3)

 END FUNCTION delgammaDavies
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
 FUNCTION gammaWATEQDH(speciesCharge, ionicStrength, AA, BB) RESULT(logGamma)  !
   !-- Incoming                                                                !
   INTEGER,     INTENT(IN)  :: speciesCharge
   DOUBLETYPE, INTENT(IN)  :: ionicStrength
   DOUBLETYPE, INTENT(IN)  :: AA
   DOUBLETYPE, INTENT(IN)  :: BB
   !-- Outgoing
   DOUBLETYPE  :: logGamma

   logGamma = (-AA * ((REAL(speciesCharge,GCHP))**2.0) * SQRT(ionicStrength)) /&
              (1.0 +  9.0*BB*SQRT(ionicStrength)) + 0.10*ionicStrength


 END FUNCTION gammaWATEQDH
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
 FUNCTION delgammaWATEQDH(speciesCharge, ionicStrength, AA, BB) RESULT(delGamma)
   !-- Incoming
   INTEGER, INTENT(IN)      :: speciesCharge
   DOUBLETYPE, INTENT(IN)  :: ionicStrength
   DOUBLETYPE, INTENT(IN)  :: AA
   DOUBLETYPE, INTENT(IN)  :: BB
   !-- Outgoing
   DOUBLETYPE  :: delGamma

   delGamma = -LOG(REAL(10.0,DP))*( (AA*(REAL(speciesCharge,GCHP))**2.0) /     &
            (2.0*SQRT(ionicStrength)*(BB*9.0*SQRT(ionicStrength)+1)**2.0) + 0.1)

  END FUNCTION delgammaWATEQDH
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! SetEqIndiciesForSpeciation                                                   !
!                                                                              !
! Routine called to assign each active component an "eq" (i.e. column) number  !
! in the inequality matrix being solved by SolveEquilibriumReactions           !
! This routine assumes NO pure phase assemblage is active in the solution      !
! matrix                                                                       !
!------------------------------------------------------------------------------!
 SUBROUTINE seteqIndiciesforSpeciation(comps,nComps2Process,componentList)     !
   !-- Incoming                                                                !
   TYPE (gcUnknowns), DIMENSION(:),  INTENT(INOUT) :: comps                    !
   INTEGER,           DIMENSION(:),  INTENT(OUT)   :: componentList            !
   INTEGER,                          INTENT(IN)    :: nComps2Process           !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "seteqIndiciesforSpeciation"         !
   INTEGER  :: cIndex, eqCounter


   eqCounter = 1
   DO cIndex = 1,nComponents
     comps(cIndex)%eqIndex = 0

     IF(comps(cIndex)%CompType == MOLEBLNCE) THEN
       comps(cIndex)%eqIndex = eqCounter
       componentList(eqCounter) = comps(cIndex)%CompIndex

       eqCounter = eqCounter + 1
     END IF
   END DO

   ionStrength%eqIndex   = eqCounter
   componentList(eqCounter) = ionStrength%CompIndex
   eqCounter = eqCounter + 1

   activityWater%eqIndex = eqCounter
   componentList(eqCounter) = activityWater%CompIndex

   IF( .NOT. pHisFixed ) THEN
     eqCounter = eqCounter + 1

     hydrogenIon%eqIndex = eqCounter
     componentList(eqCounter) = hydrogenIon%CompIndex
   END IF


   IF(eqCounter /= nComps2Process) THEN
     print *,'Problem setting eqIndex for Speciation:',eqCounter,nComps2Process
     STOP
   END IF


 END SUBROUTINE seteqIndiciesforSpeciation
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! SetEqIndiciesForBatchRxn                                                     !
!                                                                              !
! Routine called to assign each active component an "eq" (i.e. column) number  !
! in the inequality matrix being solved by SolveEquilibriumReactions           !
! This routine assumes a pure phase assemblage IS active in the solution       !
! matrix                                                                       !
!------------------------------------------------------------------------------!
 SUBROUTINE SeteqIndiciesforBatchRxn(comps,nComps2Process,componentList)       !
   !-- Incoming                                                                !
   TYPE (gcUnknowns), DIMENSION(:), INTENT(INOUT) :: comps                     !
   INTEGER,           DIMENSION(:), INTENT(OUT)   :: componentList             !
   INTEGER,                         INTENT(IN)    :: nComps2Process            !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "seteqIndiciesforSpeciation"         !
   INTEGER  :: cIndex, eqCounter, listCounter                   !


   eqCounter = 1
   DO cIndex = 1,nComponents
     comps(cIndex)%eqIndex = 0
   END DO

   DO cIndex = 1,nComponents
     IF(comps(cIndex)%CompType == PUREPHASE) THEN
       comps(cIndex)%eqIndex = eqCounter
       componentList(eqCounter) = comps(cIndex)%CompIndex

       eqCounter = eqCounter + 1
     END IF
   END DO

   DO cIndex = 1,nComponents
     IF(comps(cIndex)%CompType == MOLEBLNCE) THEN
       comps(cIndex)%eqIndex = eqCounter
       componentList(eqCounter) = comps(cIndex)%CompIndex

       eqCounter = eqCounter + 1
     END IF
   END DO

   ionStrength%eqIndex   = eqCounter
   eqCounter = eqCounter + 1

   activityWater%eqIndex = eqCounter
   eqCounter = eqCounter + 1

   hydrogenIon%eqIndex = eqCounter
   eqCounter = eqCounter + 1

   massHydrogen%eqIndex = eqCounter
   eqCounter = eqCounter + 1

   massWater%eqIndex = eqCounter


   listCounter = 1
   DO cIndex = 1,nComponents
     IF(comps(cIndex)%CompType == MOLEBLNCE) THEN
       componentList(listCounter) = comps(cIndex)%CompIndex

       listCounter = listCounter + 1
     END IF
   END DO

   componentList(listCounter) = ionStrength%CompIndex
   listCounter = listCounter + 1

   componentList(listCounter) = activityWater%CompIndex
   listCounter = listCounter + 1

   componentList(listCounter) = hydrogenIon%CompIndex
   listCounter = listCounter + 1

   componentList(listCounter) = massHydrogen%CompIndex
   listCounter = listCounter + 1

   componentList(listCounter) = massWater%CompIndex
   listCounter = listCounter + 1

   DO cIndex = 1,nComponents
     IF(comps(cIndex)%CompType == PUREPHASE) THEN
       componentList(listCounter) = comps(cIndex)%CompIndex

       listCounter = listCounter + 1
     END IF
   END DO



   IF(eqCounter /= nComps2Process .OR. (listCounter-1) /= nComps2Process) THEN
     PRINT *,'Problem setting eqIndex for Batch Rxn:',                         &
              eqCounter,nComps2Process,listCounter-1
     STOP
   END IF


 END SUBROUTINE SeteqIndiciesforBatchRxn
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! UpdateUnknownsWithdX                                                         !
!                                                                              !
! Routine to:                                                                  !
!    * Checks deltas (changes to unknowns) returned by the simplex solver to   !
!      make sure they are physically reasonable                                !
!    * Scales deltas if necessary                                              !
!    * Updates unknowns with deltas                                            !
!                                                                              !
!------------------------------------------------------------------------------!
SUBROUTINE UpdateUnknownsWithdX(deltaConc,comps,cList)                        !
   !-- Incoming                                                                !
   TYPE (gcUnknowns), DIMENSION(:), INTENT(INOUT) :: comps                     !
   DOUBLETYPE,        DIMENSION(:), INTENT(INOUT) :: deltaConc                 !
   INTEGER,           DIMENSION(:), INTENT(IN)    :: cList                     !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "updateUnknownsWithdX"               !
   INTEGER  :: nComp                                            !
   INTEGER  :: cIndex,componentIndex,miscIndex                  !
   REAL (DP)               :: highPrecBuffer                                   !
   DOUBLETYPE              :: factor,f0,up,down,step_up,sum_deltas,mu_calc,dd  !


   nComp = SIZE(cList)

   !--  Calculate interphase mass transfers
   step_up = log(100.0)
   factor  = gc_one


   !---------------------------------------------------------------------------!
   IF(SolveWithPP) THEN

     !-- Checks deltas to make sure they are physically reasonable
     DO componentIndex = 1,nComp

       ! All components have been passed through to this function
       ! even though the only a few may have been included in the
       ! inequality solver. We therefore find the actual component
       ! index using the cList map.

       cIndex = cList(componentIndex)
       IF(comps(cIndex)%CompType == PUREPHASE) THEN

         !-- Check for over-dissolution
         IF (comps(cIndex)%ppData%moles > gc_zero .AND.                        &
                    deltaConc(componentIndex) > comps(cIndex)%ppData%moles) THEN

           f0 = deltaConc(componentIndex) / comps(cIndex)%ppData%moles
           IF (f0 > factor) THEN
             IF(verbosity > 5 ) THEN
               print *,'Removing more mineral than is available:',             &
                        comps(cIndex)%CompName,deltaConc(componentIndex),      &
                        comps(cIndex)%ppData%moles,f0
             END IF
             factor = f0
           END IF

         !-- Check for over-dissolution
         ELSE IF ( deltaConc(componentIndex) > gc_zero .AND.                   &
                                    comps(cIndex)%ppData%moles <= gc_zero ) THEN
           IF(verbosity > 5 ) THEN
              print *,'Dissolving mineral with 0 mass:',                       &
                 comps(cIndex)%CompName, deltaConc(componentIndex),            &
                 comps(cIndex)%ppData%moles
           END IF
           deltaConc(componentIndex) = gc_zero

         !-- Check over precipitation
     ELSE IF ( deltaConc(componentIndex) < gc_zero  .AND.                  &
                  comps(cIndex)%ppData%moles > gc_zero  .AND.                  &
                  deltaConc(componentIndex) < -100.0 )  THEN

           f0 = -deltaConc(componentIndex) / 100.0
           IF (f0 > factor) THEN
             IF(verbosity > 5 ) THEN
               print *,'Precipitating too much mineral:',                      &
                     comps(cIndex)%CompName,deltaConc(componentIndex),         &
                     comps(cIndex)%ppData%moles, f0
             END IF
            factor = f0
           END IF

         END IF

       END IF
     END DO

   END IF


   !-- Calculate change in element concentrations due to pure phases
   IF(SolveWithPP) THEN
     DO componentIndex = 1,nComp
       ! All components have been passed through to this function
       ! even though the only a few may have been included in the
       ! inequality solver. We therefore find the actual component
       ! index using the cList map.

       cIndex = cList(componentIndex)
       comps(cIndex)%delta = gc_zero
       IF(comps(cIndex)%CompType == MOLEBLNCE   .OR.                           &
          comps(cIndex)%CompType == MASSHYDGN   .OR.                           &
          comps(cIndex)%CompType == MASSOXYGN) THEN

         DO miscIndex = 1,nComp
           IF(comps(cList(miscIndex))%CompType == PUREPHASE) THEN

              comps(cIndex)%delta = comps(cIndex)%delta +                      &
                          comps(cList(miscIndex))%ppData%stoich(cIndex)        &
                          * deltaConc(miscIndex)

           END IF
         END DO
       END IF
     END DO

     !-- Apply factor from minerals to deltas
     DO componentIndex = 1,nComp
       ! All components have been passed through to this function
       ! even though the only a few may have been included in the
       ! inequality solver. We therefore find the actual component
       ! index using the cList map.

       cIndex = cList(componentIndex)

       comps(cIndex)%delta = comps(cIndex)%delta / factor

       IF(comps(cIndex)%CompType == PUREPHASE) THEN
         deltaConc(componentIndex) = deltaConc(componentIndex) / factor
       END IF
     END DO
   END IF



   !---------------------------------------------------------------------------!
   !-- Calc factor for mass balance equations for aqueous unknowns
   factor      = gc_one
   sum_deltas  = gc_zero

   DO componentIndex = 1,nComp

     ! All components have been passed through to this function
     ! even though the only a few may have been included in the
     ! inequality solver. We therefore find the actual component
     ! index using the cList map.

     cIndex     = cList(componentIndex)
     sum_deltas = sum_deltas + ABS(deltaConc(componentIndex))

     up         = step_up
     down       = up

     IF (comps(cIndex)%CompType == MOLEBLNCE .OR. &
         comps(cIndex)%CompType == CHARGEBAL ) THEN

       up   = step_up
       down = 1.3 * up

     ELSE IF (comps(cIndex)%CompType == IONSTNGTH) THEN

       up   = 100.0 * ionStrength%Total
       down = ionStrength%Total

     ELSE IF (comps(cIndex)%CompType == ACTIVYH2O) THEN

       down = up

     ELSE IF (comps(cIndex)%CompType == MASSHYDGN) THEN

       up   = log(10.0)
       down =  1.3*up

     ELSE IF (comps(cIndex)%CompType == MASSOXYGN) THEN

       up   = log(10.)
       down = log(4.)

     ELSE IF (comps(cIndex)%CompType == PUREPHASE) THEN

       continue

     END IF


     IF (deltaConc(componentIndex) > gc_zero) THEN

       f0 = deltaConc(componentIndex) / up
       IF (f0 > factor) THEN
         factor = f0
       END IF

     ELSE

       f0 = deltaConc(componentIndex) / (-down)
       IF (f0 > factor) THEN
         factor = f0
       END IF
     END IF

   END DO

   factor = gc_one /factor

   !-- Update the deltas for the mass balance equations based on the
   !   above constraints
   DO componentIndex = 1,nComp
     ! All components have been passed through to this function
     ! even though the only a few may have been included in the
     ! inequality solver. We therefore find the actual component
     ! index using the cList map.
     cIndex = cList(componentIndex)

     IF(comps(cIndex)%CompType /= PUREPHASE) THEN
       deltaConc(componentIndex) = deltaConc(componentIndex) * factor
     END IF
   END DO



   !---------------------------------------------------------------------------!
   !-- Now apply the changes
   DO componentIndex = 1,nComp

     ! All components have been passed through to this function
     ! even though the only a few may have been included in the
     ! inequality solver. We therefore find the actual component
     ! index using the cList map.
     cIndex = cList(componentIndex)


     !-- Pure Phases
     IF(comps(cIndex)%CompType == PUREPHASE) THEN
       IF (verbosity > 5) THEN
         print *,'PUREPHASE: X=X+dX',cIndex,comps(cIndex)%ppData%moles,        &
                                     comps(cIndex)%ppData%moles -              &
                                     deltaConc(componentIndex),                &
                     deltaConc(componentIndex)
       END IF

       comps(cIndex)%ppData%moles = comps(cIndex)%ppData%moles -               &
                                                       deltaConc(componentIndex)


     ! Ionic Strength
     ELSE IF(comps(cIndex)%CompType == IONSTNGTH ) THEN
       IF (verbosity > 5) THEN
         print *,'MU: X=X+dX',cIndex,comps(cIndex)%Total, &
              comps(cIndex)%Total + deltaConc(componentIndex),&
                                          deltaConc(componentIndex)
       END IF

       mu_calc = 0.500 * comps(cIndex)%Value / Waq
       dd      = comps(cIndex)%Total + deltaConc(componentIndex)

       IF (dd < 1e-7) THEN
         deltaConc(componentIndex) = SQRT(mu_calc * comps(cIndex)%Total)       &
                                                           - comps(cIndex)%Total
         comps(cIndex)%Total = SQRT(mu_calc * comps(cIndex)%Total)
       ELSE
         comps(cIndex)%Total = comps(cIndex)%Total + deltaConc(componentIndex)
       END IF

       IF (comps(cIndex)%Total <= 1e-8) THEN
         comps(cIndex)%Total = 1e-8
       END IF


     !-- Mole Balances
     ELSE IF(comps(cIndex)%CompType == MOLEBLNCE ) THEN
       !ln(a) = ln(a) + dln(a)
       IF (verbosity > 5) THEN
         print *,'MB: X=X+dX',cIndex,comps(cIndex)%master%logActivity,         &
                              comps(cIndex)%master%logActivity +               &
                              deltaConc(componentIndex) / LOG(REAL(10.0,DP)),  &
                  deltaConc(componentIndex)
       END IF

       comps(cIndex)%master%logActivity = comps(cIndex)%master%logActivity     &
                                + deltaConc(componentIndex) / LOG(REAL(10.0,DP))


     !-- Charge Balance
     ELSE IF(comps(cIndex)%CompType == CHARGEBAL ) THEN
       !ln(a) = ln(a) + dln(a)
       IF (verbosity > 5) THEN
         print *,'CB: X=X+dX',cIndex,comps(cIndex)%master%logActivity,         &
                              comps(cIndex)%master%logActivity +               &
                              deltaConc(componentIndex) / LOG(REAL(10.0,DP)),  &
                  deltaConc(componentIndex)
       END IF

       comps(cIndex)%master%logActivity = comps(cIndex)%master%logActivity     &
                                + deltaConc(componentIndex) / LOG(REAL(10.0,DP))


     !-- Activity of Water
     ELSE IF(comps(cIndex)%CompType == ACTIVYH2O ) THEN
       IF (verbosity > 5) THEN
         PRINT *,'AH2O X=X+dX',cIndex,comps(cIndex)%master%logActivity,        &
                  Waq*EXP(comps(cIndex)%master%logActivity*LOG(REAL(10.0,DP))) &
                - Waq,deltaConc(componentIndex),deltaConc(componentIndex)      &
        / LOG(REAL(10.0,DP))
       END IF

       highPrecBuffer = REAL(comps(cIndex)%master%logActivity, DP)             &
                      + REAL(deltaConc(componentIndex), DP) / LOG(REAL(10.0,DP))
       IF(highPrecBuffer < -1.0) THEN
         highPrecBuffer = -1.0
       END IF

       comps(cIndex)%master%logActivity = REAL(highPrecBuffer,GCHP)

       comps(cIndex)%Total = Waq* (antiLOG(activityWater%master%logActivity,   &
                                                         REAL(10.0,GCHP)) - 1.0)


     !-- Mass of Hydrogen
     ELSE IF(comps(cIndex)%CompType == MASSHYDGN ) THEN
       !ln(a) = ln(a) + dln(a)
       IF (verbosity > 5) THEN
         print *,'MH: X=X+dX',cIndex,comps(cIndex)%master%logActivity,         &
                              comps(cIndex)%master%logActivity                 &
                + deltaConc(componentIndex) / LOG(REAL(10.0,DP)),  &
                              deltaConc(componentIndex)
       END IF

       comps(cIndex)%master%logActivity = comps(cIndex)%master%logActivity     &
                               +  deltaConc(componentIndex) / LOG(REAL(10.0,DP))


     !-- Mass of Water
     ELSE IF(comps(cIndex)%CompType == MASSOXYGN ) THEN
       highPrecBuffer = REAL(deltaConc(componentIndex),DP)
       Waq = Waq * REAL(EXP(highPrecBuffer),GCHP)
       IF (verbosity > 5) THEN
         print *,'MH20: X=X+dX',comps(cIndex)%Total, Waq / molWgtH2O ,         &
                                deltaConc(componentIndex)
       END IF

       activityWater%master%Moles = Waq / molWgtH2O

       IF (Waq < 1e-10) THEN
         print *, 'Mass of water is less than 1e-10 kilogram, STOPPING'
         STOP
       END IF

     END IF
   END DO



   !---------------------------------------------------------------------------!
   !--    Reset total molalities in mass balance equations

   IF(SolveWithPP) THEN
     DO componentIndex = 1,nComp
         ! All components have been passed through to this function
         ! even though the only a few may have been included in the
         ! inequality solver. We therefore find the actual component
         ! index using the cList map.

         cIndex = cList(componentIndex)

         IF(comps(cIndex)%CompType == MOLEBLNCE    .OR. &
            comps(cIndex)%CompType == MASSHYDGN    .OR. &
            comps(cIndex)%CompType == MASSOXYGN    .OR. &
            comps(cIndex)%CompType == CHARGEBAL   ) THEN

            IF (verbosity > 5) THEN
              print *,'adjust Total',comps(cIndex)%CompName,comps(cIndex)%Total&
                                ,comps(cIndex)%delta
            END IF

            comps(cIndex)%Total = comps(cIndex)%Total + comps(cIndex)%delta

         END IF

     END DO
   END IF



END SUBROUTINE UpdateUnknownsWithdX
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! MoleBalanceSums                                                              !
!                                                                              !
! Estimates the total number of moles for each component based on current      !
! speciation                                                                   !
!------------------------------------------------------------------------------!
 SUBROUTINE MoleBalanceSums(specs, comps)                                      !
   !-- Incoming                                                                !
   TYPE (gcSpecies),  DIMENSION(:), INTENT(IN)    :: specs                     !
   TYPE (gcUnknowns), DIMENSION(:), INTENT(INOUT) :: comps                     !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "MoleBalanceSums"                    !
   INTEGER  :: nSpec, nComp
   INTEGER  :: sIndex, cIndex, miscIndex


   nComp = SIZE(comps)
   nSpec = SIZE(specs)

   DO cIndex = 1,nComp
     comps(cIndex)%Value = gc_zero

     IF(comps(cIndex)%CompType == MOLEBLNCE) THEN
       ! X = sumoveri aijC
       DO sIndex = 1,nSpec
         comps(cIndex)%Value = comps(cIndex)%Value                             &
                             + specs(sIndex)%stoich(comps(cIndex)%CompIndex)   &
                             * specs(sIndex)%Moles
       END DO

     ELSE IF(comps(cIndex)%CompType == PUREPHASE) THEN
       DO miscIndex = 1,nComp
         IF(comps(miscIndex)%CompType == MOLEBLNCE .OR.                        &
            comps(miscIndex)%CompType == CHARGEBAL .OR.                        &
            comps(miscIndex)%CompType == ACTIVYH2O) THEN

           comps(cIndex)%Value = comps(cIndex)%Value                           &
                    + comps(cIndex)%ppData%stoich(comps(miscIndex)%CompIndex)  &
                    * comps(miscIndex)%master%logActivity * LOG(REAL(10.0,DP))
         END IF
       END DO

     ELSE IF(comps(cIndex)%CompType == IONSTNGTH) THEN
       comps(cIndex)%Value = updateIonicStrength(specs)

     ELSE IF(comps(cIndex)%CompType == ACTIVYH2O) THEN
       ! X = 0.017 * sumoveri Ci (but don;t include H2O!)
       DO sIndex = 1,nSpec
         IF(specs(sIndex)%Name /= TRIM('H2O') .AND. &
            specs(sIndex)%Name /= TRIM('e-')) THEN

        comps(cIndex)%Value = comps(cIndex)%Value + specs(sIndex)%Moles
         END IF
       END DO

       comps(cIndex)%Value =  - 0.01700 * comps(cIndex)%Value

     ELSE IF(comps(cIndex)%CompType == CHARGEBAL) THEN
       ! X = sumoveri ziCi
       DO sIndex = 1,nSpec
         comps(cIndex)%Value = comps(cIndex)%Value                             &
                             + (specs(sIndex)%charge * specs(sIndex)%Moles)
       END DO

     ELSE IF(comps(cIndex)%CompType == MASSHYDGN) THEN
       ! X = sumoveri aijC
       DO sIndex = 1,nSpec
         comps(cIndex)%Value = comps(cIndex)%Value                             &
                             + specs(sIndex)%Hstoich * specs(sIndex)%Moles
       END DO

     ELSE IF(comps(cIndex)%CompType == MASSOXYGN) THEN
       ! X = sumoveri aijC
       DO sIndex = 1,nSpec
         comps(cIndex)%Value = comps(cIndex)%Value                             &
                             + specs(sIndex)%H2Ostoich * specs(sIndex)%Moles
       END DO


     END IF

   END DO


 END SUBROUTINE MoleBalanceSums
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! UpdateSpeciesMoles                                                           !
!                                                                              !
! Routine to calculate concentration of each aqueous species using the         !
! current component concentrations                                             !
!------------------------------------------------------------------------------!
 SUBROUTINE UpdateSpeciesMoles(specs, comps)                                   !
   !-- Incoming                                                                !
   TYPE (gcUnknowns),  DIMENSION(:), INTENT(IN)    :: comps                    !
   TYPE (gcSpecies),   DIMENSION(:), INTENT(INOUT) :: specs                    !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "updateSpeciesMoles"                 !
   INTEGER  :: nSpec,  nComp                                    !
   INTEGER  :: sIndex, cIndex                                   !
   DOUBLETYPE             :: sum                                              !


   nSpec = SIZE(specs)
   nComp = SIZE(comps)

   DO sIndex = 1, nSpec
     sum = gc_zero
     sum = specs(sIndex)%logK - specs(sIndex)%logGamma

     DO cIndex = 1,nComp
       IF(comps(cIndex)%CompType == MOLEBLNCE .OR.                             &
          comps(cIndex)%CompType == CHARGEBAL .OR.                             &
          comps(cIndex)%CompType == MASSHYDGN .OR.                             &
          comps(cIndex)%CompType == ACTIVYH2O ) THEN

         sum = sum + specs(sIndex)%stoich(comps(cIndex)%CompIndex) *           &
                                   comps(cIndex)%master%logActivity
       END IF
     END DO


     IF(TRIM(specs(sIndex)%Name) /= "H2O" .AND. &
        TRIM(specs(sIndex)%Name) /= "e-") THEN
       specs(sIndex)%Moles = antiLOG(sum,REAL(10.0,GCHP)) / Waq
     END IF

   END DO


 END SUBROUTINE UpdateSpeciesMoles                                             !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! UpdateIonicStrength                                                          !
!                                                                              !
! Function to calculate ionic strength based on concentrations of the          !
! aqueous species and the charge                                               !
!------------------------------------------------------------------------------!
 FUNCTION UpdateIonicStrength(specs) RESULT(ionicStrength)                     !
   !-- Incoming                                                                !
   TYPE (gcSpecies), DIMENSION(:), INTENT(IN) :: specs                         !
   !-- Outgoing                                                                !
   DOUBLETYPE                                :: ionicStrength                 !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "updateIonicStrength"                !
   INTEGER  :: nSpec
   INTEGER  :: sIndex


   nSpec = SIZE(specs)
   ionicStrength = gc_zero

   DO sIndex = 1, nSpec
    IF(TRIM(specs(sIndex)%Name) /= "e-") THEN
     ionicStrength = ionicStrength + &
                     (specs(sIndex)%charge**2.0)*(specs(sIndex)%Moles/Waq)
    END IF
   END DO

   ionicStrength = 0.5000000 * ionicStrength

   IF(ionicStrength<1e-8) THEN
     ionicStrength = 1e-8
   END IF


 END FUNCTION UpdateIonicStrength                                              !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! InitialIonicStrength                                                         !
!                                                                              !
! This function roughly calculates ionic strength of the solution as an        !
! initial guess for the equilibrim calculation. It assumes all component       !
! moles are in the form of the master specie.                                  !
!------------------------------------------------------------------------------!
 FUNCTION InitialIonicStrength(comps) RESULT(ionicStrength)                    !
   !-- Incoming                                                                !
   TYPE (gcUnknowns), DIMENSION(:), INTENT(IN) :: comps                        !
   !-- Outgoing                                                                !
   DOUBLETYPE                                 :: ionicStrength                !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "InitialIonicStrength"               !
   INTEGER  :: nComp                                            !
   INTEGER  :: cIndex                                           !


   nComp = SIZE(comps)
   ionicStrength = gc_zero

   DO cIndex = 1, nComp
     IF(comps(cIndex)%CompType == MOLEBLNCE .OR. &
        comps(cIndex)%CompType == CHARGEBAL) THEN

       ionicStrength = ionicStrength + &
             (comps(cIndex)%master%charge**2.0)*(comps(cIndex)%Total/Waq)

     END IF
   END DO

   ionicStrength =  0.50 * ionicStrength

   IF(ionicStrength < 1e-8) THEN
     ionicStrength = 1e-8
   END IF


 END FUNCTION InitialIonicStrength                                             !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! BuildJacobian                                                                !
!                                                                              !
! Function to calculate Jacobian (derivative) matrix required for the Newton-  !
! Rhapson solution. It is returned to the inequality array where it is then    !
! passed to the Simplex algorithm for each iteration                           !
!------------------------------------------------------------------------------!
 FUNCTION BuildJacobian(comps, specs, nComps2Process) RESULT(jacobian)         !
   !-- Incoming                                                                !
   TYPE (gcUnknowns), DIMENSION(:),   INTENT(IN) :: comps                      !
   TYPE (gcSpecies),  DIMENSION(:),   INTENT(IN) :: specs                      !
   INTEGER, INTENT(IN)                           :: nComps2Process             !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "buildJacobian"                      !
   DOUBLETYPE, DIMENSION(nOPTEqs+nEQLeqs,nComps2Process) :: jacobian          !
   INTEGER  :: rowCompIndex,colCompIndex, speciesIndex, rowCounter              !
   DOUBLETYPE :: speciesSums                                                  !


   jacobian = gc_zero
   rowCounter = 1

   !---------------------------------------------------------------------------!
   !-- Start at the top of the inequality array (ie the 'A' matrix) and        !
   !-- successively calculate the (j,k) jacobian element :                     !
   !-- Zjk = dYj/dXk                                                           !

   DO rowCompIndex = 1, nComps2Process            ! j index

     !-- All PP components get an equation to optimise in the "A" array
     IF(comps(rowCompIndex)%CompType == PUREPHASE) THEN

       DO colCompIndex = 1, nComps2Process        ! k index
         IF(comps(colCompIndex)%CompType == MOLEBLNCE   .OR.           &   ! MB
            comps(colCompIndex)%CompType == CHARGEBAL   .OR.           &   ! H+
            comps(colCompIndex)%CompType == ACTIVYH2O ) THEN               ! H2O

           ! dfp/dln(aj) in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = &
                         comps(rowCompIndex)%ppData%stoich(colCompIndex)

         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfp/du in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = gc_zero

         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfp/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = gc_zero

         END IF
       END DO

       !--OPTIONAL CHECK:
       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter,' eqIndex:',comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1

     END IF
   END DO

   !---------------------------------------------------------------------------!
   !-- Now move onto the mole-balance equations (ie the 'B' matrix) and        !
   !-- calculate the (j,k) jacobian element :                                  !
   !-- Zjk = dYj/dXk                                                           !

   DO rowCompIndex = 1, nComps2Process            ! j index

     !--------------------------
     !-- All MOLEBLNCE components first get an equation in the "B" array!
     IF(comps(rowCompIndex)%CompType == MOLEBLNCE) THEN

       DO colCompIndex = 1, nComps2Process        ! k index
         IF(comps(colCompIndex)%CompType == MOLEBLNCE   .OR. &   ! MB eqn
            comps(colCompIndex)%CompType == MASSHYDGN   .OR. &   ! MH
            comps(colCompIndex)%CompType == CHARGEBAL   ) THEN   ! MO

           ! dfm/dln(aj) in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               specs(speciesIndex)%stoich(comps(rowCompIndex)%CompIndex) *     &
               specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *     &
               specs(speciesIndex)%Moles

           END DO

           jacobian(rowCounter,colCompIndex) = speciesSums

         ELSE IF(comps(colCompIndex)%CompType == ACTIVYH2O ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           IF(pHisFixed) THEN
             jacobian(rowCounter,colCompIndex) = 0.00 ! Not sure about that
           ELSE
             speciesSums = gc_zero
             DO speciesIndex = 1,nSpecies
               speciesSums = speciesSums +                                     &
                 specs(speciesIndex)%stoich(comps(rowCompIndex)%CompIndex) *   &
                 specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *   &
                 specs(speciesIndex)%Moles

             END DO
             jacobian(rowCounter,colCompIndex) = speciesSums
           END IF


         ELSE IF(comps(colCompIndex)%CompType == MASSOXYGN ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
               speciesSums = speciesSums +                                     &
               specs(speciesIndex)%stoich(comps(rowCompIndex)%CompIndex) *     &
               specs(speciesIndex)%Moles
           END DO
           jacobian(rowCounter,colCompIndex) = speciesSums


         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfm/du in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               specs(speciesIndex)%stoich(comps(rowCompIndex)%CompIndex) *     &
               specs(speciesIndex)%delGamma

           END DO
           jacobian(rowCounter,colCompIndex) = -speciesSums

         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfm/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) =                                 &
                         -comps(colCompIndex)%ppData%stoich(rowCompIndex)
         END IF
       END DO

       !##OPTIONAL CHECK:
       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter,' eqIndex:',comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1


     !--------------------------
     !-- Now the MU components gets an equation in the "B" array
     ELSE IF(comps(rowCompIndex)%CompType == IONSTNGTH) THEN

       DO colCompIndex = 1, nComps2Process
         IF(comps(colCompIndex)%CompType == MOLEBLNCE .OR. &     ! MB
            comps(colCompIndex)%CompType == CHARGEBAL ) THEN     ! H+
           ! dfm/dln(aj) in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
              ((specs(speciesIndex)%charge)**2.00) *                           &
               specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *     &
               specs(speciesIndex)%Moles

           END DO
           jacobian(rowCounter,colCompIndex) = 0.50000*speciesSums

         ELSE IF(comps(colCompIndex)%CompType == ACTIVYH2O ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           IF(pHisFixed) THEN
             jacobian(rowCounter,colCompIndex) = 0.00
           ELSE
             speciesSums = gc_zero
             DO speciesIndex = 1,nSpecies
               speciesSums = speciesSums +                                     &
                ((specs(speciesIndex)%charge)**2.00) *                         &
                 specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *   &
                 specs(speciesIndex)%Moles

             END DO
             jacobian(rowCounter,colCompIndex) = 0.50000*speciesSums
           END IF

         ELSE IF(comps(colCompIndex)%CompType == MASSHYDGN ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = 0.00  ! NOT

         ELSE IF(comps(colCompIndex)%CompType == MASSOXYGN ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
              ((specs(speciesIndex)%charge)**2.00) *                           &
               specs(speciesIndex)%Moles

           END DO
           jacobian(rowCounter,colCompIndex) = ionStrength%Total-0.5*speciesSums

         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfu/du in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums + specs(speciesIndex)%delGamma *        &
                (specs(speciesIndex)%charge ** 2.0)
           END DO
           jacobian(rowCounter,colCompIndex) = -Waq - 0.5000*speciesSums

         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfu/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = gc_zero
         END IF
       END DO

       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter,' eqIndex:',comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1

     !--------------------------
     !-- Now the ACTIVYH2O components gets an equation in the "B" array        !
     ELSE IF(comps(rowCompIndex)%CompType == ACTIVYH2O) THEN

       DO colCompIndex = 1, nComps2Process                       ! k index
         IF(comps(colCompIndex)%CompType == MOLEBLNCE   .OR. &   ! MB
            comps(colCompIndex)%CompType == CHARGEBAL   ) THEN   ! H+

           ! dfm/dln(aj) in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *     &
               specs(speciesIndex)%Moles
           END DO
           jacobian(rowCounter,colCompIndex) = -0.017 * speciesSums


         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfm/du in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               !specs(speciesIndex)%Moles * &
               specs(speciesIndex)%delGamma

           END DO
           jacobian(rowCounter,colCompIndex) = 0.017 * speciesSums

         ELSE IF(comps(colCompIndex)%CompType == ACTIVYH2O ) THEN
           ! dfah2o/dln(ah2o) in PHREEQC notation

           jacobian(rowCounter,colCompIndex) = -Waq  *                         &
             antiLOG(comps(colCompIndex)%master%logActivity,REAL(10.0,GCHP))

         ELSE IF(comps(colCompIndex)%CompType == MASSHYDGN ) THEN
           ! dfah2o/dln(ae-) in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *     &
               specs(speciesIndex)%Moles
           END DO

           jacobian(rowCounter,colCompIndex) =  -0.017 * speciesSums

         ELSE IF(comps(colCompIndex)%CompType == MASSOXYGN ) THEN
           ! dfah2o/dln(Waq) in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
              IF(specs(speciesIndex)%Name /= TRIM('H2O') )THEN
                speciesSums = speciesSums +                                    &
                  specs(speciesIndex)%Moles
             END IF
           END DO
           jacobian(rowCounter,colCompIndex) = -(Waq *                         &
            (antiLOG(activityWater%master%logActivity,REAL(10.0,GCHP))-gc_one) &
                                                          + 0.017 * speciesSums)


         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfm/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = gc_zero
         END IF
       END DO

       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter,' eqIndex:',comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1

     !--------------------------
     !-- Now the CB components gets an equation in the "B" array               !
     ELSE IF(comps(rowCompIndex)%CompType == CHARGEBAL) THEN

       DO colCompIndex = 1, nComps2Process        ! k index
         IF(comps(colCompIndex)%CompType == MOLEBLNCE   .OR. &   ! MB
            comps(colCompIndex)%CompType == CHARGEBAL   .OR. &   ! H+
            comps(colCompIndex)%CompType == ACTIVYH2O) THEN      ! H2O

           ! dfm/dln(aj) in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums + &
               specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) * &
               specs(speciesIndex)%Moles * specs(speciesIndex)%charge
           END DO
           jacobian(rowCounter,colCompIndex) =  speciesSums


         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfm/du in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums + &
               specs(speciesIndex)%charge * &
               specs(speciesIndex)%delGamma

           END DO
           jacobian(rowCounter,colCompIndex) = -speciesSums

         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfm/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = gc_zero


         ELSE IF(comps(colCompIndex)%CompType == MASSHYDGN ) THEN
           ! dfm/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) = gc_zero


         ELSE IF(comps(colCompIndex)%CompType == MASSOXYGN ) THEN
           ! dfm/dnp in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums + &
               specs(speciesIndex)%charge * &
               specs(speciesIndex)%Moles

           END DO
           jacobian(rowCounter,colCompIndex) = speciesSums
         END IF
       END DO

       !##OPTIONAL CHECK:
       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter, ' eqIndex:',comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1

     !--------------------------
     ELSE IF(comps(rowCompIndex)%CompType == MASSHYDGN) THEN

       DO colCompIndex = 1, nComps2Process        ! k index
         IF(comps(colCompIndex)%CompType == MOLEBLNCE   .OR. &   ! MB
            comps(colCompIndex)%CompType == CHARGEBAL   .OR. &   ! H+
            comps(colCompIndex)%CompType == MASSHYDGN   .OR. &   ! MH
            comps(colCompIndex)%CompType == ACTIVYH2O ) THEN     ! H2O

           ! dfm/dln(aj) in PHREEQC notation
           speciesSums = gc_zero
            DO speciesIndex = 1,nSpecies
              speciesSums = speciesSums +                                      &
                specs(speciesIndex)%HStoich *                                  &
                specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *    &
                specs(speciesIndex)%Moles

            END DO
           jacobian(rowCounter,colCompIndex) = speciesSums


         ELSE IF(comps(colCompIndex)%CompType == MASSOXYGN ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
               speciesSums = speciesSums +                                     &
               specs(speciesIndex)%Hstoich *                                   &
               specs(speciesIndex)%Moles
           END DO
           jacobian(rowCounter,colCompIndex) = speciesSums


         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfm/du in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               specs(speciesIndex)%Hstoich *                                   &
               specs(speciesIndex)%delGamma

           END DO
           jacobian(rowCounter,colCompIndex) = -speciesSums

         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfm/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) =                                 &
                                -comps(colCompIndex)%ppData%stoich(rowCompIndex)
         END IF
       END DO

       !##OPTIONAL CHECK:
       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter,' eqIndex:',comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1

     !--------------------------
     ELSE IF(comps(rowCompIndex)%CompType == MASSOXYGN) THEN

       DO colCompIndex = 1, nComps2Process                       ! k index
         IF(comps(colCompIndex)%CompType == MOLEBLNCE   .OR. &   ! MB
            comps(colCompIndex)%CompType == MASSHYDGN   .OR. &   ! MH
            comps(colCompIndex)%CompType == CHARGEBAL   ) THEN   ! H+

           ! dfm/dln(aj) in PHREEQC notation
           speciesSums = gc_zero
            DO speciesIndex = 1,nSpecies
              speciesSums = speciesSums +                                      &
                specs(speciesIndex)%H2OStoich *                                &
                specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *    &
                specs(speciesIndex)%Moles

            END DO
           jacobian(rowCounter,colCompIndex) = speciesSums


         ELSE IF(comps(colCompIndex)%CompType == MASSOXYGN ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums + &
               specs(speciesIndex)%Moles

           END DO
           jacobian(rowCounter,colCompIndex) = speciesSums


         ELSE IF(comps(colCompIndex)%CompType == ACTIVYH2O ) THEN
           ! dfm/dlnah2o in PHREEQC notation
           speciesSums = gc_zero
            DO speciesIndex = 1,nSpecies
              IF(specs(speciesIndex)%Name /= TRIM('H2O')) THEN
               speciesSums = speciesSums +                                     &
                specs(speciesIndex)%H2OStoich *                                &
                specs(speciesIndex)%stoich(comps(colCompIndex)%CompIndex) *    &
                specs(speciesIndex)%Moles
              END IF
            END DO
           jacobian(rowCounter,colCompIndex) = speciesSums


         ELSE IF(comps(colCompIndex)%CompType == IONSTNGTH ) THEN
           ! dfm/du in PHREEQC notation

           speciesSums = gc_zero
           DO speciesIndex = 1,nSpecies
             speciesSums = speciesSums +                                       &
               specs(speciesIndex)%H2Ostoich *                                 &
               specs(speciesIndex)%delGamma

           END DO
           jacobian(rowCounter,colCompIndex) = -speciesSums

         ELSE IF(comps(colCompIndex)%CompType == PUREPHASE ) THEN
           ! dfm/dnp in PHREEQC notation
           jacobian(rowCounter,colCompIndex) =                                 &
                         -comps(colCompIndex)%ppData%stoich(rowCompIndex)
         END IF
       END DO

       IF(rowCounter /= comps(rowCompIndex)%eqIndex) THEN
         print *,'rowCounter /= comps(rowCompIndex)%eqIndex'
         print *,'rowCounter:',rowCounter, ' eqIndex:',                        &
                                              comps(rowCompIndex)%eqIndex
         STOP
       END IF

       rowCounter = rowCounter + 1

     END IF
   END DO


   !-- Consistency check: nOPTEqs+nEQLEqs = rowCounter
   IF(rowCounter-1 /= (nOPTEqs+nEQLEqs)) THEN
     print *,'rowCounter /= nOPTEqs+nEQLEqs after A&B matrix built in ',       &
             THIS_PROC
     print *,'rowCounter = ',rowCounter-1
     print *,'nOPTEqs+nEQLEqs = ',nOPTEqs+nEQLEqs
     print *,'jacobian:'
     DO rowCompIndex = 1,nComps2Process
       print *,'-> ',rowCompIndex,jacobian(rowCompIndex,:)
     END DO

     STOP
   END IF


 END FUNCTION BuildJacobian                                                    !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! SetComponentTotalConc                                                        !
!                                                                              !
! Routine to set total molalities of components prior to solution during the   !
! equilibration. Based on the water column or sed DICHM/PICHM array values     !
!------------------------------------------------------------------------------!
 SUBROUTINE SetComponentTotalConc(dissolvedConcs, particleConcs,               &
                                  components,     concMode)                    !
   !-- Incoming                                                                !
!  DOUBLETYPE,        DIMENSION(:) :: dissolvedConcs   ! WQ3F(i,DICHM(:))      !
   AED_REAL,        DIMENSION(:) :: dissolvedConcs   ! WQ3F(i,DICHM(:))      !
!  DOUBLETYPE,        DIMENSION(:) :: particleConcs    ! WQ3F(i,PICHM(:))      !
   AED_REAL,        DIMENSION(:) :: particleConcs    ! WQ3F(i,PICHM(:))      !
   TYPE (gcUnknowns), DIMENSION(:) :: components                               !
   INTEGER,           INTENT(IN)   :: concMode                                 !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER         :: THIS_PROC = "SetComponentTotalConc"      !
   INTEGER  :: componentIndex                           !
   DOUBLETYPE                      :: conversion                               !
                                                                               !

   DO componentIndex = 1,nComponents

     !-----------------------------------------------------------
     IF(components(componentIndex)%CompType == PUREPHASE) THEN

       !-- Number of pure-phase moles
       components(componentIndex)%ppData%moles = MAX(REAL(0.0,KIND=GCHP),      &
                              particleConcs(components(componentIndex)%wqIndex))

       !-- Natural Log of log K is the target saturation value
       components(componentIndex)%Total =                                      &
                    components(componentIndex)%ppData%logK  * LOG(REAL(10.0,DP))


     !-----------------------------------------------------------
     ELSE IF(components(componentIndex)%CompType == MOLEBLNCE) THEN

       !-- Total number of moles of component calculated from conc
       !-- Depending in input units (mol/L or mg/L)

       IF(concMode == MMOLPERL) THEN
         conversion = REAL(1e3,GCHP)
       ELSE IF(concMode == MMOLPERM3) THEN
         conversion = REAL(1e6,GCHP)
       ELSE IF(concMode == MGPERL) THEN
         conversion = components(componentIndex)%MolWeight * REAL(1e3,GCHP)
       END IF

       components(componentIndex)%Total = MAX( REAL(0.0,KIND=GCHP),            &
                  dissolvedConcs(components(componentIndex)%wqIndex) /         &
                  conversion)

       !-- First guess of master species concentration assuming gamma=1
       IF (components(componentIndex)%Total > gc_zero) THEN
         components(componentIndex)%master%logActivity =                       &
                                         LOG10(components(componentIndex)%Total)
       ELSE
         components(componentIndex)%master%logActivity = minLogActivity
       END IF

     !-----------------------------------------------------------
     ELSE IF(components(componentIndex)%CompType == CHARGEBAL) THEN

       !-- pH
       hydrogenIon%master%logActivity =                                        &
                            - dissolvedConcs(components(componentIndex)%wqIndex)

       !-- Moles of H+
       hydrogenIon%master%Moles =                                              &
                         antiLOG(hydrogenIon%master%logActivity,REAL(10.0,GCHP))

       !-- Target for CB eqn
       hydrogenIon%Total = hydrogenIon%master%Moles


     !-----------------------------------------------------------
     ELSE IF(components(componentIndex)%CompType == MASSHYDGN) THEN

       !-- pe (ln[ae-] is the master unknown for MH)
       !components(componentIndex)%master%logActivity =                         &
       !                     - dissolvedConcs(components(componentIndex)%wqIndex)
       components(componentIndex)%master%logActivity = -8.0


     !-----------------------------------------------------------
     ELSE IF(components(componentIndex)%CompType == MASSOXYGN) THEN

       !-- Number of moles of water per kilogram
       components(componentIndex)%Total =  Waq / molWgtH2O


     !-----------------------------------------------------------
     ELSE IF(components(componentIndex)%CompType == ACTIVYH2O) THEN

       components(componentIndex)%Total        = gc_zero
       components(componentIndex)%master%Moles = Waq / molWgtH2O

     END IF


     IF(verbosity > 5) THEN
       print *,'TOTAL: ',components(componentIndex)%CompName,                  &
                         components(componentIndex)%Total
     END IF

   END DO


 END SUBROUTINE SetComponentTotalConc
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! UpdateWQArray                                                                !
!                                                                              !
! Routine that takes output after equilibration and reassigns the values       !
! back to the main model values, after performing any necessary conversions    !
!------------------------------------------------------------------------------!
 SUBROUTINE UpdateWQArray(dissolvedConcs, particleConcs, components, concMode) !
   !-- Incoming                                                                !
!  DOUBLETYPE,       DIMENSION(:) :: dissolvedConcs  ! WQ3F(i,DICHM(:))       !
   AED_REAL,       DIMENSION(:) :: dissolvedConcs  ! WQ3F(i,DICHM(:))       !
!  DOUBLETYPE,       DIMENSION(:) :: particleConcs   ! WQ3F(i,PICHM(:))       !
   AED_REAL,       DIMENSION(:) :: particleConcs   ! WQ3F(i,PICHM(:))       !
   TYPE (gcUnknowns), DIMENSION(:) :: components                               !
   INTEGER,           INTENT(IN)   :: concMode                                 !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER         :: THIS_PROC = "UpdateWQArray"              !
   INTEGER  :: componentIndex                           !
   DOUBLETYPE                      :: conversion                               !


   DO componentIndex = 1,nComponents

     !-- Find the chemistry variables that are transported
     IF(components(componentIndex)%wqIndex /= 0) THEN

       !-- Particulate chemistry array
       IF(components(componentIndex)%CompType == PUREPHASE) THEN
         particleConcs(components(componentIndex)%wqIndex) =                   &
                                         components(componentIndex)%ppData%moles

       !-- Normal dissolved component
       ELSE IF(components(componentIndex)%CompType == MOLEBLNCE) THEN
         IF(concMode == MMOLPERL) THEN
           conversion = REAL(1e3,GCHP)
         ELSE IF(concMode == MMOLPERM3) THEN
           conversion = REAL(1e3,GCHP)
         ELSE IF(concMode == MGPERL) THEN
           conversion = (components(componentIndex)%MolWeight*REAL(1e3,GCHP))
         END IF

         dissolvedConcs(components(componentIndex)%wqIndex) =                  &
                                   components(componentIndex)%Total * conversion


       !-- pH
       ELSE IF(components(componentIndex)%CompType == CHARGEBAL) THEN
         dissolvedConcs(components(componentIndex)%wqIndex) =                  &
                                 - components(componentIndex)%master%logActivity

       !-- pe
       ELSE IF(components(componentIndex)%CompType == MASSHYDGN) THEN
         dissolvedConcs(components(componentIndex)%wqIndex) =                  &
                                 - components(componentIndex)%master%logActivity

       END IF
     END IF
   END DO

 END SUBROUTINE UpdateWQArray
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 FUNCTION GetDissChemNames(nDICHM, comps) RESULT(listOfNames)                  !
   !-- Incoming
   INTEGER,                        INTENT(IN) :: nDICHM
   TYPE(gcUnknowns), DIMENSION(:), INTENT(IN) :: comps
   !-- Outgoing
   CHARACTER(LEN=idl), DIMENSION(1:nDICHM) :: listOfNames
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER :: THIS_PROC = "GetDissChemNames"                   !
   INTEGER  :: i

!   !-- Compulsory DICHM values
!   listOfNames(1) = "DO"
!   listOfNames(2) = "PO4"
!   listOfNames(3) = "NH4"
!   listOfNames(4) = "NO3"
!   IF(simSiO2) THEN
!     listOfNames(5) = "SiO2"
!   END IF

   !-- Extra values
   DO i = 1,SIZE(comps)
     IF(comps(i)%wqIndex>0) THEN
       listOfNames(comps(i)%wqIndex) = TRIM(comps(i)%EltName)
     END IF
   END DO
   listOfNames(chargeBalCol) = "ubalchg"


 END FUNCTION GetDissChemNames                                                 !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 FUNCTION GetPPChemNames(nPICHM, comps) RESULT(listOfNames)                  !
   !-- Incoming
   INTEGER,                        INTENT(IN) :: nPICHM
   TYPE(gcUnknowns), DIMENSION(:), INTENT(IN) :: comps
   !-- Outgoing
   CHARACTER(LEN=idl), DIMENSION(1:nPICHM) :: listOfNames
   !-- Local
   CHARACTER(*), PARAMETER :: THIS_PROC = "GetPPChemNames"
   INTEGER  :: i

   DO i = 1,SIZE(comps)
     IF(comps(i)%wqIndex>0) THEN
       listOfNames(comps(i)%wqIndex) = TRIM(comps(i)%EltName)
     END IF
   END DO

 END FUNCTION GetPPChemNames                                                 !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 FUNCTION GetDissChemIndex(Name, wqColCounter) RESULT(Index)                 !
   !-- Incoming
   CHARACTER(LEN=idl), INTENT(IN)     :: Name
   INTEGER,            INTENT(INOUT)  :: wqColCounter
   !-- Outgoing
   INTEGER  :: Index
   !-- Local
   CHARACTER(*), PARAMETER :: THIS_PROC = "GetDissChemIndex"

     wqColCounter = wqColCounter + 1
     Index    = wqColCounter

   IF(Index == 0) THEN
     print *,' Error in: ', THIS_PROC
     print *,' DICHM column = ', Index
     STOP
   END IF

 END FUNCTION GetDissChemIndex                                                   !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 FUNCTION GetPartChemIndex(wqColCounter) RESULT(Index)                 !
   !-- Incoming
   INTEGER,            INTENT(INOUT)  :: wqColCounter
   !-- Outgoing
   INTEGER  :: Index
   !-- Local
   CHARACTER(*), PARAMETER :: THIS_PROC = "GetPartChemIndex"

   wqColCounter =  wqColCounter + 1
   Index    =  wqColCounter

   IF(Index == 0) THEN
     print *,'Error in: ', THIS_PROC
     print *,'PICHM column = ', Index
     STOP
   END IF

 END FUNCTION GetPartChemIndex                                                   !
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 SUBROUTINE createSpeciesSubSetforSim(comps)
    IMPLICIT NONE
    !-- Incoming
    TYPE (gcUnknowns), DIMENSION(:) :: comps

    ! ------ Local Constants --------------------------
    CHARACTER(*), PARAMETER ::  THIS_PROC = "createSpeciesSubSetforSim"

    ! ------ Local Variables --------------------------
    INTEGER  :: i,ii,megaCompIndex     ,allocStatus   !loop counter
    INTEGER  :: nComp, nSpec
    TYPE (gcSpecies) :: dummySpecies
    CHARACTER(LEN=16) :: dummyName
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: bufferList
    INTEGER, DIMENSION(:), POINTER :: simSpeciesList
    LOGICAL, DIMENSION(:), ALLOCATABLE :: simulatedStatus
    INTEGER, DIMENSION(1) :: compLoc
    nComp = SIZE(comps)

    !--------------------------------------------------
    ALLOCATE(bufferList(NumOfSpecies()))
    bufferList = 0

    ALLOCATE(simulatedStatus(NumOfComps()))
    simulatedStatus = .FALSE.

    ALLOCATE(dummySpecies%stoich(NumOfComps()))


    DO i = 1,NumOfSpecies()

      CALL GetSpeciesName(i, dummyName)
      dummySpecies%Name = TRIM(dummyName)
      CALL GetSpeciesInfo(i, dummySpecies)


      simulatedStatus = .FALSE.
      DO megaCompIndex = 1, NumOfComps()
        IF(dummySpecies%stoich(megaCompIndex) == 0) THEN
          simulatedStatus(megaCompIndex) = .TRUE.

        ELSE

          CALL GetCompCompName(megaCompIndex, dummyName)
          !!## There exists optional 'success' 3rd argument.

          DO ii = 1,nComp
            IF( TRIM(dummyName) == TRIM(comps(ii)%CompName)) THEN

                simulatedStatus(megaCompIndex) = .TRUE.
                EXIT

            ENDIF
          END DO
        END IF
      END DO
      IF(verbosity >10) THEN
        print *,'speciesList: ',dummySpecies%Name, ANY(.NOT.simulatedStatus)
      END IF
      IF (.NOT. ANY(.NOT.simulatedStatus)) THEN
        bufferList(i) = i
      END IF

    END DO


    nSpec = COUNT(bufferList /= 0)
    ALLOCATE(allSpecies(nSpec),STAT=allocStatus)
    DO i = 1,nSpec
     ALLOCATE(allSpecies(i)%stoich(nComponents),STAT=allocStatus)
      allSpecies(i)%stoich = gc_zero
    END DO

    bufferList(1:nSpec) = PACK(bufferList, MASK = (bufferList /= 0))
    bufferList(nSpec+1:SIZE(bufferList)) =  0
    simSpeciesList => bufferList(1:nSpec)


    !-- Now we know which of the possible species are to be included
    !-- in the simulation. Now populate the local allSpecies vector
    DO i = 1,nSpec

      !-- Name of actual simulated species and extract necessary info
      CALL GetSpeciesName(simSpeciesList(i), dummyName)
      dummySpecies%Name = TRIM(dummyName)
      CALL GetSpeciesInfo(simSpeciesList(i), dummySpecies)

      !-- Set local species info
      allSpecies(i)%SpeciesIndex = i
      allSpecies(i)%Name         = dummySpecies%Name
      allSpecies(i)%Charge       = dummySpecies%Charge
      allSpecies(i)%logK         = dummySpecies%logKat25
      allSpecies(i)%logKat25     = dummySpecies%logKat25
      allSpecies(i)%deltaH       = dummySpecies%deltaH
      allSpecies(i)%logActivity  = gc_zero
      allSpecies(i)%logGamma     = gc_zero
      allSpecies(i)%delGamma     = gc_zero
      allSpecies(i)%Moles        = gc_zero
      allSpecies(i)%HStoich      = gc_zero
      allSpecies(i)%H2OStoich    = gc_zero

      !-- Populate the condensed stoich array for the simulated species
      !-- using data from the species store and the MegaComponentList
      DO megaCompIndex = 1, NumOfComps()
        IF(dummySpecies%stoich(megaCompIndex) == 0) THEN
          CYCLE
        ELSE

          CALL GetCompCompName(megaCompIndex, dummyName)  !!##GEOCH_MOD

          DO ii = 1,nComp
            IF( TRIM(dummyName) == TRIM(comps(ii)%CompName) ) THEN

                allSpecies(i)%Stoich(ii) = dummySpecies%Stoich(megaCompIndex)

                IF(TRIM(dummyName) == "H2O") THEN
                  allSpecies(i)%H2OStoich = dummySpecies%Stoich(megaCompIndex)
                ELSE IF(TRIM(dummyName) == "H+") THEN
                  allSpecies(i)%HStoich = dummySpecies%Stoich(megaCompIndex)
                END IF

                EXIT

            ENDIF
          END DO
        END IF
      END DO


      !-- This section sets the Activity Coefficient
      !-- calculation method (gflag)
      IF(allSpecies(i)%Charge == 0) THEN
        allSpecies(i)%Gflag = UNCHGD
      ELSE
        allSpecies(i)%Gflag = DAVIES
      END IF
      IF(TRIM(allSpecies(i)%Name) == "H+") THEN
        allSpecies(i)%Gflag = WATEQDH
      ELSE IF (TRIM(allSpecies(i)%Name) == "H2O") THEN
        allSpecies(i)%Gflag = 3
      ELSE IF (TRIM(allSpecies(i)%Name) == "e-")  THEN
        allSpecies(i)%Gflag = 3
      END IF

      !-- We now set the pointers between components and their
      !-- relevant free species data.
      IF(COUNT(allSpecies(i)%Stoich /= gc_zero) == 1) THEN
        compLoc = MAXLOC(ABS(allSpecies(i)%Stoich))
        ii = compLoc(1)
        IF(comps(ii)%CompType == MOLEBLNCE) THEN
          comps(ii)%master => allSpecies(i)
        ELSE IF(comps(ii)%CompType == ACTIVYH2O) THEN
          activityWater%master    => allSpecies(i)
        ELSE IF(comps(ii)%CompType == CHARGEBAL) THEN
          hydrogenIon%master      => allSpecies(i)
        ELSE IF(comps(ii)%CompType == MASSHYDGN) THEN
          massHydrogen%master     => allSpecies(i)
        END IF
      END IF

    END DO

    nSpecies = nSpec

    !-- Final check for master species pointers:
    DO ii = 1,nComp
       IF( comps(ii)%CompType == MOLEBLNCE ) THEN
         IF(.NOT. ASSOCIATED(comps(ii)%Master) ) THEN
           !-- Special case for redox master species
           IF(TRIM(comps(ii)%CompName) == "Fe+2") THEN
             DO i = 1,nSpec
               IF(TRIM(allSpecies(i)%Name) == "Fe+2") THEN
                 comps(ii)%master => allSpecies(i)
               END IF
             END DO
           ELSE IF(TRIM(comps(ii)%CompName) == "H2S") THEN
             DO i = 1,nSpec
               IF(TRIM(allSpecies(i)%Name) == "H2S") THEN
                 comps(ii)%master => allSpecies(i)
               END IF
             END DO
           ELSE IF(TRIM(comps(ii)%CompName) == "CH4") THEN
             DO i = 1,nSpec
               IF(TRIM(allSpecies(i)%Name) == "CH4") THEN
                 comps(ii)%master => allSpecies(i)
               END IF
             END DO
           ELSE IF(TRIM(comps(ii)%CompName) == "NH4") THEN
             DO i = 1,nSpec
               IF(TRIM(allSpecies(i)%Name) == "NH4") THEN
                 comps(ii)%master => allSpecies(i)
               END IF
             END DO
           ELSE
             WRITE(*,'(1X,A,/,1X,2X,"Free species conc for",A)')               &
                                             THIS_PROC, TRIM(comps(ii)%CompName)
             WRITE(*,'(1X,4X,"is not associated with its component")')
             WRITE(*,'(/,1X,2X,"Check it has not been requested twice in the con file")')
             WRITE(*,'(1X,2X,"STOPPING")')
             STOP
           END IF
         END IF
       END IF
    END DO
    !--------------------------------------------------

 END SUBROUTINE createSpeciesSubSetforSim
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 SUBROUTINE setupWCRedoxConfiguration(comps)
    !-- Incoming
    TYPE (gcUnknowns), DIMENSION(:) :: comps
    !-- Local
    INTEGER  :: nComp, cIndex

    nComp = SIZE(comps)

    !--------------------------------------------------
    !-- First search for IRON redox couple
    FEII  = -999
    FEIII = -999
    DO cIndex = 1, nComp

      !-- Search for FeII
      IF( TRIM(comps(cIndex)%EltName) == "FEII" .OR. &
          TRIM(comps(cIndex)%EltName) == "FeII" ) THEN
           FEII = comps(cIndex)%wqIndex
           CYCLE
      END IF
      !-- Search for FeIII
      IF( TRIM(comps(cIndex)%EltName) == "FEIII" .OR. &
          TRIM(comps(cIndex)%EltName) == "FeIII" ) THEN
           FEIII = comps(cIndex)%wqIndex
           EXIT
      END IF

    END DO
    IF (FEII > 0 .AND. FEIII > 0) THEN
      !-- Both redox states simulated, so include redox connection
      simIronRedox = .TRUE.
    ELSE
      !-- No redox rate transormations
      simIronRedox = .FALSE.
    END IF



    !--------------------------------------------------
    !-- Next search for MANGANESE redox couple
    MNII  = -999
    MNVII = -999
    DO cIndex = 1, nComp

      !-- Search for MnII
      IF( TRIM(comps(cIndex)%EltName) == "MNII" .OR. &
          TRIM(comps(cIndex)%EltName) == "MnII" ) THEN

           MNII = comps(cIndex)%wqIndex
           CYCLE
      END IF
      !-- Search for MnIV
      IF( TRIM(comps(cIndex)%EltName) == "MNIV" .OR. &
          TRIM(comps(cIndex)%EltName) == "MnIV" ) THEN

           MNIV = comps(cIndex)%wqIndex
           EXIT
      END IF

      !-- Search for MnVII
      IF( TRIM(comps(cIndex)%EltName) == "MNVII" .OR. &
          TRIM(comps(cIndex)%EltName) == "MnVII" ) THEN

           MNVII = comps(cIndex)%wqIndex
           EXIT
      END IF

    END DO
    IF (MNII > 0 .AND. MNIV > 0) THEN
      !-- Both redox states simulated, so include redox connection
      simManganRedox = .TRUE.
    ELSE
      !-- No redox rate transormations
      simManganRedox = .FALSE.
    END IF


    !--------------------------------------------------
    !-- Next search for SULFUR redox couple
    H2S  = -999
    SO4  = -999
    DO cIndex = 1, nComp

      !-- Search for H2S
      IF( TRIM(comps(cIndex)%EltName) == "H2S" ) THEN
           H2S = comps(cIndex)%wqIndex
           CYCLE
      END IF
      !-- Search for SO4
      IF( TRIM(comps(cIndex)%EltName) == "SO4" ) THEN
           SO4 = comps(cIndex)%wqIndex
           EXIT
      END IF

    END DO
    IF (H2S > 0 .AND. SO4 > 0) THEN
      !-- Both redox states simulated, so include redox connection
      simSulfurRedox = .TRUE.
    ELSE
      !-- No redox rate transormations
      simSulfurRedox = .FALSE.
    END IF

    !-- Next search for ALUMINIUM (not redox but DSF)
    ALIII = -999
    DO cIndex = 1, nComp

      !-- Search for Al+3
      IF( TRIM(comps(cIndex)%EltName) == "Al" .OR. &
          TRIM(comps(cIndex)%EltName) == "AL" ) THEN

           ALIII = comps(cIndex)%wqIndex
           EXIT
      END IF
    END DO


    !--------------------------------------------------
    !-- Next search for ZINC (not redox but DSF)
    ZNII = -999
    DO cIndex = 1, nComp

      !-- Search for Zn+2
      IF( TRIM(comps(cIndex)%EltName) == "Zn" .OR. &
          TRIM(comps(cIndex)%EltName) == "ZN" ) THEN

           ZNII = comps(cIndex)%wqIndex
           EXIT
      END IF
    END DO


    !--------------------------------------------------
    !-- Next search for CADMIUM (not redox but DSF)
    CDII = -999
    DO cIndex = 1, nComp

      !-- Search for Cd+2
      IF( TRIM(comps(cIndex)%EltName) == "Cd" .OR. &
          TRIM(comps(cIndex)%EltName) == "CD" ) THEN

           CDII = comps(cIndex)%wqIndex
           EXIT
      END IF
    END DO


    !--------------------------------------------------
    !-- Next search for LEAD (not redox but DSF)
    PBII = -999
    DO cIndex = 1, nComp

      !-- Search for Pb+2
      IF( TRIM(comps(cIndex)%EltName) == "Pb" .OR. &
          TRIM(comps(cIndex)%EltName) == "PB" ) THEN

           PBII = comps(cIndex)%wqIndex
           EXIT
      END IF
    END DO


    !--------------------------------------------------
    !-- Next search for NICKEL (not redox but DSF)
    NiII = -999
    DO cIndex = 1, nComp

      !-- Search for Pb+2
      IF( TRIM(comps(cIndex)%EltName) == "Ni" .OR. &
          TRIM(comps(cIndex)%EltName) == "NI" ) THEN

           NiII = comps(cIndex)%wqIndex
           EXIT
      END IF
    END DO


    !--------------------------------------------------
    !-- Next search for ARSENIC redox couple
    ASIII  = -999
    ASV    = -999
    DO cIndex = 1, nComp


      !-- Search for AsIII
      IF( TRIM(comps(cIndex)%EltName) == "ASIII" .OR. &
          TRIM(comps(cIndex)%EltName) == "AsIII" ) THEN

           ASIII = comps(cIndex)%wqIndex
           CYCLE
      END IF
      !-- Search for AsV
      IF( TRIM(comps(cIndex)%EltName) == "AsV" .OR. &
          TRIM(comps(cIndex)%EltName) == "ASV" ) THEN

           ASV = comps(cIndex)%wqIndex
      END IF

    END DO
    IF (ASIII > 0 .AND. ASV > 0) THEN
      !-- Both redox states simulated, so include redox connection
      simArsenicRedox = .TRUE.
    ELSE
      !-- No redox rate transormations
      simArsenicRedox = .FALSE.
    END IF

    !--------------------------------------------------
    !-- Next search for CARBON redox couple
    CH4  = -999
    CO3  = -999
    DO cIndex = 1, nComp


      !-- Search for CH4
      IF( TRIM(comps(cIndex)%EltName) == "CH4") THEN

           CH4 = comps(cIndex)%wqIndex
           CYCLE
      END IF
      !-- Search for CO3
      IF( TRIM(comps(cIndex)%EltName) == "CO3") THEN

           CO3 = comps(cIndex)%wqIndex
      END IF

    END DO
    IF (CO3 > 0 .AND. CH4 > 0) THEN
      !-- Both redox states simulated, so include redox connection
      simCarbonRedox = .TRUE.
    ELSE
      !-- No redox rate transormations
      simCarbonRedox = .FALSE.
    END IF

   !--------------------------------------------------



 END SUBROUTINE setupWCRedoxConfiguration
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 FUNCTION setupPhaseStoich(inComp, comps) RESULT(outComp)                      !
    IMPLICIT NONE                                                              !
    !-- Incoming                                                               !
    TYPE (gcUnknowns), DIMENSION(:) :: comps                                   !
    TYPE (gcUnknowns), INTENT(IN)   :: inComp                                  !
    !-- Outgoing                                                               !
    TYPE (gcUnknowns) :: outComp                                               !
    !-- Local Constants --------------------------
    CHARACTER(*), PARAMETER ::  THIS_PROC = "setupPhaseStoich"
    !-- Local Variables --------------------------
    INTEGER  :: ii,megaCompIndex,allocStatus,nComp
    CHARACTER(LEN=16) :: dummyName

    nComp = SIZE(comps)


    !--------------------------------------------------
    ALLOCATE(outComp%PPData,STAT=allocStatus)
    CheckAllocStatus(allocStatus,THIS_PROC,"outComp%PPData")

    ALLOCATE(outComp%PPData%stoich(nComponents),STAT=allocStatus)
    CheckAllocStatus(allocStatus,THIS_PROC,"outComp%PPData%Stoich")
    outComp%PPData%stoich = gc_zero

    outComp%EltName   =  inComp%EltName
    outComp%CompName  =  inComp%CompName
    outComp%CompType  =  inComp%CompType
    outComp%CompIndex =  inComp%CompIndex
    outComp%MolWeight =  inComp%MolWeight
    outComp%Charge    =  inComp%Charge

    outComp%PPData%logK     =  inComp%PPData%logKat25
    outComp%PPData%logKat25 =  inComp%PPData%logKat25
    outComp%PPData%deltaH   =  inComp%PPData%deltaH
    outComp%PPData%Diameter =  inComp%PPData%Diameter
    outComp%PPData%Moles    =  inComp%PPData%Moles

    !-- Populate the condensed stoich array for the simulated phases
    !-- using data from the component store and the MegaComponentList
    DO megaCompIndex = 1, NumOfComps()

      IF(inComp%PPData%Stoich(megaCompIndex) == 0) THEN
        CYCLE
      ELSE

        CALL GetCompCompName(megaCompIndex, dummyName)  !!##GEOCH_MOD

        DO ii = 1,nComp
          IF( TRIM(dummyName) == TRIM(comps(ii)%CompName) ) THEN
            !-- Set local stoich value to the appropriate value
            outComp%PPData%Stoich(ii)=inComp%PPData%stoich(megaCompIndex)

            !-- Special case for MASSHYDGN and MASSOXYGN
            IF(TRIM(dummyName) == "H2O") THEN
              outComp%PPData%Stoich(massWater%CompIndex) =              &
                                      inComp%PPData%stoich(megaCompIndex)
            ELSE IF(TRIM(dummyName) == "H+") THEN
              outComp%PPData%Stoich(massHydrogen%CompIndex) =           &
                                      inComp%PPData%stoich(megaCompIndex)
            END IF

            !-- Finally, the master species stoich arrays must be updated
            !-- based on the pure-phase composition.
            comps(ii)%Master%Stoich(inComp%CompIndex) =                 &
                                      inComp%PPData%stoich(megaCompIndex)

            EXIT

          ENDIF
        END DO
      END IF

    END DO
    !--------------------------------------------------


 END FUNCTION setupPhaseStoich
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
! reportGeochemConfig                                                          !
!                                                                              !
!   Generic subroutine for writing a brief summary of the geochemistry         !
!   configuration of the simulation. Can be used to write to screen            !
!   (set lun=6) or to any file, using its logical unit number (lun)            !
!                                                                              !
!------------------------------------------------------------------------------!
 SUBROUTINE reportGeochemConfig(lun)                                           !
    IMPLICIT NONE                                                              !
    !-- Incoming                                                               !
    INTEGER, INTENT(IN) :: lun        ! Logical Unit Number to write to        !
                                                                               !
    ! ------ Local Constants --------------------------                        !
    CHARACTER(*), PARAMETER ::  THIS_PROC = "reportGeochemConfig"              !
                                                                               !
    ! ------ Local Variables --------------------------                        !
    INTEGER  :: i        !loop counter                                          !

    WRITE(lun,"(/,'-----------------------------------------------')")
    WRITE(lun,"('GEOCHEMISTRY CONFIGURATION')")

    !--------------------------------------------------
    WRITE(lun,"(/,1X,'Prinicipal components being simulated:')")
    WRITE(lun,"(5X,15X,'Name        Charge  MolWeight(g/mol)')")
    DO i = 1, nComponents
      IF(allComponents(i)%CompType == MOLEBLNCE) THEN
      WRITE(lun,"(5X,'#',I3,':   ',A10,' (',A5,')',I8,'  ',F8.3)")i,    &
                           TRIM(allComponents(i)%CompName),             &
                           TRIM(allComponents(i)%EltName),              &
                           allComponents(i)%Charge,                     &
                           allComponents(i)%MolWeight
      END IF
    END DO
    !--------------------------------------------------
    WRITE(lun,"(/,1X,'Compulsory components being simulated:')")
    WRITE(lun,"(5X,15X,'Name')")
    DO i = 1, nComponents
      IF(allComponents(i)%CompType /= MOLEBLNCE .AND.                   &
                          allComponents(i)%CompType /= PUREPHASE) THEN
      WRITE(lun,"(5X,'#',I3,':   ',A10,' (',A10,')')")i,                &
                           TRIM(allComponents(i)%CompName),             &
                           TRIM(allComponents(i)%EltName)
      END IF
    END DO
    !--------------------------------------------------
!    IF(SIZE(PICHM)>0) THEN
      WRITE(lun,"(/,1X,'Pure phases being simulated:')")
      WRITE(lun,"(5X,15X,'Name             logK(25C) deltaH(kj/mol) settlingVel(m/s)')")
      DO i = 1, nComponents
        IF(allComponents(i)%CompType == PUREPHASE) THEN
          WRITE(lun,"(5X,'#',I3,':   ',A10,' (',A10,')','  ',3(3X,F8.3))")    &
                                         i,                                   &
                                         TRIM(allComponents(i)%CompName),     &
                                         TRIM(allComponents(i)%EltName),      &
                                         allComponents(i)%ppData%logKat25,    &
                                         allComponents(i)%ppData%deltaH,      &
                                         allComponents(i)%ppData%SettlingVel
        END IF
      END DO
!    ELSE
!      WRITE(lun,"(/,'Pure phases being simulated:')")
!      WRITE(lun,"(5X,'None')")
!    END IF
    !--------------------------------------------------
    WRITE(lun,"(/,1X,'Aqueous species being simulated (dynamically ')")
    WRITE(lun,"('      determined from available components): ')")
    WRITE(lun,"(5X,10X,'   Name Charge  logK(25C)  deltaH(kj/mol)')")
    DO i = 1, nSpecies
      WRITE(lun,"(5X,'#',I3,':','  ',A10,I8,F8.3,F10.3)")i, &
                           TRIM(allSpecies(i)%Name), &
                           allSpecies(i)%Charge,     &
                           allSpecies(i)%logKat25,   &
                           allSpecies(i)%deltaH
    END DO


 END SUBROUTINE reportGeochemConfig
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 SUBROUTINE outputComponentData(comp)
    IMPLICIT NONE
    !-- Incoming
    TYPE (gcUnknowns), DIMENSION(:) :: comp

    ! ------ Local Constants --------------------------
    CHARACTER(*), PARAMETER ::  THIS_PROC = "outputComponentData"

    ! ------ Local Variables --------------------------
    INTEGER  :: i        !loop counter
    INTEGER  :: nComp

    nComp = SIZE(comp)

    !--------------------------------------------------
    DO i = 1, nComp
      print *,'----> :     ',i
      print *,' CompName:  ',comp(i)%CompName
      print *,' EltName:   ',comp(i)%EltName
      print *,' CompType:  ',comp(i)%CompType
      print *,' CompIndex: ',comp(i)%CompIndex
      print *,' eqIndex:   ',comp(i)%eqIndex
      print *,' wqIndex:   ',comp(i)%wqIndex
      print *,' delta:     ',comp(i)%delta
      print *,' Value:     ',comp(i)%Value
      print *,' Total:     ',comp(i)%Total
      IF(comp(i)%CompType==PUREPHASE) THEN
        print *,' ppData%logK  :  ',comp(i)%ppData%logK
        print *,' ppData%moles :  ',comp(i)%ppData%moles
        print *,' ppData%stoich:  ',comp(i)%ppData%stoich
      END IF
      print *,'  --<'
    END DO
    !--------------------------------------------------

 END SUBROUTINE outputComponentData
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 SUBROUTINE outputComponentSummary(comp)
    IMPLICIT NONE
    !-- Incoming
    TYPE (gcUnknowns), DIMENSION(:) :: comp

    ! ------ Local Constants --------------------------
    CHARACTER(*), PARAMETER ::  THIS_PROC = "outputComponentData"

    ! ------ Local Variables --------------------------
    INTEGER  :: i        !loop counter
    INTEGER  :: nComp

    nComp = SIZE(comp)

    !--------------------------------------------------
    DO i = 1, nComp
      print *,'----> :     ',i
      print *,' CompName:  ',TRIM(comp(i)%CompName),'  Total:     ',comp(i)%Total
      print *,'  --<'
    END DO
    !--------------------------------------------------

 END SUBROUTINE outputComponentSummary
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 SUBROUTINE outputSpeciesData(spec)
    IMPLICIT NONE
    !-- Incoming
    TYPE (gcSpecies), DIMENSION(:) :: spec

    ! ------ Local Constants --------------------------
    CHARACTER(*), PARAMETER ::  THIS_PROC = "outputSpeciesData"

    ! ------ Local Variables --------------------------
    INTEGER  :: i,nSpec

    nSpec = SIZE(spec)
    !--------------------------------------------------
    DO i = 1, nSpec
      print *,'----> :    '
      print *,' SpeciesIndex: ',spec(i)%SpeciesIndex
      print *,' Name:     ',spec(i)%Name
      print *,' Charge:   ',spec(i)%Charge
      print *,' logK:     ',spec(i)%logK
      print *,' deltaH:   ',spec(i)%deltaH
      print *,' Gflag:    ',spec(i)%Gflag
      print *,' H2OStoich:',spec(i)%H2OStoich
      print *,' HStoich:  ',spec(i)%HStoich
      print *,' Moles:    ',spec(i)%Moles
      print *,' log(act): ',spec(i)%logActivity
      print *,' loggamma: ',spec(i)%logGamma
      print *,' delgamma: ',spec(i)%delGamma
      print *,' stoich:   ',INT(spec(i)%stoich)
      print *,'  --<'
    END DO
    !--------------------------------------------------

 END SUBROUTINE outputSpeciesData
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! SimplexAlgorithm                                                             !
!                                                                              !
! THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX                           !
! METHOD OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION                     !
! TO A K BY N SYSTEM OF LINEAR EQUATIONS                                       !
!             AX=B                                                             !
! SUBJECT TO L LINEAR EQUALITY CONSTRAINTS                                     !
!             CX=D                                                             !
! AND M LINEAR INEQUALITY CONSTRAINTS                                          !
!             EX.LE.F.                                                         !
! DESCRIPTION OF PARAMETERS                                                    !
! K         NUMBER OF ROWS OF THE MATRIX A (K.GE.1).                           !
! L         NUMBER OF ROWS OF THE MATRIX C (L.GE.0).                           !
! M         NUMBER OF ROWS OF THE MATRIX E (M.GE.0).                           !
! N         NUMBER OF COLUMNS OF THE MATRICES A,C,E (N.GE.1).                  !
! KLMD      SET TO AT LEAST K+L+M FOR ADJUSTABLE DIMENSIONS.                   !
! KLM2D     SET TO AT LEAST K+L+M+2 FOR ADJUSTABLE DIMENSIONS.                 !
! NKLMD     SET TO AT LEAST N+K+L+M FOR ADJUSTABLE DIMENSIONS.                 !
! N2D       SET TO AT LEAST N+2 FOR ADJUSTABLE DIMENSIONS                      !
! Q         TWO DIMENSIONAL REAL ARRAY WITH KLM2D ROWS AND                     !
!           AT LEAST N2D COLUMNS.                                              !
!           ON ENTRY THE MATRICES A,C AND E, AND THE VECTORS                   !
!           B,D AND F MUST BE STORED IN THE FIRST K+L+M ROWS                   !
!           AND N+1 COLUMNS OF Q AS FOLLOWS                                    !
!                A B                                                           !
!            Q = C D                                                           !
!                E F                                                           !
!           THESE VALUES ARE DESTROYED BY THE SUBROUTINE.                      !
! KODE      A CODE USED ON ENTRY TO, AND EXIT                                  !
!           FROM, THE SUBROUTINE.                                              !
!           ON ENTRY, THIS SHOULD NORMALLY BE SET TO 0.                        !
!           HOWEVER, IF CERTAIN NONNEGATIVITY CONSTRAINTS                      !
!           ARE TO BE INCLUDED IMPLICITLY, RATHER THAN                         !
!           EXPLICITLY IN THE CONSTRAINTS EX.LE.F, THEN KODE                   !
!           SHOULD BE SET TO 1, AND THE NONNEGATIVITY                          !
!           CONSTRAINTS INCLUDED IN THE ARRAYS X AND                           !
!           RES (SEE BELOW).                                                   !
!           ON EXIT, KODE HAS ONE OF THE                                       !
!           FOLLOWING VALUES                                                   !
!                0- OPTIMAL SOLUTION FOUND,                                    !
!                1- NO FEASIBLE SOLUTION TO THE                                !
!                   CONSTRAINTS,                                               !
!                2- CALCULATIONS TERMINATED                                    !
!                   PREMATURELY DUE TO ROUNDING ERRORS,                        !
!                3- MAXIMUM NUMBER OF ITERATIONS REACHED.                      !
! TOLER     A SMALL POSITIVE TOLERANCE. EMPIRICAL                              !
!           EVIDENCE SUGGESTS TOLER = 10**(-D*2/3),                            !
!           WHERE D REPRESENTS THE NUMBER OF DECIMAL                           !
!           DIGITS OF ACCURACY AVAILABLE. ESSENTIALLY,                         !
!           THE SUBROUTINE CANNOT DISTINGUISH BETWEEN ZERO                     !
!           AND ANY QUANTITY WHOSE MAGNITUDE DOES NOT EXCEED                   !
!           TOLER. IN PARTICULAR, IT WILL NOT PIVOT ON ANY                     !
!           NUMBER WHOSE MAGNITUDE DOES NOT EXCEED TOLER.                      !
! ITER      ON ENTRY ITER MUST CONTAIN AN UPPER BOUND ON                       !
!           THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.                          !
!           A SUGGESTED VALUE IS 10*(K+L+M). ON EXIT ITER                      !
!           GIVES THE NUMBER OF SIMPLEX ITERATIONS.                            !
! X         ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST N2D.                   !
!           ON EXIT THIS ARRAY CONTAINS A                                      !
!           SOLUTION TO THE L1 PROBLEM. IF KODE=1                              !
!           ON ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE                       !
!           SIMPLE NONNEGATIVITY CONSTRAINTS ON THE                            !
!           VARIABLES. THE VALUES -1, 0, OR 1                                  !
!           FOR X(J) INDICATE THAT THE J-TH VARIABLE                           !
!           IS RESTRICTED TO BE .LE.0, UNRESTRICTED,                           !
!           OR .GE.0 RESPECTIVELY.                                             !
! RES       ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST KLMD.                  !
!           ON EXIT THIS CONTAINS THE RESIDUALS B-AX                           !
!           IN THE FIRST K COMPONENTS, D-CX IN THE                             !
!           NEXT L COMPONENTS (THESE WILL BE =0),AND                           !
!           F-EX IN THE NEXT M COMPONENTS. IF KODE=1 ON                        !
!           ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE SIMPLE                   !
!           NONNEGATIVITY CONSTRAINTS ON THE RESIDUALS                         !
!           B-AX. THE VALUES -1, 0, OR 1 FOR RES(I)                            !
!           INDICATE THAT THE I-TH RESIDUAL (1.LE.I.LE.K) IS                   !
!           RESTRICTED TO BE .LE.0, UNRESTRICTED, OR .GE.0                     !
!           RESPECTIVELY.                                                      !
! ERROR     ON EXIT, THIS GIVES THE MINIMUM SUM OF                             !
!           ABSOLUTE VALUES OF THE RESIDUALS.                                  !
! CU        A TWO DIMENSIONAL REAL ARRAY WITH TWO ROWS AND                     !
!           AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.                         !
! IU        A TWO DIMENSIONAL INTEGER ARRAY WITH TWO ROWS AND                  !
!           AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.                         !
! S         INTEGER ARRAY OF SIZE AT LEAST KLMD, USED FOR                      !
!           WORKSPACE.                                                         !
!------------------------------------------------------------------------------!
  SUBROUTINE simplexAlgorithm(K, L, M, N,                                      &
                              KLMD, NKLMD,                                     &
                              Q,                                               &
                              KODE, TOLER,                                     &
                              ITER, X, RES, ERROR)                             !
                                                                               !
    !-- Incoming                                                               !
    INTEGER, INTENT(IN)    :: K, L, M, N                                       !
    INTEGER, INTENT(IN)    :: KLMD, NKLMD                                      !
    INTEGER, INTENT(INOUT) :: KODE, ITER                                       !
    DOUBLETYPE            :: TOLER                                             !
    DOUBLETYPE            :: ERROR                                             !
    DOUBLETYPE, DIMENSION(:,:), INTENT(INOUT) :: Q                             !
    DOUBLETYPE, DIMENSION(:),   INTENT(INOUT) :: X                             !
    DOUBLETYPE, DIMENSION(:),   INTENT(INOUT) :: RES                           !
                                                                               !
    !-- Local                                                                  !
    DOUBLETYPE :: Z, SN, ZU, ZV, CUV, XMAX, XMIN, PIVOT, TPIVOT, ABS           !
    REAL (DP)   :: SUM, DBLE                                                   !
    DOUBLETYPE :: CU(2,NKLMD)                                                  !
                                                                               !
    INTEGER  :: IU(2,NKLMD), S(KLMD)                                        !
    INTEGER  :: I, J,IA, II, IN, JS, KK,                                    &
                   NK, N1, N2, JMN, JPN, KLM, NKL, NK1, IIMN,                  &
                   IOUT, KLM1, KLM2, NKLM, NKL1,                               &
                   MAXIT, IPHASE, KFORCE, IINEG, IABS                          !


!
! INITIALIZATION.
!
      iout=0; in=0
      MAXIT = ITER
      N1 = N + 1
      N2 = N + 2
      NK = N + K
      NK1 = NK + 1
      NKL = NK + L
      NKL1 = NKL + 1
      KLM = K + L + M
      KLM1 = KLM + 1
      KLM2 = KLM + 2
      NKLM = N + KLM
      KFORCE = 1
      ITER = 0
      JS = 1
      IA = 0

      ! Set up labels in Q.
      DO J=1,N
         Q(KLM2,J) = J
      END DO
      DO I=1,KLM
         Q(I,N2) = N + I
         IF (Q(I,N1).LT.0.) THEN
            DO J=1,N2
               Q(I,J) = -Q(I,J)
            END DO
         END IF
      END DO

      ! Set up phase 1 costs.
      IPHASE = 2
      DO J=1,NKLM
         CU(1,J) = 0.
         CU(2,J) = 0.
         IU(1,J) = 0
         IU(2,J) = 0
      END DO
      IF (L.GT.0) THEN
         DO J=NK1,NKL
            CU(1,J) = 1.
            CU(2,J) = 1.
            IU(1,J) = 1
            IU(2,J) = 1
         END DO
         IPHASE = 1
      END IF
      IF (M.NE.0) THEN
         DO J=NKL1,NKLM
            CU(2,J) = 1.
            IU(2,J) = 1
            JMN = J - N
            IF (Q(JMN,N2).LT.0.) IPHASE = 1
         END DO
      END IF
      IF (KODE.EQ.0) GO TO 150
      DO 110 J=1,N
         IF (X(J)) 90, 110, 100
   90    CU(1,J) = 1.
         IU(1,J) = 1
         GO TO 110
  100    CU(2,J) = 1.
         IU(2,J) = 1
  110 CONTINUE
      DO 140 J=1,K
         JPN = J + N
         IF (RES(J)) 120, 140, 130
  120    CU(1,JPN) = 1.
         IU(1,JPN) = 1
         IF (Q(J,N2).GT.0.0) IPHASE = 1
         GO TO 140
  130    CU(2,JPN) = 1.
         IU(2,JPN) = 1
         IF (Q(J,N2).LT.0.0) IPHASE = 1
  140 CONTINUE
  150 IF (IPHASE.EQ.2) GO TO 500
! COMPUTE THE MARGINAL COSTS.
  160 DO J=JS,N1
         SUM = 0.D0
         DO I=1,KLM
            II = Q(I,N2)
            IF (II.LT.0) GO TO 170
            Z = CU(1,II)
            GO TO 180
  170       IINEG = -II
            Z = CU(2,IINEG)
  180       SUM = SUM + DBLE(Q(I,J))*DBLE(Z)
         END DO
         Q(KLM1,J) = SUM
      END DO
      DO J=JS,N
         II = Q(KLM2,J)
         IF (II.LT.0) GO TO 210
         Z = CU(1,II)
         GO TO 220
  210    IINEG = -II
         Z = CU(2,IINEG)
  220    Q(KLM1,J) = Q(KLM1,J) - Z
      END DO
! DETERMINE THE VECTOR TO ENTER THE BASIS.
  240 XMAX = 0.
      IF (JS.GT.N) GO TO 490
      DO 280 J=JS,N
         ZU = Q(KLM1,J)
         II = Q(KLM2,J)
         IF (II.GT.0) GO TO 250
         II = -II
         ZV = ZU
         ZU = -ZU - CU(1,II) - CU(2,II)
         GO TO 260
  250    ZV = -ZU - CU(1,II) - CU(2,II)
  260    IF (KFORCE.EQ.1 .AND. II.GT.N) GO TO 280
         IF (IU(1,II).EQ.1) GO TO 270
         IF (ZU.LE.XMAX) GO TO 270
         XMAX = ZU
         IN = J
  270    IF (IU(2,II).EQ.1) GO TO 280
         IF (ZV.LE.XMAX) GO TO 280
         XMAX = ZV
         IN = J
  280 CONTINUE
      IF (XMAX.LE.TOLER) GO TO 490
      IF (Q(KLM1,IN).EQ.XMAX) GO TO 300
      DO 290 I=1,KLM2
         Q(I,IN) = -Q(I,IN)
  290 CONTINUE
      Q(KLM1,IN) = XMAX
! DETERMINE THE VECTOR TO LEAVE THE BASIS.
  300 IF (IPHASE.EQ.1 .OR. IA.EQ.0) GO TO 330
      XMAX = 0.
      DO 310 I=1,IA
         Z = ABS(Q(I,IN))
         IF (Z.LE.XMAX) GO TO 310
         XMAX = Z
         IOUT = I
  310 CONTINUE
      IF (XMAX.LE.TOLER) GO TO 330
      DO 320 J=1,N2
         Z = Q(IA,J)
         Q(IA,J) = Q(IOUT,J)
         Q(IOUT,J) = Z
  320 CONTINUE
      IOUT = IA
      IA = IA - 1
      PIVOT = Q(IOUT,IN)
      GO TO 420
  330 KK = 0
      DO 340 I=1,KLM
         Z = Q(I,IN)
         IF (Z.LE.TOLER) GO TO 340
         KK = KK + 1
         RES(KK) = Q(I,N1)/Z
         S(KK) = I
  340 CONTINUE
  350 IF (KK.GT.0) GO TO 360
      KODE = 2
      GO TO 590
  360 XMIN = RES(1)
      IOUT = S(1)
      J = 1
      IF (KK.EQ.1) GO TO 380
      DO 370 I=2,KK
         IF (RES(I).GE.XMIN) GO TO 370
         J = I
         XMIN = RES(I)
         IOUT = S(I)
  370 CONTINUE
      RES(J) = RES(KK)
      S(J) = S(KK)
  380 KK = KK - 1
      PIVOT = Q(IOUT,IN)
      II = Q(IOUT,N2)
      IF (IPHASE.EQ.1) GO TO 400
      IF (II.LT.0) GO TO 390
      IF (IU(2,II).EQ.1) GO TO 420
      GO TO 400
  390 IINEG = -II
      IF (IU(1,IINEG).EQ.1) GO TO 420
  400 II = IABS(II)
      CUV = CU(1,II) + CU(2,II)
      IF (Q(KLM1,IN)-PIVOT*CUV.LE.TOLER) GO TO 420
! BYPASS INTERMEDIATE VERTICES.
      DO 410 J=JS,N1
         Z = Q(IOUT,J)
         Q(KLM1,J) = Q(KLM1,J) - Z*CUV
         Q(IOUT,J) = -Z
  410 CONTINUE
      Q(IOUT,N2) = -Q(IOUT,N2)
      GO TO 350
! GAUSS-JORDAN ELIMINATION.
  420 IF (ITER.LT.MAXIT) GO TO 430
      KODE = 3
      GO TO 590
  430 ITER = ITER + 1
      DO 440 J=JS,N1
         IF (J.NE.IN) Q(IOUT,J) = Q(IOUT,J)/PIVOT
  440 CONTINUE
! IF PERMITTED, USE SUBROUTINE COLVEC OF THE DESCRIPTION
! SECTION AND REPLACE THE FOLLOWING SEVEN STATEMENTS DOWN
! TO AND INCLUDING STATEMENT NUMBER 460 BY..

     DO 460 J=JS,N1
        IF(J .EQ. IN) GO TO 460
        Z = -Q(IOUT,J)
        CALL COLVEC(Q(:,J), Q(:,IN), Z, IOUT, KLM1)
 460 CONTINUE

! IF COLVEC IS A PROBLEM THEN USE THIS INSTEAD OF ABOVE:
!      DO 460 J=JS,N1
!         IF (J.EQ.IN) GO TO 460
!         Z = -Q(IOUT,J)
!         DO 450 I=1,KLM1
!            IF (I.NE.IOUT) Q(I,J) = Q(I,J) + Z*Q(I,IN)
!  450    CONTINUE
!  460 CONTINUE

      TPIVOT = -PIVOT
      DO 470 I=1,KLM1
         IF (I.NE.IOUT) Q(I,IN) = Q(I,IN)/TPIVOT
  470 CONTINUE
      Q(IOUT,IN) = 1./PIVOT
      Z = Q(IOUT,N2)
      Q(IOUT,N2) = Q(KLM2,IN)
      Q(KLM2,IN) = Z
      II = ABS(Z)
      IF (IU(1,II).EQ.0 .OR. IU(2,II).EQ.0) GO TO 240
      DO 480 I=1,KLM2
         Z = Q(I,IN)
         Q(I,IN) = Q(I,JS)
         Q(I,JS) = Z
  480 CONTINUE
      JS = JS + 1
      GO TO 240
! TEST FOR OPTIMALITY.
  490 IF (KFORCE.EQ.0) GO TO 580
      IF (IPHASE.EQ.1 .AND. Q(KLM1,N1).LE.TOLER) GO TO 500
      KFORCE = 0
      GO TO 240
! SET UP PHASE 2 COSTS.
  500 IPHASE = 2
      DO J=1,NKLM
         CU(1,J) = 0.
         CU(2,J) = 0.
      END DO
      DO J=N1,NK
         CU(1,J) = 1.
         CU(2,J) = 1.
      END DO
      DO 560 I=1,KLM
         II = Q(I,N2)
         IF (II.GT.0) GO TO 530
         II = -II
         IF (IU(2,II).EQ.0) GO TO 560
         CU(2,II) = 0.
         GO TO 540
  530    IF (IU(1,II).EQ.0) GO TO 560
         CU(1,II) = 0.
  540    IA = IA + 1
         DO 550 J=1,N2
            Z = Q(IA,J)
            Q(IA,J) = Q(I,J)
            Q(I,J) = Z
  550    CONTINUE
  560 CONTINUE
      GO TO 160
  570 IF (Q(KLM1,N1).LE.TOLER) GO TO 500
      KODE = 1
      GO TO 590
  580 IF (IPHASE.EQ.1) GO TO 570
! PREPARE OUTPUT.
      KODE = 0
  590 SUM = 0.D0
      DO 600 J=1,N
         X(J) = 0.
  600 CONTINUE
      DO 610 I=1,KLM
         RES(I) = 0.
  610 CONTINUE
      DO 640 I=1,KLM
         II = Q(I,N2)
         SN = 1.
         IF (II.GT.0) GO TO 620
         II = -II
         SN = -1.
  620    IF (II.GT.N) GO TO 630
         X(II) = SN*Q(I,N1)
         GO TO 640
  630    IIMN = II - N
         RES(IIMN) = SN*Q(I,N1)
         IF (II.GE.N1 .AND. II.LE.NK) SUM = SUM + &
         DBLE(Q(I,N1))
  640 CONTINUE
      ERROR = SUM
      RETURN

 END SUBROUTINE simplexAlgorithm
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! IF YOUR FORTRAN COMPILER PERMITS A SINGLE COLUMN OF A TWO                    !
! DIMENSIONAL ARRAY TO BE PASSED TO A ONE DIMENSIONAL ARRAY                    !
! THROUGH A SUBROUTINE CALL, CONSIDERABLE SAVINGS IN                           !
! EXECUTION TIME MAY BE ACHIEVED THROUGH THE USE OF THE                        !
! FOLLOWING SUBROUTINE, WHICH OPERATES ON COLUMN VECTORS.                      !
! SEE COMMENTS FOLLOWING STATEMENT LABELLED 440 FOR                            !
! INSTRUCTIONS ON THE IMPLEMENTATION OF THIS MODIFICATION.                     !
!------------------------------------------------------------------------------!
SUBROUTINE COLVEC(V1, V2, XMLT, NOTROW, K)
  ! THIS SUBROUTINE ADDS TO THE VECTOR V1 A MULTIPLE OF THE
  ! VECTOR V2 (ELEMENTS 1 THROUGH K EXCLUDING NOTROW).
     INTEGER K, NOTROW, KEND, KSTART,I
     DOUBLETYPE, DIMENSION(K) :: V1, V2
     DOUBLETYPE :: XMLT
!BEGIN
     KEND = NOTROW - 1
     KSTART = NOTROW + 1
     IF (KEND .GE. 1) THEN
        DO I=1,KEND
           V1(I) = V1(I) + XMLT*V2(I)
        END DO
        IF(KSTART .GT. K) RETURN
     END IF
     DO I=KSTART,K
       V1(I) = V1(I) + XMLT*V2(I)
     END DO

END SUBROUTINE COLVEC
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
!-- Below Here Relates to Derived Output Only
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! GetListOfGeochemDiagnostics                                                         !
!                                                                              !
! Routine called during setup/initialisation to determine necessary variables  !
! that need to be included in the geochem derived variable array               !
!------------------------------------------------------------------------------!
 SUBROUTINE GetListOfGeochemDiagnostics(speciesOutput,derivedNames)
    IMPLICIT NONE

    CHARACTER(LEN=64), DIMENSION(10), INTENT(IN) :: speciesOutput
    CHARACTER(LEN=64), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: derivedNames

    CHARACTER(LEN=64), DIMENSION(:), ALLOCATABLE              :: tmpNames
    INTEGER  :: i, sIndex, gcCntr
    INTEGER  :: start,finish
    LOGICAL :: found

      start  = 1
      finish = SIZE(speciesOutput)


    !-- Allocate a character array to store the names of all non
    !-- state variable variables that must be remembered. This array is
    !-- purposely allocated too large
    ALLOCATE(tmpNames(finish+10))

    !--------------------------------------------------------------------------!
    !-- First we must count the number of derived variables that are relevant
    !-- to this geochemistry module
    !-- Loop through the derived/total variables the user selected to
    !-- output, and determine if it is relevant for the geochem module
    nGCDerivedVars = 0


    !-- Aqueous Species
    DO i = start,finish
      DO sIndex = 1, nSpecies
        IF(TRIM(speciesOutput(i)) == TRIM(allSpecies(sIndex)%Name)) THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = TRIM(allSpecies(sIndex)%Name)
          EXIT
        END IF
      END DO
    END DO

    !-- If Iron Redox is being simulated then AED_geochem needs us
    !-- to remember the some FeII & FeIII species since it uses these for
    !-- calculating the FeII oxidation, and FeIII photoreduction
    IF (simIronRedox) THEN
      !-- Fe+2 (FeII)
      found = .FALSE.
      DO i = start,finish
        IF(TRIM(speciesOutput(i)) == "Fe+2") THEN
           found = .TRUE.
           EXIT
        END IF
      END DO
      IF(.NOT.found) THEN
        nGCDerivedVars = nGCDerivedVars + 1
        tmpNames(nGCDerivedVars) = "Fe+2      "
      END IF
      !-- FeOH+ (FeII)
      found = .FALSE.
      DO i = start,finish
        IF(TRIM(speciesOutput(i)) == "FeOH+") THEN
           found = .TRUE.
           EXIT
        END IF
      END DO
      IF(.NOT.found) THEN
        nGCDerivedVars = nGCDerivedVars + 1
        tmpNames(nGCDerivedVars) = "FeOH+     "
      END IF
      !-- Fe(OH)2 (FeII)
      found = .FALSE.
      DO i = start,finish
        IF(TRIM(speciesOutput(i)) == "Fe(OH)2") THEN
           found = .TRUE.
           EXIT
        END IF
      END DO
      IF(.NOT.found) THEN
        nGCDerivedVars = nGCDerivedVars + 1
        tmpNames(nGCDerivedVars) = "Fe(OH)2   "
      END IF
      !-- FeOH+2 (FeIII)
      found = .FALSE.
      DO i = start,finish
        IF(TRIM(speciesOutput(i)) == "FeOH+2") THEN
           found = .TRUE.
           EXIT
        END IF
      END DO
      IF(.NOT.found) THEN
        nGCDerivedVars = nGCDerivedVars + 1
        tmpNames(nGCDerivedVars) = "FeOH+2    "
      END IF
    END IF

    !-- Non convergence
    DO i = start,finish
      IF(TRIM(speciesOutput(i)) == "NONCON") THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = "NONCON    "
      END IF
    END DO

    !-- Saturation Indicies for upto five mineral phases
    DO i = start,finish
      IF(TRIM(speciesOutput(i)) == "SI_PP1") THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = "SI_PP1    "
      ELSE IF(TRIM(speciesOutput(i)) == "SI_PP2") THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = "SI_PP2    "
      ELSE IF(TRIM(speciesOutput(i)) == "SI_PP3") THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = "SI_PP3    "
      ELSE IF(TRIM(speciesOutput(i)) == "SI_PP4") THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = "SI_PP4    "
      ELSE IF(TRIM(speciesOutput(i)) == "SI_PP5") THEN
          nGCDerivedVars = nGCDerivedVars + 1
          tmpNames(nGCDerivedVars) = "SI_PP5    "
      END IF
    END DO

    !-- If DIC is being simulated then we need to add one for pCO2
    IF (simC_DIC) THEN
      nGCDerivedVars = nGCDerivedVars + 1
      tmpNames(nGCDerivedVars) = "pCO2      "
    END IF

    WRITE(*,"(5X,'nGCDerivedVars: ',I3,/)")nGCDerivedVars
    DO i = 1,nGCDerivedVars
      WRITE(*,"(7X,'Var: ',A10,/)")tmpNames(i)
    END DO


    !--------------------------------------------------------------------------!
    !-- Now allocate a vector which maps to the relevant species               !
    ALLOCATE(derivedGCList(nGCDerivedVars))

    !-- Populate the map
    derivedGCList = 0
    gcCntr = 1
    DO i = 1,nGCDerivedVars
      DO sIndex = 1, nSpecies
        IF(TRIM(tmpNames(i)) == TRIM(allSpecies(sIndex)%Name)) THEN
          derivedGCList(gcCntr) = allSpecies(sIndex)%SpeciesIndex
          gcCntr = gcCntr + 1
        END IF
      END DO
    END DO

    DO i = 1,nGCDerivedVars
      IF(TRIM(tmpNames(i)) == "SI_PP1") THEN
          derivedGCList(gcCntr) = -1
          gcCntr = gcCntr + 1
      ELSE IF(TRIM(tmpNames(i)) == "SI_PP2") THEN
          derivedGCList(gcCntr) = -2
          gcCntr = gcCntr + 1
      ELSE IF(TRIM(tmpNames(i)) == "SI_PP3") THEN
          derivedGCList(gcCntr) = -3
          gcCntr = gcCntr + 1
      ELSE IF(TRIM(tmpNames(i)) == "SI_PP4") THEN
          derivedGCList(gcCntr) = -4
          gcCntr = gcCntr + 1
      ELSE IF(TRIM(tmpNames(i)) == "SI_PP5") THEN
          derivedGCList(gcCntr) = -5
          gcCntr = gcCntr + 1
      END IF
    END DO

    DO i = 1,nGCDerivedVars
      IF(TRIM(tmpNames(i)) == "NONCON") THEN
          derivedGCList(gcCntr) = -998
          gcCntr = gcCntr + 1
      END IF
    END DO

    !-- If DIC is being simulated then we need to add one for pCO2
    IF (simC_DIC) THEN
      derivedGCList(gcCntr) = -999
    END IF

   !-- Allocate array which will store the derived values
   ALLOCATE(derivedGCVals(nGCDerivedVars))
   derivedGCVals = gc_zero


   ALLOCATE(derivedNames(nGCDerivedVars))
   DO i = 1,nGCDerivedVars
          derivedNames(i) = tmpNames(i)
   END DO
   DEALLOCATE(tmpNames)

 END SUBROUTINE GetListOfGeochemDiagnostics
!------------------------------------------------------------------------------!








!------------------------------------------------------------------------------!
! CalcpCO2                                                                     !
!                                                                              !
! Function that returns the partial pressure of CO2 in the water when supplied !
! with the current species molalities                                          !
!------------------------------------------------------------------------------!
 FUNCTION calcpCO2(specs, temp, salinity) RESULT(pCO2)                         !
   !-- Incoming                                                                !
   TYPE (gcSpecies), DIMENSION(:), INTENT(IN) :: specs                         !
   DOUBLETYPE,                    INTENT(IN) :: temp                          !
   DOUBLETYPE,                    INTENT(IN) :: salinity                      !
   !-- Outgoing                                                                !
   DOUBLETYPE                                :: pCO2                          !
   !-- Local                                                                   !
   CHARACTER(*), PARAMETER                    :: THIS_PROC = "calcpCO2"        !
   DOUBLETYPE                                :: d1,Solubility,CO2,H2CO3       !
   INTEGER  :: nSpec, sIndex                 !


   nSpec = SIZE(specs)
   pCO2  = 350e-6

   CO2    = gc_zero
   H2CO3  = gc_zero

   DO sIndex = 1, nSpec
     IF(TRIM(specs(sIndex)%Name) == "CO2") THEN
       CO2 = specs(sIndex)%Moles
       CYCLE
     END IF
     IF(TRIM(specs(sIndex)%Name) == "H2CO3") THEN
       H2CO3 = specs(sIndex)%Moles
       CYCLE
     END IF
   END DO

   ! From Stumm and Morgan:
   ! pCO2 = H2CO3* / Ko
   ! [atm] = [mol/L] / ([mol/L/atm])
   ! where:
   ! H2CO3* = CO2(aq) + H2CO3

   ! Temperature (K)
   d1 = temp + WaterFreezingInKelvins
   ! Solubility, Ko (mol/L/atm)
   Solubility = -58.0931 + 90.5069*(100.0/d1) + 22.294*LOG(d1/100.0)           &
              +  0.027766*salinity - 0.025888*salinity*(d1/100.0)

   Solubility = Solubility + 0.0050578*salinity*(d1/100.0)*(d1/100.0)

   Solubility = exp(Solubility)

   pCO2 = (CO2 + H2CO3) / Solubility

   END FUNCTION calcpCO2
!------------------------------------------------------------------------------!





!------------------------------------------------------------------------------!
! CalcSIforPP                                                                  !
!                                                                              !
! This is a derived function that calculates a minerals saturation index (SI)  !
! based on the simulated results for output to the NetCDF files                !
!------------------------------------------------------------------------------!
 FUNCTION calcSIforPP(PICHMIndex) RESULT(SI_PPX)                               !
   !-- Incoming                                                                !
   INTEGER  :: PICHMIndex                                                   !
   !-- Outgoing                                                                !
   DOUBLETYPE  :: SI_PPX                                                       !
   !-- Local                                                                   !
   INTEGER  :: cIndex                                                       !


   SI_PPX = -999
   DO cIndex = 1,nComponents
     IF(allComponents(cIndex)%CompType == PUREPHASE) THEN
       IF(allComponents(cIndex)%eqIndex == PICHMIndex) THEN
         ! Saturation Index, SI = log(Q)-log(Ksp) , which here is corrected for
         ! the fact Value and Total are ntaural logs for the solution scheme.
         SI_PPX = (allComponents(cIndex)%Value -                               &
                               allComponents(cIndex)%Total) / LOG(REAL(10.0,DP))
       END IF
     END IF
   END DO


 END FUNCTION calcSIforPP
!------------------------------------------------------------------------------!




!------------------------------------------------------------------------------!
! CalcIAPforPP                                                                 !
!                                                                              !
! This is a function that calculates a minerals saturation index (SI)          !
! based on the non-simulated results for availability to REDOX calcs           !
!------------------------------------------------------------------------------!
FUNCTION CalcIAPforPP(PICHMIndex) RESULT(IAPonK)                               !
   !-- Incoming                                                                !
   INTEGER  :: PICHMIndex                                                   !
   !-- Outgoing                                                                !
   DOUBLETYPE  :: IAPonK                                                       !
   !-- Local
   DOUBLETYPE  :: Q,K

   IAPonK = -999.

   IF(allComponents(PICHMIndex)%CompType == PUREPHASE) THEN

     ! For CANDI-AED we need "IAP/Ksp" to limit rate laws (note in Bethke the
     ! term "Q" (Eq 3.35,p47) is used as ion activity product
     ! In this module, %Value = log(Q)*ln(10) &
     !                 %Total = log(Ksp)*ln(10
     Q = antiLOG( allComponents(PICHMIndex)%Value / LOG(REAL(10.0,DP)), REAL(10.0,DP) )
     K = antiLOG( allComponents(PICHMIndex)%ppData%logK, REAL(10.0,DP) )
     !Q = allComponents(PICHMIndex)%Value / LOG(REAL(10.0,DP))
     !K = allComponents(PICHMIndex)%ppData%logK

     IAPonK = Q/K

   ELSE
     print *,'CalcIAP error - CalcIAPforPP called on a non-PurePhase component'
   END IF


END FUNCTION CalcIAPforPP
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! CalcIAPforPP2                                                                 !
!                                                                              !
! This is a function that is exactly the same as the one above,                !
! but it produces K as a result, for use in equilibrium calculations. -Dan     !
!------------------------------------------------------------------------------!
FUNCTION CalcIAPforPP2(PICHMIndex) RESULT(KIAP)                                 !
   !-- Incoming                                                                !
   INTEGER  :: PICHMIndex                                                   !
   !-- Outgoing                                                                !
   DOUBLETYPE  :: KIAP                                                       !
   !-- Local
   DOUBLETYPE  :: Q,K

   KIAP = -999.

   IF(allComponents(PICHMIndex)%CompType == PUREPHASE) THEN

     ! For CANDI-AED we need "IAP/Ksp" to limit rate laws (note in Bethke the
     ! term "Q" (Eq 3.35,p47) is used as ion activity product
     ! In this module, %Value = log(Q)*ln(10) &
     !                 %Total = log(Ksp)*ln(10
     Q = antiLOG( allComponents(PICHMIndex)%Value / LOG(REAL(10.0,DP)), REAL(10.0,DP) )
     K = antiLOG( allComponents(PICHMIndex)%ppData%logK, REAL(10.0,DP) )
     !Q = allComponents(PICHMIndex)%Value / LOG(REAL(10.0,DP))
     !K = allComponents(PICHMIndex)%ppData%logK

     KIAP = K

   ELSE
     print *,'CalcIAP error - CalcIAPforPP2 called on a non-PurePhase component'
   END IF


END FUNCTION CalcIAPforPP2
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! CalcIAPforPP3                                                                 !
!                                                                              !
! This is a function that is exactly the same as the one above,                !
! but it produces Q as a result, for use in equilibrium calculations.
! Q here is the IAP that you might see in the references, for example
! Tufano 2009 p 1008.   -Dan     !
!------------------------------------------------------------------------------!
FUNCTION CalcIAPforPP3(PICHMIndex) RESULT(QIAP)                                 !
   !-- Incoming                                                                !
   INTEGER  :: PICHMIndex                                                   !
   !-- Outgoing                                                                !
   DOUBLETYPE  :: QIAP                                                       !
   !-- Local
   DOUBLETYPE  :: Q,K

   QIAP = -999.

   IF(allComponents(PICHMIndex)%CompType == PUREPHASE) THEN

     ! For CANDI-AED we need "IAP/Ksp" to limit rate laws (note in Bethke the
     ! term "Q" (Eq 3.35,p47) is used as ion activity product
     ! In this module, %Value = log(Q)*ln(10) &
     !                 %Total = log(Ksp)*ln(10
     Q = antiLOG( allComponents(PICHMIndex)%Value / LOG(REAL(10.0,DP)), REAL(10.0,DP) )
     K = antiLOG( allComponents(PICHMIndex)%ppData%logK, REAL(10.0,DP) )
     !Q = allComponents(PICHMIndex)%Value / LOG(REAL(10.0,DP))
     !K = allComponents(PICHMIndex)%ppData%logK

     QIAP = Q

   ELSE
     print *,'CalcIAP error - CalcIAPforPP3 called on a non-PurePhase component'
   END IF


END FUNCTION CalcIAPforPP3
!------------------------------------------------------------------------------!



!------------------------------------------------------------------------------!
END MODULE aed2_sedeqsolver
