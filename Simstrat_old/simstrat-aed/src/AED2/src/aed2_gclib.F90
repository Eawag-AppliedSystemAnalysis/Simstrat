!###############################################################################
!#                                                                             #
!# File: aed2_gclib.F90                                                        #
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
!#                                                                             #
!###############################################################################

#include "aed2.h"
#include "aed2_debug.h"

MODULE aed2_gclib

   USE aed2_gctypes

   IMPLICIT NONE

   PRIVATE

   TYPE(AEDChmConstType) :: AEDConst

   CHARACTER(512)        :: geochemFile = ''

   PUBLIC antiLOG,NumOfComps,NumOfSpecies,AED_GC_Input, &
          GetCompCompName,GetSpeciesName,GetCompInfo,   &
          GetSpeciesInfo, GetPhaseInfo, printTables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define IsSign(ch) ((ch) == PLUS .or. (ch) == MINUS)
#define IsDigit(ch) ((ch) >= '0' .and. (ch) <= '9')
#define IsHorizWhiteSpace(ch) ((ch) == ' ' .or. (ch) == '\t')
#define FoundSpeciesName(speciesName) (Locn_Species(speciesName) /= NOT_FOUND)
#define FoundPhaseName(phaseName)     (Locn_Phase(phaseName) /= NOT_FOUND)
!#define    HALT .true.
!#define NO_HALT .false.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   INTEGER, PARAMETER :: MAX_NUM_SPECIES = 500
   INTEGER, PARAMETER :: MAX_NUM_PHASES  = 500
   INTEGER, PARAMETER :: MAX_NUM_COMPS   = 200

   INTEGER, PARAMETER :: NOT_FOUND = -666

!#define NAME_TYPE 1
#define ELT_NAME  2

!  INTEGER(NAME_TYPE), PARAMETER ::  COMP_NAME = 1

   TYPE DATA_STORE_gcUnknowns
      TYPE(gcUnknowns), DIMENSION(:), POINTER :: Data_V => NULL()
      INTEGER  :: Last   =  0
      LOGICAL                                 :: Exists = .FALSE.
   END TYPE DATA_STORE_gcUnknowns


   TYPE DATA_STORE_gcSpecies
      TYPE(gcSpecies), DIMENSION(:), POINTER :: Data_V => NULL()
      INTEGER  :: Last   =  0
      LOGICAL                                :: Exists = .FALSE.
   END TYPE DATA_STORE_gcSpecies


   TYPE DATA_STORE_PP
      TYPE(gcUnknowns), DIMENSION(:), POINTER :: Data_V => NULL()
      INTEGER  :: Last   =  0
      LOGICAL                                 :: Exists = .FALSE.
   END TYPE DATA_STORE_PP



   ! ------ Module-Level Variables -----------

   ! -- The main data stores for geochem info:
   TYPE(DATA_STORE_gcUnknowns), SAVE :: compList        !Components Data Store
   TYPE(DATA_STORE_gcSpecies),  SAVE :: speciesList     !Species Data Store
   TYPE(DATA_STORE_PP),         SAVE :: phaseList       !Pure Phases Data Store


!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE AED_GC_Input(geoChemFile)
!-------------------------------------------------------------------------------
   CHARACTER(*),INTENT(in) :: geoChemFile
   CHARACTER(512)  :: t_word1
   CHARACTER(512)  :: symb
   INTEGER  :: data_count
   INTEGER  :: status
   INTEGER  :: infile
!
!-------------------------------------------------------------------------------
!BEGIN
   AEDConst%GCH_DATA = .FALSE.
   AEDConst%CNP_PARS = .FALSE.
   AEDConst%IRN_PARS = .FALSE.
   AEDConst%MAN_PARS = .FALSE.
   AEDConst%ASN_PARS = .FALSE.
   AEDConst%ALM_PARS = .FALSE.
   AEDConst%DOX_PARS = .FALSE.
   AEDConst%CH4_PARS = .FALSE.
   AEDConst%TRA_PARS = .FALSE.
   AEDConst%BAC_PARS = .FALSE.

   infile = f_get_lun()
   OPEN(UNIT=infile,FILE=geochemFile,STATUS = "OLD",ACTION="READ",IOSTAT=status)

   data_count = 1
   DO WHILE (data_count /= 0)
      t_word1 = ' '

      read(infile,*) t_word1

      IF (.not. (t_word1(1:1) == '#' .or. t_word1(1:1) == '!') ) THEN
         IF (t_word1 /= 'GEOCHEMISTRY' ) CYCLE

         AEDConst%GCH_DATA = .TRUE.

         ! Skip leading comment line of file
         READ(INFILE,*,IOSTAT=status)

         !-------
         ! -- Process the Components part of the chem input file:
         CALL CreateCompStore()
         CALL ProcessComponentsBlock(infile)

         !-------
         ! -- Process the Species part of the chem input file:
         CALL CreateSpeciesStore()
         CALL ProcessSpeciesBlock(infile, symb)

         !-------
         IF (symb == 'PHASES') THEN
            ! -- Process the Phases part of the chem input file if it exists:
            CALL CreatePhasesStore()
            CALL ProcessPhasesBlock(infile, symb)
         ELSE IF (symb == 'ENDFILE') THEN
            EXIT
         ENDIF

         EXIT
      ENDIF
   ENDDO

   close(infile)
END SUBROUTINE AED_GC_Input
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ProcessComponentsBlock(infile)
!-------------------------------------------------------------------------------
! Reads the Components block of the input file and extracts
! appropriate parameters for later use.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,intent(in) :: infile
!
!LOCALS
   INTEGER  :: i

   CHARACTER(INP_LINE_LEN) :: inpLine
   CHARACTER(INP_LINE_LEN) :: symb

   TYPE(gcUnknowns) :: oneChemComponent

   ! ------ Non-Local Variables ------------
   !  geochemFile
   !----------------------------------------

!
!-------------------------------------------------------------------------------
!BEGIN
   !----------
   ! -- 1. Skip leading blank and comment lines:
   DO
      CALL GetNextSignificantLine(infile,inpLine)
      CALL ReadLeadingSymbol(inpLine,symb)
      IF (symb == 'COMPONENTS') EXIT
   ENDDO

   !----------
   ! -- 2. Process COMPONENTS block:
   i = 0
   CALL GetNextSignificantLine(infile,inpLine)

   DO WHILE (symb /= 'SPECIES')
      i = i + 1
      CALL ProcessComponentsRecord(inpLine,i,oneChemComponent)
      CALL AddComp(oneChemComponent)
      CALL GetNextSignificantLine(infile,inpLine)
      CALL ReadLeadingSymbol(inpLine,symb)
   ENDDO
END SUBROUTINE ProcessComponentsBlock
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ProcessComponentsRecord(inpLine,i,gcUnknown)
!-------------------------------------------------------------------------------
! ProcessComponentsRecord
! -----------------------
!
! Reads the Components block of the input file and extracts
! appropriate parameters for later use.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(INP_LINE_LEN),INTENT(IN)  :: inpLine     !the input line
   INTEGER,INTENT(IN)                  :: i           !count
   TYPE(gcUnknowns),INTENT(OUT)        :: gcUnknown   !element of vector
!
!LOCALS
   INTEGER  :: charge
   INTEGER  :: status

   DOUBLETYPE :: molWt

   CHARACTER(STR_LEN) :: element
   CHARACTER(STR_LEN) :: component

!
!-------------------------------------------------------------------------------
!BEGIN
   READ(inpLine,*,IOSTAT=status) element, component, charge, molWt
   _CheckFileIOStatus(status)

   gcUnknown%EltName   = TRIM(element)
   gcUnknown%CompName  = TRIM(component)
   gcUnknown%CompIndex = i
   gcUnknown%Charge    = charge
   gcUnknown%MolWeight = molWt
END SUBROUTINE ProcessComponentsRecord
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ProcessSpeciesBlock(infile, symb)
!-------------------------------------------------------------------------------
!  ProcessSpeciesBlock
!  -------------------
!
!  Reads the Species block of the input file and extracts
!  appropriate parameters for later use.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                  :: infile
   CHARACTER(INP_LINE_LEN),INTENT(OUT) :: symb
!
!LOCALS
   INTEGER  :: status     !dynamic memory allocation status
   INTEGER  :: i               !loop counter

   CHARACTER(INP_LINE_LEN) :: inpLine

   TYPE(gcSpecies) :: oneChemSpecies  !species structured variable

   ! ------ Non-Local Variables ------------
   !  geochemFile
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(oneChemSpecies%stoich(1:NumOfComps()),STAT=status)
   _CheckAllocStatus(status)

   !----------
   ! -- Process SPECIES block:
   i = 0
   CALL GetNextSignificantLine(infile,inpLine)
   CALL ReadLeadingSymbol(inpLine,symb)
   DO WHILE (symb /= 'PHASES' .and. symb /= 'ENDFILE')
      i = i + 1
      oneChemSpecies%stoich(:) = 0     !initialisation
      CALL ProcessSpeciesRecord(inpLine,i,oneChemSpecies)
      CALL AddSpecies(oneChemSpecies)
      CALL GetNextSignificantLine(infile,inpLine)
      CALL ReadLeadingSymbol(inpLine,symb)
   ENDDO

   DEALLOCATE(oneChemSpecies%stoich,STAT=status)
   _CheckAllocStatus(status)
END SUBROUTINE ProcessSpeciesBlock
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE ProcessSpeciesRecord(inpLine,i,oneSpecies)
!-------------------------------------------------------------------------------
!  ProcessSpeciesRecord
!  --------------------
!
!  Reads the Species block of the input file and extracts
!  appropriate parameters for later use.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(INP_LINE_LEN),INTENT(in) :: inpLine     !the input line
   INTEGER,INTENT(in)                 :: i           !count
   TYPE(gcSpecies),INTENT(inout)      :: oneSpecies  !element of vector
!
!LOCALS
   ! ------ Local Variables ----------------
   INTEGER  :: charge
   INTEGER  :: startLocn
   INTEGER  :: endLocn

   DOUBLETYPE :: deltaH
   DOUBLETYPE :: logKat25

   CHARACTER(INP_LINE_LEN) :: speciesRec       !species info record
   CHARACTER(INP_LINE_LEN) :: chemParamsStr    !chem. params for the species
   CHARACTER(INP_LINE_LEN) :: chemReactionStr  !chem. reaction for the species
   CHARACTER(STR_LEN)      :: speciesName      !no charge included in name
   CHARACTER(STR_LEN)      :: speciesCompName  !component name includes charge
!
!-------------------------------------------------------------------------------
!BEGIN
   !----------
   speciesRec = ADJUSTL(inpLine)  !make a copy of input line to play with

   ! -- Get species name:
   CALL GetNextSymbol(speciesRec,speciesCompName)
   CALL ExtractChargeAndChemName(speciesCompName,speciesName,charge)
   !!$WRITE(*,'(1X,3X,"speciesName = ",A,3X,"charge = ",SP,I0)')   &
   !!$                                    TRIM(speciesName), charge

   oneSpecies%SpeciesIndex = i
   !!$WRITE(*,'(1X,3X,"speciesCompName = ",A)')  TRIM(speciesCompName)

   startLocn = SCAN(speciesRec,REACTION_START_CH)
   !!$ EXPAND
   IF (startLocn == 0) STOP "ProcessSpeciesRecord:  no chem. reaction start symbol."
   endLocn   = SCAN(speciesRec,REACTION_END_CH)
   !!$  EXPAND
   IF (endLocn == 0) STOP "ProcessSpeciesRecord:  no chem. reaction end symbol."

   ! -- Only pass the string of chars between "[" and "]" inclusive; i.e., the
   ! -- chemical reaction string:
   chemReactionStr = speciesRec(startLocn:endLocn)
   chemParamsStr   = speciesRec(endLocn+1: )

   CALL ParseChemReaction(speciesCompName,chemReactionStr,oneSpecies%Stoich)
   CALL ReadDataTrailingChemReaction(chemParamsStr,logKat25,deltaH)

   oneSpecies%Name     = speciesCompName
   oneSpecies%logKat25 = logKat25
   oneSpecies%deltaH   = deltaH
   oneSpecies%charge   = ChargeFromName(speciesCompName)
END SUBROUTINE ProcessSpeciesRecord
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE ProcessPhasesBlock(infile, symb)
!-------------------------------------------------------------------------------
!  ProcessPhasesBlock
!  ------------------
!
!  Reads the Phases block of the input file and extracts
!  appropriate parameters for later use.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                    :: infile
   CHARACTER(INP_LINE_LEN),INTENT(inout) :: symb    !input symbol
!
!LOCALS
   INTEGER  :: status  ! memory allocation status
   INTEGER  :: i       ! line counter

   CHARACTER(INP_LINE_LEN) :: inpLine !the input line (as a string)

   TYPE(gcUnknowns) :: oneChemPhase   !species structured variable

   ! ------ Non-Local Variables ------------
   !  geochemFile
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(oneChemPhase%PPData,STAT=status)
   _CheckAllocStatus(status)

   ALLOCATE(oneChemPhase%PPData%Stoich(1:NumOfComps()),STAT=status)
   _CheckAllocStatus(status)

   !----------
   ! -- Process PHASES block:
   i = 0
   CALL GetNextSignificantLine(infile,inpLine)
   CALL ReadLeadingSymbol(inpLine,symb)

   DO WHILE (inpLine /= "ENDFILE")
      i = i + 1
      oneChemPhase%PPData%stoich(:) = 0     !initialisation
      CALL ProcessPhaseRecord(inpLine,i,oneChemPhase)
      CALL AddPhase(oneChemPhase)
      CALL GetNextSignificantLine(infile,inpLine)
      CALL ReadLeadingSymbol(inpLine,symb)
   ENDDO
END SUBROUTINE ProcessPhasesBlock
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE ProcessPhaseRecord(inpLine,i,onePhase)
!-------------------------------------------------------------------------------
!  ProcessPhaseRecord
!  ------------------
!
!  Reads the Phases block of the input file and extracts
!  appropriate parameters for later use.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(INP_LINE_LEN),INTENT(in) :: inpLine   !the input line
   INTEGER,INTENT(in)                 :: i         !count
   TYPE(gcUnknowns),INTENT(inout)     :: onePhase  !element of vector
!
!LOCALS
   INTEGER  :: startLocn
   INTEGER  :: endLocn

   DOUBLETYPE :: deltaH
   DOUBLETYPE :: logKat25
   DOUBLETYPE :: settlingVel

   CHARACTER(INP_LINE_LEN) :: phaseRec         !phase info record
   CHARACTER(INP_LINE_LEN) :: chemParamsStr    !chem. paramaters for the phase
   CHARACTER(INP_LINE_LEN) :: chemReactionStr  !chem. reaction for the phase
   CHARACTER(STR_LEN)      :: phaseCompName    !phase component-name
   CHARACTER(STR_LEN)      :: phaseEltName     !phase element-name

!
!-------------------------------------------------------------------------------
!BEGIN
   phaseRec = ADJUSTL(inpLine)  !make a copy of input line to play with

   ! -- Get phase element-name:
   CALL GetNextSymbol(phaseRec,phaseEltName)
   !!$WRITE(*,'(1X,3X,"phaseEltName = ",A)')  TRIM(phaseEltName)

   onePhase%CompIndex = i

   startLocn = SCAN(phaseRec,REACTION_START_CH)
   !!$ EXPAND
   IF (startLocn == 0) STOP "ProcessPhaseRecord:  no chem. reaction start symbol."
   endLocn   = SCAN(phaseRec,REACTION_END_CH)
   !!$  EXPAND
   IF (endLocn == 0) STOP "ProcessPhaseRecord:  no chem. reaction end symbol."

   ! -- Only pass the string of chars between "[" and "]" inclusive; i.e.,
   ! -- th chemical reaction string:
   chemReactionStr = phaseRec(startLocn:endLocn)
   chemParamsStr   = phaseRec(endLocn+1: )

   ! -- The phase component-name is the string on the LHS of the chemical
   ! -- equation:
   READ(phaseRec(startLocn+1:),*)  phaseCompName

   CALL ParseChemReaction(phaseEltName,chemReactionStr,onePhase%PPData%Stoich)

   CALL ReadDataTrailingChemReaction(chemParamsStr,logKat25,deltaH,settlingVel)

   onePhase%EltName  = phaseEltName
   onePhase%CompName = phaseCompName
   onePhase%Charge   = ChargeFromName(phaseCompName)

   onePhase%PPData%logKat25    = logKat25
   onePhase%PPData%deltaH      = deltaH
   onePhase%PPData%SettlingVel = settlingVel
END SUBROUTINE ProcessPhaseRecord
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ParseChemReaction(chemName,chemReaction,stoichCoeffs_V)
!-------------------------------------------------------------------------------
!  ParseChemReaction
!  -----------------
!
!  Parses the chemical equation of a species and returns
!  the stoichiometric coefficients vector for the species.
!  The chemical reaction string is returned with the initial
!  leading chemical component removed.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(STR_LEN),INTENT(in)         :: chemName       !chem. name (for ID)
   CHARACTER(INP_LINE_LEN),INTENT(inout) :: chemReaction   !the chem eqn
   REAL,DIMENSION(:),INTENT(inout)       :: stoichCoeffs_V !stoichiometric coeffs..
!                                                          !..for this reaction
!LOCALS
   INTEGER  :: charge       ! charge of the chemical component
   INTEGER  :: coeff        ! unsigned stoichiometric coefficient
   INTEGER  :: signedCoeff  ! signed stoichiometric coefficient
   INTEGER  :: mult         ! unitary multiplier to get signed stoich. coeff.

   CHARACTER(1)            :: ch
   CHARACTER(INP_LINE_LEN) :: reactionStr
   CHARACTER(STR_LEN)      :: tokenRead

   CHARACTER(INP_LINE_LEN) :: theFollowSet
!
!-------------------------------------------------------------------------------
!BEGIN
   reactionStr = ADJUSTL(chemReaction)

   !----------

   ! -- Read "[":
   CALL Expect(reactionStr,REACTION_START_CH)

   CALL ReadChemComp(reactionStr,tokenRead,coeff,charge)

   CALL Expect(reactionStr,EQUALS)

   ! -- Leading symbol on RHS of reaction may be a sign:
   mult = 1   !initialisation
   ch = reactionStr(1:1)
   IF (IsSign(ch) ) THEN
      IF (ch == MINUS) THEN
         mult = -1
      ENDIF
      reactionStr = reactionStr(2:)    !eat the sign character
   ELSE IF (.NOT. IsAlphaNumeric(ch)) THEN
      WRITE(*,'(1X,A,":  error parsing reactionStr for ",A,".")') &
                                             TRIM(chemName)
      WRITE(*,'(1X,3X,"reactionStr = ''",A,"''")') reactionStr
      WRITE(*,'(1X,3X,"char read = ",A)') TRIM(ch)
      WRITE(*,'(1X,3X,"Expected a sign or alphanumeric character.")')
      STOP "ParseChemReaction:  syntax error in reactionStr"
   ENDIF

   !-------------------------------
   ! -- Repeat-Until loop:
   ParseChemReactionString:  &
   DO
      CALL ReadChemComp(reactionStr,tokenRead,coeff,charge);
      signedCoeff = mult * coeff

      !!##EXPAND
      !!##Add chem component and stoich coeff into local stoich vector. Pass
      CALL AddToStoichiometry(tokenRead,signedCoeff,stoichCoeffs_V)

      CALL GetNextSymbol(reactionStr,tokenRead)

      IF (.NOT. (tokenRead == PLUS .or. tokenRead == MINUS .or. &
                                          tokenRead == REACTION_END_CH)) THEN
!     IF (.NOT. IsFollowSymbol(FOLLOW_SET_6,tokenRead)) THEN
!        theFollowSet = " "
!        DO i = 1,SIZE(FOLLOW_SET_6)
!           theFollowSet = TRIM(theFollowSet) // ", " // TRIM(FOLLOW_SET_6(i))
!        ENDDO
         WRITE(*,'(1X,3X,"Expected one of the following:  ",A)')              &
                                           TRIM(theFollowSet)
         STOP "ParseChemReaction:  syntax error in reactionStr"
      ENDIF

      ! -- Remember the sign in front of the stoichiometric coeff.:
      IF (tokenRead == PLUS) THEN
         mult = 1
      ELSE IF (tokenRead == MINUS) THEN
         mult = -1
      ENDIF

      ! -- IF token read is end token THEN exit the loop (this emulates a
      ! -- REPEAT...UNTIL loop):
      IF (TRIM(tokenRead) == REACTION_END_CH) EXIT
   ENDDO ParseChemReactionString
   !-------------------------------

   chemReaction = reactionStr
END SUBROUTINE ParseChemReaction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION ChargeFromName(compName)
!-------------------------------------------------------------------------------
!  ChargeFromName
!  --------------
!
!  Returns the charge of the given chemical component name.
!  charge is derived from the component name.
!  EBNF defn.:
!    ComponentName = ChemName[Charge]
!         ChemName = Alpha{AlphaNum}
!           Charge = Sign | SignedInteger
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(IN) :: compName
!
!LOCALS
   INTEGER  :: charge
   CHARACTER(INP_LINE_LEN) :: dummyStr
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL ExtractChargeAndChemName(compName,dummyStr,charge)
   ChargeFromName = charge
END FUNCTION ChargeFromName
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ExtractChargeAndChemName(compName,chemName,charge)
!-------------------------------------------------------------------------------
!  ExtractChargeAndChemName
!  ------------------------
!
!  Given the component name extract the chemical name and
!  charge from this component name.
!  The input string ('compName') is NOT modified by this PROCEDURE.
!  EBNF:
!  ComponentName = ChemName[Charge]
!       ChemName = Alpha{AlphaNum}
!         Charge = Sign | SignedInteger
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(IN)  :: compName  !component name
   CHARACTER(*),INTENT(OUT) :: chemName  !chemical name
   INTEGER,INTENT(OUT)      :: charge    !the charge of the chem. component
!
!LOCALS
   INTEGER  :: i
   INTEGER  :: ioStatus

   INTEGER  :: len
   INTEGER  :: nameEndLocn

   LOGICAL :: isChemNameChar  !T <==> char is alphanumeric or a parenthesis

   CHARACTER(1) :: ch
!
!-------------------------------------------------------------------------------
!BEGIN
   ! -- Initialisations:
   charge = 0
   len    = LEN_TRIM(compName)

   ! -- Leading char. must be alphabetic:
   !!##EXPAND
   IF (.NOT. IsAlphaNumeric(compName(1:1))) THEN
      WRITE(*,'(1X,A,":  chem component name must start with a letter &
                                                   &or digit")')
      WRITE(*,'(1X,3X,"chem. name = ''",A,"''")') TRIM(compName)
      STOP "ExtractChargeAndChemName:  chem comp. must start with alphanum char"
   ENDIF

   isChemNameChar = .TRUE.
   i = 0
   DO WHILE (isChemNameChar  .AND.  i < len)
      i = i + 1
      isChemNameChar = IsAlphaNumeric(compName(i:i)) .OR. IsParen(compName(i:i))
   ENDDO
   ! -- ASSERT: i >= len  OR  compName(i:i) is not an alphanumeric AND not
   !            a left or right parenthesis.

   ! -- Note: Can build in a semantic check wrt parentheses.

   !!##EXPAND
   nameEndLocn = i-1   !initialisation
   ch = compName(i:i)
   IF (i >= len) THEN
      IF (IsAlphaNumeric(ch)  .OR.  ch == R_PAREN) THEN
         charge = 0
         nameEndLocn = i
      ELSE IF (ch == PLUS) THEN
         charge = +1
      ELSE IF (ch == MINUS) THEN
         charge = -1
      ELSE
         WRITE(*,'(1X,3X,"comp. name = ''",A,"''")') TRIM(compName)
         WRITE(*,'(1X,3X,"End character of chem. comp. name = ''",A,"''")') ch
         STOP "ExtractChargeAndChemName:  illegal character at end of chem. comp. name"
      ENDIF
   ELSE
      READ(compName(i: ),*,IOSTAT=ioStatus) charge  !read the charge
      IF (ioStatus /= 0) THEN
         WRITE(*,'(1X,A,":  error on read of charge from compName")')
         WRITE(*,'(1X,3X,"compName = ",A)') TRIM(compName)
         WRITE(*,'(1X,"compName(",I0,":",I0,") = ''",A,"''")') i, len, compName(i:len)
         STOP "ExtractChargeAndChemName:  error on read of charge from compName"
      ENDIF
   ENDIF

   READ(compName(1:nameEndLocn),*,IOSTAT=ioStatus) chemName   !read the chemical name
   IF (ioStatus /= 0) THEN
      WRITE(*,'(1X,A,":  error reading chemical name from component name")')
      STOP "ExtractChargeAndChemName:  error reading chemName from compName"
   ENDIF
END SUBROUTINE ExtractChargeAndChemName
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ProcessChemLibrary

END SUBROUTINE ProcessChemLibrary
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION IsLetter(ch)
   CHARACTER,INTENT(in) :: ch
!BEGIN
   IsLetter = (ch .ge. 'a' .and. ch .le. 'z') .or. (ch .ge. 'A' .and. ch .le. 'Z')
END FUNCTION IsLetter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION IsAlphaNumeric(ch)
   CHARACTER,INTENT(in) :: ch
!BEGIN
   IsAlphaNumeric = (ch .ge. 'a' .and. ch .le. 'z') .or. &
                    (ch .ge. 'A' .and. ch .le. 'Z') .or. &
                    (ch .ge. '0' .and. ch .le. '9')
END FUNCTION IsAlphaNumeric
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION IsParen(ch)
   CHARACTER,intent(in) :: ch

   IsParen = ( ch .eq. '(' .or. ch .eq. ')' )
END FUNCTION IsParen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
DOUBLETYPE FUNCTION antiLOG(loggedNumber, base)
   DOUBLETYPE,INTENT(in)  :: loggedNumber
   DOUBLETYPE,INTENT(in)  :: base
!BEGIN
   IF (loggedNumber < -40.0) THEN
     antiLOG = 0.000
     RETURN
   ENDIF

   antiLOG = base**loggedNumber
END FUNCTION antiLOG
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION f_get_lun( )
!-------------------------------------------------------------------------------
! Find the first free logical unit number
!-------------------------------------------------------------------------------
!LOCALS
    INTEGER  :: lun
    LOGICAL :: opened
!BEGIN
    DO lun = 10,99
      inquire(unit=lun, opened=opened)
      IF ( .NOT. opened ) THEN
        f_get_lun = lun
        RETURN
      ENDIF
    ENDDO
!
    f_get_lun = lun
END FUNCTION f_get_lun
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE AddComp(component)
!-------------------------------------------------------------------------------
!  AddComp
!  -------
!
!  Inserts a new chemical component into the data structure. The key for
!  insertion is the component name. If an element with this key already
!  exists the the program stops with an error message.
!
!  Search Key:  component name, i.e., component%CompName
!--------------------
!    The component must NOT already exist in this data structure.
!-------------------------------------------------------------------------------
   TYPE(gcUnknowns), INTENT(IN) :: component

   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ASSERT(compList%Exists)

   IF (compList%Last >= MAX_NUM_COMPS) THEN
      WRITE(*,'(1X,A,":  chem components data store is full.")')
      WRITE(*,'(1X,3X,"maximum number allowed = ",I0)') MAX_NUM_COMPS
      STOP "AddComp:  data store for chem components is full."
   !ELSE IF (Locn_Comp(component%CompName,COMP_NAME) /= NOT_FOUND) THEN
   ELSE IF (Locn_Comp(component%CompName) /= NOT_FOUND) THEN
      WRITE(*,'(1X,A,":  component named ''",A,"'' already in components &
                           &data store.")') TRIM(component%CompName)
      WRITE(*,'(1X,3X,"Check that you have not multiply defined this component")')
      WRITE(*,'(1X,3X,"in your chemical data input file.")')
      STOP "AddComp:  component to add is already in data store."
   ENDIF

   compList%Last                  = compList%Last + 1  !increment
   compList%Data_V(compList%Last) = component

   ! -- Assign correct updated component index:
   compList%Data_V(compList%Last)%CompIndex = compList%Last
END SUBROUTINE AddComp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE AddPhase(phase)
!-------------------------------------------------------------------------------
!  AddPhase
!  --------
!
!  Inserts a new chemical pure phase into the data structure. The key for
!  insertion is the phase element-name. If an element of the data structure
!  with this key already exists then the program stops with an error message.
!    The pure phase must NOT already exist in this data structure.
!-------------------------------------------------------------------------------
    TYPE(gcUnknowns),INTENT(inout) :: phase
!
!LOCALS
    INTEGER  :: status

    ! ------ Non-Local Variables ------------
    !  phaseList
    !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ASSERT(phaseList%Exists)

    ! -- Validation:
   IF (phaseList%Last >= MAX_NUM_PHASES) THEN
      WRITE(*,'(1X,A,":  chem pure phase data store is full.")')
      WRITE(*,'(1X,3X,"maximum number allowed = ",I0)') MAX_NUM_PHASES
      STOP "AddPhase:  data store for chem  is full."
   ELSE IF (FoundPhaseName(phase%EltName)) THEN
      WRITE(*,'(1X,A,":  phase named ''",A,"'' already in  &
                               &data store.")') TRIM(phase%EltName)
      WRITE(*,'(1X,3X,"Check that you have not multiply defined this phase")')
      WRITE(*,'(1X,3X,"in your chemical data input file.")')
      STOP "AddPhase:  phase to add is already in data store."
   ELSE IF (NumOfComps() < 0) THEN
      WRITE(*,'(1X,A,":  programming error: the PhaseStore module &
                                            &variable ''phaseList%Last''")')
      WRITE(*,'(1X,3X,"is negative.")')
      STOP "AddPhase:  phaseList%Last < 0"
   ENDIF


   phaseList%Last                   = phaseList%Last + 1   !update counter
   phaseList%Data_V(phaseList%Last) = phase

   ! -- Assign phase index:
   phaseList%Data_V(phaseList%Last)%CompIndex = phaseList%Last


   ! -- Allocate space for pure phase info:
   ALLOCATE(phaseList%Data_V(phaseList%Last)%PPData,STAT=status)
   _CheckAllocStatus(status)

   ! -- Assign the data
   phaseList%Data_V(phaseList%Last)%PPdata = phase%PPData

   ASSERT(NumOfComps() == SIZE(phase%PPData%Stoich(:)))

   ! -- Allocate space to store stoichiometric coefficients for phase:
   ALLOCATE(phaseList%Data_V(phaseList%Last)%PPData%Stoich(1:NumOfComps()),   &
                                                STAT=status)
   _CheckAllocStatus(status)

   phaseList%Data_V(phaseList%Last)%PPdata%Stoich(:) = phase%PPData%Stoich(:)
END SUBROUTINE AddPhase
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE AddSpecies(species)
!-------------------------------------------------------------------------------
!  AddSpecies
!  ----------
!
!  Inserts a new chemical species into the data structure. The key for
!  insertion is the species name. If an element with this key already
!  exists the the program stops with an error message.
!--------------------
!    The species must NOT already exist in this data structure.
!-------------------------------------------------------------------------------
   TYPE(gcSpecies), INTENT(IN OUT) :: species
!
!LOCAL
   INTEGER  :: status

   ! ------ Non-Local Variables ------------
   !  speciesList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ASSERT(speciesList%Exists)

   ! -- Validation:
   IF (speciesList%Last >= MAX_NUM_SPECIES) THEN
      WRITE(*,'(1X,A,":  chem species data store is full.")')
      WRITE(*,'(1X,3X,"maximum number allowed = ",I0)')          &
                                                     MAX_NUM_SPECIES
      STOP "AddSpecies:  data store for chem  is full."
   ELSE IF (FoundSpeciesName(species%Name)) THEN
      WRITE(*,'(1X,A,":  species named ''",A,"'' already in  &
                               &data store.")') TRIM(species%Name)
      WRITE(*,'(1X,3X,"Check that you have not multiply defined &
                                                            &this species")')
      WRITE(*,'(1X,3X,"in your chemical data input file.")')
      STOP "AddSpecies:  species to add is already in data store."
   ELSE IF (NumOfComps() < 0) THEN
      WRITE(*,'(1X,A,":  programming error: the ComponentsStore &
                                    &module variable ''speciesList%Last''")')
      WRITE(*,'(1X,3X,"is negative.")')
      STOP "AddSpecies:  speciesList%Last < 0"
   ENDIF

   speciesList%Last = speciesList%Last + 1
   speciesList%Data_V(speciesList%Last) = species

   ! -- Assign species index:
   speciesList%Data_V(speciesList%Last)%SpeciesIndex = speciesList%Last

   ASSERT(NumOfComps() == SIZE(species%Stoich(:)))

   ! -- Allocate space to store stoichiometric coefficients:
   ALLOCATE(speciesList%Data_V(speciesList%Last)%Stoich(1:NumOfComps()),      &
                                                        STAT=status)
   _CheckAllocStatus(status)

   speciesList%Data_V(speciesList%Last)%Stoich(:) = species%Stoich(:)   !assign stoich. coeffs.
END SUBROUTINE AddSpecies
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE AddToStoichiometry(compName,signedCoeff,stoichCoeffs_V)
!ARGUMENTS
   CHARACTER(*),INTENT(in)         :: compName       !chem comp name
   INTEGER,INTENT(in)              :: signedCoeff    !signed stoich coeff
   REAL,DIMENSION(:),INTENT(inout) :: stoichCoeffs_V !vector of stoich coeffs
!
!LOCALS
   ! ------ Local Variables ----------------
   INTEGER  :: compIndex    !index of this component in stoich vector
!
!-------------------------------------------------------------------------------
!BEGIN
   compIndex = Locn_Comp(compName)
   IF ( compIndex /= NOT_FOUND ) compIndex = compList%Data_V(compIndex)%CompIndex

   !-------
   ! -- Verify that component name is in the components data store:
   ASSERT( compIndex /= NOT_FOUND)
   !-------
   ! -- Ensure that the index value returned is within bounds:
   ASSERT( compIndex <= SIZE(stoichCoeffs_V))

   stoichCoeffs_V(compIndex) = stoichCoeffs_V(compIndex) + signedCoeff
END SUBROUTINE AddToStoichiometry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Expect(inpLine,expectedToken,token)
   CHARACTER(*),INTENT(inout)        :: inpLine         !line to read
   CHARACTER(*),INTENT(in)           :: expectedToken   !expected token
   CHARACTER(*),OPTIONAL,INTENT(out) :: token           !token read
!
!LOCALS
    CHARACTER(STR_LEN) :: symbol
!
!-------------------------------------------------------------------------------
!BEGIN
   ! -- Get leading token:
   CALL GetNextSymbol(inpLine,symbol)   !NB: this proc. removes leading symbol

   IF (TRIM(symbol) /= TRIM(expectedToken)) THEN
      WRITE(*,'(1X,A,":  read token does not match expected token")')
      WRITE(*,'(1X,3X,"expected token = ",A)') TRIM(expectedToken)
      WRITE(*,'(1X,3X,"    token read = ",A)') TRIM(symbol)

      STOP "Expect: expected token /= actual token read"
   ENDIF

   IF (PRESENT(token)) token = symbol
END SUBROUTINE Expect
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetNextSignificantLine(lun,inpLine)
!-------------------------------------------------------------------------------
!  Returns the next line which is not a comment line nor a blank line.
   INTEGER,INTENT(in)                  :: lun
   CHARACTER(INP_LINE_LEN),INTENT(out) :: inpLine
!
!LOCALS
   INTEGER  :: status

   CHARACTER(INP_LINE_LEN) :: leftJustLine
!
!-------------------------------------------------------------------------------
!BEGIN
   DO
      READ(lun,'(A)',IOSTAT=status) inpLine
      _CheckFileIOStatus(status)
      leftJustLine = ADJUSTL(inpLine)
      IF ( .NOT. (leftJustLine(1:1)=='#' .OR. &
                  leftJustLine(1:1)=='!' .OR. &
                  NumberOfFields(inpLine)==0) ) EXIT
   ENDDO
END SUBROUTINE GetNextSignificantLine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetNextSymbol(inpLine,nextSymb)
!-------------------------------------------------------------------------------
!  GetNextSymbol
!  -------------
!
!  Get the next symbol in the given input line. The input
!  line is then modified by removing the symbol and shifting
!  all non-white space characters as far to the left as possible
!  (i.e., shifting them all so that the left-most character is
!  in location 1). This prepares the input line for the next
!  symbol read.
!  The symbol read is returned left-justified.
!-------------------------------------------------------------------------------
   CHARACTER(*), INTENT(IN OUT) :: inpLine      !input line
   CHARACTER(*), INTENT(OUT)    :: nextSymb     !next symbol

   ! ------ Local Variables ----------------
   INTEGER  :: len             !num of characters in the symbol
   INTEGER  :: status

   !----------------------------------------

   inpLine = ADJUSTL(inpLine)
   len     = LastLocnOfField(inpLine,1)

   ! -- Read the next symbol:
   READ(inpLine,*,IOSTAT=status) nextSymb
   _CheckFileIOStatus(status)

   ! -- Left-justify line by overwriting left-most symbol:
   inpLine = inpLine(len+1:)      !left-justify
   inpLine = ADJUSTL(inpLine)     !remove leading blanks
END SUBROUTINE GetNextSymbol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ReadChemComp(reactionStr,compName,coeff,charge)
!-------------------------------------------------------------------------------
    CHARACTER(INP_LINE_LEN),INTENT(inout) :: reactionStr  !the chem eqn
   CHARACTER(STR_LEN),INTENT(out)         :: compName     !the chem component name
   INTEGER,INTENT(out)                    :: coeff        !stoichiometric coeff
   INTEGER,INTENT(out)                    :: charge       !charge of chem. comp.
!
!LOCALS
   INTEGER  :: startChemName   !posn. of 1st char of chem. name

   CHARACTER(1)       :: ch
    CHARACTER(STR_LEN) :: chemName
!
!-------------------------------------------------------------------------------
!BEGIN
   coeff = IVOID
   reactionStr = ADJUSTL(reactionStr)

   CALL GetNextSymbol(reactionStr,compName)

   CALL DetermineStoichCoeff(compName,coeff,startChemName)

   CALL ShiftStringLeft(compName,startChemName)

   ! -- Ensure that first character is a letter:
   ch = compName(1:1)   !read 1st letter of name
   IF (.NOT. IsLetter(ch)) THEN
      WRITE(*,'(1X,A,":  a letter of the alphabet was expected as the 1st &
                                              &character of ")')
      WRITE(*,'(1X,3X,"a chemical component''s name.")')
      IF (IsHorizWhiteSpace(ch)) THEN
         WRITE(*,'(1X,3X,"It may be that the stoichiometric coefficient &
                                           &is separated from the chemical")')
         WRITE(*,'(1X,3X,"name by a space. This is not allowed.")')
      ENDIF
      !!$WRITE(*,'(1X,3X,"      compName = ''",A,"''")')  TRIM(compName)
      WRITE(*,'(1X,3X,"character read = ''",A,"''")') ch
      STOP "ReadChemComp:  leading character of chemical component name &
                                                      &must be a letter."
   ENDIF

   CALL ExtractChargeAndChemName(compName,chemName,charge)
END SUBROUTINE ReadChemComp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ReadDataTrailingChemReaction(inpString,logKat25,deltaH,settlingVel)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),         INTENT(IN)  :: inpString
   DOUBLETYPE,           INTENT(OUT) :: logKat25
   DOUBLETYPE,           INTENT(OUT) :: deltaH
   DOUBLETYPE, OPTIONAL, INTENT(OUT) :: settlingVel
!
!LOCALS
   INTEGER  :: status
!
!-------------------------------------------------------------------------------
!BEGIN
   logKat25 = VOID
   deltaH   = VOID

   IF (PRESENT(settlingVel)) THEN
      settlingVel = VOID
      READ(inpString,*,IOSTAT=status) logKat25, deltaH, settlingVel
      IF (status /= 0) STOP "ReadDataTrailingChemReaction:  error on &
                                 &read of logKat25, deltaH and settlingVel"
   ELSE
      READ(inpString,*,IOSTAT=status) logKat25, deltaH
      IF (status /= 0) STOP "ReadDataTrailingChemReaction:  error on &
                                              &read of logKat25 and deltaH"
   ENDIF
END SUBROUTINE ReadDataTrailingChemReaction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ReadLeadingSymbol(inpLine,symb)
!-------------------------------------------------------------------------------
!ARGUMENTS
    CHARACTER(*), INTENT(IN OUT) :: inpLine      !input line
    CHARACTER(*), INTENT(OUT)    :: symb         !next symbol
!
!LOCALS
    INTEGER  :: status
!
!-------------------------------------------------------------------------------
!BEGIN
    ! -- Read the next symbol:
    READ(inpLine,*,IOSTAT=status) symb
    _CheckFileIOStatus(status)
END SUBROUTINE ReadLeadingSymbol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION LastLocnOfField(line,startLocn)
!-------------------------------------------------------------------------------
!  LastLocnOfField
!  ---------------
!
!  Returns the position of the last character of a given
!  field  The field is defined by its start location
!  ('startLocn'), i.e., the index value of the start
!  character in its input line.
!--------------------
!       Creation Date:  2005-02-09
!         Last Change:
!              Author:  Alan Imerito
!--------------------
!  Modifications:
!  <date>;<change_author>:          <version>
!  - <description of change>
!-------------------------------------------------------------------------------
   CHARACTER(*), INTENT(IN) :: line        !the line
   INTEGER,      INTENT(IN) :: startLocn   !locn. of 1st char of field
!
!LOCALS
   INTEGER  :: i            !character loop counter
   INTEGER  :: strLen       !number of characters in the string

   LOGICAL :: found        !working Boolean variable (because Fortran
                           !..does not do short-circuit evaluations)
!
!-------------------------------------------------------------------------------
!BEGIN
   strLen = LEN(line)

   IF (strLen < startLocn) THEN
      WRITE(*,'(1X,3X,A,":  given start location > string length")')
      WRITE(*,'(1X,3X,2X,"strLen = ",I0,3X,"startLocn = ",I0)') strLen, startLocn
      STOP "LastLocnOfField:  start locn > str len"
   ELSE IF (startLocn <= 0) THEN
      WRITE(*,'(1X,3X,A,":  given start location <= 0")')
      WRITE(*,'(1X,3X,2X,"startLocn = ",I0)') startLocn
      STOP "LastLocnOfField:  start locn <= 0"
   ENDIF

   i = startLocn
   found = .FALSE.
   DO WHILE  (i < strLen  .AND.  .NOT. found)
      i = i + 1
      found = IsHorizWhiteSpace(line(i:i))
   ENDDO

   IF (found) THEN
      LastLocnOfField = i - 1
   ELSE
      LastLocnOfField = i
   ENDIF
END FUNCTION LastLocnOfField
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE INTEGER FUNCTION NumberOfFields(inpLine)
!-------------------------------------------------------------------------------
!  NumberOfFields
!  --------------
!
!  Returns the number of fields in a line of text. Blanks
!  and tabs are considered to be the field delimeters.
!
!--------------------
!       Creation Date:  2004-02-29
!         Last Change:
!              Author:  Alan Imerito
!--------------------
!  Modifications:
!  <date>;<change_author>:          <version>
!  - <description of change>
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*), INTENT(IN) :: inpLine     !the line
!
!LOCALS
   INTEGER  :: i            !character loop counter
   INTEGER  :: strLen       !number of characters in the string
   INTEGER  :: numFields    !number of fields in the input string

   LOGICAL :: inString     !T <==> if we have just been in a string
!
!-------------------------------------------------------------------------------
!BEGIN
   strLen = LEN(inpLine)

   numFields = 0
   inString = .FALSE.
   DO i = 1, strLen
      IF ( IsHorizWhiteSpace(inpLine(i:i)) .AND. inString ) THEN
         inString = .FALSE.
      ELSE IF (.NOT. IsHorizWhiteSpace(inpLine(i:i)) .AND. .NOT. inString) THEN
         inString = .TRUE.
         numFields = numFields + 1
      ENDIF
   ENDDO

   NumberOfFields = numFields
END FUNCTION NumberOfFields
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE DetermineStoichCoeff(compName,coeff,startChemName)
!-------------------------------------------------------------------------------
!  DetermineStoichCoeff
!  --------------------
!
!  Read the leading stoichiometric coefficient of a chemical
!  species. If there is no leading integer then 1 will be
!  returned as the coefficient.
!  NB: the coefficient returned is always +ve.
!
!    compName is left-justified.
!-------------------------------------------------------------------------------
   CHARACTER(*), INTENT(IN)  :: compName      !chemical component
   INTEGER,      INTENT(OUT) :: coeff         !stoichiometric coeff.
   INTEGER,      INTENT(OUT) :: startChemName !locn. of 1st char of chemical
!                                             !..chem. compound name
!LOCALS
   INTEGER  :: status
   INTEGER  :: i               !char loop counter
   INTEGER  :: coeffEnd        !last posn. of coefficient in string
   INTEGER  :: len             !num of characters of the chem component

   LOGICAL :: chIsDigit       !TRUE <==> ch is a digit

   CHARACTER(INP_LINE_LEN) :: chemStr   !chemical name + coeff
!
!-------------------------------------------------------------------------------
!BEGIN
   chemStr = compName
   len     = LEN_TRIM(chemStr)

   ! -- Check leading character:
   IF (.NOT. IsAlphaNumeric(chemStr(1:1))) THEN
      STOP "DetermineStoichCoeff:  alpha-numeric character expected."
   ENDIF

   !-------
   i = 1
   DO WHILE ( IsDigit(chemStr(i:i))  .AND. i < len)
      i = i + 1
   ENDDO
   chIsDigit = IsDigit(chemStr(i:i))

   !-------
   IF (.NOT. chIsDigit  .AND.  i == 1) THEN
      ! -- i.e., there is no leading stoich. coeff., so we interpret
      ! -- this as meaning that the stoich. coeff. is 1:
      coeffEnd      = 0
      startChemName = 1
      coeff         = 1
   ELSE
      IF (.NOT. chIsDigit) THEN
         coeffEnd = i-1
         startChemName = i
      ELSE
         ! ** ASSERT: i == len
         coeffEnd = i
         startChemName = i+1
      ENDIF
      READ(chemStr(1:coeffEnd),*,IOSTAT=status) coeff
      _CheckFileIOStatus(status)
   ENDIF
END SUBROUTINE DetermineStoichCoeff
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ShiftStringLeft(string,startLocn)
!-------------------------------------------------------------------------------
!  ShiftStringLeft
!  ---------------
!
!  Shifts the string starting at locn. 'startLocn' to the
!  left-most position of the string variable. Any characters
!  pre-existing characters are removed.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(INOUT) :: string    !string var. to be operated on
   INTEGER,INTENT(IN)         :: startLocn !posn. of 1st non-white space
!                                             !..char in string.
!-------------------------------------------------------------------------------
!BEGIN
    ! -- Left-justify string starting at 'startLocn':
    string = string(startLocn:)
END SUBROUTINE ShiftStringLeft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!# The Data Store Routines                                                     #
!###############################################################################

!###############################################################################
SUBROUTINE CreateCompStore()
!-------------------------------------------------------------------------------
!  CreateCompStore
!  ---------------
!
!  Creates the data structure for storing the chemical components info.
!-------------------------------------------------------------------------------
!LOCALS
   INTEGER  :: status

   ! ------ Non-Local Variables ------------
   !  compList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
    ASSERT(.NOT. compList%Exists)

    ALLOCATE(compList%Data_V(1:MAX_NUM_COMPS),STAT=status)
    _CheckAllocStatus(status)

    compList%Exists = .TRUE.
END SUBROUTINE CreateCompStore
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE CreatePhasesStore()
!-------------------------------------------------------------------------------
!  CreatePhasesStore
!  ------------------
!
!  Creates the data structure for storing the chemical pure phases info.
!
!-------------------------------------------------------------------------------
!LOCALS
   INTEGER  :: status

   ! ------ Non-Local Variables ------------
   !  phaseList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
    ASSERT(.NOT. phaseList%Exists)

    ALLOCATE(phaseList%Data_V(1:MAX_NUM_PHASES),STAT=status)
    _CheckAllocStatus(status)

    phaseList%Exists = .TRUE.
END SUBROUTINE CreatePhasesStore
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE CreateSpeciesStore()
!-------------------------------------------------------------------------------
!  CreateSpeciesStore
!  ------------------
!
!  Creates the data structure for storing the chemical species info.
!-------------------------------------------------------------------------------
!LOCALS
   INTEGER  :: status

   ! ------ Non-Local Variables ------------
   !  speciesList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ASSERT(.NOT. speciesList%Exists)

   ALLOCATE(speciesList%Data_V(1:MAX_NUM_SPECIES),STAT=status)
   _CheckAllocStatus(status)

   speciesList%Exists = .TRUE.
END SUBROUTINE CreateSpeciesStore
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION Locn_Species(speciesName)
!-------------------------------------------------------------------------------
!  Locn_Species
!  ------------
!
!  Returns the index location in speciesList%Data_V of the species named
!  <speciesName> if it is in the data store; ELSE returns value <NOT_FOUND>.
!
!-------------------------------------------------------------------------------
!ARGUMENTS
    CHARACTER(*), INTENT(IN)  :: speciesName
!
!LOCALS
   INTEGER  :: i
   LOGICAL :: found

   ! ------ Non-Local Variables ------------
   !  speciesList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ! -- Pre-condition: data structure must exist.
   ASSERT(speciesList%Exists)

   !-------
   ! -- Search for key until found or all of data store checked:
   ! -- (Linear search on unsorted data for now. Replace with better later.
   ! --  This is an O(n) algorithm.)
   i = 0
   found = .FALSE.
   DO WHILE (i < speciesList%Last  .AND.  .NOT. found)
      i = i+1
      found = speciesName == speciesList%Data_V(i)%Name
   ENDDO

   !-------
   ! -- Return appropriate value for location:
   IF (found) THEN
      Locn_Species = i
   ELSE
      Locn_Species = NOT_FOUND
   ENDIF

 END FUNCTION Locn_Species
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION Locn_Phase(phaseName)
!-------------------------------------------------------------------------------
!  Locn_Phase
!  ----------
!
!  Returns the index location in phaseList%Data_V of the phase which has the
!  element name <phaseName> if it is in the data store; ELSE returns value
!  <NOT_FOUND>.
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*), INTENT(IN)  :: phaseName
!
!LOCALS
   INTEGER  :: i
   LOGICAL :: found

   ! ------ Non-Local Variables ------------
   !  phaseList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ! -- Pre-condition: data structure must exist.
   ASSERT(phaseList%Exists)

   !-------
   ! -- Search for key until found or all of data store checked:
   ! -- (Linear search on unsorted data for now. Replace with better later.
   ! --  This is an O(n) algorithm.)
   i = 0
   found = .FALSE.
   DO WHILE (i < phaseList%Last  .AND.  .NOT. found)
      i = i+1
      found = phaseName == phaseList%Data_V(i)%EltName
   ENDDO

   !-------
   ! -- Return appropriate value for location:
   IF (found) THEN
      Locn_Phase = i
   ELSE
      Locn_Phase = NOT_FOUND
   ENDIF
 END FUNCTION Locn_Phase
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION Locn_Comp(name)!,key_type)
!-------------------------------------------------------------------------------
!  Locn_Comp
!  ---------
!
!  Returns the index location in compList%Data_V of the component named
!  <compName> if it is in the data store; ELSE returns value <NOT_FOUND>.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),       INTENT(IN) :: name       !name to search for
 ! INTEGER(NAME_TYPE), INTENT(IN) :: key_type   !search key (name) name type
!
!LOCALS
   INTEGER  :: i
   LOGICAL :: found

   ! ------ Non-Local Variables ------------
   !  compList
   !----------------------------------------
!
!-------------------------------------------------------------------------------
!BEGIN
   ! -- Pre-condition: data structure must exist.
   ASSERT(compList%Exists)

   !-------
   ! -- Search for key until found or all of data store checked:
   ! -- (Linear search on unsorted data for now. Replace with better later.
   ! --  This is an O(n) algorithm.)
   i = 0
   found = .FALSE.
!# key_type was only ever COMP_NAME
!  IF (key_type == COMP_NAME) THEN
      DO WHILE (i < compList%Last  .AND.  .NOT. found)
         i = i+1
         found = (name == compList%Data_V(i)%CompName)
      ENDDO

#if 0
   ELSE IF (key_type == ELT_NAME) THEN
      DO WHILE (i < compList%Last  .AND.  .NOT. found)
         i = i+1
         found = (name == compList%Data_V(i)%EltName)
      ENDDO

   ELSE
      WRITE(*,'(1X, A,":  programming error. Search key name &
                                     &type is unrecognized.")')
      WRITE(*,'(1X,3X,"key_type = ",I0)')  key_type
      STOP "Locn_Comp:  search key name type is unrecognized."
   ENDIF
#endif

   !-------
   ! -- Return appropriate value for location:
   IF (found) THEN
      Locn_Comp = i
   ELSE
      Locn_Comp = NOT_FOUND
   ENDIF
  END FUNCTION Locn_Comp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if 0
!###############################################################################
!========================================================================
!  GetCompIndex
!  ------------
!
!  Returns the component index value corresponding to the given component
!  name.
!========================================================================
SUBROUTINE GetCompIndex(compName,compIndex)
    CHARACTER(*), INTENT(IN)  :: compName
    INTEGER,      INTENT(OUT) :: compIndex

    ! ------ Local Variables ----------------
    INTEGER  :: i

    i = Locn_Comp(compName,COMP_NAME)

    IF (i == NOT_FOUND) THEN  !component name NOT found in data store
      compIndex = NOT_FOUND
    ELSE                      !component name found in data store
      compIndex = compList%Data_V(i)%CompIndex
    ENDIF

END SUBROUTINE GetCompIndex
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


!###############################################################################
INTEGER FUNCTION NumOfComps()
   NumOfComps = (compList%Last)
END FUNCTION NumOfComps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION NumOfSpecies()
   NumOfSpecies = (speciesList%Last)
END FUNCTION NumOfSpecies
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetCompCompName(index, compName)
   INTEGER,      INTENT(IN)  :: index
   CHARACTER(*), INTENT(OUT) :: compName
!
!-------------------------------------------------------------------------------
!BEGIN
    compName = compList%Data_V(index)%CompName
END SUBROUTINE GetCompCompName
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetSpeciesName(index, speciesName)
   INTEGER,      INTENT(IN)  :: index
   CHARACTER(*), INTENT(OUT) :: speciesName
!
!-------------------------------------------------------------------------------
!BEGIN
   speciesName = speciesList%Data_V(index)%Name
END SUBROUTINE GetSpeciesName
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetCompInfo(comp)
   TYPE(gcUnknowns),INTENT(IN OUT) :: comp    !component data
   INTEGER  :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i = 1,compList%Last
      IF ( compList%Data_V(i)%EltName == comp%EltName ) THEN
         comp = compList%Data_V(i)
         RETURN
      ENDIF
   ENDDO
END SUBROUTINE GetCompInfo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetSpeciesInfo(index,species)
   INTEGER,         INTENT(IN)     :: index   !component index search key
   TYPE(gcSpecies), INTENT(IN OUT) :: species !component data
   species = speciesList%Data_V(index)
END SUBROUTINE GetSpeciesInfo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GetPhaseInfo(phase)
!-------------------------------------------------------------------------------
   TYPE(gcUnknowns),INTENT(IN OUT) :: phase   !component data
   INTEGER  :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i = 1,phaseList%Last
      IF ( phaseList%Data_V(i)%EltName == phase%EltName ) THEN
          phase = phaseList%Data_V(i)
          RETURN
      ENDIF
   ENDDO
END SUBROUTINE GetPhaseInfo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE printTables
!-------------------------------------------------------------------------------
   INTEGER  :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"GEOCHEMISTRY data:"
   print *,"components list:"
   DO i = 1,compList%Last
!     print *, ' ',compList%Data_V(i)%EltName,compList%Data_V(i)%CompName,compList%Data_V(i)%Charge
      print *,'----> :     ',i
      print *,' CompName:  ',compList%Data_V(i)%CompName
      print *,' EltName:   ',compList%Data_V(i)%EltName
      print *,' CompType:  ',compList%Data_V(i)%CompType
      print *,' CompIndex: ',compList%Data_V(i)%CompIndex
      print *,' eqIndex:   ',compList%Data_V(i)%eqIndex
      print *,' wqIndex:   ',compList%Data_V(i)%wqIndex
      print *,' delta:     ',compList%Data_V(i)%delta
      print *,' Value:     ',compList%Data_V(i)%Value
      print *,' Total:     ',compList%Data_V(i)%Total
 !    IF(compList%Data_V(i)%CompType==PUREPHASE) THEN
      IF(compList%Data_V(i)%CompType==1) THEN
        print *,' ppData%logK  :  ',compList%Data_V(i)%ppData%logK
        print *,' ppData%moles :  ',compList%Data_V(i)%ppData%moles
        print *,' ppData%stoich:  ',compList%Data_V(i)%ppData%stoich
      ENDIF
      print *,'  --<'
   ENDDO

   print *,"species list:"
   DO i = 1,speciesList%Last
!     print *, " ",speciesList%Data_V(i)%Name
      print *,'----> :    '
      print *,' SpeciesIndex: ',speciesList%Data_V(i)%SpeciesIndex
      print *,' Name:     ',speciesList%Data_V(i)%Name
      print *,' Charge:   ',speciesList%Data_V(i)%Charge
      print *,' logKat25  ',speciesList%Data_V(i)%logKat25
      print *,' deltaH:   ',speciesList%Data_V(i)%deltaH
      print *,' Gflag:    ',speciesList%Data_V(i)%Gflag
      print *,' H2OStoich:',speciesList%Data_V(i)%H2OStoich
      print *,' HStoich:  ',speciesList%Data_V(i)%HStoich
      print *,' Moles:    ',speciesList%Data_V(i)%Moles
      print *,' log(act): ',speciesList%Data_V(i)%logActivity
      print *,' loggamma: ',speciesList%Data_V(i)%logGamma
      print *,' delgamma: ',speciesList%Data_V(i)%delGamma
      print *,' stoich:   ',INT(speciesList%Data_V(i)%stoich)
      print *,'  --<'
   ENDDO

   print *,"phase list:"
   DO i = 1,phaseList%Last
      print *, " ",phaseList%Data_V(i)%EltName,phaseList%Data_V(i)%CompName
   ENDDO
END SUBROUTINE printTables
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_gclib
!===============================================================================
