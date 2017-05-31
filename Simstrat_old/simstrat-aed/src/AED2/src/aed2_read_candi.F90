MODULE aed2_read_candi

   USE aed2_gctypes

   IMPLICIT NONE

   PUBLIC read_sed_candi_params

CONTAINS


!###############################################################################
SUBROUTINE read_sed_candi_params(init_values_file, scp, count)
!-------------------------------------------------------------------------------
   USE aed2_csv_reader
!
!ARGUMENTS
   CHARACTER(*) :: init_values_file
   TYPE(aed2_sed_candi_param_t),DIMENSION(:),ALLOCATABLE,INTENT(inout) :: scp
   INTEGER, INTENT(in) :: count
!-------------------------------------------------------------------------------
!LOCALS
   INTEGER :: unit, nccols, ccol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0, numv = 0, numd = 0, zon
   LOGICAL :: meh
!
!BEGIN
!-------------------------------------------------------------------------------
   unit = aed_csv_read_header(init_values_file, csvnames, nccols)
   print *,'benthic vars initialised from file : ', csvnames
   IF (unit <= 0) RETURN !# No file found
   DO ccol=1,nccols
      IF ( csvnames(ccol) == "zone" ) THEN
         idx_col = ccol
         EXIT
      ENDIF
   ENDDO

   ALLOCATE(scp(count))
   ALLOCATE(values(nccols))
   CALL default_sed_parm(scp)

   zon = 0
   IF ( idx_col > 0 ) THEN
      DO WHILE ( aed_csv_read_row(unit, values) )
         zon = zon + 1
         DO ccol=1,nccols
            SELECT CASE (csvnames(ccol))
                     !-- Sediment Physical Properties --!
               CASE ('db0')                ; scp(zon)%db0                = extract_double(values(ccol))
               CASE ('imix')               ; scp(zon)%imix               = extract_integer(values(ccol))
               CASE ('xs')                 ; scp(zon)%xs                 = extract_double(values(ccol))
               CASE ('x1')                 ; scp(zon)%x1                 = extract_double(values(ccol))
               CASE ('x2')                 ; scp(zon)%x2                 = extract_double(values(ccol))
               CASE ('irrg')               ; scp(zon)%irrg               = extract_integer(values(ccol))
               CASE ('alpha0')             ; scp(zon)%alpha0             = extract_double(values(ccol))
               CASE ('xirrig')             ; scp(zon)%xirrig             = extract_double(values(ccol))
               CASE ('ventflow')           ; scp(zon)%ventflow           = extract_double(values(ccol))
               CASE ('w00')                ; scp(zon)%w00                = extract_double(values(ccol))
               CASE ('p0')                 ; scp(zon)%p0                 = extract_double(values(ccol))
               CASE ('p00')                ; scp(zon)%p00                = extract_double(values(ccol))
               CASE ('bp')                 ; scp(zon)%bp                 = extract_double(values(ccol))
               CASE ('torteq')             ; scp(zon)%torteq             = extract_integer(values(ccol))
               CASE ('an')                 ; scp(zon)%an                 = extract_double(values(ccol))
               CASE ('aa')                 ; scp(zon)%aa                 = extract_double(values(ccol))
               CASE ('ab')                 ; scp(zon)%ab                 = extract_double(values(ccol))
               CASE ('xl')                 ; scp(zon)%xl                 = extract_double(values(ccol))
               CASE ('maxnpts')            ; scp(zon)%maxnpts            = extract_integer(values(ccol))

                     !-- Biogeochemical Configuration --!
               CASE ('numOM')              ; scp(zon)%numOM              = extract_integer(values(ccol))
               CASE ('simDOM')             ; scp(zon)%simDOM             = extract_logical(values(ccol))

               CASE ('OMapproach')         ; scp(zon)%OMapproach         = extract_integer(values(ccol))
               CASE ('OMModel')            ; scp(zon)%OMModel            = extract_integer(values(ccol))  ! Dan added
               CASE ('FTemswitch')         ; scp(zon)%FTemswitch         = extract_integer(values(ccol))  ! Dan added
               CASE ('FTswitch')           ; scp(zon)%FTswitch           = extract_integer(values(ccol))  ! Dan added
               CASE ('FBIOswitch')         ; scp(zon)%FBIOswitch         = extract_integer(values(ccol))  ! Dan added
               CASE ('FINswitch')          ; scp(zon)%FINswitch          = extract_integer(values(ccol))  ! Dan added
               CASE ('FInO2OnlySwitch')    ; scp(zon)%FInO2OnlySwitch    = extract_integer(values(ccol))
               CASE ('FOMswitch')          ; scp(zon)%FOMswitch          = extract_integer(values(ccol))  ! Dan added
               CASE ('Bsolidswitch')       ; scp(zon)%Bsolidswitch       = extract_integer(values(ccol))  ! Dan added

               CASE ('simMnFe')            ; scp(zon)%simMnFe            = extract_logical(values(ccol))
               CASE ('simFeS')             ; scp(zon)%simFeS             = extract_logical(values(ccol))
               CASE ('simX')               ; scp(zon)%simX               = extract_logical(values(ccol))
               CASE ('simCaCO3')           ; scp(zon)%simCaCO3           = extract_logical(values(ccol))
               CASE ('simFeCO3')           ; scp(zon)%simFeCO3           = extract_logical(values(ccol))
               CASE ('simMnCO3')           ; scp(zon)%simMnCO3           = extract_logical(values(ccol))
               CASE ('simNPAds')           ; scp(zon)%simNPAds           = extract_logical(values(ccol))

                     !-- Reaction Rates --!                              !Rate constants - Dan
                     !-- Organic matter decomposition
                     !-- OMModel 1
               CASE ('poml2dic')           ; scp(zon)%poml2dic           = extract_double(values(ccol))   ! Dan added
               CASE ('pomr2dic')           ; scp(zon)%pomr2dic           = extract_double(values(ccol))   ! Dan added
               CASE ('pomspecial2dic')     ; scp(zon)%pomspecial2dic     = extract_double(values(ccol))   ! Dan added

                     !-- OMModel 2
               CASE ('docl2dic')           ; scp(zon)%docl2dic           = extract_double(values(ccol))
               CASE ('donl2din')           ; scp(zon)%donl2din           = extract_double(values(ccol))
               CASE ('dopl2dip')           ; scp(zon)%dopl2dip           = extract_double(values(ccol))
               CASE ('pocl2docl')          ; scp(zon)%pocl2docl          = extract_double(values(ccol))
               CASE ('ponl2donl')          ; scp(zon)%ponl2donl          = extract_double(values(ccol))
               CASE ('popl2dopl')          ; scp(zon)%popl2dopl          = extract_double(values(ccol))
               CASE ('docr2docl')          ; scp(zon)%docr2docl          = extract_double(values(ccol))
               CASE ('donr2donl')          ; scp(zon)%donr2donl          = extract_double(values(ccol))
               CASE ('dopr2dopl')          ; scp(zon)%dopr2dopl          = extract_double(values(ccol))
               CASE ('pocr2docr')          ; scp(zon)%pocr2docr          = extract_double(values(ccol))
               CASE ('ponr2donr')          ; scp(zon)%ponr2donr          = extract_double(values(ccol))
               CASE ('popr2dopr')          ; scp(zon)%popr2dopr          = extract_double(values(ccol))
               CASE ('pocvr2docr')         ; scp(zon)%pocvr2docr         = extract_double(values(ccol))
               CASE ('ponvr2donr')         ; scp(zon)%ponvr2donr         = extract_double(values(ccol))
               CASE ('popvr2dopr')         ; scp(zon)%popvr2dopr         = extract_double(values(ccol))

                     !-- OMModel 3
 !##           CASE ('dHyd2dic')           ; scp(zon)%dHyd2dic           = extract_double(values(ccol))   ! Dan added
 !##           CASE ('dFer2dic             ; scp(zon)%dFer2dic           = extract_double(values(ccol))   ! Dan added
 !##           CASE ('poml2dic             ; scp(zon)%poml2dic           = extract_double(values(ccol))   ! Dan added
               CASE ('domr2dic')           ; scp(zon)%domr2dic           = extract_double(values(ccol))   ! Dan added
 !##           CASE ('dHyd2dFer            ; scp(zon)%dHyd2dFer          = extract_double(values(ccol))   ! Dan added
 !##           CASE ('dHyd2domr')          ; scp(zon)%dHyd2domr          = extract_double(values(ccol))   ! Dan added
               CASE ('domr2pomr')          ; scp(zon)%domr2pomr          = extract_double(values(ccol))   ! Dan added
               CASE ('poml2doml')          ; scp(zon)%poml2doml          = extract_double(values(ccol))   ! Dan added
 !##           CASE ('dFer2domr            ; scp(zon)%poml2doml          = extract_double(values(ccol))   ! Dan added

                     ! Stoichiometric coefficients
                     ! OM model 1

               CASE ('FTR')                ; scp(zon)%FTR                = extract_double(values(ccol))   ! Dan added
               CASE ('FTT')                ; scp(zon)%FTT                = extract_double(values(ccol))   ! Dan added
               CASE ('deltaGATP')          ; scp(zon)%deltaGATP          = extract_double(values(ccol))   ! Dan added
               CASE ('dG0FerDHyd')         ; scp(zon)%dG0FerDHyd         = extract_double(values(ccol))
               CASE ('dG0AerDHy')          ; scp(zon)%dG0AerDHy          = extract_double(values(ccol))
               CASE ('dG0AerOAc')          ; scp(zon)%dG0AerOAc          = extract_double(values(ccol))
               CASE ('dG0DenDHy')          ; scp(zon)%dG0DenDHy          = extract_double(values(ccol))
               CASE ('dG0DenOAc')          ; scp(zon)%dG0DenOAc          = extract_double(values(ccol))
               CASE ('dG0DenH2')           ; scp(zon)%dG0DenH2           = extract_double(values(ccol))
               CASE ('dG0ManOAc')          ; scp(zon)%dG0ManOAc          = extract_double(values(ccol))
               CASE ('dG0IroOAc')          ; scp(zon)%dG0IroOAc          = extract_double(values(ccol))
               CASE ('dG0IroH2')           ; scp(zon)%dG0IroH2           = extract_double(values(ccol))
               CASE ('dG0SulOAc')          ; scp(zon)%dG0SulOAc          = extract_double(values(ccol))
               CASE ('dG0SulH2')           ; scp(zon)%dG0SulH2           = extract_double(values(ccol))
               CASE ('dG0MetOAc')          ; scp(zon)%dG0MetOAc          = extract_double(values(ccol))
               CASE ('dG0MetH2')           ; scp(zon)%dG0MetH2           = extract_double(values(ccol))

               CASE ('kgrowthFer')         ; scp(zon)%kgrowthFer         = extract_double(values(ccol))
               CASE ('kgrowthAer')         ; scp(zon)%kgrowthAer         = extract_double(values(ccol))
               CASE ('kgrowthDen')         ; scp(zon)%kgrowthDen         = extract_double(values(ccol))
               CASE ('kgrowthMan')         ; scp(zon)%kgrowthMan         = extract_double(values(ccol))
               CASE ('kgrowthIro')         ; scp(zon)%kgrowthIro         = extract_double(values(ccol))
               CASE ('kgrowthSul')         ; scp(zon)%kgrowthSul         = extract_double(values(ccol))
               CASE ('kgrowthMet')         ; scp(zon)%kgrowthMet         = extract_double(values(ccol))
               CASE ('kdeathFer')          ; scp(zon)%kdeathFer          = extract_double(values(ccol))
               CASE ('kdeathAer')          ; scp(zon)%kdeathAer          = extract_double(values(ccol))
               CASE ('kdeathDen')          ; scp(zon)%kdeathDen          = extract_double(values(ccol))
               CASE ('kdeathMan')          ; scp(zon)%kdeathMan          = extract_double(values(ccol))
               CASE ('kdeathIro')          ; scp(zon)%kdeathIro          = extract_double(values(ccol))
               CASE ('kdeathSul')          ; scp(zon)%kdeathSul          = extract_double(values(ccol))
               CASE ('kdeathMet')          ; scp(zon)%kdeathMet          = extract_double(values(ccol))

               CASE ('kHyd1')              ; scp(zon)%kHyd1              = extract_double(values(ccol))
               CASE ('kHyd2')              ; scp(zon)%kHyd2              = extract_double(values(ccol))
               CASE ('kHyd3')              ; scp(zon)%kHyd3              = extract_double(values(ccol))
               CASE ('kHyd4')              ; scp(zon)%kHyd4              = extract_double(values(ccol))
               CASE ('kHydN')              ; scp(zon)%kHydN              = extract_double(values(ccol))
               CASE ('CellWeight')         ; scp(zon)%CellWeight         = extract_double(values(ccol))
               CASE ('fuse')               ; scp(zon)%fuse               = extract_double(values(ccol))
               CASE ('Tiny')               ; scp(zon)%Tiny               = extract_double(values(ccol))
               CASE ('Temporary_proton')   ; scp(zon)%Temporary_proton   = extract_double(values(ccol))
               CASE ('BMax')               ; scp(zon)%BMax               = extract_double(values(ccol))
               CASE ('e')                  ; scp(zon)%e                  = extract_double(values(ccol))
               CASE ('F')                  ; scp(zon)%F                  = extract_double(values(ccol))
               CASE ('n')                  ; scp(zon)%n                  = extract_double(values(ccol))
               CASE ('dPsi')               ; scp(zon)%dPsi               = extract_double(values(ccol))
               CASE ('KDHyd')              ; scp(zon)%KDHyd              = extract_double(values(ccol))
               CASE ('KOAc')               ; scp(zon)%KOAc               = extract_double(values(ccol))
               CASE ('KH2')                ; scp(zon)%KH2                = extract_double(values(ccol))
               CASE ('YDHyAer')            ; scp(zon)%YDHyAer            = extract_double(values(ccol))
               CASE ('YDHyFer')            ; scp(zon)%YDHyFer            = extract_double(values(ccol))
               CASE ('YDenDHy')            ; scp(zon)%YDenDHy            = extract_double(values(ccol))
               CASE ('YAerOAc')            ; scp(zon)%YAerOAc            = extract_double(values(ccol))
               CASE ('YDenOAc')            ; scp(zon)%YDenOAc            = extract_double(values(ccol))
               CASE ('YDenH2')             ; scp(zon)%YDenH2             = extract_double(values(ccol))
               CASE ('YManOAc')            ; scp(zon)%YManOAc            = extract_double(values(ccol))
               CASE ('YIroOAc')            ; scp(zon)%YIroOAc            = extract_double(values(ccol))
               CASE ('YIroH2')             ; scp(zon)%YIroH2             = extract_double(values(ccol))
               CASE ('YSulOAc')            ; scp(zon)%YSulOAc            = extract_double(values(ccol))
               CASE ('YSulH2')             ; scp(zon)%YSulH2             = extract_double(values(ccol))
               CASE ('YMetOAc')            ; scp(zon)%YMetOAc            = extract_double(values(ccol))
               CASE ('YMetH2')             ; scp(zon)%YMetH2             = extract_double(values(ccol))

               CASE ('knh4ox')             ; scp(zon)%knh4ox             = extract_double(values(ccol))
               CASE ('ktsno3')             ; scp(zon)%ktsno3             = extract_double(values(ccol))
               CASE ('ktsox')              ; scp(zon)%ktsox              = extract_double(values(ccol))
               CASE ('kmnox')              ; scp(zon)%kmnox              = extract_double(values(ccol))
               CASE ('kmnno3')             ; scp(zon)%kmnno3             = extract_double(values(ccol))
               CASE ('kfeox')              ; scp(zon)%kfeox              = extract_double(values(ccol))
               CASE ('kfesox')             ; scp(zon)%kfesox             = extract_double(values(ccol))
               CASE ('kfes2ox')            ; scp(zon)%kfes2ox            = extract_double(values(ccol))
 !##           CASE ('kfe1ox')             ; scp(zon)%kfe1ox             = extract_double(values(ccol))
 !##           CASE ('kfe2ox')             ; scp(zon)%kfe2ox             = extract_double(values(ccol))
               CASE ('kfeAge')             ; scp(zon)%kfeAge             = extract_double(values(ccol))
               CASE ('kmnAge')             ; scp(zon)%kmnAge             = extract_double(values(ccol))
               CASE ('kfeno3')             ; scp(zon)%kfeno3             = extract_double(values(ccol))
 !##           CASE ('kfe2no3')            ; scp(zon)%kfe2no3            = extract_double(values(ccol))
               CASE ('kfesmn')             ; scp(zon)%kfesmn             = extract_double(values(ccol))
               CASE ('kfesfe')             ; scp(zon)%kfesfe             = extract_double(values(ccol))
               CASE ('kfesppt')            ; scp(zon)%kfesppt            = extract_double(values(ccol))
               CASE ('ktsfe')              ; scp(zon)%ktsfe              = extract_double(values(ccol))
               CASE ('kpyrite')            ; scp(zon)%kpyrite            = extract_double(values(ccol))
               CASE ('kmnfe')              ; scp(zon)%kmnfe              = extract_double(values(ccol))
               CASE ('ktsmn')              ; scp(zon)%ktsmn              = extract_double(values(ccol))
               CASE ('kch4ox')             ; scp(zon)%kch4ox             = extract_double(values(ccol))
               CASE ('kch4so4')            ; scp(zon)%kch4so4            = extract_double(values(ccol))
               CASE ('ko2')                ; scp(zon)%ko2                = extract_double(values(ccol))
               CASE ('kpo2')               ; scp(zon)%kpo2               = extract_double(values(ccol))
               CASE ('kno3')               ; scp(zon)%kno3               = extract_double(values(ccol))
               CASE ('kpno3')              ; scp(zon)%kpno3              = extract_double(values(ccol))
               CASE ('kmn')                ; scp(zon)%kmn                = extract_double(values(ccol))
               CASE ('kpmn')               ; scp(zon)%kpmn               = extract_double(values(ccol))
               CASE ('kfe')                ; scp(zon)%kfe                = extract_double(values(ccol))
               CASE ('kpfe')               ; scp(zon)%kpfe               = extract_double(values(ccol))
               CASE ('kso4')               ; scp(zon)%kso4               = extract_double(values(ccol))
               CASE ('kpso4')              ; scp(zon)%kpso4              = extract_double(values(ccol))
               CASE ('lo2')                ; scp(zon)%lo2                = extract_double(values(ccol))
               CASE ('lpo2')               ; scp(zon)%lpo2               = extract_double(values(ccol))
               CASE ('lno3')               ; scp(zon)%lno3               = extract_double(values(ccol))
               CASE ('lpno3')              ; scp(zon)%lpno3              = extract_double(values(ccol))
               CASE ('lmn')                ; scp(zon)%lmn                = extract_double(values(ccol))
               CASE ('lpmn')               ; scp(zon)%lpmn               = extract_double(values(ccol))
               CASE ('lfe')                ; scp(zon)%lfe                = extract_double(values(ccol))
               CASE ('lpfe')               ; scp(zon)%lpfe               = extract_double(values(ccol))
               CASE ('lso4')               ; scp(zon)%lso4               = extract_double(values(ccol))
               CASE ('lpso4')              ; scp(zon)%lpso4              = extract_double(values(ccol))
               CASE ('kanh4')              ; scp(zon)%kanh4              = extract_double(values(ccol))
               CASE ('kapo4')              ; scp(zon)%kapo4              = extract_double(values(ccol))

               CASE ('kSidppt')            ; scp(zon)%kSidppt            = extract_double(values(ccol))
               CASE ('kCalppt')            ; scp(zon)%kCalppt            = extract_double(values(ccol))
               CASE ('kRodppt')            ; scp(zon)%kRodppt            = extract_double(values(ccol))
               CASE ('kMnO2Appt')          ; scp(zon)%kMnO2Appt          = extract_double(values(ccol))
               CASE ('kMnO2Bppt')          ; scp(zon)%kMnO2Bppt          = extract_double(values(ccol))
               CASE ('kFeOHAppt')          ; scp(zon)%kFeOHAppt          = extract_double(values(ccol))
               CASE ('kFeOHBppt')          ; scp(zon)%kFeOHBppt          = extract_double(values(ccol))

               CASE ('Xname')              ; scp(zon)%Xname              = extract_string(values(ccol))
               CASE ('Xmk')                ; scp(zon)%Xmk                = extract_double(values(ccol))
               CASE ('Xfl')                ; scp(zon)%Xfl                = extract_double(values(ccol))
               CASE ('Xfm')                ; scp(zon)%Xfm                = extract_double(values(ccol))
               CASE ('kXSppt')             ; scp(zon)%kXSppt             = extract_double(values(ccol))
               CASE ('rxn_mode')           ; scp(zon)%rxn_mode           = extract_integer(values(ccol))
               CASE ('PO4AdsorptionModel') ; scp(zon)%PO4AdsorptionModel = extract_integer(values(ccol))
               CASE ('KPO4p')              ; scp(zon)%KPO4p              = extract_double(values(ccol))
               CASE ('Kadsratio')          ; scp(zon)%Kadsratio          = extract_double(values(ccol))
               CASE ('Qmax')               ; scp(zon)%Qmax               = extract_double(values(ccol))
               CASE ('ads_use_pH')         ; scp(zon)%ads_use_pH         = extract_logical(values(ccol))
               CASE DEFAULT ; print *, 'Unknown column ', TRIM(csvnames(ccol))
            END SELECT
         ENDDO
      ENDDO
   ENDIF

   meh = aed_csv_close(unit)
   !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CONTAINS

   !############################################################################
   CHARACTER FUNCTION tolower(c)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      CHARACTER, INTENT(in) :: c
   !LOCALS
      INTEGER :: ic
   !BEGIN
   !----------------------------------------------------------------------------
      ic = ichar(c)
      if (ic >= 65 .and. ic < 90) ic = (ic+32)
      tolower = char(ic)
   END FUNCTION tolower
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   FUNCTION same_str_icase(a, b) RESULT(res)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      CHARACTER(len=*), INTENT(in) :: a,b
   !LOCALS
      INTEGER :: len, i
      LOGICAL :: res
   !
   !BEGIN
   !----------------------------------------------------------------------------
      res = .FALSE.
      len = LEN_TRIM(a)
      IF ( len /= LEN_TRIM(b) ) RETURN
      DO i=1, len
         if (tolower(a(i:i)) /= tolower(b(i:i)) ) RETURN
      ENDDO
      res = .TRUE.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   END FUNCTION same_str_icase
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE read_sed_candi_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE default_sed_parm(sed_parm)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed2_sed_candi_param_t),DIMENSION(:),INTENT(out) :: sed_parm
!
!-------------------------------------------------------------------------------
!BEGIN
      !-- Sediment Physical Properties --!
   sed_parm%db0                 = 20.00
   sed_parm%imix                = 1
   sed_parm%xs                  = 20.0E+0
   sed_parm%x1                  = 15.0e-0
   sed_parm%x2                  = 35.0e-0
   sed_parm%irrg                = 1
   sed_parm%alpha0              = 30.00
   sed_parm%xirrig              = 7
   sed_parm%ventflow            = -15.0E-05 !
   sed_parm%w00                 = 0.15
   sed_parm%p0                  = 0.9
   sed_parm%p00                 = 0.6
   sed_parm%density             = 2.5
   sed_parm%bp                  = 0.6
   sed_parm%torteq              = 2
   sed_parm%an                  = 2.14
   sed_parm%aa                  = 3.79
   sed_parm%ab                  = 2.02
   sed_parm%xl                  = 40.0
   sed_parm%maxnpts             = 100
   sed_parm%job                 = 1

      !-- Biogeochemical Configuration --!
   sed_parm%numOM               = 3
   sed_parm%simDOM              = .true.
   sed_parm%OMapproach          = 1
   sed_parm%OMModel             = 1
   sed_parm%FTemswitch          = 2
   sed_parm%FTswitch            = 1
   sed_parm%FBIOswitch          = 1
   sed_parm%FINswitch           = 2
   sed_parm%FInO2Onlyswitch     = 2
   sed_parm%FOMswitch           = 1
   sed_parm%Bsolidswitch        = 1

   sed_parm%VCW                 = .false.
   sed_parm%simMnFe             = .true.
   sed_parm%simFeS              = .true.
   sed_parm%simX                = .false.
   sed_parm%simCaCO3            = .false.
   sed_parm%simFeCO3            = .false.
   sed_parm%simMnCO3            = .false.
   sed_parm%simNPAds            = .false.
   sed_parm%rxn_mode            = 2
   sed_parm%ads_use_pH          = .true.

      !-- Reaction rate constants --!
      !-- Primary reactions --!
      !-- OMModel 1
   sed_parm%poml2dic            = 0.5
   sed_parm%pomr2dic            = 0.0003
   sed_parm%pomspecial2dic      = 2.00
   sed_parm%R0                  = 135.0
   sed_parm%VCWBeta             = 0.208
      !-- OMModel 2
   sed_parm%docl2dic            = 0.0
   sed_parm%donl2din            = 0.0
   sed_parm%dopl2dip            = 0.0
   sed_parm%pocl2docl           = 0.0
   sed_parm%ponl2donl           = 0.0
   sed_parm%popl2dopl           = 0.0
   sed_parm%docr2docl           = 0.0
   sed_parm%donr2donl           = 0.0
   sed_parm%dopr2dopl           = 0.0
   sed_parm%pocr2docr           = 0.0
   sed_parm%ponr2donr           = 0.0
   sed_parm%popr2dopr           = 0.0
   sed_parm%pocvr2docr          = 0.0
   sed_parm%ponvr2donr          = 0.0
   sed_parm%popvr2dopr          = 0.0
      !-- OMModel 3
   sed_parm%kgrowthFer          = 5.0E1
   sed_parm%kgrowthAer          = 5.0E1
   sed_parm%kgrowthDen          = 5.0E1
   sed_parm%kgrowthMan          = 5.0E1
   sed_parm%kgrowthIro          = 5.0E1
   sed_parm%kgrowthSul          = 5.0E1
   sed_parm%kgrowthMet          = 5.0E1

   sed_parm%kdeathFer           = 5.0E-1
   sed_parm%kdeathAer           = 5.0E-1
   sed_parm%kdeathDen           = 5.0E-1
   sed_parm%kdeathMan           = 5.0E-1
   sed_parm%kdeathIro           = 5.0E-1
   sed_parm%kdeathSul           = 5.0E-1
   sed_parm%kdeathMet           = 5.0E-1

   sed_parm%kHyd1               = 0.95
   sed_parm%kHyd2               = 1.
   sed_parm%kHyd3               = 0.1
   sed_parm%kHyd4               = 0.01
   sed_parm%kHydN               = 0.

   sed_parm%domr2dic            = 0.01
   sed_parm%domr2pomr           = 0.001
   sed_parm%poml2doml           = 1.0

   sed_parm%BMax                = 5.50E-00

   sed_parm%BMin                = 0.050E-00
   sed_parm%fuse                = 0.1
   sed_parm%Cellweight          = 113.
   sed_parm%Tiny                = 1.0000000E-08
   sed_parm%Temporary_proton    = 123456789

      !-- Stoichiometric constants for OMModel 1 --!
      !--C--!  For the VCW or Canavan stoichiometry, this is like "x"
   sed_parm%xlab                = 140
   sed_parm%xref                = 140
   sed_parm%xspecial            = 140!9.09
      !--N--! For the VCW or Canavan stoichiometry, this is like "y"
   sed_parm%ylab                = 10.04545
   sed_parm%yref                = 10.04545
   sed_parm%yspecial            = 10 !0.76
      !--P--!   For the VCW or Canavan stoichiometry, this is like "z"
   sed_parm%zlab                = 0.5
   sed_parm%zref                = 0.5
   sed_parm%zspecial            = 0.5
      !--XMetal--!   For the VCW or Canavan stoichiometry, this is like "z"
   sed_parm%XMetal_lab          = 0.01
   sed_parm%XMetal_ref          = 0.01
   sed_parm%XMetal_special      = 1

      !-- Stoichiometric constants for OMModel 3 --!
   sed_parm%xPOM1               = 108.
   sed_parm%yPOM1               = 16.
   sed_parm%zPOM1               = 1.
   sed_parm%XMetal_POM1         = 1
   sed_parm%xPOM2               = 108.
   sed_parm%yPOM2               = 16.
   sed_parm%zPOM2               = 1.
   sed_parm%XMetal_POM2         = 1
   sed_parm%xPOM3               = 108.
   sed_parm%yPOM3               = 16.
   sed_parm%zPOM3               = 1.
   sed_parm%XMetal_POM3         = 1
   sed_parm%xPOM4               = 108.
   sed_parm%yPOM4               = 16.
   sed_parm%zPOM4               = 1.
   sed_parm%XMetal_POM4         = 1
   sed_parm%xDHyd               = 6.
   sed_parm%yDHyd               = 0.
   sed_parm%zDHyd               = 0.
   sed_parm%XMetal_DHyd         = 1
   sed_parm%xOAc                = 2.
   sed_parm%yOAc                = 0.
   sed_parm%zOAc                = 0.
   sed_parm%xH2                 = 1.
   sed_parm%yH2                 = 0.
   sed_parm%zH2                 = 0.

      !-- Secondary reactions --!
   sed_parm%kNH4OX              = 5.0e3
   sed_parm%kTSNO3              = 0.0
   sed_parm%kTSOX               = 1.6e2
   sed_parm%kMnOX               = 0.0
   sed_parm%kMnadsOx            = 5.0e3
   sed_parm%kMnNO3              = 0.0
   sed_parm%kFeOX               = 1.4e5
   sed_parm%kFeadsOX            = 5.0e4
   sed_parm%kFeSOX              = 2.2E+2
   sed_parm%kFeS2OX             = 1.0
   sed_parm%kMnAge              = 0.0
   sed_parm%kFeAge              = 0.00
   sed_parm%kFeNO3              = 0.0
   sed_parm%kFeSMn              = 0.0
   sed_parm%kFeSFe              = 0.0
   sed_parm%kFeSppt             = 1.5e-1
   sed_parm%kFeSdis             = 1.0e-3
   sed_parm%kTSFe               = 8.0e0
   sed_parm%kPyrite             = 2.3e-4
   sed_parm%kMnFe               = 3.0e3
   sed_parm%kTSMn               = 2.0e1
   sed_parm%kCH4OX              = 1.0e7
   sed_parm%kCH4SO4             = 1.0e+1
   sed_parm%kSidppt             = 4.5e+2
   sed_parm%kSiddis             = 2.5e-1
   sed_parm%kRodppt             = 1e+2
   sed_parm%kRoddis             = 2.5e-1
   sed_parm%kCalppt             = 0.0
   sed_parm%kMnO2Appt           = 0.0
   sed_parm%kMnO2Bppt           = 0.0e-2
   sed_parm%kFeOHAppt           = 0.0
   sed_parm%kFeOHBppt           = 0.0e-2

      !------------------------------------------------------------------------
      !#########!           FT setup  v        ################################
      !------------------------------------------------------------------------
   sed_parm%FTR                 = 8.31e-3
   sed_parm%FTT                 = 298.0
   sed_parm%deltaGATP           = 45.0
   sed_parm%e                   = 2.71828182845904523536028747135266249775724709369995
   sed_parm%F                   = 96.48534
   sed_parm%n                   = 1.
   sed_parm%dPsi                = 0.120

      !-- FT energies kJ mol substrate^-1
   sed_parm%dG0FerDHyd          = -30.
   sed_parm%dG0AerDHy           = -2883.
   sed_parm%dG0AerOAc           = -847.0

   sed_parm%dG0DenDHy           = -2774.
   sed_parm%dG0DenOAc           = -813.0
   sed_parm%dG0DenH2            = -226.0
   sed_parm%dG0ManOAc           = -625.0

   sed_parm%dG0IroOAc           = -736.6
   sed_parm%dG0IroH2            = -230.7
   sed_parm%dG0SulOAc           = -64.7
   sed_parm%dG0SulH2            = -38.8
   sed_parm%dG0MetOAc           = -31.7
   sed_parm%dG0MetH2            = -31.3

   sed_parm%YDHyAer             = 68.41
   sed_parm%YDHyFer             = 11.95
   sed_parm%YDenDHy             = 67.10
   sed_parm%YAerOAc             = 17.25
   sed_parm%YDenOAc             = 16.82
   sed_parm%YDenH2              = 04.49
   sed_parm%YManOAc             = 14.14
   sed_parm%YIroOAc             = 18.70
   sed_parm%YIroH2              = 04.54
   sed_parm%YSulOAc             = 02.03
   sed_parm%YSulH2              = 01.15
   sed_parm%YMetOAc             = 01.02
   sed_parm%YMetH2              = 00.94

      !------------------------------------------------------------------------
      !##########!             FT setup ^           ###########################
      !------------------------------------------------------------------------
      !-- FTEA and FIN constants
      !-- Approach 1
   sed_parm%kO2                 = 20.0E-3
   sed_parm%kpO2                = 20.0E-3
   sed_parm%kNO3                = 3.0E-3
   sed_parm%kpNO3               = 3.0E-3
   sed_parm%kMnO2               = 16.0E-0
   sed_parm%kpMnO2              = 16.0E-0
   sed_parm%kFeOH               = 100.0E-0
   sed_parm%kpFeOH              = 100.0E-0
   sed_parm%kSO4                = 1600.0E-3
   sed_parm%kpSO4               = 1600.0E-3

   sed_parm%kpO2NO3             = 3E-3
   sed_parm%kpO2MnO2            = 200E-0
   sed_parm%kpO2FeOH            = 200E-0
   sed_parm%kpO2SO4             = 2E-3
   sed_parm%kpO2CH4             = 2E-3

      !-- Approach 2
   sed_parm%lO2                 = 20e-3
   sed_parm%lpO2                = 20e-3
   sed_parm%lNO3                = 5e-3
   sed_parm%lpNO3               = 5e-3

   sed_parm%lMnO2               = 16e-0
   sed_parm%lpMnO2              = 16e-0
   sed_parm%lFeOH               = 100e-0
   sed_parm%lpFeOH              = 100e-0

   sed_parm%lSO4                = 1600e-3
   sed_parm%lpSO4               = 1600e-3

   sed_parm%lpo2no3             = 20.0e-3
   sed_parm%lpo2mno2            = 20.0e-3
   sed_parm%lpo2feoh            = 20.0e-3
   sed_parm%lpo2so4             = 20.0e-3
   sed_parm%lpo2ch4             = 20.0e-3

      !-- FOM
   sed_parm%KDHyd               = 1.0E-00
   sed_parm%KOAc                = 1.0E-2
   sed_parm%KH2                 = 1.0E-2
   sed_parm%fracOAc             = 0.666687
   sed_parm%fracH2              = 0.333333

      !-- Geochemcial setup
   sed_parm%kNH4Ads             = 1.4
   sed_parm%VCWSb               = 30
   sed_parm%GammaS              = 0.11
   sed_parm%KSMnads             = 3.5
   sed_parm%KSFeads             = 3.7
   sed_parm%kPO4Ads             = 3.0
   sed_parm%PO4AdsorptionModel  = 1
   sed_parm%KPO4p               = 1.05
   sed_parm%Kadsratio           = 1.05
   sed_parm%Qmax                = 1.05
   sed_parm%KPSid               = -8.4
   sed_parm%KPRod               = -8.5
   sed_parm%KPFeS               = 2.2
   sed_parm%Xname               = "Zn"
   sed_parm%Xmk                 = 0.06
   sed_parm%Xfl                 = 0.065
   sed_parm%Xfm                 = 0.0001
   sed_parm%kXSppt              = 0.0

   sed_parm%kmn                 = 0.0
   sed_parm%kpmn                = 0.0
   sed_parm%kfe                 = 0.0
   sed_parm%kpfe                = 0.0
   sed_parm%lmn                 = 0.0
   sed_parm%lpmn                = 0.0
   sed_parm%lfe                 = 0.0
   sed_parm%lpfe                = 0.0
   sed_parm%kanh4               = 0.0
   sed_parm%kapo4               = 0.0
END SUBROUTINE default_sed_parm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_read_candi
