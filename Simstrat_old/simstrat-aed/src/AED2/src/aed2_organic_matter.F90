!###############################################################################
!#                                                                             #
!# aed2_organic_matter.F90                                                     #
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
!# Created 6 June 2011                                                         #
!# Updated 11 May 2015 : Add refractory groups, based on Swan estuary model    #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_organic_matter
!-------------------------------------------------------------------------------
! aed2_organic_matter --- organic matter biogeochemical model
!
! The Organic Matter (OM) module contains equations for mineralisation
! of particulate and dissolved organic matter pools.
!
! Users can configure the number of pools
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE  ! By default make everything private
!
   PUBLIC aed2_organic_matter_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_organic_matter_data_t

      !# Variable identifiers
      INTEGER  :: id_poc,id_doc           ! particulate & dissolved organic carbon
      INTEGER  :: id_pon,id_don           ! particulate & dissolved organic nitrogen
      INTEGER  :: id_pop,id_dop           ! particulate & dissolved organic phosphorus
      INTEGER  :: id_docr,id_donr,id_dopr ! dissolved refractory OM pools
      INTEGER  :: id_cpom                 ! coarse particulate OM (bulk) pool
      INTEGER  :: id_photolysis           ! photolysis

      INTEGER  :: id_oxy,id_amm,id_frp,id_dic
      INTEGER  :: id_Fsed_pon, id_Fsed_don ! sed. rate organic nitrogen
      INTEGER  :: id_Fsed_pop, id_Fsed_dop ! sed. rate organic phosphorus
      INTEGER  :: id_Fsed_poc, id_Fsed_doc ! sed. rate organic carbon
      INTEGER  :: id_Psed_poc, id_Psed_pon, id_Psed_pop ! sedimentation rates
      INTEGER  :: id_temp, id_salt, id_vis, id_uva, id_uvb
      INTEGER  :: id_pon_miner, id_don_miner
      INTEGER  :: id_pop_miner, id_dop_miner
      INTEGER  :: id_poc_miner, id_doc_miner
      INTEGER  :: id_sed_pon, id_sed_don
      INTEGER  :: id_sed_pop, id_sed_dop
      INTEGER  :: id_sed_poc, id_sed_doc
      INTEGER  :: id_bod, id_cdom

      !# Model parameters
      AED_REAL :: w_pon,Rpon_miner,Rdon_miner,Fsed_pon,Fsed_don,           &
                          Kpon_miner, Kdon_miner, Ksed_don,                &
                          theta_pon_miner, theta_don_miner, theta_sed_don
      AED_REAL :: w_pop,Rpop_miner,Rdop_miner,Fsed_pop,Fsed_dop,           &
                          Kpop_miner, Kdop_miner, Ksed_dop,                &
                          theta_pop_miner, theta_dop_miner, theta_sed_dop
      AED_REAL :: w_poc,Rpoc_miner,Rdoc_miner,Fsed_poc,Fsed_doc,           &
                          Kpoc_miner, Kdoc_miner, Ksed_doc,                &
                          theta_poc_miner, theta_doc_miner, theta_sed_doc, &
                          KeDOM, KePOM
      AED_REAL ::    docr_initial,donr_initial,dopr_initial,cpom_initial,  &
                          Rdocr_miner,Rdonr_miner,Rdopr_miner,Rcpom_bdown, &
                          w_cpom,X_cpom_n,X_cpom_p,KeDOMR,KeCPOM

      AED_REAL :: photo_fmin

      LOGICAL  :: use_oxy, use_amm, use_frp, use_dic, use_Fsed_model, use_Psed_model
      LOGICAL  :: simRPools, extra_diag, simphotolysis

     CONTAINS
         PROCEDURE :: define            => aed2_define_organic_matter
         PROCEDURE :: calculate         => aed2_calculate_organic_matter
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_organic_matter
!        PROCEDURE :: mobility          => aed2_mobility_organic_matter
         PROCEDURE :: light_extinction  => aed2_light_extinction_organic_matter
!        PROCEDURE :: delete            => aed2_delete_organic_matter

   END TYPE

   AED_REAL :: c

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_organic_matter(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                              :: namlst
   CLASS (aed2_organic_matter_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status

   AED_REAL                  :: poc_initial = 4.5
   AED_REAL                  :: doc_initial = 4.5
   AED_REAL                  :: poc_min=0.
   AED_REAL                  :: poc_max=nan_
   AED_REAL                  :: doc_min=0.
   AED_REAL                  :: doc_max=nan_
   AED_REAL                  :: w_poc       = 0.0
   AED_REAL                  :: Rpoc_miner  = 0.01
   AED_REAL                  :: Rdoc_miner  = 0.01
   AED_REAL                  :: Fsed_poc    =  0.0
   AED_REAL                  :: Fsed_doc    = 30.0
   AED_REAL                  :: Kpoc_miner  = 30.0
   AED_REAL                  :: Kdoc_miner  = 30.0
   AED_REAL                  :: Ksed_doc    = 4.5
   AED_REAL                  :: theta_poc_miner = 1.0
   AED_REAL                  :: theta_doc_miner = 1.0
   AED_REAL                  :: theta_sed_doc   = 1.0
   AED_REAL                  :: KeDOM = 0.0001
   AED_REAL                  :: KePOM = 0.0001
   CHARACTER(len=64)         :: doc_miner_product_variable=''
   CHARACTER(len=64)         :: doc_miner_reactant_variable=''
   CHARACTER(len=64)         :: Fsed_poc_variable=''
   CHARACTER(len=64)         :: Fsed_doc_variable=''

   AED_REAL                  :: pon_initial = 4.5
   AED_REAL                  :: don_initial = 4.5
   AED_REAL                  :: pon_min=0.
   AED_REAL                  :: pon_max=nan_
   AED_REAL                  :: don_min=0.
   AED_REAL                  :: don_max=nan_
   AED_REAL                  :: w_pon       = 0.0
   AED_REAL                  :: Rpon_miner  = 0.01
   AED_REAL                  :: Rdon_miner  = 0.01
   AED_REAL                  :: Fsed_pon    =  0.0
   AED_REAL                  :: Fsed_don    = 30.0
   AED_REAL                  :: Kpon_miner  = 30.0
   AED_REAL                  :: Kdon_miner  = 30.0
   AED_REAL                  :: Ksed_don    = 4.5
   AED_REAL                  :: theta_pon_miner = 1.0
   AED_REAL                  :: theta_don_miner = 1.0
   AED_REAL                  :: theta_sed_don   = 1.0
   CHARACTER(len=64)         :: don_miner_product_variable=''
   CHARACTER(len=64)         :: Fsed_pon_variable=''
   CHARACTER(len=64)         :: Fsed_don_variable=''

   AED_REAL                  :: pop_initial = 4.5
   AED_REAL                  :: dop_initial = 4.5
   AED_REAL                  :: pop_min=0.
   AED_REAL                  :: pop_max=nan_
   AED_REAL                  :: dop_min=0.
   AED_REAL                  :: dop_max=nan_
   AED_REAL                  :: w_pop       = 0.0
   AED_REAL                  :: Rpop_miner  = 0.01
   AED_REAL                  :: Rdop_miner  = 0.01
   AED_REAL                  :: Fsed_pop    =  0.0
   AED_REAL                  :: Fsed_dop    = 30.0
   AED_REAL                  :: Kpop_miner  = 30.0
   AED_REAL                  :: Kdop_miner  = 30.0
   AED_REAL                  :: Ksed_dop    = 4.5
   AED_REAL                  :: theta_pop_miner = 1.0
   AED_REAL                  :: theta_dop_miner = 1.0
   AED_REAL                  :: theta_sed_dop   = 1.0
   CHARACTER(len=64)         :: dop_miner_product_variable=''
   CHARACTER(len=64)         :: Fsed_pop_variable=''
   CHARACTER(len=64)         :: Fsed_dop_variable=''

   CHARACTER(len=64)         :: Psed_poc_variable=''
   CHARACTER(len=64)         :: Psed_pon_variable=''
   CHARACTER(len=64)         :: Psed_pop_variable=''

   LOGICAL                   :: simRPools
   AED_REAL                  :: docr_initial = 4.5
   AED_REAL                  :: donr_initial = 4.5
   AED_REAL                  :: dopr_initial = 4.5
   AED_REAL                  :: cpom_initial = 4.5
   AED_REAL                  :: Rdocr_miner  = 0.01
   AED_REAL                  :: Rdonr_miner  = 0.01
   AED_REAL                  :: Rdopr_miner  = 0.01
   AED_REAL                  :: Rcpom_bdown  = 0.01
   AED_REAL                  :: w_cpom       = 0.0
   AED_REAL                  :: X_cpom_n     = 0.0
   AED_REAL                  :: X_cpom_p     = 0.0
   AED_REAL                  :: KeDOMR = 0.0001
   AED_REAL                  :: KeCPOM = 0.0001
   CHARACTER(len=64)         :: docr_miner_product_variable='OGM_doc'
   CHARACTER(len=64)         :: donr_miner_product_variable='OGM_don'
   CHARACTER(len=64)         :: dopr_miner_product_variable='OGM_dop'

   LOGICAL                   :: extra_diag = .FALSE.
   LOGICAL                   :: simphotolysis = .FALSE.

   AED_REAL                  :: photo_fmin = 0.9, photo_c = 7.52


   NAMELIST /aed2_organic_matter/ &
             poc_min, poc_max, doc_min, doc_max, &
             pon_min, pon_max, don_min, don_max, &
             pop_min, pop_max, dop_min, dop_max, &
             poc_initial, doc_initial, w_poc, Rpoc_miner, Rdoc_miner, Fsed_poc, Fsed_doc, &
             Kpoc_miner, Kdoc_miner, Ksed_doc,                &
             theta_poc_miner, theta_doc_miner, theta_sed_doc, KeDOM, KePOM, &
             doc_miner_reactant_variable, doc_miner_product_variable, &
             Fsed_poc_variable, Fsed_doc_variable, &
             pon_initial, don_initial, w_pon, Rpon_miner, Rdon_miner, Fsed_pon, Fsed_don, &
             Kpon_miner, Kdon_miner, Ksed_don,                &
             theta_pon_miner, theta_don_miner, theta_sed_don, &
             don_miner_product_variable,                      &
             Fsed_pon_variable, Fsed_don_variable,            &
             pop_initial, dop_initial, w_pop, Rpop_miner, Rdop_miner, Fsed_pop, Fsed_dop, &
             Kpop_miner, Kdop_miner, Ksed_dop,                &
             theta_pop_miner, theta_dop_miner, theta_sed_dop, &
             dop_miner_product_variable,                      &
             Fsed_pop_variable, Fsed_dop_variable,            &
             Psed_poc_variable, Psed_pon_variable, Psed_pop_variable, &
             simRPools, docr_initial, donr_initial, dopr_initial, cpom_initial, &
             Rdocr_miner, Rdonr_miner, Rdopr_miner, Rcpom_bdown, &
             w_cpom, X_cpom_n, X_cpom_p, KeDOMR, KeCPOM, &
           ! docr_miner_product_variable, donr_miner_product_variable,dopr_miner_product_variable, &
             extra_diag, simphotolysis, photo_fmin, photo_c


!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_organic_matter,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_organic_matter'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   ! Organic Carbon
   data%w_poc           = w_poc/secs_per_day
   data%Rpoc_miner      = Rpoc_miner/secs_per_day
   data%Rdoc_miner      = Rdoc_miner/secs_per_day
 ! data%Fsed_poc        = Fsed_poc/secs_per_day
   data%Fsed_poc        = 0.0
   data%Fsed_doc        = Fsed_doc/secs_per_day
   data%Kpoc_miner      = Kpoc_miner
   data%Kdoc_miner      = Kdoc_miner
   data%Ksed_doc        = Ksed_doc
   data%theta_poc_miner = theta_poc_miner
   data%theta_doc_miner = theta_doc_miner
   data%theta_sed_doc   = theta_sed_doc
   data%KeDOM           = KeDOM
   data%KePOM           = KePOM
   ! Organic Nitrogen
   data%w_pon           = w_pon/secs_per_day
   data%Rpon_miner      = Rpon_miner/secs_per_day
   data%Rdon_miner      = Rdon_miner/secs_per_day
 ! data%Fsed_pon        = Fsed_pon/secs_per_day
   data%Fsed_pon        = 0.0
   data%Fsed_don        = Fsed_don/secs_per_day
   data%Kpon_miner      = Kpon_miner
   data%Kdon_miner      = Kdon_miner
   data%Ksed_don        = Ksed_don
   data%theta_pon_miner = theta_pon_miner
   data%theta_don_miner = theta_don_miner
   data%theta_sed_don   = theta_sed_don
   ! Organic Phosphorus
   data%w_pop           = w_pop/secs_per_day
   data%Rpop_miner      = Rpop_miner/secs_per_day
   data%Rdop_miner      = Rdop_miner/secs_per_day
 ! data%Fsed_pop        = Fsed_pop/secs_per_day
   data%Fsed_pop        = 0.0
   data%Fsed_dop        = Fsed_dop/secs_per_day
   data%Kpop_miner      = Kpop_miner
   data%Kdop_miner      = Kdop_miner
   data%Ksed_dop        = Ksed_dop
   data%theta_pop_miner = theta_pop_miner
   data%theta_dop_miner = theta_dop_miner
   data%theta_sed_dop   = theta_sed_dop
   ! Refractory Organic Matter (Optional)
   data%simRPools       = simRPools
   data%docr_initial    = docr_initial
   data%donr_initial    = donr_initial
   data%dopr_initial    = dopr_initial
   data%cpom_initial    = cpom_initial
   data%Rdocr_miner     = Rdocr_miner/secs_per_day
   data%Rdonr_miner     = Rdonr_miner/secs_per_day
   data%Rdopr_miner     = Rdopr_miner/secs_per_day
   data%Rcpom_bdown     = Rcpom_bdown/secs_per_day
   data%w_cpom          = w_cpom
   data%X_cpom_n        = X_cpom_n
   data%X_cpom_p        = X_cpom_p
   data%KeDOMR          = KeDOMR
   data%KeCPOM          = KeCPOM

   data%extra_diag      = extra_diag
   data%simphotolysis  = simphotolysis
   data%photo_fmin      = photo_fmin
   c = photo_c

   ! Register state variables
   data%id_doc = aed2_define_variable('doc','mmol/m**3','dissolved organic carbon',       &
                                    doc_initial,minimum=doc_min,maximum=doc_max)
   data%id_poc = aed2_define_variable('poc','mmol/m**3','particulate organic carbon',     &
                                    poc_initial,minimum=poc_min,maximum=poc_max,mobility=data%w_poc)
   data%id_don = aed2_define_variable('don','mmol/m**3','dissolved organic nitrogen',     &
                                    don_initial,minimum=don_min,maximum=don_max)
   data%id_pon = aed2_define_variable('pon','mmol/m**3','particulate organic nitrogen',   &
                                    pon_initial,minimum=pon_min,maximum=pon_max,mobility=data%w_pon)
   data%id_dop = aed2_define_variable('dop','mmol/m**3','dissolved organic phosphorus',   &
                                    dop_initial,minimum=dop_min,maximum=dop_max)
   data%id_pop = aed2_define_variable('pop','mmol/m**3','particulate organic phosphorus', &
                                    pop_initial,minimum=pop_min,maximum=pop_max,mobility=data%w_pop)

   IF (simRPools) THEN
     data%id_docr = aed2_define_variable('docr','mmol/m**3','dissolved organic carbon R', &
                                    docr_initial,minimum=zero_)
     data%id_donr = aed2_define_variable('donr','mmol/m**3','dissolved organic nitrogen R', &
                                    donr_initial,minimum=zero_)
     data%id_dopr = aed2_define_variable('dopr','mmol/m**3','dissolved organic phosphorus R', &
                                    dopr_initial,minimum=zero_)
     data%id_cpom = aed2_define_variable('cpom','mmol/m**3','coarse particulate matter',  &
                                    cpom_initial,minimum=zero_,mobility=data%w_cpom)
     data%id_vis= aed2_locate_global('par')
     data%id_uva= aed2_locate_global('uva')
     data%id_uvb= aed2_locate_global('uvb')
   ENDIF



   ! Register external state variable dependencies (carbon)
   data%use_oxy = doc_miner_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) &
     data%id_oxy = aed2_locate_variable(doc_miner_reactant_variable)

   data%use_dic = doc_miner_product_variable .NE. '' !This means carbon module switched on
   IF (data%use_dic) &
     data%id_dic = aed2_locate_variable(doc_miner_product_variable)

   ! Register external state variable dependencies (nitrogen)
   data%use_amm = don_miner_product_variable .NE. '' !This means nitrogen module switched on
   IF (data%use_amm) &
     data%id_amm = aed2_locate_variable(don_miner_product_variable)

   ! Register external state variable dependencies (phosphorous)
   data%use_frp = dop_miner_product_variable .NE. '' !This means phosphorus module switched on
   IF (data%use_frp) &
     data%id_frp = aed2_locate_variable(dop_miner_product_variable)

   data%id_Fsed_pon = -1 ; data%id_Fsed_pop = -1 ; data%id_Fsed_poc = -1
   data%id_Fsed_don = -1 ; data%id_Fsed_dop = -1 ; data%id_Fsed_doc = -1

   data%use_Fsed_model = Fsed_pon_variable .NE. ''
   IF (data%use_Fsed_model) THEN
     data%id_Fsed_pon = aed2_locate_global_sheet(Fsed_pon_variable)
     IF (Fsed_don_variable .NE. '') &
        data%id_Fsed_don = aed2_locate_global_sheet(Fsed_don_variable)
     IF (Fsed_pop_variable .NE. '') &
        data%id_Fsed_pop = aed2_locate_global_sheet(Fsed_pop_variable)
     IF (Fsed_dop_variable .NE. '') &
        data%id_Fsed_dop = aed2_locate_global_sheet(Fsed_dop_variable)
     IF (Fsed_poc_variable .NE. '') &
        data%id_Fsed_poc = aed2_locate_global_sheet(Fsed_poc_variable)
     IF (Fsed_doc_variable .NE. '') &
        data%id_Fsed_doc = aed2_locate_global_sheet(Fsed_doc_variable)
   ENDIF

   data%id_Psed_poc = -1 ; data%id_Psed_pon = -1 ; data%id_Psed_pop = -1

   data%use_Psed_model = Psed_poc_variable .NE. ''
   IF (data%use_Psed_model) THEN
     data%id_Psed_poc = aed2_locate_global_sheet(Psed_poc_variable)
     IF (Psed_pon_variable .NE. '') &
        data%id_Psed_pon = aed2_locate_global_sheet(Psed_pon_variable)
     IF (Psed_pop_variable .NE. '') &
        data%id_Psed_pop = aed2_locate_global_sheet(Psed_pop_variable)
   ENDIF


   !# Register diagnostic variables

   data%id_cdom = aed2_define_diag_variable('CDOM','/m',  'Chromophoric DOM (CDOM)')

   IF (extra_diag) THEN
     data%id_pon_miner = aed2_define_diag_variable('pon_miner','mmol/m**3/d',  'PON mineralisation')
     data%id_don_miner = aed2_define_diag_variable('don_miner','mmol/m**3/d',  'DON mineralisation')
     data%id_sed_pon = aed2_define_sheet_diag_variable('sed_pon','mmol/m**2/d',  'PON sediment flux')
     data%id_sed_don = aed2_define_sheet_diag_variable('sed_don','mmol/m**2/d',  'DON sediment flux')

     data%id_pop_miner = aed2_define_diag_variable('pop_miner','mmol/m**3/d',  'POP mineralisation')
     data%id_dop_miner = aed2_define_diag_variable('dop_miner','mmol/m**3/d',  'DOP mineralisation')
     data%id_sed_pop = aed2_define_sheet_diag_variable('sed_pop','mmol/m**2/d',  'POP sediment flux')
     data%id_sed_dop = aed2_define_sheet_diag_variable('sed_dop','mmol/m**2/d',  'DOP sediment flux')

     data%id_poc_miner = aed2_define_diag_variable('poc_miner','mmol/m**3/d',  'POC mineralisation')
     data%id_doc_miner = aed2_define_diag_variable('doc_miner','mmol/m**3/d',  'DOC mineralisation')
     data%id_sed_poc = aed2_define_sheet_diag_variable('sed_poc','mmol/m**2/d',  'POC sediment flux')
     data%id_sed_doc = aed2_define_sheet_diag_variable('sed_doc','mmol/m**2/d',  'DOC sediment flux')

     data%id_bod = aed2_define_diag_variable('BOD5','mg O2/L',  'Biochemical Oxygen Demand (BOD)')
     IF ( simphotolysis ) &
        data%id_photolysis = aed2_define_diag_variable('photolysis','mmol C/m3/d',  'photolysis rate of breakdown of DOC')
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')

END SUBROUTINE aed2_define_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION photo(Q, cdom, band)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: Q          ! [W/m2] incoming light intensity
   AED_REAL,INTENT(in) :: cdom       ! [/m] current cdom absorbance (based on doc)
   INTEGER, INTENT(in) :: band       ! [-] bandwidth identifier (VIS, UVA, UVB)
!LOCALS
!  AED_REAL,PARAMETER :: c = 7.52    ! now in NML as photo_c and module global c
   AED_REAL,PARAMETER :: d = 0.0122  ! ref
   AED_REAL,PARAMETER :: S = 0.0188  ! ref
   AED_REAL,PARAMETER :: x = 440.    ! [nm] wavelength at which cdom is measured
   AED_REAL,PARAMETER :: lamda(3) = (/440., 358., 298./) ! [nm] wavelengths
! vis = 390-700 nanometers
! uva = 315-400 nanometers
! uvb = 280-315 nanometers
!
   AED_REAL :: phi, alpha, QQ
!
!-------------------------------------------------------------------------------
!BEGIN
   ! convert incoming light to quantum units: [mol photons /m2 /hr]
   QQ = Q * 4.53 * 1e-6 * 3.6e3

   ! apparent quantum yield [mol produced /(mol absorbed photons) /nm-1]
   phi = c * 10.**(-(d*lamda(band)))

   ! absorption coefficient (/m)
   alpha = cdom*exp(S*(x-lamda(band)))

   ! return the rate at which DOC is broken down [mmol C /m3 /s]
   photo = phi * QQ * alpha * 1e3 / 3.6e3
!  print *,'ph',Q,phi,cdom,alpha,photo
END FUNCTION photo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_organic_matter(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_organic_matter model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: pon,don,amm,oxy,temp !State variables
   AED_REAL :: pon_mineralisation, don_mineralisation
   AED_REAL :: pop,dop,frp !State variables
   AED_REAL :: pop_mineralisation, dop_mineralisation
   AED_REAL :: poc,doc,dic !State variables
   AED_REAL :: poc_mineralisation, doc_mineralisation
   AED_REAL :: docr,donr,dopr,cpom,cdom
   AED_REAL :: docr_mineralisation, donr_mineralisation
   AED_REAL :: dopr_mineralisation, cpom_breakdown
   AED_REAL :: photolysis
   AED_REAL :: vis, uva, uvb, photo_fmin

   AED_REAL, PARAMETER :: Yoxy_doc_miner = 32./12. !ratio of oxygen to carbon utilised during doc mineralisation

!-----------------------------------------------------------------------
!BEGIN
   !CALL log_message('model aed2_organic_matter enter do loop successfully.')

   photolysis = zero_
   photo_fmin = data%photo_fmin

   ! Retrieve current (local) state variable values.
   pon = _STATE_VAR_(data%id_pon)! particulate organic nitrogen
   don = _STATE_VAR_(data%id_don)! dissolved organic nitrogen
   pop = _STATE_VAR_(data%id_pop)! particulate organic phosphorus
   dop = _STATE_VAR_(data%id_dop)! dissolved organic phosphorus
   poc = _STATE_VAR_(data%id_poc)! particulate organic carbon
   doc = _STATE_VAR_(data%id_doc)! dissolved organic carbon

   docr = zero_
   IF(data%simRPools) THEN
      docr = _STATE_VAR_(data%id_docr)
      donr = _STATE_VAR_(data%id_donr)
      dopr = _STATE_VAR_(data%id_dopr)
      cpom = _STATE_VAR_(data%id_cpom)
      vis = _STATE_VAR_(data%id_vis)
      uva = _STATE_VAR_(data%id_uva)
      uvb = _STATE_VAR_(data%id_uvb)
      cdom = _DIAG_VAR_(data%id_cdom)
   ENDIF


   IF (data%use_oxy) THEN ! & use_oxy
      oxy = _STATE_VAR_(data%id_oxy)! oxygen
   ELSE
      oxy = 300.0
   ENDIF
   IF (data%use_dic) THEN ! & use_dic
      dic = _STATE_VAR_(data%id_dic)! disolved inorganic carbon
   ELSE
      dic = 100.0
   ENDIF
   IF (data%use_amm) THEN ! & use_amm
      amm = _STATE_VAR_(data%id_amm)! ammonium
   ELSE
      amm = 0.0
   ENDIF
   IF (data%use_frp) THEN ! & use_frp
      frp = _STATE_VAR_(data%id_frp)! phosphate
   ELSE
      frp = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp) ! temperature

   ! Define some intermediate quantities units mmol N/m3/day
   pon_mineralisation = fpon_miner(data%use_oxy,data%Rpon_miner,data%Kpon_miner,data%theta_pon_miner,oxy,temp)
   don_mineralisation = fdon_miner(data%use_oxy,data%Rdon_miner,data%Kdon_miner,data%theta_don_miner,oxy,temp)
   pop_mineralisation = fpop_miner(data%use_oxy,data%Rpop_miner,data%Kpop_miner,data%theta_pop_miner,oxy,temp)
   dop_mineralisation = fdop_miner(data%use_oxy,data%Rdop_miner,data%Kdop_miner,data%theta_dop_miner,oxy,temp)
   poc_mineralisation = fpoc_miner(data%use_oxy,data%Rpoc_miner,data%Kpoc_miner,data%theta_poc_miner,oxy,temp)
   doc_mineralisation = fdoc_miner(data%use_oxy,data%Rdoc_miner,data%Kdoc_miner,data%theta_doc_miner,oxy,temp)
   IF(data%simRPools) THEN
      docr_mineralisation = fdoc_miner(data%use_oxy,data%Rdocr_miner,data%Kdoc_miner,data%theta_doc_miner,oxy,temp)
      donr_mineralisation = fdoc_miner(data%use_oxy,data%Rdonr_miner,data%Kdoc_miner,data%theta_doc_miner,oxy,temp)
      dopr_mineralisation = fdoc_miner(data%use_oxy,data%Rdopr_miner,data%Kdoc_miner,data%theta_doc_miner,oxy,temp)
      cpom_breakdown = fpoc_miner(data%use_oxy,data%Rcpom_bdown,data%Kpoc_miner,data%theta_poc_miner,oxy,temp)
      IF ( data%simphotolysis ) THEN
         photolysis = photo(vis,cdom,1) + photo(uva,cdom,2) + photo(uvb,cdom,3)
         !# Limit photolysis to 90% of doc pool within 1 hour
         IF(photolysis > 0.9*docr/3.6e3) photolysis = 0.9*docr/3.6e3
      ENDIF
   ENDIF
!print *,'aa',photolysis,donr,docr

   ! Set temporal derivatives
   _FLUX_VAR_(data%id_pon)  = _FLUX_VAR_(data%id_pon) + (-pon*pon_mineralisation)
   _FLUX_VAR_(data%id_don)  = _FLUX_VAR_(data%id_don) + (pon*pon_mineralisation-don*don_mineralisation)
   _FLUX_VAR_(data%id_pop)  = _FLUX_VAR_(data%id_pop) + (-pop*pop_mineralisation)
   _FLUX_VAR_(data%id_dop)  = _FLUX_VAR_(data%id_dop) + (pop*pop_mineralisation-dop*dop_mineralisation)
   _FLUX_VAR_(data%id_poc)  = _FLUX_VAR_(data%id_poc) + (-poc*poc_mineralisation)
   _FLUX_VAR_(data%id_doc)  = _FLUX_VAR_(data%id_doc) + (poc*poc_mineralisation-doc*doc_mineralisation)
   IF(data%simRPools) THEN
      docr = MAX(1.e-3,docr)
      _FLUX_VAR_(data%id_docr) = _FLUX_VAR_(data%id_docr) - docr*docr_mineralisation - photolysis
      _FLUX_VAR_(data%id_doc)  = _FLUX_VAR_(data%id_doc)  + docr*docr_mineralisation + photolysis*photo_fmin
      _FLUX_VAR_(data%id_donr) = _FLUX_VAR_(data%id_donr) - donr*donr_mineralisation - photolysis*(donr/docr)
      _FLUX_VAR_(data%id_don)  = _FLUX_VAR_(data%id_don)  + donr*donr_mineralisation + photolysis*photo_fmin*(donr/docr)
      _FLUX_VAR_(data%id_dopr) = _FLUX_VAR_(data%id_dopr) - dopr*dopr_mineralisation - photolysis*(dopr/docr)
      _FLUX_VAR_(data%id_dop)  = _FLUX_VAR_(data%id_dop)  + dopr*dopr_mineralisation + photolysis*photo_fmin*(dopr/docr)
      _FLUX_VAR_(data%id_cpom) = 0.00000 !_FLUX_VAR_(data%id_cpom) - cpom*cpom_breakdown
      _FLUX_VAR_(data%id_poc)  = _FLUX_VAR_(data%id_poc) + (cpom*cpom_breakdown)
      _FLUX_VAR_(data%id_pon)  = _FLUX_VAR_(data%id_pon) + (cpom*cpom_breakdown)*data%X_cpom_n
      _FLUX_VAR_(data%id_pop)  = _FLUX_VAR_(data%id_pop) + (cpom*cpom_breakdown)*data%X_cpom_p
   ENDIF


   ! If an externally maintained oxygen pool is present, take mineralisation from it
   IF (data%use_oxy) THEN
      !_FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-Yoxy_doc_miner*doc*doc_mineralisation)
      _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-doc*doc_mineralisation)
   ENDIF
   if (data%use_dic) THEN
      _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (doc*doc_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + photolysis*(1.-photo_fmin)
   ENDIF
   IF (data%use_amm) THEN
      _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (don*don_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + photolysis*(1.-photo_fmin)*(donr/docr)
    ! IF ( data%simRPools ) _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + photolysis*(1.-photo_fmin)*MAX((donr/docr),1e-3)
   ENDIF
   IF (data%use_frp) THEN
      _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + (dop*dop_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + photolysis*(1.-photo_fmin)*(dopr/docr)
   ENDIF

   !# Export diagnostic variables
   ! CDOM computed as a function of DOC amount, as empirically defined by
   ! Kostoglidis et al 2005 for Swan-Canning

   IF ( data%simRPools ) THEN
      _DIAG_VAR_(data%id_cdom) = 0.35*exp(0.1922*(doc+docr)*(12./1e3))
   ELSE
      _DIAG_VAR_(data%id_cdom) = 0.35*exp(0.1922*(doc)*(12./1e3))
   ENDIF

   IF (data%extra_diag) THEN
     _DIAG_VAR_(data%id_pon_miner) = pon_mineralisation*pon*secs_per_day
     _DIAG_VAR_(data%id_don_miner) = don_mineralisation*don*secs_per_day
     _DIAG_VAR_(data%id_pop_miner) = pop_mineralisation*pop*secs_per_day
     _DIAG_VAR_(data%id_dop_miner) = dop_mineralisation*dop*secs_per_day
     _DIAG_VAR_(data%id_poc_miner) = poc_mineralisation*poc*secs_per_day
     _DIAG_VAR_(data%id_doc_miner) = doc_mineralisation*doc*secs_per_day
     IF ( data%simphotolysis ) &
        _DIAG_VAR_(data%id_photolysis)= photolysis*secs_per_day
     ! BOD5 is computed as the amount of oxygen consumed over a 5-day period (mgO2/L)
     _DIAG_VAR_(data%id_bod) = Yoxy_doc_miner*(doc*doc_mineralisation*secs_per_day*12./1e3)*5.0
   ENDIF
END SUBROUTINE aed2_calculate_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_organic_matter(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: pon,don
   AED_REAL :: pop,dop
   AED_REAL :: poc,doc

   ! Temporary variables
   AED_REAL :: pon_flux,don_flux
   AED_REAL :: pop_flux,dop_flux
   AED_REAL :: poc_flux,doc_flux

   AED_REAL :: Fsed_pon,Fsed_don
   AED_REAL :: Fsed_pop,Fsed_dop
   AED_REAL :: Fsed_poc,Fsed_doc
   AED_REAL :: Psed_poc, Psed_pon, Psed_pop

!-------------------------------------------------------------------------------
!BEGIN


   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   ! Retrieve current (local) state variable values.
   pon = _STATE_VAR_(data%id_pon)! particulate organic matter
   don = _STATE_VAR_(data%id_don)! particulate organic matter
   pop = _STATE_VAR_(data%id_pop)! particulate organic matter
   dop = _STATE_VAR_(data%id_dop)! particulate organic matter
   poc = _STATE_VAR_(data%id_poc)! particulate organic matter
   doc = _STATE_VAR_(data%id_doc)! particulate organic matter

   IF (data%use_Fsed_model) THEN
      Fsed_pon = _STATE_VAR_S_(data%id_Fsed_pon)
      IF ( data%id_Fsed_don .NE. -1 ) Fsed_don = _STATE_VAR_S_(data%id_Fsed_don)
      IF ( data%id_Fsed_pop .NE. -1 ) Fsed_pop = _STATE_VAR_S_(data%id_Fsed_pop)
      IF ( data%id_Fsed_dop .NE. -1 ) Fsed_dop = _STATE_VAR_S_(data%id_Fsed_dop)
      IF ( data%id_Fsed_poc .NE. -1 ) Fsed_poc = _STATE_VAR_S_(data%id_Fsed_poc)
      IF ( data%id_Fsed_doc .NE. -1 ) Fsed_doc = _STATE_VAR_S_(data%id_Fsed_doc)
   ELSE
      Fsed_pon = data%Fsed_pon
      Fsed_don = data%Fsed_don * data%Ksed_don/(data%Ksed_don+don) * (data%theta_sed_don**(temp-20.0))
      Fsed_pop = data%Fsed_pop
      Fsed_dop = data%Fsed_dop * data%Ksed_dop/(data%Ksed_dop+dop) * (data%theta_sed_dop**(temp-20.0))
      Fsed_poc = data%Fsed_poc
      Fsed_doc = data%Fsed_doc * data%Ksed_doc/(data%Ksed_doc+doc) * (data%theta_sed_doc**(temp-20.0))
   ENDIF

   ! Calculate sedimentation flux (mmmol/m2/s) loss from benthos.
   IF (data%use_Psed_model) THEN
       Psed_poc = data%w_poc * max(zero_,poc)
       IF ( data%id_Psed_pon .NE. -1 ) Psed_pon = data%w_pon * max(zero_,pon)
       IF ( data%id_Psed_pop .NE. -1 ) Psed_pop = data%w_pop * max(zero_,pop)
   ELSE
       Psed_poc = zero_
       Psed_pon = zero_
       Psed_pop = zero_
   ENDIF

   pon_flux = Fsed_pon + Psed_pon
   don_flux = Fsed_don
   pop_flux = Fsed_pop + Psed_pop
   dop_flux = Fsed_dop
   poc_flux = Fsed_poc + Psed_poc
   doc_flux = Fsed_doc

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   _FLUX_VAR_(data%id_pon) = _FLUX_VAR_(data%id_pon) + (pon_flux)
   _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + (don_flux)
   _FLUX_VAR_(data%id_pop) = _FLUX_VAR_(data%id_pop) + (pop_flux)
   _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) + (dop_flux)
   _FLUX_VAR_(data%id_poc) = _FLUX_VAR_(data%id_poc) + (poc_flux)
   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + (doc_flux)

   ! Set sedimentation flux (mmmol/m2) as calculated by organic matter.
   IF (data%use_Psed_model) THEN
      _STATE_VAR_S_(data%id_Psed_poc) = Psed_poc
      _STATE_VAR_S_(data%id_Psed_pon) = Psed_pon
      _STATE_VAR_S_(data%id_Psed_pop) = Psed_pop
   ENDIF

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_amm) = _FLUX_VAR_B_(data%id_ben_amm) + (-amm_flux/secs_per_day)

   ! Also store sediment flux as diagnostic variable.
   IF (data%extra_diag) THEN
      _DIAG_VAR_S_(data%id_sed_pon) = -pon_flux
      _DIAG_VAR_S_(data%id_sed_don) = -don_flux
      _DIAG_VAR_S_(data%id_sed_pop) = -pop_flux
      _DIAG_VAR_S_(data%id_sed_dop) = -dop_flux
      _DIAG_VAR_S_(data%id_sed_poc) = -poc_flux
      _DIAG_VAR_S_(data%id_sed_doc) = -doc_flux
   ENDIF
END SUBROUTINE aed2_calculate_benthic_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_organic_matter(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: doc,poc, cpom, cdom
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
   doc = _STATE_VAR_(data%id_doc)
   poc = _STATE_VAR_(data%id_poc)

   ! Self-shading with explicit contribution from background OM concentration.
   extinction = extinction + (data%KeDOM*doc +data%KePOM*poc)

   IF (data%simRPools) THEN
     cdom = _DIAG_VAR_(data%id_cdom)    ! CDOM is in "/m" units here
     cpom = _STATE_VAR_(data%id_cpom)   ! CPOM is in "mmol/m3" units here

     extinction = extinction + (data%KeDOMR*cdom)
     extinction = extinction + (data%KeCPOM*cpom)
   ENDIF
END SUBROUTINE aed2_light_extinction_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fpon_miner(use_oxy,Rpon_miner,Kpon_miner,theta_pon_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Nitrogen
!
! Michaelis-Menten formulation for mineralisation
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpon_miner,Kpon_miner,theta_pon_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpon_miner = Rpon_miner * oxy/(Kpon_miner+oxy) * (theta_pon_miner**(temp-20.0))
   ELSE
      fpon_miner = Rpon_miner * (theta_pon_miner**(temp-20.0))
   ENDIF

END FUNCTION fpon_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdon_miner(use_oxy,Rdon_miner,Kdon_miner,theta_don_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdon_miner,Kdon_miner,theta_don_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdon_miner = Rdon_miner * oxy/(Kdon_miner+oxy) * (theta_don_miner**(temp-20.0))
   ELSE
      fdon_miner = Rdon_miner * (theta_don_miner**(temp-20.0))
   ENDIF

END FUNCTION fdon_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION fpop_miner(use_oxy,Rpop_miner,Kpop_miner,theta_pop_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Phosphorus
!
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpop_miner,Kpop_miner,theta_pop_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpop_miner = Rpop_miner * oxy/(Kpop_miner+oxy) * (theta_pop_miner**(temp-20.0))
   ELSE
      fpop_miner = Rpop_miner * (theta_pop_miner**(temp-20.0))
   ENDIF

END FUNCTION fpop_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdop_miner(use_oxy,Rdop_miner,Kdop_miner,theta_dop_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdop_miner,Kdop_miner,theta_dop_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdop_miner = Rdop_miner * oxy/(Kdop_miner+oxy) * (theta_dop_miner**(temp-20.0))
   ELSE
      fdop_miner = Rdop_miner * (theta_dop_miner**(temp-20.0))
   ENDIF

END FUNCTION fdop_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION fpoc_miner(use_oxy,Rpoc_miner,Kpoc_miner,theta_poc_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Carbon
!
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpoc_miner,Kpoc_miner,theta_poc_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpoc_miner = Rpoc_miner * oxy/(Kpoc_miner+oxy) * (theta_poc_miner**(temp-20.0))
   ELSE
      fpoc_miner = Rpoc_miner * (theta_poc_miner**(temp-20.0))
   ENDIF

END FUNCTION fpoc_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdoc_miner(use_oxy,Rdoc_miner,Kdoc_miner,theta_doc_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdoc_miner,Kdoc_miner,theta_doc_miner,oxy,temp
!
!-----------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdoc_miner = Rdoc_miner * oxy/(Kdoc_miner+oxy) * (theta_doc_miner**(temp-20.0))
   ELSE
      fdoc_miner = Rdoc_miner * (theta_doc_miner**(temp-20.0))
   ENDIF

END FUNCTION fdoc_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_organic_matter
