!################################################################################
!#                                                                              #
!# aed2_isotope.F90                                                             #
!#                                                                              #
!# Developed by :                                                               #
!#     AquaticEcoDynamics (AED) Group                                           #
!#     School of Earth & Environment                                            #
!# (C) The University of Western Australia                                      #
!#                                                                              #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org   #
!#                                                                              #
!#  -------------------------------------------------------------------------   #
!#                                                                              #
!# Created 9 May 2013                                                           #
!################################################################################


#include "aed2.h"
#define PURE
!
MODULE aed2_isotope
!----------------------------------------------------------------------------------
! aed2_isotope --- isotope biogeochemical model
!
! Created for Caboolture estuary to test O, N and C isotope submodel
! Sri Adiyanti & Matt Hipsey
! Modifications:
! (1) 07Nov2013 Inclusion of diagnostic variables for fluxes between pools (SA)
! (2) 05Feb2014 Inclusion of POC pool (SA)
! (3) 12Dec2014 Inclusion of Nitrogen and Oxygen pool (SA)
! (4) 01Jul2015 Inclusion of diagnostic for all C,N,O isotopes
!-----------------------------------------------------------------------------------
    USE aed2_core
    USE aed2_util,  ONLY: aed2_gas_piston_velocity, aed2_oxygen_sat

    IMPLICIT NONE

    PRIVATE    ! By default make everything private
!
    PUBLIC aed2_isotope_data_t
!
    TYPE,extends(aed2_model_data_t) :: aed2_isotope_data_t
      !# Variable identifiers
      INTEGER  :: id_oxy, id_oxy18 !18O in nox (nitrogen oxides)
      INTEGER  :: id_nit, id_nit15 !15N in nox (nitrogen oxides)
      INTEGER  :: id_no3, id_no315 !15N in no3 (nitrite)
      INTEGER  :: id_amm, id_amm15!15N in nh4 (ammonia)
      INTEGER  :: id_don, id_don15 !
      INTEGER  :: id_doc, id_doc13, id_toc
      INTEGER  :: id_dic, id_dic13
      INTEGER  :: id_pH, id_po4
      INTEGER  :: id_Fsed_amm,id_Fsed_nit,id_Fsed_oxy,id_sed_oxy
      INTEGER  :: id_ta
      INTEGER  :: id_ret  ! retention time
      INTEGER  :: id_chl  ! chlorophyll-a
      INTEGER  :: id_pCO2m ! pCO2 measured

      INTEGER  :: id_temp,id_salt,id_wind,id_I_0,id_par,id_h

      INTEGER  :: id_del13CDOC,id_del13CDIC
      INTEGER  :: id_del15NNOx, id_del18ONOx, id_del15NAmm

      ! IN/OUT Flux ================================================
      INTEGER  :: id_docsed2wat, id_dicsed2wat
      INTEGER  :: id_doc13sed2wat, id_dic13sed2wat
      INTEGER  :: id_dicwat2atm, id_dic13wat2atm
      INTEGER  :: id_dicatm2wat, id_dic13atm2wat

      ! Carbon
      INTEGER  :: id_dic2POC, id_dic132POC13      ! Primary Production
      INTEGER  :: id_poc2DIC, id_poc132DIC13      ! Respiration
      INTEGER  :: id_doc2DIC, id_doc132DIC13      ! Mineralisation
      INTEGER  :: id_phymor2doc, id_phymor2doc13  ! PhytoMortality
      INTEGER  :: id_phyex2doc, id_phyex2doc13    ! Exudation
      INTEGER  :: id_poc2docsed, id_poc132doc13sed! Sedimentation
      INTEGER  :: id_resp2dic, id_resp2dic13      ! Respiration to DIC
      INTEGER  :: id_det2doc, id_det2doc13        ! Decomposition Det to DOC
      INTEGER  :: id_ppdic, id_ppdic13            ! PP DIC uptake
      ! Nitrogen
      INTEGER  :: id_donsed2wat, id_don15sed2wat, id_dinsed2wat,id_amm15sed2wat,id_oxysed2wat
      INTEGER  :: id_atm_oxy_exch, id_oxy_sat, id_amm2n2, id_n22amm,id_oxy18sed2wat
      INTEGER  :: id_pon, id_totN
      INTEGER  :: id_ammsed2wat
      INTEGER  :: id_nitsed2wat, id_nit15sed2wat
      INTEGER  :: id_don2din, id_don2din15
      INTEGER  :: id_amm2no3, id_amm2no315       !nitrification
      INTEGER  :: id_no32n2, id_no32n215         !denitrification

      INTEGER  :: id_nitrif,id_denit,id_sed_amm,id_sed_nit

      !# Model parameters
      AED_REAL :: atmco2,atmn2,ionic,photoEff,                          &
        Rnitrif,Rdenit,Rommin,Rmort,Rdecomp,Rresp,Rresus,Rmain,         &
        fgrow,fdoc,fdetr,Cupt,CphyPOC,maxU,fexud,omega,Knitrif,Kdenit,  &
        Knit,Ktoc,Ko2,Kinhoxy,Kdon_miner,Rdon_miner,                    &
        Fsed_oxy,Fsed_amm,Fsed_nit,Fsed_don,Fsed_doc,Fsed_dic,          &
        Ksed_oxy,Ksed_amm,Ksed_nit,                                     &
        theta_nitrif,theta_nitrif_sed,theta_denit,theta_denit_sed,      &
        theta_don_miner,theta_ommin,theta_ommin_sed,                    &
        theta_oxy_sed,theta_resp,                                       &
        theta_don_exud,alpha_doc_sed,alpha_dic_sed,                     &
        alpha_doc_min,alpha_dic_phy,alpha_dic_atm,                      &
        alpha_doc_decomp,alpha_dic_resp,alpha_doc_mort,alpha_doc_exud,  &
        alpha_oxy_sed,alpha_nox_sed,alpha_don_sed,alpha_amm_sed,        &
        alpha_don_decomp,alpha_don_min,alpha_don_exud,alpha_don_mort,   &
        alpha_nit_nit,alpha_nit_denit, alpha_nit_ammon,                 &
        alpha_n2_pp,alpha_resp,alpha_cons

      LOGICAL  :: use_sed_model

      CONTAINS
         PROCEDURE :: define             => aed2_define_isotope
         PROCEDURE :: calculate_surface  => aed2_calculate_surface_isotope
         PROCEDURE :: calculate          => aed2_calculate_isotope
         PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_isotope
!        PROCEDURE :: mobility           => aed2_mobility_isotope
!        PROCEDURE :: light_extinction   => aed2_light_extinction_isotope
!        PROCEDURE :: delete             => aed2_delete_isotope
    END TYPE

!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE aed2_define_isotope(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS

    CLASS (aed2_isotope_data_t),INTENT(inout) :: data
    INTEGER,INTENT(in)                       :: namlst

!LOCALS

    AED_REAL    :: oxy_initial=4.5, oxy18_initial=4.5
    AED_REAL    :: doc_initial=4.5, doc13_initial=4.5
    AED_REAL    :: dic_initial=4.5, dic13_initial=4.5
    AED_REAL    :: nit_initial=4.5, nit15_initial=4.5
    AED_REAL    :: dop_initial=4.5, po4_initial=5
    AED_REAL    :: amm_initial=4.5, amm15_initial=4.5
    AED_REAL    :: don_initial=4.5, don15_initial=4.5
    AED_REAL    :: toc_initial=4, ret_initial=0

    AED_REAL    :: ph_initial=7.5
    AED_REAL    :: alkalinity_initial=1000
    AED_REAL    :: chl_initial=3
    AED_REAL    :: pCO2m_initial=500
    AED_REAL    :: retention_initial=0.0
    AED_REAL    :: atmco2   = 395.
    AED_REAL    :: atmn2    = 78.
    AED_REAL    :: ionic    = 1.5
    AED_REAL    :: photoEff = 5.8e-7 !m2s/uEs
    AED_REAL    :: Fsed_don = 3.5
    AED_REAL    :: Fsed_doc = 3.5
    AED_REAL    :: Fsed_dic = 3.5
    AED_REAL    :: Fsed_oxy = 3.5
    AED_REAL    :: Ksed_oxy = 30.0
    AED_REAL    :: Rnitrif  = 0.01
    AED_REAL    :: Rdenit = 0.01
    AED_REAL    :: Rommin = 0.01
    AED_REAL    :: Rmort = 0.01
    AED_REAL    :: Rdecomp = 0.01
    AED_REAL    :: Rdon_miner = 0.005
    AED_REAL    :: Rresus = 0.01
    AED_REAL    :: Rresp = 0.01
    AED_REAL    :: Rmain = 0.01
    AED_REAL    :: fgrow = 0.3
    AED_REAL    :: fdoc = 0.2
    AED_REAL    :: fdetr = 0.2
    AED_REAL    :: Cupt = 25
    AED_REAL    :: CphyPOC = 0.13
    AED_REAL    :: maxU = 0.2
    AED_REAL    :: fexud = 0.05
    AED_REAL    :: omega = 2.8e-6
    AED_REAL    :: Fsed_amm = 3.5
    AED_REAL    :: Fsed_nit = 3.5
    AED_REAL    :: Kdon_miner = 31.25
    AED_REAL    :: Knitrif = 100.0  !KNH4 Michaelis-Menten constant for NH4
    AED_REAL    :: Kdenit = 45.0    !KNO3 Michaelis-Menten constant for NO3
    AED_REAL    :: Knit = 5.0       !KN Michaelis-Menten constant for dissolved nitrogen
    AED_REAL    :: Ko2 = 15.0       !KO2 aerobic degradation, KO2
    AED_REAL    :: Ktoc = 60.0      !KTOC for organic matter
    AED_REAL    :: Kinhoxy = 50.0   !Inhibition tern for denitrification
    AED_REAL    :: Ksed_amm = 30.0
    AED_REAL    :: Ksed_nit = 30.0
    AED_REAL    :: alpha_photo = 5.8e-7 ! photosynthetic efficiency (m2s/(uEs))
    AED_REAL    :: Irrad = 200      ! irradiance at surface (watt/s)
    !AED_REAL   :: Kd = 0.3         ! light extinction coef.
    AED_REAL    :: theta_nitrif = 1.0
    AED_REAL    :: theta_nitrif_sed = 1.0
    AED_REAL    :: theta_denit = 1.0
    AED_REAL    :: theta_denit_sed = 1.0
    AED_REAL    :: theta_don_miner = 1.08
    AED_REAL    :: theta_ommin = 1.0
    AED_REAL    :: theta_ommin_sed = 1.04
    AED_REAL    :: theta_don_exud = 1.0
    AED_REAL    :: theta_oxy_sed = 1.08
    AED_REAL    :: theta_resp = 1.04
    CHARACTER(len=64)   :: nitrif_reactant_variable=''
    CHARACTER(len=64)   :: denit_product_variable=''
    CHARACTER(len=64)   :: Fsed_amm_variable=''
    CHARACTER(len=64)   :: Fsed_nit_variable=''
    CHARACTER(len=64)   :: Fsed_oxy_variable=''
    AED_REAL    :: alpha_doc_sed = 1.03
    AED_REAL    :: alpha_dic_sed = 1.03
    AED_REAL    :: alpha_amm_sed = 1.01 !1.013
    AED_REAL    :: alpha_oxy_sed = 1.01
    AED_REAL    :: alpha_nox_sed = 1.01
    AED_REAL    :: alpha_don_sed = 1.01

    AED_REAL    :: alpha_dic_phy = 1.04
    AED_REAL    :: alpha_dic_resp  = 0.97
    AED_REAL    :: alpha_dic_atm = 1.04
    AED_REAL    :: alpha_doc_decomp = 1.01
    AED_REAL    :: alpha_doc_mort = 1.03
    AED_REAL    :: alpha_doc_exud = 1.03
    AED_REAL    :: alpha_doc_min = 1.04
    AED_REAL    :: alpha_don_mort = 1.013
    AED_REAL    :: alpha_don_decomp = 1.10
    AED_REAL    :: alpha_don_min = 1.10
    AED_REAL    :: alpha_don_exud = 1.013
    AED_REAL    :: alpha_nit_nit = 1.017
    AED_REAL    :: alpha_nit_denit = 1.025
    AED_REAL    :: alpha_nit_ammon = 1.03
    AED_REAL    :: alpha_nit_assim = 1.005
    AED_REAL    :: alpha_oxy_assim = 1.0

    INTEGER     :: status

    NAMELIST /aed2_isotope/ oxy_initial, oxy18_initial,     &
                        nit_initial, nit15_initial,         &
                        amm_initial, amm15_initial,         &
                        don_initial, don15_initial,         &
                        po4_initial,                        &
                        doc_initial, doc13_initial,         &
                        dic_initial, dic13_initial,         &
                        toc_initial,                        &
                        ph_initial, alkalinity_initial,     &
                        ret_initial,                        &
                        chl_initial, pCO2m_initial,         &
                        atmco2,atmn2,ionic,                 &
                        Rnitrif,Rdenit,Rommin,              &
                        Rmort,Rdecomp,Rdon_miner,Rresp,     &
                        Rresus,Rmain,                       &
                        fgrow,fdoc,fdetr,Cupt,CphyPOC,      &
                        maxU,fexud,omega,Kdon_miner,        &
                        Knitrif,Kdenit,Knit,Ktoc,           &
                        Ko2,Kinhoxy,                        &
                        Fsed_oxy,Fsed_amm,Fsed_nit,         &
                        Fsed_don,Fsed_doc,Fsed_dic,         &
                        Ksed_oxy,Ksed_amm,Ksed_nit,         &
                        theta_nitrif,theta_nitrif_sed,      &
                        theta_denit,theta_denit_sed,        &
                        theta_don_miner,                    &
                        theta_ommin,theta_ommin_sed,        &
                        theta_oxy_sed,                      &
                        theta_resp,theta_don_exud,          &
                        alpha_doc_sed,alpha_dic_sed,        &
                        alpha_doc_min,                      &
                        alpha_dic_phy,alpha_dic_atm,        &
                        alpha_doc_decomp,alpha_dic_resp,    &
                        alpha_doc_mort,alpha_doc_exud,      &
                        alpha_oxy_sed,alpha_nox_sed,alpha_don_sed,  &
                        alpha_amm_sed,alpha_don_decomp,     &
                        alpha_don_min,alpha_don_exud,       &
                        alpha_don_mort,alpha_nit_nit,       &
                        alpha_nit_denit,alpha_nit_assim,    &
                        alpha_oxy_assim


!--------------------------------------------------------------
!BEGIN

    ! Read the namelist
    read(namlst,nml=aed2_isotope,iostat=status)
    IF (status /= 0) THEN
        print *,'Error reading namelist aed2_isotope'
        STOP
    ENDIF

    ! Store parameter values in our own derived type
    ! NB: all rates must be provided in values per day,
    ! and are converted here to values per second.

    ! OXYGEN
    data%Fsed_oxy      = Fsed_oxy/secs_per_day
    data%Ksed_oxy      = Ksed_oxy
    data%theta_oxy_sed =  theta_oxy_sed
    !data%use_sed_model = Fsed_oxy_variable .NE. ''

    ! NITROGEN
    data%Rnitrif  = Rnitrif/secs_per_day
    data%Rdenit   = Rdenit/secs_per_day
    data%Fsed_amm = Fsed_amm/secs_per_day
    data%Fsed_nit = Fsed_nit/secs_per_day
    data%Fsed_don = Fsed_don/secs_per_day
    data%Knitrif  = Knitrif
    data%Kdenit   = Kdenit
    data%Knit  = Knit
    data%Ktoc = Ktoc
    data%Ko2 = Ko2
    data%Kinhoxy = Kinhoxy
    data%Ksed_amm  = Ksed_amm
    data%Ksed_nit  = Ksed_nit
    data%theta_nitrif = theta_nitrif
    data%theta_nitrif_sed = theta_nitrif_sed
    data%theta_denit  = theta_denit
    data%theta_denit_sed = theta_denit_sed
    data%theta_ommin  = theta_ommin
    data%theta_ommin_sed  = theta_ommin_Sed
    data%theta_oxy_sed  =  theta_oxy_sed
    data%theta_resp = theta_resp

    ! CARBON, NITROGEN and Phyto
    data%Fsed_doc = Fsed_doc/secs_per_day
    data%Fsed_dic = Fsed_dic/secs_per_day
    data%Rommin   = Rommin/secs_per_day
    data%Rmort    = Rmort/secs_per_day
    data%Rdecomp  = Rdecomp/secs_per_day
    data%Rresus   = Rresus/secs_per_day
    data%Rresp    = Rresp/secs_per_day
    data%Rmain    = Rmain/secs_per_day
    data%Cupt = Cupt/secs_per_day
    data%maxU = maxU/secs_per_day
    data%atmco2 = atmco2
    data%atmn2 = atmn2
    data%ionic = ionic
    data%fdoc = fdoc
    data%fdetr = fdetr
    data%CphyPOC = CphyPOC
    data%fexud = fexud
    data%alpha_doc_sed = alpha_doc_sed
    data%alpha_dic_sed = alpha_dic_sed
    data%alpha_doc_min = alpha_doc_min
    data%alpha_don_min = alpha_don_min
    data%alpha_dic_phy = alpha_dic_phy
    data%alpha_dic_atm = alpha_dic_atm
    data%alpha_doc_decomp = alpha_doc_decomp
    data%alpha_dic_resp = alpha_dic_resp
    data%alpha_doc_mort = alpha_doc_mort
    data%alpha_don_mort = alpha_don_mort
    data%alpha_doc_exud = alpha_doc_exud
    data%alpha_don_exud = alpha_don_exud
    data%alpha_amm_sed = alpha_amm_sed

    ! Register state variables
    ! Oxygen
    data%id_oxy = aed2_define_variable( 'oxy','mmol/m**3','oxygen',     &
        oxy_initial,minimum=zero_)
    !data%id_oxy = aed2_define_variable('oxy','mmol/m**3','oxygen',   &
    !   oxy_initial,minimum=oxy_min,maximum=oxy_max)
    data%id_oxy18 = aed2_define_variable( 'oxy18','mmol/m**3','oxygen', &
        oxy18_initial,minimum=zero_)            !14N_18O_16O2

    ! Nitrogen
    data%id_nit = aed2_define_variable( 'nit','mmol/m**3','nitrogenoxides',     &
        nit_initial,minimum=zero_)  !14N_16O3
    data%id_nit15 = aed2_define_variable( 'nit15','mmol/m**3','nitrogenoxides15',&
        nit15_initial,minimum=zero_)  !15N_16O3
    data%id_amm = aed2_define_variable( 'amm','mmol/m**3','ammonium',            &
        amm_initial,minimum=zero_)
    data%id_amm15 = aed2_define_variable( 'amm15','mmol/m**3','ammonium',        &
        amm15_initial,minimum=zero_)
    data%id_don = aed2_define_variable( 'don','mmol/m**3','don',                 &
        don_initial,minimum=zero_)
    data%id_don15 = aed2_define_variable( 'don15','mmol/m**3','don15',           &
        don15_initial,minimum=zero_)

    ! Phosphorus
    data%id_po4 = aed2_define_variable( 'po4','mmol/m**3','FRP',                &
        po4_initial,minimum=zero_)  !

    ! Carbon
    data%id_doc = aed2_define_variable( 'doc','mmol/m**3','doc',                &
        doc_initial,minimum=zero_)
    data%id_doc13 = aed2_define_variable( 'doc13','mmol/m**3','doc',            &
        doc13_initial,minimum=zero_)
    data%id_dic = aed2_define_variable( 'dic','mmol/m**3','dic',                &
        dic_initial,minimum=zero_)
    data%id_dic13 = aed2_define_variable( 'dic13','mmol/m**3','dic',            &
        dic13_initial,minimum=zero_)
    data%id_toc = aed2_define_variable( 'toc','mmol/m**3','toc',                &
        toc_initial,minimum=zero_)

    ! Others
    data%id_ph = aed2_define_variable( 'pH','pH unit','pH',                 &
        ph_initial,minimum=zero_)
    data%id_ta = aed2_define_variable( 'ta',' mmol/m**3','alkalinity',      &
        alkalinity_initial,minimum=zero_)
    data%id_ret = aed2_define_variable( 'ret','time','water_retention',     &
        ret_initial,minimum=zero_)
    data%id_chl = aed2_define_variable( 'chl','ug/L','chlorophyll-a',       &
        chl_initial,minimum=zero_)
    data%id_pCO2m = aed2_define_variable( 'pCO2m','uAtm','pCO2 measured',   &
        pCO2m_initial,minimum=zero_)

    ! Register diagnostic variables

    data%id_oxy_sat = aed2_define_sheet_diag_variable('oxySat','mmol/m2/d','Oxygen saturation')
    data%id_atm_oxy_exch = aed2_define_sheet_diag_variable('atm_oxy_exch','mmol/m2/d',  &
        'Oxygen exchange across atm/water interface')
    data%id_sed_oxy = aed2_define_sheet_diag_variable('sed_oxy', 'mmol/m**2/d', &
        'Oxygen sediment flux')

    data%id_del13CDOC = aed2_define_diag_variable( 'del13CDOC','o/oo','del ratio of 13C-DOC')
    data%id_del13CDIC = aed2_define_diag_variable( 'del13CDIC','o/oo','del ratio of 13C-DIC')
    data%id_del18ONOx = aed2_define_diag_variable( 'del18ONOx','o/oo','del ratio of 18O-NOx')
    data%id_del15NAmm = aed2_define_diag_variable( 'del15NAmm','o/oo','del ratio of 15N-Amm')
    data%id_del15NNOx = aed2_define_diag_variable( 'del15NNOx','o/oo','del ratio of 15N-NOx')

! Fluxes
   ! Carbon

    data%id_docsed2wat = aed2_define_diag_variable( 'DOCsed2wat','mmol/m**3',   &
        'DOC release')
    data%id_doc13sed2wat = aed2_define_diag_variable( 'DOC13sed2wat','mmol/m**3',  &
        'DOC13 release')
    data%id_dicsed2wat = aed2_define_diag_variable( 'DICsed2wat','mmol/m**3',   &
        'DIC release')
    data%id_dic13sed2wat = aed2_define_diag_variable( 'DIC13sed2wat','mmol/m**3',  &
        'DIC13 release')
    data%id_dicwat2atm = aed2_define_diag_variable( 'DICwat2atm','mmol/m**3',  &
        'DIC outflux')
    data%id_dicatm2wat = aed2_define_diag_variable( 'DICatm2wat','mmol/m**3',  &
        'DIC influx')
    data%id_dic13wat2atm = aed2_define_diag_variable( 'DIC13wat2atm','mmol/m**3',  &
        'DIC13 outflux')
    data%id_dic13atm2wat = aed2_define_diag_variable( 'DIC13atm2wat','mmol/m**3',  &
        'DIC13 influx')
    data%id_doc2dic = aed2_define_diag_variable( 'DOC2DIC','mmol/m**3',  &
        'DOC convert to DIC')
    data%id_doc132dic13 = aed2_define_diag_variable( 'DOC132DIC13','mmol/m**3',  &
        'DOC13 convert to DIC13')
    data%id_phymor2doc = aed2_define_diag_variable( 'PhyMor2DOC','mmol/m**3',  &
        'POC Phyto Mortality to DOC')
    data%id_phyMor2doc13 = aed2_define_diag_variable( 'PhyMor2DOC13','mmol/m**3',  &
        'POC Phyto Mortality to DOC13')
    data%id_phyEx2doc = aed2_define_diag_variable( 'PhyEx2DOC','mmol/m**3',  &
        'POC Phyto Excudates to DOC')
    data%id_phyEx2doc13 = aed2_define_diag_variable( 'PhyEx2DOC13','mmol/m**3', &
        'POC Phyto Excudates to DOC13')
    data%id_det2doc = aed2_define_diag_variable( 'Det2DOC','mmol/m**3', &
        'Detritus to DOC')
    data%id_det2doc13 = aed2_define_diag_variable( 'Det2DOC13','mmol/m**3', &
        'Detritus to DOC13')
    data%id_resp2dic = aed2_define_diag_variable( 'Resp2DIC','mmol/m**3', &
        'Resp to DIC')
    data%id_resp2DIC13 = aed2_define_diag_variable( 'Resp2DIC13','mmol/m**3', &
        'Resp to DIC13')
    data%id_ppdic = aed2_define_diag_variable( 'PPDIC','mmol/m**3',         &
        'PP to DIC')
    data%id_ppdic13 = aed2_define_diag_variable( 'PPDIC13','mmol/m**3', &
        'PP to DIC13')

    ! Nitrogen ==================================================================

    data%id_donsed2wat = aed2_define_diag_variable( 'DONsed2wat','mmol/m**3',  &
        'DON release')
    data%id_don15sed2wat = aed2_define_diag_variable('DON15sed2wat','mmol/m**2/d', &
        'DON15 sediment flux')
    !data%id_dinsed2wat = aed2_define_diag_variable( 'DINsed2wat','mmol/m**3',  &
    !   'DIN release')
    data%id_oxysed2wat = aed2_define_diag_variable( 'OXYsed2wat','mmol/m**3',  &
        'OXY release')
    data%id_oxy18sed2wat = aed2_define_diag_variable( 'OXY18sed2wat','mmol/m**3',  &
        'OXY18 release')
    data%id_totn = aed2_define_diag_variable( 'totN','mmol/m**3','total N')
    data%id_ammsed2wat = aed2_define_diag_variable('ammsed2wat','mmol/m**2/d', &
        'Ammonium sediment flux')
    data%id_amm15sed2wat = aed2_define_diag_variable('amm15sed2wat','mmol/m**2/d', &
        'Amm15 sediment flux')
    data%id_nitsed2wat = aed2_define_diag_variable('nitsed2wat','mmol/m**2/d', &
        'NOx sediment flux')
    data%id_nit15sed2wat = aed2_define_diag_variable('nit15sed2wat','mmol/m**2/d', &
        'NOx15 sediment flux')

    !Denitrification
    data%id_no32n2 = aed2_define_diag_variable( 'NO32N2','mmol/m**3',   &
        'NO3 convert to N2')
    data%id_no32n215 = aed2_define_diag_variable( 'NO32N215','mmol/m**3',   &
        'NO315 convert to N215')
    data%id_denit = aed2_define_diag_variable('denit','mmol/m**3/d',     &
        'De-nitrification rate')

    !N2 fixation (N2 to ammonium)
    data%id_n22amm = aed2_define_diag_variable( 'N22AMM','mmol/m**3',   &
        'N2 convert to ammonium')

    !Nitrification (NH4 to NO2 followed straight to NO3)
    data%id_amm2no3 = aed2_define_diag_variable( 'AMM2NO3','mmol/m**3',     &
        'NH4 convert to NO3')
    data%id_amm2no315 = aed2_define_diag_variable( 'AMM2NO315','mmol/m**3',     &
        'NH415 convert to NO315')
    data%id_nitrif = aed2_define_diag_variable('nitrif','mmol/m**3/d',  &
        'Nitrification rate')

    !Mineralisation
    data%id_don2din = aed2_define_diag_variable('DON2DIN','mmol/m**3',  &
        'DON convert to DIN')
    data%id_don2din15 = aed2_define_diag_variable('DON2DIN15','mmol/m**3',  &
        'DON15 convert to DIN15') !NH4
    !data%id_mineral = aed2_define_diag_variable('mineralisation','mmol/m**3/d',    &
    !    'mineralisation rate')

! Register external state variable dependencies
!   data%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
!   IF (data%use_oxy) THEN
!     data%id_oxy = aed2_locate_variable(nitrif_reactant_variable)
!   ENDIF
!   data%use_no2 = denit_product_variable .NE. '' !This means n2 module switched on
!   IF (data%use_no2) data%id_denit_product = aed2_locate_variable(denit_product_variable)

!   data%use_sed_model = Fsed_amm_variable .NE. ''
!   IF (data%use_sed_model) THEN
!     data%id_Fsed_amm = aed2_locate_global_sheet(Fsed_amm_variable)
!     data%id_Fsed_nit = aed2_locate_global_sheet(Fsed_nit_variable)
!   ENDIF

   ! Register environmental dependencies

    data%id_temp = aed2_locate_global( 'temperature') ! Temperature (degrees Celsius)
    data%id_salt = aed2_locate_global( 'salinity')    ! Salinity (psu)
    data%id_wind = aed2_locate_global_sheet('wind_speed')
    data%id_par = aed2_locate_global('par')
    data%id_I_0 = aed2_locate_global_sheet('par_sf')
    data%id_h = aed2_locate_global('layer_ht')
    !data%id_sed = aed2_locate_global('sed_zone')

END SUBROUTINE aed2_define_isotope
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed2_calculate_isotope(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_isotope model
!-------------------------------------------------------------------------------
!ARGUMENTS
    CLASS (aed2_isotope_data_t),INTENT(in) :: data
    TYPE (aed2_column_t),INTENT(inout) :: column(:)
    INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
    AED_REAL    :: temp,salt,wind,Io,h,par   !State variables (from global)
    AED_REAL    :: oxy,amm,don,doc,dic !State variables
    AED_REAL    :: nit, dop, toc, po4, N2Fix, NPhy, omegaFix
    AED_REAL    :: doc13,dic13 !State variables
    AED_REAL    :: nit15, amm15, don15, oxy18
    AED_REAL    :: alkalinity
    AED_REAL    :: retention
    AED_REAL    :: chl,pCO2m,GPP,NPP,PBmax,nlim
    AED_REAL    :: Tabs
    AED_REAL    :: pCO2a,pCO2d,Ko,kCO2,pH,HCO3,CO3,windHt
    AED_REAL    :: denitrification,nitrification,mineralisation
    AED_REAL    :: resuspension,respiration,don_mineralisation
    AED_REAL    :: respDIC
    AED_REAL    :: DOC2DIC,DOC132DIC13,CPhy,CMort,PhyMor2DOC,PhyMor2DOC13
    AED_REAL    :: PhyEx2DOC,PhyMor2DON,PhyEx2DOC13,Det2DOC,Det2DOC13
    AED_REAL    :: PhyEx2DON,PhyMor2DON15,PhyEx2DON15,Det2DON,Det2DON15
    AED_REAL    :: PPN2,PPN215,fNH4,NO3Cons,NO315Cons,AMMCons,AMM15Cons
    AED_REAL    :: OXYCons,OXY18Cons,oxyNit,oxy18Nit
    AED_REAL    :: OXYConsAMM,OXYConsAMM_Rat,OXY18ConsAMM,OXYConsNO3,OXY18ConsNO3
    AED_REAL    :: DOC132DIC13_Rat,Det2DOC13_Rat,Resp2DIC13_Rat,PPDIC13_Rat
    AED_REAL    :: DON152DIN15_Rat,Det2DON15_Rat,Resp2DIN15_Rat,PPDIN15_Rat
    AED_REAL    :: PPN215_Rat,NO3Cons_Rat,OXYConsno3_Rat,respAMM_Rat
    AED_REAL    :: DON152DIN15,respAMM,respAMM15
    AED_REAL    :: PhyEx2DOC13_Rat,PhyMor2DOC13_Rat,AMMCons_Rat
    AED_REAL    :: PhyEx2DON15_Rat,PhyMor2DON15_Rat,respNH4_Rat
    AED_REAL    :: Cdetr,Resp2DIC,Resp2DIC13,PPDIC,PPDIC13,alphaDICpp
    AED_REAL    :: AMM2NO3, NO32N2, DON2DIN
    AED_REAL    :: AMM2NO315,NO32N215,OXY18NO3,DON2DIN15
    AED_REAL    :: del13CDIC,del13CDOC,del15NNOx,del18ONOx,del15NDON,del15NAmm
    AED_REAL    :: sigmaE,diff_oxy
    AED_REAL,PARAMETER :: Yoxy_ommin  = 1.  !ratio of oxygen to carbon utilised during mineralisation
    AED_REAL,PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to nitrogen utilised during nitrification from NH4
    AED_REAL,PARAMETER :: photoEff = 5.8e-7 !m2s/uEs
    AED_REAL,PARAMETER :: CChl = 50.        !C/Chl ratio grC/gChl
    AED_REAL,PARAMETER :: redfN = 0.1509433962  !16/106!molN /molC in phytoplankton (redfield ratio)
    AED_REAL,PARAMETER :: redfP = 0.0094339623  !1/106 !molP /molC in phytoplankton (redfield ratio)
    AED_REAL,PARAMETER :: omegaNol = 2.8e-6 !metabolism constant (s)
    AED_REAL,PARAMETER :: stoiO2AMM = 2.    !4/2 !ratio of oxygen/ammonium during nitrification 2NH4+4O2-->2NO3+2H2O+2H+
    AED_REAL,PARAMETER :: stoiO2TOC = 1.    !ratio of oxygen/carbon during mineralisation
    AED_REAL,PARAMETER :: stoiNO3POC = 0.8905660377 !94.4/106   !Stoichio nitrate/POC during denitrification
    AED_REAL,PARAMETER :: stoiO2POC = 1.3018867925  !138/106    !Stoichio O2/POC during consumption of nitrate
    AED_REAL,PARAMETER :: stoiN2POC = 0.6666666667  !2/3    !Stoichio N2/POC during N2 fixation
    AED_REAL,PARAMETER :: del13CPOC = -23   ! median d13CPOC during CR2
    AED_REAL,PARAMETER :: del13CPhy = -34   ! Hamilton & Lewis (1992)
    AED_REAL,PARAMETER :: del13CDet = -28   ! d13C Detritus (Otero et al, 2000)
    AED_REAL,PARAMETER :: del15NPON = 7     ! median d15N Cloern (2002)
    AED_REAL,PARAMETER :: del15NN2 = -1     ! d15N at N2 atmosphere
    AED_REAL,PARAMETER :: del15NDet = 0     !
    AED_REAL,PARAMETER :: Rf_13C = 0.011180 !(VPDB), Fry(2006)Appendix 1, p. 285
    AED_REAL,PARAMETER :: Rf_15N = 0.0036765 !(VPDB)
    AED_REAL,PARAMETER :: Rf_18O = 0.0020052 !(SMOW)
    AED_REAL,PARAMETER :: Rf_17O = 0.0003799 !(SMOW)
    AED_REAL,PARAMETER :: ret_dt = 1.0      !per time step
    AED_REAL,PARAMETER :: Kd = 0.3          ! light extinction coef.
    AED_REAL,PARAMETER :: alpha_n2_pp = 0.9978      ! fractionation of N2 fixing
    AED_REAL,PARAMETER :: alpha_resp = 1.025        ! fractionation of respiration
    AED_REAL,PARAMETER :: alpha_cons = 1.01         ! fractionation of consumption

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
    oxy = _STATE_VAR_(data%id_oxy)! oxygen
    nit = _STATE_VAR_(data%id_nit)! nitrogen oxides (NO2+NO3)
    nit15 = _STATE_VAR_(data%id_nit15)! nitrogen15
    oxy18 = _STATE_VAR_(data%id_oxy18)! oxygen18
    amm = _STATE_VAR_(data%id_amm)! ammonium
    amm15 = _STATE_VAR_(data%id_amm15)! ammonium
    don = _STATE_VAR_(data%id_don)! dissolved organic nitrogen
    don15 = _STATE_VAR_(data%id_don15)!
    po4 = _STATE_VAR_(data%id_po4)!
    doc = _STATE_VAR_(data%id_doc)! dissolved organic carbon
    doc13 = _STATE_VAR_(data%id_doc13)! dissolved organic carbon C13
    dic = _STATE_VAR_(data%id_dic)! dissolved organic carbon
    dic13 = _STATE_VAR_(data%id_dic13)! dissolved organic carbon C13
    toc = _STATE_VAR_(data%id_toc)! total organic carbon
    pH = _STATE_VAR_(data%id_pH)            ! pH
    alkalinity = _STATE_VAR_(data%id_ta)    ! alkalinity
    retention = _STATE_VAR_(data%id_ret)    ! retention time
    chl = _STATE_VAR_(data%id_chl)          ! chlorophyll-a
    pCO2m = _STATE_VAR_(data%id_pCO2m)      ! measured pCO2

    ! Retrieve current environmental conditions.
    temp = _STATE_VAR_(data%id_temp)    ! temperature [C]
    salt = _STATE_VAR_(data%id_salt)    ! salinity    [psu]
    wind = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10 m above surface (m/s)
    par = _STATE_VAR_(data%id_par)      ! local photosynthetically active radiation
    Io = _STATE_VAR_S_(data%id_I_0)     ! surface short wave radiation
    h = _STATE_VAR_(data%id_h)          ! dz

    del13CDIC=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
    del13CDOC=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
    del15NNOx=(((nit15/(nit-nit15))/Rf_15N)-1.)*1e3
    del18ONOx=(((oxy18/((3*nit)-oxy18))/Rf_18O)-1.)*1e3
    del15NAmm=(((amm15/(amm-amm15))/Rf_15N)-1.)*1e3
    del15NDON=(((don15/(don-don15))/Rf_15N)-1.)*1e3

    ! Define some intermediate rates units [/s]
    nitrification   = fnitrif(data,amm,oxy,temp)
    denitrification = fdenit(data,toc,nit,oxy,temp)
    mineralisation  = fommin(data,toc,oxy,temp)
    don_mineralisation = fdon_miner(data,oxy,temp)
    resuspension = fthetased(data,oxy,temp)
    respiration = fresp(data,temp)

    !print *, '========== A1. ISO CALC, at the beginning ......'
    !print *, 'del15N_nox=',del15NNOx,'del18O_nox=',del18ONOx

    ! GPP & NPP ===================
    nlim = (nit+amm)/(nit+amm+data%Knit)    !limiting nutrient
    Pbmax = 1/CChl*exp(0.33+0.102*temp)     !Volta, et al 2014, CChl=gC/gChl
    CPhy = CChl*chl/12  !uM C
    !assumed 1 m
    h = 1
    par = 100
    GPP = Pbmax * nlim * CPhy * (1 - exp(-photoEff/Pbmax*par*exp(-Kd*h)))  !uM
    NPP = (GPP/h * (1-data%fexud) * (1-data%fgrow)) - data%Rmain*CPhy  !uM
    !print *, 'NPP,nlim,Pbmax,CPhy,GPP=', NPP,nlim,Pbmax,CPhy,GPP

    ! MINERALISATION
    DON2DIN = don*don_mineralisation
    DON2DIN15 = don15*don_mineralisation
    DON152DIN15_Rat = (del15NDON/1000+1)*Rf_15N/data%alpha_don_min
    DON152DIN15 = DON2DIN-(1/(1+DON152DIN15_Rat)*DON2DIN)

    !oxyMin = oxy*mineralisation
    oxyNit =-Yoxy_nitrif*amm*nitrification !Oxy used for nitrification

    ! DENITRIFICATION
    NO32N2 = stoiNO3POC*nit*denitrification
    NO32N215 = stoiNO3POC*nit15*denitrification
    !NO32N2 = nit*denitrification
    !NO32N215 = nit15*denitrification

    ! NITRIFICATION
    AMM2NO3 = amm*nitrification
    AMM2NO315 = amm15*nitrification
    !OXY18NO3 = oxy18*nitrification*3        !3*AMM2NO315

    NO3Cons = 0.0
    NO315Cons = 0.0
    !NO3Cons = redfN*(1-fNH4)*NPP
    !NO3Cons_Rat = (del15NNOx/1000+1)*Rf_15N/alpha_cons
    !NO315Cons = NO3Cons-(1/(1+NO3Cons_Rat)*NO3Cons)

    !OXY consumed during resp NH4 and NO3
    OXYConsAMM = fNH4*NPP
    OXYConsNO3 = stoiO2POC*(1-fNH4)*NPP

    ! Phytoplankton mortality ===============
    ! Phyto Death PON-->DON
    CPhy = CChl*chl/12  !uM C
    CMort = data%Rmort*CPhy !uM
    PhyMor2DOC = data%fdoc*CMort !uM
    PhyMor2DON = redfN*PhyMor2DOC !uM
    PhyMor2DON15_Rat = (del15NPON/1000+1)*Rf_15N/data%alpha_don_mort
    PhyMor2DON15 = PhyMor2DON-(1/(1+PhyMor2DON15_Rat)*PhyMor2DON)

    ! Phytoplankton exudation ===============
    ! PON phyto-->DON
    PhyEx2DON = redfN*data%fexud*(data%maxU*Cphy)
    PhyEx2DON15_Rat = (del15NPON/1000+1)*Rf_13C/data%alpha_don_exud
    PhyEx2DON15 = PhyEx2DON-(1/(1+PhyEx2DON15_Rat)*PhyEx2DON)

    !stoiO2AMM 4/2 !ratio of oxygen/ammonium during nitrification 2NH4+4O2-->2NO3+2H2O+2H+

    ! Set temporal derivatives
    !NIT, NIT15, OXY18
    _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + AMM2NO3 - NO32N2 - NO3Cons

    _FLUX_VAR_(data%id_nit15) = _FLUX_VAR_(data%id_nit15) + AMM2NO315 - NO32N215 - NO315Cons

    _FLUX_VAR_(data%id_oxy18) = _FLUX_VAR_(data%id_oxy18) +         &
        (-NO32N215)*oxy18/nit15 + (-NO315Cons)*oxy18/nit15 + AMM2NO315*oxy18/nit15
    !same fractionation for nit15 and oxy18

    !AMM, AMM15
    _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (-AMM2NO3) + DON2DIN
    _FLUX_VAR_(data%id_amm15) = _FLUX_VAR_(data%id_amm15) + (-AMM2NO315) + DON2DIN15

    ! OXY
    diff_oxy = 0.
    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (diff_oxy) + (-oxyNit) - &
        (OXYConsAMM + OXYConsNO3)

    ! DON
    _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + PhyMor2DON + PhyEx2DON - DON2DIN !+ N2Fix
    _FLUX_VAR_(data%id_don15) = _FLUX_VAR_(data%id_don15) + PhyMor2DON15 + PhyEx2DON15 - DON2DIN15

    ! PO4
    _FLUX_VAR_(data%id_po4) = _FLUX_VAR_(data%id_po4) -  redfP*(oxyNit+denitrification-NPP)

   ! OXY
    !-AerDeg + NH4phy - 2Nit + (O2atm)
    !_FLUX_VAR_(data%id_oxy18) = _FLUX_VAR_(data%id_oxy18) -        &
    !   oxy18Min + (OXY18ConsAMM + OXY18ConsNO3) - stoiO2AMM*AMM2NO315

    !Isotopes signals for diagnostic
    _DIAG_VAR_(data%id_del13CDIC) =  del13CDIC
    _DIAG_VAR_(data%id_del13CDOC) =  del13CDOC
    _DIAG_VAR_(data%id_del15NNOx) =  del15NNOx
    _DIAG_VAR_(data%id_del15NAmm) =  del15NAmm
    _DIAG_VAR_(data%id_del18ONOx) =  del18ONOx

    !Nitrogen Fluxes, useful for diagnostic
    _DIAG_VAR_(data%id_nitrif) =  nitrification
    _DIAG_VAR_(data%id_denit) =  denitrification
    _DIAG_VAR_(data%id_don2din) = DON2DIN
    _DIAG_VAR_(data%id_amm2no3) = AMM2NO3   !amm*nitrification
    _DIAG_VAR_(data%id_no32n2) = NO32N2     !nox*denitrification

    if ((del15NNOx<-5) .OR. (del15NNOx>60)) then
        print *,  'error.. nox=',nit, 'nox15=', nit15, 'del15N_nox=',del15NNOx, 'del18O_nox=',del18ONOx
        STOP
    endif
    if ((del18ONOx<-5) .OR. (del18ONOx>60)) then
        print *,  'error.. nox=',nit, 'nox15=', nit15, 'del15N_nox=',del15NNOx, 'del18O_nox=',del18ONOx
        STOP
    endif

  RETURN

  !===================================================================================================
  !===================================================================================================

    ! Retrieve current environmental conditions.
    temp = _STATE_VAR_(data%id_temp)    ! temperature [C]
    salt = _STATE_VAR_(data%id_salt)    ! salinity    [psu]
    wind = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10 m above surface (m/s)
    par = _STATE_VAR_(data%id_par)      ! local photosynthetically active radiation
    Io = _STATE_VAR_S_(data%id_I_0)     ! surface short wave radiation
    h = _STATE_VAR_(data%id_h)          ! dz

    del13CDIC=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
    del13CDOC=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
    del15NNOx=(((nit15/(nit-nit15))/Rf_15N)-1.)*1e3
    del18ONOx=(((oxy18/((3*nit)-oxy18))/Rf_18O)-1.)*1e3
    del15NAmm=(((amm15/(amm-amm15))/Rf_15N)-1.)*1e3
    del15NDON=(((don15/(don-don15))/Rf_15N)-1.)*1e3

    print *, '========== A1. ISO CALC, at the beginning ......'
    print *, 'del15N_nox=',del15NNOx,'del18O_nox=',del18ONOx, 'del13CDIC=',del13CDIC, 'del13CDOC=',del13CDOC
    print *, 'oxy=', oxy, 'amm=',amm, 'don=',don, 'don15=',don15,'po4=',po4,'toc=',toc, 'alkalinity =',alkalinity, 'pH=',pH
    print *, '============================='

    !correct sim temp (overheating about 1 to 2 degC)
    ! temp=temp-1

    ! Define some intermediate rates units [/day]
    nitrification   = fnitrif(data,amm,oxy,temp)
    denitrification = fdenit(data,toc,nit,oxy,temp)
    mineralisation  = fommin(data,toc,oxy,temp)
    resuspension = fthetased(data,oxy,temp)
    respiration = fresp(data,temp)

    !print *, 'ISOTOPE CALC'
    !print *, 'nit,den,min,res,resp=', nitrification, denitrification,mineralisation,resuspension,respiration

    ! CO2 uptake by Phyto ------------------
    ! 1st calculate Solubility, Ko (mol/L/atm) (Weiss, 1974)
    Tabs = temp + 273.15
    Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
         + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
    Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
    Ko = exp(Ko)

    ! 2nd calculate gas piston velocity --------
    windHt = 10.  !reference wind measurement height
    !windHt = 48-15.  ! wind measurement height at: -26.96  Lon: 152.96 Height: 48.0 m NSL: 15m
    kCO2 = aed2_gas_piston_velocity(windHt,wind,temp,salt)

     ! 3rd dissolved pCO2
    pCO2a = aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3) / Ko !pCO2 in atm
    pCO2d =  kCO2 * (Ko*1e6) * pCO2a ! dissolved CO2 in uM

    ! 4th calculate C uptake, etc, note C uptake is in mgCO2/mgChla/day
    PPDIC = chl*data%Cupt/12 !uM
    alphaDICpp = 1.0007+(0.0283*0.69*(pCO2d-5.6)/pCO2d) ! (Yoshioka,1997)
    PPDIC13_Rat = (del13CDIC/1000+1)*Rf_13C/alphaDICpp
    PPDIC13 = PPDIC-(1/(1+PPDIC13_Rat)*PPDIC)

    ! N2 fixation ----------------- N2-->NH4
    PPN2 = stoiN2POC*PPDIC
    PPN215_Rat = (del15NN2/1000+1)*Rf_15N/data%alpha_N2_pp
    PPN215 = PPN2-(1/(1+PPN215_Rat)*PPN2)

    !====================================================
    ! GPP & NPP ===================
    nlim = (nit+amm)/(nit+amm+data%Knit)    !limiting nutrient
    Pbmax = 1/CChl*exp(0.33+0.102*temp)     !Volta, et al 2014, CChl=gC/gChl
    CPhy = CChl*chl/12  !uM C
    !assumed 1 m
    h = 1
    par = 100
    GPP = Pbmax * nlim * CPhy * (1 - exp(-data%photoEff/Pbmax*par*exp(-Kd*h)))  !uM
    NPP = (GPP/h * (1-data%fexud) * (1-data%fgrow)) - data%Rmain*CPhy  !uM
    !print *, 'NPP,nlim,Pbmax,CPhy,GPP=', NPP,nlim,Pbmax,CPhy,GPP

    ! N2 fixation ===============
    NPhy = redfN*CPhy
    omegaFix = omegaNol * (0.022+(1/(0.25+exp(3/(temp-12)-0.5)+exp(-(500/(temp-12)-25)))))
    sigmaE = 1-(0.1-1)**20
    N2Fix = omegaFix * sigmaE * NPhy

    ! Phytoplankton mortality ===============
    ! Phyto Death POC phyto-->DOC, PON-->DON
    CMort = data%Rmort*CPhy !uM
    PhyMor2DOC = data%fdoc*CMort !uM
    PhyMor2DOC13_Rat = (del13CPOC/1000+1)*Rf_13C/data%alpha_doc_mort
    PhyMor2DOC13 = PhyMor2DOC-(1/(1+PhyMor2DOC13_Rat)*PhyMor2DOC)

    PhyMor2DON = redfN*PhyMor2DOC !uM
    PhyMor2DON15_Rat = (del15NPON/1000+1)*Rf_15N/data%alpha_don_mort
    PhyMor2DON15 = PhyMor2DON-(1/(1+PhyMor2DON15_Rat)*PhyMor2DON)

    !print *, 'alpha_don_mort=', data%alpha_don_mort

    ! Exudation POC phyto-->DOC, PON phyto-->DON
    PhyEx2DOC = data%fexud*(data%maxU*Cphy)
    PhyEx2DOC13_Rat = (del13CPOC/1000+1)*Rf_13C/data%alpha_doc_exud
    PhyEx2DOC13 = PhyEx2DOC-(1/(1+PhyEx2DOC13_Rat)*PhyEx2DOC)

    PhyEx2DON = data%fexud*(data%maxU*Cphy)*redfN
    PhyEx2DON15_Rat = (del15NPON/1000+1)*Rf_13C/data%alpha_don_exud
    PhyEx2DON15 = PhyEx2DON-(1/(1+PhyEx2DON15_Rat)*PhyEx2DON)

    !print *, 'alpha_don_exud=', data%alpha_don_exud

    ! Mineralisation DOC-->DIC, DON-->DIN, and OXYconsumption
    DOC2DIC = doc*mineralisation
    !del13CDOC = (((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
    DOC132DIC13_Rat = (del13CDOC/1000+1)*Rf_13C/data%alpha_doc_min
    DOC132DIC13 = DOC2DIC-(1/(1+DOC132DIC13_Rat)*DOC2DIC)

    DON2DIN = don*mineralisation
    DON2DIN15 = don15*mineralisation
    DON152DIN15_Rat = (del15NDON/1000+1)*Rf_15N/data%alpha_don_min
    DON152DIN15 = DON2DIN-(1/(1+DON152DIN15_Rat)*DON2DIN)

    oxyNit = Yoxy_nitrif*amm*nitrification

    ! Decomposition of POC Detritus -->DOC, PON -->DON
    Cdetr = (1-data%CphyPOC)*Cphy/data%CphyPOC ! C in detritus = (Det/POC)*Cphyto/(Cphyto/POC)
    Det2DOC = data%Rdecomp*Cdetr  !uM
    Det2DOC13_Rat = (del13CDet/1000+1)*Rf_13C/data%alpha_doc_decomp;
    Det2DOC13 = Det2DOC-(1/(1+Det2DOC13_Rat)*Det2DOC);

    Det2DON = data%Rdecomp*Cdetr*redfN  !uM
    Det2DON15_Rat = (del15NDet/1000+1)*Rf_15N/data%alpha_don_decomp;
    Det2DON15 = Det2DON-(1/(1+Det2DON15_Rat)*Det2DON);

    ! RESPIRATION =========================

    ! DIC due to dead phyto respiration (heterotrophic respiration)
    Resp2DIC = (1-data%fDetr-data%fdoc)*CMort
    Resp2DIC13_Rat = (del13CPOC/1000+1)*Rf_13C/data%alpha_dic_resp
    Resp2DIC13 = Resp2DIC-(1/(1+Resp2DIC13_Rat)*Resp2DIC);

    ! phyto respiration (aerobic degradation)
    respDIC = CPhy*respiration

    !NH4 produced during respiration aerobic and heterethropic
    respAMM = redfN*CPhy*respiration
    respAMM_Rat = (del15NAMM/1000+1)*Rf_15N/alpha_resp
    respAMM15 = respAMM-(1/(1+respAMM_Rat)*respAMM)

    !NH4 consumed during resp
    fNH4 = amm/(10+amm) ! swith between NH45 and NO3 utilisation
    AMMCons = redfN*fNH4*NPP
    AMMCons_Rat = (del15NAmm/1000+1)*Rf_15N/alpha_cons
    AMM15Cons = AMMCons-(1/(1+AMMCons_Rat)*AMMCons)

    !NO3 consumed during resp
    NO3Cons = stoiO2POC *(1-fNH4)*NPP
    NO3Cons_Rat = (del15NNOx/1000+1)*Rf_15N/alpha_cons
    NO315Cons = NO3Cons-(1/(1+NO3Cons_Rat)*NO3Cons)

    !OXY consumed during resp NH4
    OXYConsAMM = fNH4*NPP
    OXYConsAMM_Rat = (del18ONOx/1000+1)*Rf_18O/alpha_cons
    OXY18ConsAMM = OXYConsAMM-(1/(1+OXYConsAMM_Rat)*OXYConsAMM)

    !print *, 'alpha_cons=', data%alpha_cons

    !OXY consumed during resp NO3
    OXYConsNO3 = stoiO2AMM*(1-fNH4)*NPP
    OXYConsNO3_Rat = (del18ONOx/1000+1)*Rf_18O/alpha_cons
    OXY18ConsNO3 = OXYConsNO3 -(1/(1+OXYConsNO3_Rat)*OXYConsNO3 )

    ! NITRIFICATION
    AMM2NO3 = amm*nitrification
    AMM2NO315 = amm15*nitrification

    ! DENITRIFICATION
    NO32N2= stoiNO3POC*nit*denitrification
    NO32N215 = stoiNO3POC*nit15*denitrification

    ! Resuspension DOC from sediment
    !  resusDOC = docSed*resuspenstion
    !  resusDOC13 = (del13DOCSed+1000)/(del13CDOCSed+1000+(1000/Rf_13C))*resusDOC*data%alpha_doc_sed

    !========================================
    ! Set temporal derivatives
    ! PHYTO (chl-a)

    print *, ' A2. .....before FLUX_VAR....'
    print *, 'del15N_nox=',del15NNOx,'del18O_nox=',del18ONOx, 'del13CDIC=',del13CDIC, 'del13CDOC=',del13CDOC
    print *, 'oxy=', oxy, 'amm=',amm, 'don=',don, 'don15=',don15,'po4=',po4,'toc=',toc, 'alkalinity =',alkalinity, 'pH=',pH

    _FLUX_VAR_(data%id_chl) = _FLUX_VAR_(data%id_chl) + (NPP-CMort)*12/CChl

    !  DOC
    _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) -         &
        DOC2DIC +(PhyMor2DOC+PhyEx2DOC) + Det2DOC
        !-aeroDeg + (phyMor+phyEx+det)
    _FLUX_VAR_(data%id_doc13) = _FLUX_VAR_(data%id_doc13) -     &
        DOC132DIC13 + (PhyMor2DOC13+PhyEx2DOC13) + Det2DOC13

    !  TOC
    _FLUX_VAR_(data%id_toc) = _FLUX_VAR_(data%id_toc) -             &
        DOC2DIC +(PhyMor2DOC+PhyEx2DOC) + Det2DOC + (CMort+Cdetr)   &
        - denitrification
        ! -AerDeg + PhyMor - Denit

    !  DIC
    _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) +         &
        DOC2DIC+Resp2DIC-PPDIC
        ! +aeroDeg + Resp - uptake
    _FLUX_VAR_(data%id_dic13) = _FLUX_VAR_(data%id_dic13) +     &
        DOC132DIC13+Resp2DIC13-PPDIC13

    ! OXY
     diff_oxy = 0.
    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (diff_oxy) + (-oxyNit)

    !_FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) -            &
    !   oxyMin + (OXYConsAMM + OXYConsNO3) - stoiO2AMM*AMM2NO3
        !-AerDeg + NH4phy - 2Nit + (O2atm)
    _FLUX_VAR_(data%id_oxy18) = _FLUX_VAR_(data%id_oxy18) -     &
        oxy18Nit + (OXY18ConsAMM + OXY18ConsNO3) - stoiO2AMM*AMM2NO315

    ! Set temporal derivatives

    ! AMM
    _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) +         &
        respAMM - AMMCons - AMM2NO3 + PPN2
    _FLUX_VAR_(data%id_amm15) = _FLUX_VAR_(data%id_amm15) + &
        respAMM15 - AMM15Cons - AMM2NO315 + PPN215

   ! NOx
    _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) -         &
        NO32N2 - NO3Cons + AMM2NO3

   _FLUX_VAR_(data%id_nit15) = _FLUX_VAR_(data%id_nit15) -      &
        NO32N215 - NO315Cons + AMM2NO315

    ! DON
    _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) - DON2DIN + N2Fix
    _FLUX_VAR_(data%id_don15) = _FLUX_VAR_(data%id_don15) - DON2DIN15

   ! PO4
    _FLUX_VAR_(data%id_po4) = _FLUX_VAR_(data%id_po4) -  &
        redfP*(oxyNit+denitrification-NPP)

!-------------------------------------------------------

   del13CDIC=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
   del13CDOC=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
   del15NNOx=(((nit15/(nit-nit15))/Rf_15N)-1.)*1e3
   del18ONOx=(((oxy18/((3*nit)-oxy18))/Rf_18O)-1.)*1e3
   del15NAmm=(((amm15/(amm-amm15))/Rf_15N)-1.)*1e3
   del15NDON=(((don15/(don-don15))/Rf_15N)-1.)*1e3
!-----------------------------------------------------

   _DIAG_VAR_(data%id_del13CDOC) =  del13CDOC
   _DIAG_VAR_(data%id_del13CDIC) =  del13CDIC
   _DIAG_VAR_(data%id_del15NNOx) =  del15NNOx
   _DIAG_VAR_(data%id_del15NAmm) =  del15NAmm
   _DIAG_VAR_(data%id_del18ONOx) =  del18ONOx

   !===================================

   ! Export diagnostic variables
   !Carbon Fluxes, useful for diagnostic

   _DIAG_VAR_(data%id_doc2dic) =  DOC2DIC
   _DIAG_VAR_(data%id_doc132dic13) =  DOC132DIC13
   _DIAG_VAR_(data%id_phymor2doc) =  PhyMor2DOC
   _DIAG_VAR_(data%id_phymor2doc13) =  PhyMor2DOC13
   _DIAG_VAR_(data%id_phyex2doc) =  PhyEx2DOC
   _DIAG_VAR_(data%id_phyex2doc13) =  PhyEx2DOC13
   _DIAG_VAR_(data%id_det2doc) =  Det2DOC
   _DIAG_VAR_(data%id_det2doc13) =  Det2DOC13
   _DIAG_VAR_(data%id_resp2dic) =  Resp2DIC
   _DIAG_VAR_(data%id_resp2dic13) =  Resp2DIC13
   _DIAG_VAR_(data%id_ppdic) =  PPDIC
   _DIAG_VAR_(data%id_ppdic13) =  PPDIC13

   !Nitrogen Fluxes, useful for diagnostic
   _DIAG_VAR_(data%id_nitrif) =  nitrification
   _DIAG_VAR_(data%id_denit) =  denitrification
   _DIAG_VAR_(data%id_don2din) = DON2DIN
   _DIAG_VAR_(data%id_amm2no3) = AMM2NO3    !amm*nitrification
   _DIAG_VAR_(data%id_no32n2) = NO32N2      !nox*denitrification

 !  print *, 'toc=', toc, 'nitrification=',nitrification, 'denitrification=',denitrification
 !  print *, 'DON2DIN=',DON2DIN, 'AMM2NO3=',AMM2NO3, 'NO32AMM=',AMM2NO3

print *, '=== A3. isotope_calculate ending ==='
print *, 'del15N_nox=',del15NNOx,'del18O_nox=',del18ONOx, 'del13CDIC=',del13CDIC, 'del13CDOC=',del13CDOC
print *, '================================'

if (del15NNOx<-10) then
 print *,  'error.. nox=',nit, 'nox15=', nit15, 'del15N_nox=',del15NNOx, 'del18O_nox=',del18ONOx
 STOP
endif
if (del18ONOx<-10) then
 print *,  'error.. nox=',nit, 'nox15=', nit15, 'del15N_nox=',del15NNOx, 'del18O_nox=',del18ONOx
 STOP
endif

if (del15NNOx>70) then
 print *,  'error.. nox=',nit, 'nox15=', nit15, 'del15N_nox=',del15NNOx, 'del18O_nox=',del18ONOx
 STOP
endif

if (del18ONOx>70) then
 print *,  'error.. nox=',nit, 'nox15=', nit15, 'del15N_nox=',del15NNOx, 'del18O_nox=',del18ONOx
 STOP
endif

END SUBROUTINE aed2_calculate_isotope
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed2_calculate_benthic_isotope(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED isotope.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: amm,amm15,nit,nit15,oxy,oxy18,don,don15,doc,dic
   AED_REAL :: doc13,dic13
   AED_REAL :: pH ! for testing only

   ! Temporary variables
   AED_REAL :: oxy_flux,oxy18_flux,Fsed_oxy
   AED_REAL :: amm_flux,amm15_flux,nit_flux,nit15_flux,don_flux,don15_flux
   AED_REAL :: Fsed_amm,Fsed_nit,Fsed_don
   AED_REAL :: doc_flux,doc13_flux,dic_flux,dic13_flux
   AED_REAL :: Fsed_doc,Fsed_dic
   AED_REAL :: del13CDIC,del13CDOC
   AED_REAL :: del15NNOx, del15NDON, del15Namm, del18ONOx
   AED_REAL :: doc13_flux_Rat,dic13_flux_Rat
   AED_REAL :: nit15_flux_Rat,don15_flux_Rat,amm15_flux_Rat
   AED_REAL :: oxy18_flux_Rat

   ! Parameters
   AED_REAL,PARAMETER :: d13DOCSed = -28
   AED_REAL,PARAMETER :: Rf_13C = 0.011180 !(VPDB), Fry(2006)Appendix 1, p. 285
   AED_REAL,PARAMETER :: Rf_15N = 0.0036765; !(AIR)
   AED_REAL,PARAMETER :: Rf_18O = 0.0020052; !(SMOW)
!
!-------------------------------------------------------------------------------
!BEGIN

RETURN
   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen
   oxy18 = _STATE_VAR_(data%id_oxy18)! oxygen
   amm = _STATE_VAR_(data%id_amm)! ammonium
   amm15 = _STATE_VAR_(data%id_amm15)! ammonium15
   nit = _STATE_VAR_(data%id_nit)! nitrate and nitrite
   nit15 = _STATE_VAR_(data%id_nit15)!
   don = _STATE_VAR_(data%id_don)! dissolved organic nitrogen
   don15 = _STATE_VAR_(data%id_don15)!
   doc = _STATE_VAR_(data%id_doc)! dissolved organic carbon
   dic = _STATE_VAR_(data%id_dic)! dissolved inorganic carbon
   doc13 = _STATE_VAR_(data%id_doc13)! doc13
   dic13 = _STATE_VAR_(data%id_dic13)! dic13
   !toc =  _STATE_VAR_(data%id_toc)!
   pH = _STATE_VAR_(data%id_pH)! pH for testing only

   !-----------------------
   !correct sim temp (overheating about 1 to 2 degC)
   ! temp=temp-1
   !------------------------
  !calculate delC before fluxes
   del13Cdic=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
   del13Cdoc=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
   del15NNox=(((nit15/(nit-nit15))/Rf_15N)-1.)*1e3
   del18ONOx=(((oxy18/((3*nit)-oxy18))/Rf_18O)-1.)*1e3
   del15Namm=(((amm15/(amm-amm15))/Rf_15N)-1.)*1e3

   !print *, 'B1. Do Benthos Start =================='
   !print *, 'dic=', dic, 'dic13=', dic13, 'del13C=', (((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
   !print *, 'nit=', nit, 'nit15=', nit15,'del15N=', (((nit15/(nit-nit15))/Rf_15N)-1.)*1e3
   !print *, 'del15NNOx=',del15NNox,'del18ONOx=',del18ONOx

   !TEMP FIX ===================================
   !IF (del15NNOx>40) THEN
   !   del15NNOx=40
   !   nit15=nit*(del15NNOx+1000)/(del15NNOx+1000+(1000/Rf_15N))
   !ENDIF

   !IF (del18ONOx>30) THEN
    !   del18ONOx=30
    !   oxy18=3*nit*(del18ONOx+1000)/(del18ONOx+1000+(1000/Rf_18O))
    !   print *, 'fixing, del18ONOx=', del18ONOx
   !ENDIF
   !print *, 'AFTER FIXING, del15NNOx=', (((nit15/(nit-nit15))/Rf_15N)-1.)*1e3, 'del18ONOx=', (((oxy18/((3*nit)-oxy18))/Rf_18O)-1.)*1e3
   !------
   Fsed_oxy = data%Fsed_oxy
   Fsed_amm = data%Fsed_amm
   Fsed_nit = data%Fsed_nit
   Fsed_don = data%Fsed_don
   Fsed_doc = data%Fsed_doc
   Fsed_dic = data%Fsed_dic

   ! Sediment flux dependent on oxygen and temperature

   !oxy_flux = Fsed_oxy * oxy/(data%Ksed_oxy+oxy) * (data%theta_sed_oxy**(temp-20.0))
   oxy_flux = Fsed_oxy * oxy/(data%Ksed_oxy+oxy) * (data%theta_oxy_sed**(temp-20.0))
   amm_flux = Fsed_amm * data%Ksed_amm/(data%Ksed_amm+oxy) * (data%theta_ommin_sed**(temp-20.0))
   nit_flux = Fsed_nit * oxy/(data%Ksed_nit+oxy) * (data%theta_ommin_sed**(temp-20.0))
   don_flux = Fsed_don * (data%theta_ommin_sed**(temp-20.0))
   doc_flux = Fsed_doc * (data%theta_ommin_sed**(temp-20.0))
   !dic_flux = -oxy_flux
   dic_flux = -0.8*oxy_flux

    !Calculate C13, N15, O18 Ratios after fractionation
    doc13_flux_Rat=(del13CDOC/1000+1)*Rf_13C/data%alpha_doc_sed
    dic13_flux_Rat=(del13CDIC/1000+1)*Rf_13C/data%alpha_dic_sed
    nit15_flux_Rat=(del15NNOx/1000+1)*Rf_15N/data%alpha_nox_sed
    amm15_flux_Rat=(del15Namm/1000+1)*Rf_15N/data%alpha_amm_sed
    don15_flux_Rat=(del15NDON/1000+1)*Rf_15N/data%alpha_don_sed
    oxy18_flux_Rat=(del18ONOx/1000+1)*Rf_18O/data%alpha_oxy_sed

    doc13_flux=(1-1/(1+doc13_flux_Rat))*doc_flux
    dic13_flux=(1-1/(1+dic13_flux_Rat))*dic_flux
    nit15_flux=(1-1/(1+nit15_flux_Rat))*nit_flux
    amm15_flux=(1-1/(1+amm15_flux_Rat))*amm_flux
    don15_flux=(1-1/(1+don15_flux_Rat))*don_flux
    oxy18_flux=(1-1/(1+oxy18_flux_Rat))*nit_flux*3


   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (oxy_flux)
    _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (amm_flux)
    _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + (nit_flux)
    _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + (don_flux)
    _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + (doc_flux)
    _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (dic_flux)

   _FLUX_VAR_(data%id_doc13) = _FLUX_VAR_(data%id_doc13) + (doc13_flux)
   _FLUX_VAR_(data%id_dic13) = _FLUX_VAR_(data%id_dic13) + (dic13_flux)
   _FLUX_VAR_(data%id_nit15) = _FLUX_VAR_(data%id_nit15) + (nit15_flux)
   _FLUX_VAR_(data%id_amm15) = _FLUX_VAR_(data%id_amm15) + (amm15_flux)
   _FLUX_VAR_(data%id_don15) = _FLUX_VAR_(data%id_don15) + (don15_flux)
   _FLUX_VAR_(data%id_oxy18) = _FLUX_VAR_(data%id_oxy18) + (oxy18_flux)

   !_FLUX_VAR_(data%id_nit15) = _FLUX_VAR_(data%id_nit15) + &
    ! (nit_flux-(1/(1+nit15_flux_Rat)*nit_flux))
   !_FLUX_VAR_(data%id_oxy18) = _FLUX_VAR_(data%id_oxy18) + &
   ! 3*(nit_flux-(1/(1+oxy18_flux_Rat)*nit_flux))
   !_FLUX_VAR_(data%id_amm15) = _FLUX_VAR_(data%id_amm15) + &
    ! (amm_flux-(1/(1+amm15_flux_Rat)*amm_flux))
   !_FLUX_VAR_(data%id_don15) = _FLUX_VAR_(data%id_don15) + &
    !(don_flux-(1/(1+don15_flux_Rat)*don_flux))

 ! o.m is degraded, release isotopically light C, lowering delC13 of dic in lower layer
 ! o.m is degraded, release isotopically light N, lowering delN15 of NOx in lower layer
   _DIAG_VAR_(data%id_DICsed2wat) =  dic_flux
   _DIAG_VAR_(data%id_DOCsed2wat) =  doc_flux
   _DIAG_VAR_(data%id_nitsed2wat) =  nit_flux
   _DIAG_VAR_(data%id_ammsed2wat) =  amm_flux
   _DIAG_VAR_(data%id_donsed2wat) =  don_flux
   _DIAG_VAR_(data%id_oxysed2wat) =  oxy_flux

   _DIAG_VAR_(data%id_DIC13sed2wat) =  dic13_flux
   _DIAG_VAR_(data%id_DOC13sed2wat) =  doc13_flux
   _DIAG_VAR_(data%id_nit15sed2wat) =  nit15_flux
   _DIAG_VAR_(data%id_amm15sed2wat) =  amm15_flux
   _DIAG_VAR_(data%id_don15sed2wat) =  don15_flux
   _DIAG_VAR_(data%id_oxy18sed2wat) =  oxy18_flux

   !_DIAG_VAR_(data%id_DIC13sed2wat) =  dic_flux-(1/(1+dic13_flux_Rat)*dic_flux)
   !_DIAG_VAR_(data%id_DOC13sed2wat) =  doc_flux-(1/(1+doc13_flux_Rat)*doc_flux)
   !_DIAG_VAR_(data%id_nit15sed2wat) =  nit_flux-(1/(1+nit15_flux_Rat)*nit_flux)
   !_DIAG_VAR_(data%id_amm15sed2wat) =  amm_flux-(1/(1+amm15_flux_Rat)*amm_flux)
   !_DIAG_VAR_(data%id_don15sed2wat) =  don_flux-(1/(1+don15_flux_Rat)*don_flux)
   !_DIAG_VAR_(data%id_oxy18sed2wat) =  3*(nit_flux-(1/(1+oxy18_flux_Rat)*nit_flux))

 ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_oxy) = oxy_flux
   _DIAG_VAR_S_(data%id_sed_amm) = amm_flux
   _DIAG_VAR_S_(data%id_sed_nit) = nit_flux

   !calculate delC after fluxes
   !del13Cdic=(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
   !del13Cdoc=(((doc13/(doc-doc13))/Rf_13C)-1.)*1e3
   !del15NNox=(((nit15/(nit-nit15))/Rf_15N)-1.)*1e3
   !del18ONOx=(((oxy18/((3*nit)-oxy18))/Rf_18O)-1.)*1e3
   !del15Namm=(((amm15/(amm-amm15))/Rf_15N)-1.)*1e3

!print *, '== B2. set bottom exchange ===='
!print *, 'doc=', doc, 'doc13=', doc13, 'del13CDOC=',del13CDOC
!print *, 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',del13CDIC
!print *, 'del15NNOx=',del15Nnox, 'del18ONOx=',del18ONOx
!print *, 'del13CDOC=',del13CDOC, 'del13CDIC=',del13CDIC
!print *, '================================'

END SUBROUTINE aed2_calculate_benthic_isotope
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed2_calculate_surface_isotope(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed oxygen model
!
! see aed2_oxygen_do_surface_exchange
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, pCO2m

   ! State ! Temporary variables
   AED_REAL :: oxy,dic,pH,alkalinity,dic13,atm13C ! Temporary variables
   AED_REAL :: nit,nit15,oxy18
   AED_REAL :: oxy_atm_flux = zero_
   AED_REAL :: Coxy_air = zero_ !Dissolved oxygen in the air phase
   AED_REAL :: koxy_trans = zero_
   AED_REAL :: f_pres  = 1.0      ! Pressure correction function only applicable at high altitudes
   AED_REAL :: pCO2,FCO2
   AED_REAL :: Ko, KCO2
   AED_REAL :: Tabs,windHt
   AED_REAL :: Faw, Fwa, FC13
   AED_REAL :: HCO3, CO3, CO2
   AED_REAL :: Rdic, RCO2a
   AED_REAL :: eta_dic_g, eta_HCO3_g, eta_CO3_g, eta_aq_g, eta_as
   AED_REAL :: alpha_dic_g, alpha_HCO3, alpha_CO3, alpha_aq_g, alpha_k
   AED_REAL :: FC13wat, FC13air, diff_oxy
   AED_REAL :: del15NNox, del18ONOx !temporary

   AED_REAL, PARAMETER :: Rf_13C = 0.011180 !(VPDB), Fry(2006)Appendix 1, p. 285
   AED_REAL, PARAMETER :: Rf_15N=0.0036765; !(AIR)
   AED_REAL, PARAMETER :: Rf_18O=0.0020052; !(SMOW)
   AED_REAL, PARAMETER :: delC13_atm = -7.8 !o/oo, delC13 of atmospheric CO2
   AED_REAL, PARAMETER :: delN15_atm = -1   !o/oo, delN15 of atmospheric N2
   AED_REAL, PARAMETER :: eta_k = -0.81 !o/oo, kinetic isotope fract during CO2 gas transfer  (Zhang, 1995)
   !AED_REAL, PARAMETER :: eta_aq_g = -1.19 !o/oo, equilibrium fractination during CO2 dissolution (Zhang, 1995)
!
!-------------------------------------------------------------------------------
!BEGIN

   !---- First, define environment context
   !Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind = _STATE_VAR_S_(data%id_wind) ! Wind speed at 10 m above surface (m/s)
   windHt = 10.

    !---- Second, flux oxygen
    ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! Concentration of oxygen in surface layer

   koxy_trans = aed2_gas_piston_velocity(windHt,wind,temp,salt)

   !print *, 'koxy_trans=', koxy_trans, 'wind=', wind, 'temp=', temp, 'salt=', salt

   ! First get the oxygen concentration in the air phase at interface
   ! Taken from Riley and Skirrow (1974)
   f_pres = 1.0
   Coxy_air = f_pres * aed2_oxygen_sat(salt,temp)

   ! Get the oxygen flux
   oxy_atm_flux = koxy_trans * (Coxy_air - oxy)

   !print *, 'Coxy_air=',Coxy_air, 'oxy_atm_flux=',oxy_atm_flux,'oxy=',oxy

   !---- Third, flux CO2 --------
   ! Retrieve current (local) state variable values.
   dic = _STATE_VAR_(data%id_dic)          ! Concentration of inorganic carbon in surface layer
   dic13 = _STATE_VAR_(data%id_dic13)      ! Concentration of inorganic carbon13 in surface layer
   pH = _STATE_VAR_(data%id_pH)            ! pH in surface layer
   alkalinity = _STATE_VAR_(data%id_ta)    ! alkalinity in surface layer

   !-----------------------------

   !print *,'======================================'
!print *, '== A. Start surface exchange ', 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',((dic13/(dic-dic13)/Rf_13C)-1.)*1e3
!print *, 'pH=',pH, 'TA=', alkalinity

   kCO2 = aed2_gas_piston_velocity(windHt,wind,temp,salt)

   ! Solubility, Ko (mol/L/atm) (Weiss, 1974)
   Tabs = temp + 273.15
   Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
        + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
   Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
   Ko = exp(Ko)

   ! pCO2 in surface water layer

   !print *, '== A. surface exchange ', 'dic=', dic, 'dic13=', dic13, 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3, 'pH=', pH

   pCO2 = aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3) / Ko

   !print *, '== aed2_carbon_co2 , pH=',pH,' alkalinity=',alkalinity,' HCO3=',HCO3,' CO3=',CO3,' pCO2=',pCO2,' Ko=',Ko

   ! mmol/m2/s = m/s * mmol/L/atm * atm
   !FCO2 = kCO2 * Ko*1e6 * (pCO2 - data%atmco2)

   !------ Yanti correction (20/5/2013)
   ! pCO2 is actually in uatm (=ppm)
     FCO2 = - kCO2 * Ko*1e6 * (pCO2 - data%atmco2) *1e-6 ! dCO2/dt

   ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
   !_FLUX_VAR_T_(data%id_dic) = FCO2 later in the end..

!------------------------------------------------------------------------
!
!   !Lynch-Stieglitz etal (1995)
!   Faw=kCO2 * Ko*1e6 * data%atmco2 *1e-6        !atm to surface
!   Fwa=kCO2 * Ko*1e6 * pCO2 *1e-6               !water to atm
!
!   atm13C = (delC13_atm+1000)/((delC13_atm+1000+(1000/Rf_13C))) !note this is fraction of 13C_atm,13C/C
!
!   FC13=((-0.373/Tabs+1.00019)*data%alpha_dic_atm*Faw*atm13C)-((-9.866/Tabs+1.02412)*data%alpha_dic_atm*Fwa*dic13/dic)
!---------------------------------------------------------------------------
!
!   Based on Zhang (1995)
!   Note, Zhang definition: alpha=k13/k (not k12/k13), R=13C/12C
!
!   First, calculate the individual temperature-dependent fractination of each CO2, HCO3, CO3

    eta_HCO3_g = -(0.1141*temp)+10.78     ! fractination between bicarb and CO2
    eta_CO3_g = -(0.052*temp)+7.22        ! fractination between carbonate and CO2
    eta_aq_g = -(0.0049*temp)-1.31        ! fractination during gass dissolution CO2

!   Second, calculate kinetic isotope fractionation during CO2 exhange (alpha_k) and
!   equilibrium fractionation during CO2(g) dissolution (alpha_aq_g)

    alpha_k = (eta_k/1000)+1
    alpha_aq_g = (eta_aq_g/1000)+1

!   Third, calculate total fractination during CO2 expansion
    eta_as = (alpha_k*alpha_aq_g-1)*1000

!   Fourth, calculate the equilibrium fractination factor between DIC and CO2(g)
    CO2=aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3)
    eta_dic_g = (CO2/dic*eta_as)+(HCO3/dic*eta_HCO3_g)+(CO3/dic*eta_CO3_g)
    alpha_dic_g = (eta_dic_g/1000)+1

!   Last, calculate the d13CO2/dt
    Rdic = dic13/(dic-dic13) ! 13C/12C ratio of dic in water
    RCO2a = (delC13_atm/1000+1)*Rf_13C  ! 13C/12C ratio of CO2 in atmosphere

    !print *,'Rdic=',Rdic,'RCO2a=',Rco2a,'kCO2=',kco2

    FC13 = - kCO2*alpha_k*alpha_aq_g*(Ko*1e6)*((pCO2*1e-6)*Rdic/alpha_dic_g-(data%atmco2*1e-6)*RCO2a)

    FC13wat=-kCO2*alpha_k*alpha_aq_g*(Ko*1e6)*((pCO2*1e-6)*Rdic/alpha_dic_g)
    FC13air=-kCO2*alpha_k*alpha_aq_g*(Ko*1e6)*((data%atmco2*1e-6)*RCO2a)

!   Diagnosis

    _DIAG_VAR_(data%id_DICwat2atm) =  - kCO2 * Ko*1e6 * pCO2*1e-6
    _DIAG_VAR_(data%id_DICatm2wat) =  - kCO2 * Ko*1e6 * (-data%atmco2)*1e-6

    _DIAG_VAR_(data%id_DIC13wat2atm) =  FC13wat
    _DIAG_VAR_(data%id_DIC13atm2wat) =  FC13air

!print *, 'pCO2=', pCO2, 'atmco2=', data%atmco2,'(pCO2*1e-6)*Rdic/alpha_dic_g=', &
!   (pCO2*1e-6)*Rdic/alpha_dic_g, '(data%atmco2*1e-6)*RCO2a=', (data%atmco2*1e-6)*RCO2a
!print *, 'pCO2=', pCO2, 'atmco2=', data%atmco2

  !_FLUX_VAR_T_(data%id_dic13) = FC13 later in the end

!print *, '== A. End surface exchange ', 'dic=', dic, 'dic13=', dic13, &
! 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3, 'FC13=',FC13

!----------------------
   IF(FCO2>0.) THEN
!     ! pCO2 is high and flux is from the water to air
     _FLUX_VAR_T_(data%id_dic13) = FCO2*(dic13/dic)/data%alpha_dic_atm !losing CO2

   ELSE !FCO2<0., gaining
    ! pCO2 is low and flux is from air to water
     !atm13C = ((delC13_atm/1000.)+1)* Rf_13C !note, this is Rs
      atm13C = (delC13_atm+1000)/((delC13_atm+1000+(1000/Rf_13C))) !note this is fraction of 13C_atm, 13C/C
     _FLUX_VAR_T_(data%id_dic13) = FCO2*atm13C

   END IF
!----------------------

!   if (FCO2>0) then
!      print *, 'influx,wtr=', FC13wat,'air=',FC13air,'FC13=',FC13, &
!   'Rs=', dic13/(dic-dic13),'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
!   else
!      print *, 'outflux,wtr=', FC13wat,'air=',FC13air,'FC13=',FC13, &
!   'Rs=', dic13/(dic-dic13), 'del13CDIC=',(((dic13/(dic-dic13))/Rf_13C)-1.)*1e3
!   endif

   !---- Fourth, set diagnostics
   ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).

   ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
   _FLUX_VAR_T_(data%id_oxy) = oxy_atm_flux
   _FLUX_VAR_T_(data%id_dic) = FCO2
   _FLUX_VAR_T_(data%id_dic13) = FC13

   _DIAG_VAR_S_(data%id_atm_oxy_exch) = oxy_atm_flux
   _DIAG_VAR_S_(data%id_oxy_sat) =  Coxy_air

END SUBROUTINE aed2_calculate_surface_isotope
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION fnitrif(data,amm,oxy,temp)

!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in)   :: data
   AED_REAL,INTENT(in)                      :: amm,oxy,temp
   AED_REAL                                 :: fnit
!-------------------------------------------------------------------------------
!BEGIN
   fnit = data%theta_nitrif**(temp-20.0)
   fnitrif = data%Rnitrif * fnit * (amm/(data%Knitrif+amm)) * (oxy/(data%Ko2+oxy)) !Volta
   !fnitrif = data%Rnitrif * fnit * oxy/(data%Knitrif+oxy)
   !fnitrif = Rnitrif * oxy/(Knitrif+oxy) * (theta_nitrif**(temp-20.0)) aed2 N module

END FUNCTION fnitrif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdenit(data,toc,nox,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                    :: toc,nox,oxy,temp
   AED_REAL                               :: fden
!
!-------------------------------------------------------------------------------
!BEGIN
   fden = data%theta_denit**(temp-20.0)
   fdenit = data%Rdenit * fden * (toc/(data%Ktoc+toc)) *        &
    (nox/data%Kdenit+nox) * data%Kinhoxy/(data%Kinhoxy+oxy)     ! Volta

   !fdenit = data%Rdenit * data%Kdenit/(data%Kdenit+oxy) * fden  !aed N module


END FUNCTION fdenit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION fdon_miner(data,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                    :: oxy,temp
   AED_REAL                               :: fdonminer
!
!-------------------------------------------------------------------------------
!BEGIN
   fdonminer = data%theta_don_miner**(temp-20.0)

   fdon_miner = data%Rdon_miner*fdonminer*oxy/(data%Kdon_miner+oxy)

END FUNCTION fdon_miner
!###############################################################################
PURE AED_REAL FUNCTION fommin(data,toc,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for organic matter mineralisation
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated. Ko2 = K_O2 for aerobic degradation
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in)   :: data
   AED_REAL,INTENT(in)                      :: toc,oxy,temp
   AED_REAL                                 :: fom
!
!-------------------------------------------------------------------------------
!BEGIN

   fom = data%theta_ommin**(temp-20.0)
   fommin = data%Rommin * fom *  (toc/(data%Ktoc+toc)) * oxy/(data%Ko2+oxy)

END FUNCTION fommin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION fresp(data,temp)
!-------------------------------------------------------------------------------
! Temperature function for resuspension
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                   :: temp
!
!-------------------------------------------------------------------------------
!BEGIN

   fresp = data%Rresp * (data%theta_resp**(temp-20.0))

END FUNCTION fresp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION fthetased(data,oxy,temp)
!-------------------------------------------------------------------------------
! Temperature function for resuspension
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                   :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN

   fthetased = data%Rresus * oxy/(data%Ko2+oxy) * (data%theta_oxy_sed**(temp-20.0))

END FUNCTION fthetased
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
!PURE AED_REAL FUNCTION aed2_carbon_co2(data,temp,dic,pH,alkalinity)
AED_REAL FUNCTION aed2_carbon_co2(data,temp,dic,pH,alkalinity,HCO3,CO3)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_isotope_data_t),INTENT(in)  :: data
   AED_REAL, INTENT(IN)                   :: dic, temp, alkalinity
   AED_REAL, INTENT(OUT)                  :: HCO3, CO3
   AED_REAL, INTENT(INOUT)                :: pH

!
!LOCALS
   ! Temporary variables
   AED_REAL :: pKh, pKa1, pKa2, pKw, TA
   AED_REAL :: Kh, Kw, Ka1, Ka2, is_f
   AED_REAL :: alf_CO2, alf_HCO3, alf_CO3
   AED_REAL :: H, CO2
!-------------------------------------------------------------------------------
!BEGIN

   ! Acidity constants temperature dependence (mol/kg atm)at ionic strength=0
   ! Note (Y13/4/13), have checked the values are within the range of values provided
   ! in Butler 1998 Table 10.1 (selected values of equilibrium constant extrapolated to I=0)

   RETURN

   pKh = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
   pKa1 = 0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
   pKa2 = 0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
   pKw = 0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909

   ! Ionic strength dependence

   ! 1st calculate ionic strength function is_f
   is_f = (((SQRT(data%ionic)) / (1+SQRT(data%ionic))) - 0.20*data%ionic) * &
             (298.0/(temp+273.))**0.666667

   ! pKh = pKh(0) + bI
   ! b = 0.105 (Butler, 1982)
   pKh = pKh + 0.105*data%ionic

   ! pKw = pKw(0) - f
   pKw = pKw - is_f

   ! pKa1 = pKa1(0) - f - bI
   pKa1 = pKa1 - is_f - 0.105*data%ionic

   !pKa2 = pKa2(0) - 2f
   !pKa2 = pKa2 + 2.0*is_f !YA
   pKa2 = pKa2 - 2.0*is_f

   ! Convert from pK etc to Kh, Kw, Ka1, Ka2
   Kh  = 10.**(-pKh)
   Ka1 = 10.**(-pKa1)
   Ka2 = 10.**(-pKa2)
   Kw  = 10.**(-pKw)

   ! Calculate the speciation to know the molar mass of DIC                                                             !
   H    = 10.**(-pH)
   alf_CO2  = (H*H)/(H*H + Ka1*H + Ka1*Ka2)  !alfa0, this is (H2CO3+CO2)
   alf_HCO3 = (Ka1*H)/(H*H + Ka1*H + Ka1*Ka2) !alfa1
   alf_CO3  = (Ka1*Ka2)/(H*H + Ka1*H + Ka1*Ka2) !alfa2

   ! and update speciation (mol C/L),fractination
   CO2  = dic*alf_CO2
   HCO3 = dic*alf_HCO3
   CO3  = dic*alf_CO3

   ! calculate TA
   TA = dic * (Ka1*H + 2.0*Ka1*Ka2) / (H*H + Ka1*H + Ka1*Ka2) + (Kw/H) - H

   aed2_carbon_co2 = CO2

   !print *, 'End aed2_carbon_co2,pH=', pH, '% diff |meas-calc| alk=', (alkalinity-TA)/alkalinity*100

END FUNCTION aed2_carbon_co2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed2_isotope
