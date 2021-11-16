!###############################################################################
!#                                                                             #
!# aed_organic_matter.F90                                                      #
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
!#  Created June 2011                                                          #
!#  Track changes on GitHub @ https://github.com/AquaticEcoDynamics/libaed-water
!#                                                                             #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ____     | || |    ______    | || | ____    ____ | |         !
!         | |   .'    `.   | || |  .' ___  |   | || ||_   \  /   _|| |         !
!         | |  /  .--.  \  | || | / .'   \_|   | || |  |   \/   |  | |         !
!         | |  | |    | |  | || | | |    ____  | || |  | |\  /| |  | |         !
!         | |  \  `--'  /  | || | \ `.___]  _| | || | _| |_\/_| |_ | |         !
!         | |   `.____.'   | || |  `._____.'   | || ||_____||_____|| |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

MODULE aed_organic_matter
!-------------------------------------------------------------------------------
! aed_organic_matter --- organic matter biogeochemical model
!
! The Organic Matter (OM) module contains equations for mineralisation
! and breakdown of particulate and dissolved organic matter pools, as well as
! settling, resuspension and other processes
!
! Users can optionally configure refractory pools
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,ONLY : water_viscosity

   IMPLICIT NONE

   PRIVATE  ! By default make everything private
!
   PUBLIC aed_organic_matter_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_organic_matter_data_t

      !# Variable identifiers (OGM)
      INTEGER  :: id_poc, id_doc          ! particulate & dissolved organic carbon
      INTEGER  :: id_pon, id_don          ! particulate & dissolved organic nitrogen
      INTEGER  :: id_pop, id_dop          ! particulate & dissolved organic phosphorus
      INTEGER  :: id_docr,id_donr,id_dopr ! dissolved refractory OM pools
      INTEGER  :: id_cpom                 ! coarse particulate OM (bulk) pool
      INTEGER  :: id_sed_pon, id_sed_don  ! sediment ON pool
      INTEGER  :: id_sed_pop, id_sed_dop  ! sediment OP pool
      INTEGER  :: id_sed_poc, id_sed_doc  ! sediment OC pool
      !# Variable identifiers (links)
      INTEGER  :: id_oxy,  id_dic
      INTEGER  :: id_nit,  id_no2,  id_n2o, id_amm
      INTEGER  :: id_frp,  id_fe,   id_ch4, id_so4
      INTEGER  :: id_temp, id_salt, id_vis, id_uva, id_uvb, id_extc, id_rho, id_dz
      INTEGER  :: id_Psed_poc, id_Psed_pon, id_Psed_pop, id_Psed_cpom ! sedimentation rates
      INTEGER  :: id_Fsed_pon, id_Fsed_don ! sed. rate organic nitrogen
      INTEGER  :: id_Fsed_pop, id_Fsed_dop ! sed. rate organic phosphorus
      INTEGER  :: id_Fsed_poc, id_Fsed_doc ! sed. rate organic carbon
      !# Variable identifiers (diags)
      INTEGER  :: id_pon_miner,  id_don_miner
      INTEGER  :: id_pop_miner,  id_dop_miner
      INTEGER  :: id_poc_miner,  id_doc_miner
      INTEGER  :: id_swi_pon,    id_swi_don
      INTEGER  :: id_swi_pop,    id_swi_dop
      INTEGER  :: id_swi_poc,    id_swi_doc
      INTEGER  :: id_docr_miner, id_donr_miner,  id_dopr_miner
      INTEGER  :: id_sed_toc,    id_sed_ton,     id_sed_top
      INTEGER  :: id_l_resus,    id_denit,       id_anaerobic
      INTEGER  :: id_pom_vvel,   id_cpom_vvel
      INTEGER  :: id_bod,        id_cdom
      INTEGER  :: id_photolysis

      !# Model options
      INTEGER  :: resuspension, settling
      INTEGER  :: simDenitrification,simFeReduction,simSO4Reduction,simMethanogenesis,simSedimentOM
      LOGICAL  :: simRPools,simPhotolysis
      LOGICAL  :: use_oxy,use_amm,use_frp,use_dic,use_nit,use_no2,use_n2o,use_fe3,use_ch4,use_so4
      LOGICAL  :: use_Fsed_link_don, use_Fsed_link_dop, use_Fsed_link_doc
      LOGICAL  :: use_Fsed_link_pon, use_Fsed_link_pop, use_Fsed_link_poc
      LOGICAL  :: use_Psed_link_pon, use_Psed_link_pop, use_Psed_link_poc

      !# Model parameters
      AED_REAL :: Rpoc_hydrol,  Rpon_hydrol,  Rpop_hydrol,                     &
                  Rdoc_minerl,  Rdon_minerl,  Rdop_minerl, Rdom_minerl,        &
                  theta_hydrol, theta_minerl, Kpom_hydrol, Kdom_minerl, f_an
      AED_REAL :: Rdomr_minerl, Rcpom_bdown,  X_cpom_n,    X_cpom_p
      AED_REAL :: Ko2_0,Kin_denitrat,Kin_denitrit,Kin_denitrous,Klim_denitrous,Klim_denitrit,Kpart_denitrit
      AED_REAL :: K_nit,K_so4,K_fe
      AED_REAL :: KeDOM, KePOM, KeDOMR, KeCPOM, photo_fmin
      AED_REAL :: w_pom, d_pom, rho_pom, w_cpom, d_cpom, rho_cpom
      AED_REAL :: sedimentOMfrac, Xsc, Xsn, Xsp
      AED_REAL :: Ksed_dom, theta_sed_dom, Fsed_doc, Fsed_don, Fsed_dop

     CONTAINS
         PROCEDURE :: define            => aed_define_organic_matter
         PROCEDURE :: calculate         => aed_calculate_organic_matter
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_organic_matter
         PROCEDURE :: mobility          => aed_mobility_organic_matter
         PROCEDURE :: light_extinction  => aed_light_extinction_organic_matter

   END TYPE

! MODULE GLOBALS
   LOGICAL  :: extra_diag = .false.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

   AED_REAL :: c

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_define_organic_matter(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED module
!
!  Here, the aed namelist is read and the variables required
!  by the module are registered for AED simulation.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                              :: namlst
   CLASS (aed_organic_matter_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER                   :: status
   AED_REAL                  :: initial_sedom

!  %% NAMELIST    %%  /aed_organic_matter/
!  %% Last Checked 20/08/2021
   !-- Default initial values and limits
   AED_REAL                  :: om_min         = zero_
   AED_REAL                  :: om_max         =   1e8
   AED_REAL                  :: poc_initial    = 100.0
   AED_REAL                  :: doc_initial    = 100.0
   AED_REAL                  :: pon_initial    =  16.0
   AED_REAL                  :: don_initial    =  16.0
   AED_REAL                  :: pop_initial    =   1.0
   AED_REAL                  :: dop_initial    =   1.0
   AED_REAL                  :: docr_initial   = zero_
   AED_REAL                  :: donr_initial   = zero_
   AED_REAL                  :: dopr_initial   = zero_
   AED_REAL                  :: cpom_initial   = zero_
   !-- Breakdown and mineralisation (basic pool)
   INTEGER                   :: simDenitrification = 0
   INTEGER                   :: simFeReduction     = 0
   INTEGER                   :: simSO4Reduction    = 0
   INTEGER                   :: simMethanogenesis  = 0
   AED_REAL                  :: Rpom_hydrol    = zero_
   AED_REAL                  :: Rpoc_hydrol    = zero_
   AED_REAL                  :: Rpon_hydrol    = zero_
   AED_REAL                  :: Rpop_hydrol    = zero_
   AED_REAL                  :: Rdom_minerl    = zero_
   AED_REAL                  :: Rdoc_minerl    = zero_
   AED_REAL                  :: Rdon_minerl    = zero_
   AED_REAL                  :: Rdop_minerl    = zero_
   AED_REAL                  :: theta_hydrol   =  1.0
   AED_REAL                  :: theta_minerl   =  1.0
   AED_REAL                  :: Kpom_hydrol    = 20.0
   AED_REAL                  :: Kdom_minerl    = 20.0
   AED_REAL                  :: f_an           =  1.0
   AED_REAL                  :: Kin_denitrat   = 20.0
   AED_REAL                  :: Kin_denitrit   =  0.297
   AED_REAL                  :: Kin_denitrous  =  0.205
   AED_REAL                  :: K_nit          =  1.
   AED_REAL                  :: Klim_denitrous =  1.
   AED_REAL                  :: Klim_denitrit  =  1.
   AED_REAL                  :: Kpart_denitrit =  1.
   CHARACTER(len=64)         :: denit_link     = 'NIT_denit'
   CHARACTER(len=64)         :: dom_miner_oxy_reactant_var=''
   CHARACTER(len=64)         :: dom_miner_nit_reactant_var=''
   CHARACTER(len=64)         :: dom_miner_no2_reactant_var=''
   CHARACTER(len=64)         :: dom_miner_n2o_reactant_var=''
   CHARACTER(len=64)         :: dom_miner_fe3_reactant_var=''
   CHARACTER(len=64)         :: dom_miner_so4_reactant_var=''
   CHARACTER(len=64)         :: dom_miner_ch4_reactant_var=''
   CHARACTER(len=64)         :: doc_miner_product_variable=''
   CHARACTER(len=64)         :: don_miner_product_variable=''
   CHARACTER(len=64)         :: dop_miner_product_variable=''

   !-- Refractory organic matter (optional)
   LOGICAL                   :: simRPools     = .false.
   AED_REAL                  :: Rdomr_minerl  =   0.001
   AED_REAL                  :: Rcpom_bdown   =   0.001
   AED_REAL                  :: X_cpom_n      =  16./106.
   AED_REAL                  :: X_cpom_p      =   1./106.
   CHARACTER(len=64)         :: docr_miner_product_variable='OGM_doc'
   CHARACTER(len=64)         :: donr_miner_product_variable='OGM_don'
   CHARACTER(len=64)         :: dopr_miner_product_variable='OGM_dop'

   !-- Light related parameters
   AED_REAL                  :: KeDOM         = 0.0001
   AED_REAL                  :: KePOM         = 0.0001
   AED_REAL                  :: KeDOMR        = 0.0001
   AED_REAL                  :: KeCPOM        = 0.0001
   AED_REAL                  :: photo_fmin    = 0.9
   AED_REAL                  :: photo_c       = 7.52
   LOGICAL                   :: simPhotolysis = .FALSE.

   !-- Particle settling parameters
   INTEGER                   :: settling      = 0
   AED_REAL                  :: w_pom         = zero_
   AED_REAL                  :: d_pom         = zero_
   AED_REAL                  :: rho_pom       = zero_
   AED_REAL                  :: w_cpom        = zero_
   AED_REAL                  :: d_cpom        = zero_
   AED_REAL                  :: rho_cpom      = zero_
   CHARACTER(len=64)         :: Psed_poc_variable=''
   CHARACTER(len=64)         :: Psed_pon_variable=''
   CHARACTER(len=64)         :: Psed_pop_variable=''
   CHARACTER(len=64)         :: Psed_cpom_variable=''

   !-- Sediment interaction parameters (basic model)
   INTEGER                   :: resuspension     = 0
   INTEGER                   :: simSedimentOM    = 0
   AED_REAL                  :: sedimentOMfrac   = zero_
   AED_REAL                  :: sedimentBulkDens = 1.3e3
   AED_REAL                  :: sedimentDepth    = 0.5
   AED_REAL                  :: Xsc              = zero_
   AED_REAL                  :: Xsn              = zero_
   AED_REAL                  :: Xsp              = zero_
   AED_REAL                  :: Fsed_doc         = zero_
   AED_REAL                  :: Fsed_don         = zero_
   AED_REAL                  :: Fsed_dop         = zero_
   AED_REAL                  :: Ksed_dom         = zero_
   AED_REAL                  :: theta_sed_dom    = 1.0
   CHARACTER(len=64)         :: resus_link       =''
   CHARACTER(len=64)         :: Fsed_poc_variable=''
   CHARACTER(len=64)         :: Fsed_doc_variable=''
   CHARACTER(len=64)         :: Fsed_pon_variable=''
   CHARACTER(len=64)         :: Fsed_don_variable=''
   CHARACTER(len=64)         :: Fsed_pop_variable=''
   CHARACTER(len=64)         :: Fsed_dop_variable=''

   !-- From Module Globals
!  LOGICAL  :: extra_diag = .FALSE.      !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST    %%  /aed_organic_matter/

   !-- Parameters to read in
   NAMELIST /aed_organic_matter/   om_min, om_max,                             &
                      poc_initial,pon_initial,pop_initial,                     &
                      doc_initial,don_initial,dop_initial,                     &
                      docr_initial, donr_initial, dopr_initial, cpom_initial,  &
                      Rpoc_hydrol,Rpon_hydrol,Rpop_hydrol,Rpom_hydrol,         &
                      Rdoc_minerl,Rdon_minerl,Rdop_minerl,Rdom_minerl,         &
                      theta_hydrol,theta_minerl,Kpom_hydrol,Kdom_minerl,f_an,  &
                      dom_miner_oxy_reactant_var, dom_miner_nit_reactant_var,  &
                      dom_miner_no2_reactant_var, dom_miner_n2o_reactant_var,  &
                      dom_miner_fe3_reactant_var, dom_miner_so4_reactant_var,  &
                      dom_miner_ch4_reactant_var, doc_miner_product_variable,  &
                      don_miner_product_variable, dop_miner_product_variable,  &
                      simRPools,Rdomr_minerl,                                  &
                      Rcpom_bdown,X_cpom_n,X_cpom_p,                           &
                      Kin_denitrat,Kin_denitrit,Kin_denitrous,K_nit,           &
                      Klim_denitrous,Klim_denitrit,Kpart_denitrit,             &
                      KeDOM, KePOM, KeDOMR, KeCPOM,                            &
                      simPhotolysis, photo_fmin, photo_c,                      &
                      settling,w_pom,d_pom,rho_pom,w_cpom,d_cpom,rho_cpom,     &
                      Psed_poc_variable, Psed_pon_variable,                    &
                      Psed_pop_variable,Psed_cpom_variable,                    &
                      resuspension, resus_link,                                &
                      sedimentOMfrac,sedimentBulkDens,sedimentDepth,           &
                      Xsc,Xsn,Xsp,                                             &
                      Fsed_doc,Fsed_don,Fsed_dop,Ksed_dom,theta_sed_dom,       &
                      Fsed_poc_variable, Fsed_doc_variable,                    &
                      Fsed_pop_variable, Fsed_dop_variable,                    &
                      Fsed_pon_variable, Fsed_don_variable,                    &
                      simDenitrification, simFeReduction,                      &
                      simSO4Reduction, simMethanogenesis, simSedimentOM,       &
                      extra_diag, diag_level

!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_organic_matter initialization"

   ! Default parameters set above at declaration

   ! Read the namelist
   read(namlst,nml=aed_organic_matter,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_organic_matter'

   IF( extra_diag )   diag_level = 10           ! legacy use of extra_diag

   ! Store parameter values in the modules data structure
   !   Note: all rates are provided via input as values per day,
   !       and are converted here to values per second.

   !-- Breakdown and mineralisation (basic pool)
   data%Rpoc_hydrol         = Rpoc_hydrol/secs_per_day
   data%Rpon_hydrol         = Rpon_hydrol/secs_per_day
   data%Rpop_hydrol         = Rpop_hydrol/secs_per_day
   data%Rdoc_minerl         = Rdoc_minerl/secs_per_day
   data%Rdon_minerl         = Rdon_minerl/secs_per_day
   data%Rdop_minerl         = Rdop_minerl/secs_per_day
   data%Rdom_minerl         = Rdom_minerl/secs_per_day
   data%theta_hydrol        = theta_hydrol
   data%theta_minerl        = theta_minerl
   data%Kpom_hydrol         = Kpom_hydrol
   data%Kdom_minerl         = Kdom_minerl
   data%f_an                = f_an
   data%simDenitrification  = simDenitrification
   data%simFeReduction      = simFeReduction
   data%simSO4Reduction     = simSO4Reduction
   data%simMethanogenesis   = simMethanogenesis
   data%K_nit               = K_nit
   data%Kin_denitrat        = Kin_denitrat
   data%Kin_denitrit        = Kin_denitrit
   data%Kin_denitrous       = Kin_denitrous
   data%Klim_denitrit       = Klim_denitrit
   data%Klim_denitrous      = Klim_denitrous
   data%Kpart_denitrit      = Kpart_denitrit
   !-- Refractory organic matter (optional)
   data%simRPools           = simRPools
   data%Rdomr_minerl        = Rdomr_minerl/secs_per_day
   data%Rcpom_bdown         = Rcpom_bdown/secs_per_day
   data%X_cpom_n            = X_cpom_n
   data%X_cpom_p            = X_cpom_p
   !-- Light related parameters
   data%KeDOM               = KeDOM
   data%KePOM               = KePOM
   data%KeDOMR              = KeDOMR
   data%KeCPOM              = KeCPOM
   data%simPhotolysis       = simPhotolysis
   data%photo_fmin          = photo_fmin
   c = photo_c
   !-- Particle settling parameters
   data%settling            = settling
   data%w_pom               = w_pom/secs_per_day
   data%d_pom               = d_pom
   data%rho_pom             = rho_pom
   data%w_cpom              = w_cpom/secs_per_day
   data%d_cpom              = d_cpom
   data%rho_cpom            = rho_cpom
   IF( settling==0 ) THEN
     w_pom=zero_; w_cpom=zero_
   ENDIF
   !-- Sediment interaction parameters (basic model)
   data%resuspension        = resuspension
   data%simSedimentOM       = simSedimentOM
   data%sedimentOMfrac      = sedimentOMfrac
   data%Xsc                 = Xsc
   data%Xsn                 = Xsn
   data%Xsp                 = Xsp
   data%Fsed_doc            = Fsed_doc/secs_per_day
   data%Fsed_don            = Fsed_don/secs_per_day
   data%Fsed_dop            = Fsed_dop/secs_per_day
   data%Ksed_dom            = Ksed_dom
   data%theta_sed_dom       = theta_sed_dom

   ! Register state variables

   data%id_doc = aed_define_variable('doc','mmol/m**3','dissolved organic carbon', &
                                    doc_initial,minimum=om_min,maximum=om_max)
   data%id_poc = aed_define_variable('poc','mmol/m**3','particulate organic carbon', &
                                    poc_initial,minimum=om_min,maximum=om_max,mobility=data%w_pom)
   data%id_don = aed_define_variable('don','mmol/m**3','dissolved organic nitrogen', &
                                    don_initial,minimum=om_min,maximum=om_max)
   data%id_pon = aed_define_variable('pon','mmol/m**3','particulate organic nitrogen', &
                                    pon_initial,minimum=om_min,maximum=om_max,mobility=data%w_pom)
   data%id_dop = aed_define_variable('dop','mmol/m**3','dissolved organic phosphorus', &
                                    dop_initial,minimum=om_min,maximum=om_max)
   data%id_pop = aed_define_variable('pop','mmol/m**3','particulate organic phosphorus', &
                                    pop_initial,minimum=om_min,maximum=om_max,mobility=data%w_pom)

   data%id_docr = 0 ; data%id_donr = 0
   data%id_dopr = 0 ; data%id_cpom = 0
   IF (simRPools) THEN
     data%id_docr = aed_define_variable('docr','mmol/m**3',&
                                         'refractory dissolved organic carbon',&
                                         docr_initial,minimum=zero_)
     data%id_donr = aed_define_variable('donr','mmol/m**3',&
                                         'refractory dissolved organic nitrogen',&
                                         donr_initial,minimum=zero_)
     data%id_dopr = aed_define_variable('dopr','mmol/m**3',&
                                         'refractory dissolved organic phosphorus',&
                                         dopr_initial,minimum=zero_)
     data%id_cpom = aed_define_variable('cpom','mmol/m**3',&
                                         'coarse particulate matter',&
                                         cpom_initial,minimum=zero_,mobility=data%w_cpom)
   ENDIF

   IF (simSedimentOM>0) THEN
     ! Initialise OM pools here.
     ! simSedimentOM = 1 - initialise to zero and look at relative evolution
     ! simSedimentOM = 2 - initialise mass based %OM and bulk density (assume 50:50 POM:DOM)
     initial_sedom = zero_
     IF (simSedimentOM>1) &
       initial_sedom = sedimentBulkDens * sedimentDepth * sedimentOMfrac * Xsc * 0.5 / (12.0/1e3)
     data%id_sed_doc = aed_define_sheet_variable( 'sed_doc',                   &
                                                  'mmol/m**2',                 &
                                                  'sediment doc mass',         &
                                                  initial_sedom,               &
                                                  minimum=-1e10)
     IF (simSedimentOM>1) &
       initial_sedom = sedimentBulkDens * sedimentDepth * sedimentOMfrac * Xsn * 0.5 / (14.0/1e3)
     data%id_sed_don = aed_define_sheet_variable( 'sed_don',                   &
                                                  'mmol/m**2',                 &
                                                  'sediment don mass',         &
                                                  initial_sedom,               &
                                                  minimum=-1e10)
     IF (simSedimentOM>1) &
       initial_sedom = sedimentBulkDens * sedimentDepth * sedimentOMfrac * Xsp * 0.5 / (30.9/1e3)
     data%id_sed_dop = aed_define_sheet_variable( 'sed_dop',                   &
                                                  'mmol/m**2',                 &
                                                  'sediment dop mass',         &
                                                  initial_sedom,               &
                                                  minimum=-1e10)
     IF (simSedimentOM>1) &
       initial_sedom = sedimentBulkDens * sedimentDepth * sedimentOMfrac * Xsc * 0.5 / (12.0/1e3)
     data%id_sed_poc = aed_define_sheet_variable( 'sed_poc',                   &
                                                  'mmol/m**2',                 &
                                                  'sediment poc mass',         &
                                                  initial_sedom,               &
                                                  minimum=-1e10)
     IF (simSedimentOM>1) &
       initial_sedom = sedimentBulkDens * sedimentDepth * sedimentOMfrac * Xsn * 0.5 / (14.0/1e3)
     data%id_sed_pon = aed_define_sheet_variable( 'sed_pon',                   &
                                                  'mmol/m**2',                 &
                                                  'sediment pon mass',         &
                                                  initial_sedom,               &
                                                  minimum=-1e10)
     IF (simSedimentOM>1) &
       initial_sedom = sedimentBulkDens * sedimentDepth * sedimentOMfrac * Xsp * 0.5 / (30.9/1e3)
     data%id_sed_pop = aed_define_sheet_variable( 'sed_pop',                   &
                                                  'mmol/m**2',                 &
                                                  'sediment pop mass',         &
                                                  initial_sedom,               &
                                                  minimum=-1e10)
   ENDIF

   ! Register external (linked) state variable dependencies
   !-- oxygen & oxidants
   data%use_nit = .FALSE. ; data%use_no2 = .FALSE. ; data%use_n2o = .FALSE.
   data%use_oxy = dom_miner_oxy_reactant_var .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) &
     data%id_oxy = aed_locate_variable(dom_miner_oxy_reactant_var)
   IF( simDenitrification==1 ) THEN
     data%use_nit = dom_miner_nit_reactant_var .NE. ''
     IF (data%use_nit) THEN
       data%id_nit = aed_locate_variable(dom_miner_nit_reactant_var)
     ENDIF
   ELSEIF( simDenitrification==2 ) THEN
       data%use_nit = dom_miner_nit_reactant_var .NE. ''
       IF (data%use_nit) data%id_nit = aed_locate_variable(dom_miner_nit_reactant_var)
       data%use_no2 = dom_miner_no2_reactant_var .NE. ''
       IF (data%use_no2) data%id_no2 = aed_locate_variable(dom_miner_no2_reactant_var)
       data%use_n2o = dom_miner_n2o_reactant_var .NE. ''
       IF (data%use_n2o) data%id_n2o = aed_locate_variable(dom_miner_n2o_reactant_var)
   ENDIF
   IF( simFeReduction==1 ) THEN
     data%use_fe3 = dom_miner_fe3_reactant_var .NE. ''
     IF (data%use_fe3) THEN
       data%id_fe = aed_locate_variable(dom_miner_fe3_reactant_var)
     ELSE
       simFeReduction = 0
     ENDIF
   ELSE
     data%use_fe3 = .FALSE.
   ENDIF
   IF( simSO4Reduction==1 ) THEN
     data%use_so4 = dom_miner_so4_reactant_var .NE. ''
     IF (data%use_so4) THEN
       data%id_so4 = aed_locate_variable(dom_miner_so4_reactant_var)
     ELSE
       simSO4Reduction = 0
     ENDIF
   ELSE
     data%use_so4 = .FALSE.
   ENDIF
   IF( simMethanogenesis==1 ) THEN
     data%use_ch4 = dom_miner_ch4_reactant_var .NE. ''
     IF (data%use_ch4) THEN
       data%id_ch4 = aed_locate_variable(dom_miner_ch4_reactant_var)
     ELSE
       simMethanogenesis = 0
     ENDIF
   ELSE
     data%use_ch4 = .FALSE.
   ENDIF

   !-- carbon
   data%use_dic = doc_miner_product_variable .NE. '' !This means carbon module switched on
   IF (data%use_dic) &
     data%id_dic = aed_locate_variable(doc_miner_product_variable)

   !-- nitrogen
   data%use_amm = don_miner_product_variable .NE. '' !This means nitrogen module switched on
   IF (data%use_amm) &
     data%id_amm = aed_locate_variable(don_miner_product_variable)

   !-- phosphorus
   data%use_frp = dop_miner_product_variable .NE. '' !This means phosphorus module switched on
   IF (data%use_frp) &
     data%id_frp = aed_locate_variable(dop_miner_product_variable)

   !-- sediment flux link variables
   data%id_Fsed_pon = -1 ; data%id_Fsed_pop = -1 ; data%id_Fsed_poc = -1
   data%id_Fsed_don = -1 ; data%id_Fsed_dop = -1 ; data%id_Fsed_doc = -1

   data%use_Fsed_link_doc = Fsed_doc_variable .NE. ''
   IF (data%use_Fsed_link_doc) data%id_Fsed_doc = aed_locate_sheet_variable(Fsed_doc_variable)
   data%use_Fsed_link_don = Fsed_don_variable .NE. ''
   IF (data%use_Fsed_link_don) data%id_Fsed_don = aed_locate_sheet_variable(Fsed_don_variable)
   data%use_Fsed_link_dop = Fsed_dop_variable .NE. ''
   IF (data%use_Fsed_link_dop) data%id_Fsed_dop = aed_locate_sheet_variable(Fsed_dop_variable)
   data%use_Fsed_link_poc = Fsed_poc_variable .NE. ''
   IF (data%use_Fsed_link_poc) data%id_Fsed_poc = aed_locate_sheet_variable(Fsed_poc_variable)
   data%use_Fsed_link_pon = Fsed_pon_variable .NE. ''
   IF (data%use_Fsed_link_pon) data%id_Fsed_pon = aed_locate_sheet_variable(Fsed_pon_variable)
   data%use_Fsed_link_pop = Fsed_pop_variable .NE. ''
   IF (data%use_Fsed_link_pop) data%id_Fsed_pop = aed_locate_sheet_variable(Fsed_pop_variable)

   !-- sedimentation link variables
   data%id_Psed_poc = -1 ; data%id_Psed_pon = -1 ; data%id_Psed_pop = -1

   data%use_Psed_link_poc = Psed_poc_variable .NE. ''
   IF (data%use_Psed_link_poc) THEN
     data%id_Psed_poc = aed_locate_sheet_variable(Psed_poc_variable)
   ELSE
     data%id_Psed_poc = aed_define_diag_variable('Psed_poc','mmol/m**2/s',  'POC sedimentation')
     IF (simRPools) &
       data%id_Psed_cpom = aed_define_diag_variable('Psed_cpom','mmol/m**2/s',  'CPOM sedimentation')
   ENDIF
   data%use_Psed_link_pon = Psed_pon_variable .NE. ''
   IF (data%use_Psed_link_pon) THEN
     data%id_Psed_pon = aed_locate_sheet_variable(Psed_pon_variable)
   ELSE
     data%id_Psed_pon = aed_define_diag_variable('Psed_pon','mmol/m**2/s',  'PON sedimentation')
   ENDIF
   data%use_Psed_link_pop = Psed_pop_variable .NE. ''
   IF (data%use_Psed_link_pop) THEN
     data%id_Psed_pop = aed_locate_sheet_variable(Psed_pop_variable)
   ELSE
     data%id_Psed_pop = aed_define_diag_variable('Psed_pop','mmol/m**2/s',  'POP sedimentation')
   ENDIF


   !-- resuspension link variable
   IF ( resuspension>0 .AND. .NOT.resus_link .EQ. '' ) THEN
      data%id_l_resus  = aed_locate_sheet_variable(TRIM(resus_link)) ! ('TRC_resus')
   ELSE
      data%id_l_resus = 0
      data%resuspension = 0
   ENDIF

   !-- light
   IF (simRPools) THEN
     data%id_vis = aed_locate_global('par')
     data%id_uva = aed_locate_global('uva')
     data%id_uvb = aed_locate_global('uvb')
     data%id_extc= aed_locate_global('extc_coef')
     data%id_dz  = aed_locate_global('layer_ht')
   ENDIF

   ! Register diagnostic variables
   data%id_swi_poc    = 0
   data%id_swi_doc    = 0
   data%id_swi_pon    = 0
   data%id_swi_don    = 0
   data%id_swi_pop    = 0
   data%id_swi_dop    = 0
   data%id_poc_miner  = 0
   data%id_doc_miner  = 0
   data%id_pon_miner  = 0
   data%id_don_miner  = 0
   data%id_pop_miner  = 0
   data%id_dop_miner  = 0
   data%id_denit      = 0
   data%id_bod        = 0
   data%id_photolysis = 0
   data%id_docr_miner = 0
   data%id_donr_miner = 0
   data%id_dopr_miner = 0

    IF ( diag_level>0 ) THEN
      data%id_cdom = aed_define_diag_variable('CDOM','/m',  'Chromophoric DOM (CDOM)')
      data%id_sed_toc  = aed_define_sheet_diag_variable('sed_toc','mmol/m**2',  'Sediment TOC mass')
      data%id_sed_ton  = aed_define_sheet_diag_variable('sed_ton','mmol/m**2',  'Sediment TON mass')
      data%id_sed_top  = aed_define_sheet_diag_variable('sed_top','mmol/m**2',  'Sediment TOP mass')
    ENDIF

   IF ( diag_level>1 ) THEN
     data%id_swi_poc  = aed_define_sheet_diag_variable('swi_poc','mmol/m**2/d',  'Net POC flux @ the SWI')
     data%id_swi_doc  = aed_define_sheet_diag_variable('swi_doc','mmol/m**2/d',  'Net DOC flux @ the SWI')
     data%id_swi_pon  = aed_define_sheet_diag_variable('swi_pon','mmol/m**2/d',  'Net PON flux @ the SWI')
     data%id_swi_don  = aed_define_sheet_diag_variable('swi_don','mmol/m**2/d',  'Net DON flux @ the SWI')
     data%id_swi_pop  = aed_define_sheet_diag_variable('swi_pop','mmol/m**2/d',  'Net POP flux @ the SWI')
     data%id_swi_dop  = aed_define_sheet_diag_variable('swi_dop','mmol/m**2/d',  'Net DOP flux @ the SWI')
     data%id_poc_miner= aed_define_diag_variable('poc_hydrol','mmol/m**3/d','POC hydrolosis')
     data%id_doc_miner= aed_define_diag_variable('doc_miner' ,'mmol/m**3/d','DOC mineralisation')
     data%id_pon_miner= aed_define_diag_variable('pon_hydrol','mmol/m**3/d','PON hydrolosis')
     data%id_don_miner= aed_define_diag_variable('don_miner' ,'mmol/m**3/d','DON mineralisation')
     data%id_pop_miner= aed_define_diag_variable('pop_hydrol','mmol/m**3/d','POP hydrolysis')
     data%id_dop_miner= aed_define_diag_variable('dop_miner' ,'mmol/m**3/d','DOP mineralisation')
     data%id_anaerobic= aed_define_diag_variable('anaerobic' ,'mmol/m**3/d','anaerobic metabolism')
     data%id_denit    = aed_define_diag_variable('denit'     ,'mmol/m**3/d','denitrification')

     IF (simRPools) THEN
        data%id_docr_miner = aed_define_diag_variable('docr_miner','mmol/m**3/d','DOCR mineralisation')
        data%id_donr_miner = aed_define_diag_variable('donr_miner','mmol/m**3/d','DONR mineralisation')
        data%id_dopr_miner = aed_define_diag_variable('dopr_miner','mmol/m**3/d','DOPR mineralisation')
     ENDIF

     data%id_bod = aed_define_diag_variable('BOD5','mg O2/L',  'Biochemical Oxygen Demand (BOD)')

     IF ( simPhotolysis .and. simRpools  ) data%id_photolysis = &
       aed_define_diag_variable('photolysis','mmol C/m3/d',  'photolysis rate of breakdown of DOC')

     data%id_pom_vvel  = aed_define_diag_variable('pom_vvel','m/d','POM vertical velocity')
     data%id_cpom_vvel = aed_define_diag_variable('cpom_vvel','m/d','CPOM vertical velocity')

   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed_locate_global('temperature')
   data%id_salt = aed_locate_global('salinity')
   data%id_rho  = aed_locate_global('density')

END SUBROUTINE aed_define_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_organic_matter(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_organic_matter model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: poc,pon,pop
   AED_REAL :: doc,don,dop
   AED_REAL :: docr,donr,dopr,cpom,cdom
   AED_REAL :: temp,dic,amm,oxy,nit,no2,n2o,frp,feiii,so4,ch4
   AED_REAL :: pom_hydrolysis, dom_mineralisation
   AED_REAL :: pon_hydrolysis, don_mineralisation
   AED_REAL :: pop_hydrolysis, dop_mineralisation
   AED_REAL :: poc_hydrolysis, doc_mineralisation
   AED_REAL :: docr_mineralisation = 0., donr_mineralisation = 0.
   AED_REAL :: dopr_mineralisation = 0., cpom_breakdown = 0.
   AED_REAL :: denitrification, denitratation = 0., denitritation = 0.
   AED_REAL :: denitrousation = 0., dnra, nitrous_denitritation, ammonium_release
   AED_REAL :: photolysis, vis, uva, uvb, photo_fmin, cdoc, att
   AED_REAL :: doc_min_aerobic, doc_min_anaerobic
   AED_REAL :: fereduction = 0., so4reduction = 0., methanogenesis = 0., fso4, fnit, ffe

!-----------------------------------------------------------------------
!BEGIN
   !CALL log_message('model aed_organic_matter enter do loop successfully.')

   !---------------------------------------------------------------------------+
   !# Setup
   !---------------------------------------------------------------------------+

   !# Initialisations
   nit = zero_
   docr = zero_
   photolysis = zero_
   denitrification = zero_
   photo_fmin = data%photo_fmin

   !# Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp)  ! temperature

   !# Retrieve current (local) state variable values.
   poc = _STATE_VAR_(data%id_poc)    ! particulate organic carbon
   pon = _STATE_VAR_(data%id_pon)    ! particulate organic nitrogen
   pop = _STATE_VAR_(data%id_pop)    ! particulate organic phosphorus
   doc = _STATE_VAR_(data%id_doc)    ! dissolved organic carbon
   don = _STATE_VAR_(data%id_don)    ! dissolved organic nitrogen
   dop = _STATE_VAR_(data%id_dop)    ! dissolved organic phosphorus

   IF(data%simRPools) THEN
      docr = _STATE_VAR_(data%id_docr)
      donr = _STATE_VAR_(data%id_donr)
      dopr = _STATE_VAR_(data%id_dopr)
      cpom = _STATE_VAR_(data%id_cpom)

      att = EXP( -_STATE_VAR_(data%id_extc) * _STATE_VAR_(data%id_dz)/2. )
      vis = _STATE_VAR_(data%id_vis) * att
      uva = _STATE_VAR_(data%id_uva) * att
      uvb = _STATE_VAR_(data%id_uvb) * att
      cdom = _DIAG_VAR_(data%id_cdom)
   ENDIF

!  IF (data%id_n_den>0) denitrification = _DIAG_VAR_(data%id_n_den) / secs_per_day
   IF (data%id_denit>0) denitrification = _DIAG_VAR_(data%id_denit) / secs_per_day

   IF (data%use_oxy) THEN
      oxy = _STATE_VAR_(data%id_oxy) ! oxygen
   ELSE
      oxy = 300.0
   ENDIF
   IF (data%use_nit) THEN
      nit = _STATE_VAR_(data%id_nit) ! nitrate
   ELSE
      nit = 300.0
   ENDIF
   IF (data%use_no2) THEN
      no2 = _STATE_VAR_(data%id_no2) ! nitrite
   ELSE
      no2 = 300.0
   ENDIF
   IF (data%use_n2o) THEN
      n2o = _STATE_VAR_(data%id_n2o) ! nitrous oxide
   ELSE
      n2o = 300.0
   ENDIF
   IF (data%use_dic) THEN
      dic = _STATE_VAR_(data%id_dic) ! dissolved inorganic carbon
   ENDIF
   IF (data%use_amm) THEN
      amm = _STATE_VAR_(data%id_amm) ! ammonium
   ENDIF
   IF (data%use_frp) THEN
      frp = _STATE_VAR_(data%id_frp) ! phosphate
   ENDIF
   IF (data%use_fe3) THEN
      feiii = _STATE_VAR_(data%id_fe)
   ELSE
      feiii = 0.
   ENDIF
   IF (data%use_so4) THEN
      so4 = _STATE_VAR_(data%id_so4)
   ELSE
      so4 = 100.
   ENDIF
   IF (data%use_ch4) THEN
      ch4 = _STATE_VAR_(data%id_ch4)
   ELSE
      ch4 = 0.
   ENDIF


   !---------------------------------------------------------------------------+
   !# Reaction rates
   !---------------------------------------------------------------------------+

   ! Compute hydroloysis and mineralisation (units mmol N/m3/day)
   poc_hydrolysis = poc*fpoc_miner(data%use_oxy,data%Rpoc_hydrol,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)
   pon_hydrolysis = pon*fpon_miner(data%use_oxy,data%Rpon_hydrol,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)
   pop_hydrolysis = pop*fpop_miner(data%use_oxy,data%Rpop_hydrol,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)

   dom_mineralisation = R_dom(data%Rdom_minerl,data%Kdom_minerl,data%f_an,data%theta_minerl,oxy,temp)
   doc_mineralisation = doc*dom_mineralisation
   don_mineralisation = don*dom_mineralisation
   dop_mineralisation = dop*dom_mineralisation

   !doc_mineralisation = fdoc_miner(data%use_oxy,data%Rdoc_minerl,data%Kdom_minerl,data%theta_minerl,oxy,temp)
   !don_mineralisation = fdon_miner(data%use_oxy,data%Rdon_minerl,data%Kdom_minerl,data%theta_minerl,oxy,temp)
   !dop_mineralisation = fdop_miner(data%use_oxy,data%Rdop_minerl,data%Kdom_minerl,data%theta_minerl,oxy,temp)
   IF (data%simRPools) THEN
      cpom_breakdown = cpom*fpoc_miner(data%use_oxy,data%Rcpom_bdown,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)

      dom_mineralisation  = R_dom(data%Rdomr_minerl,data%Kdom_minerl,data%f_an,data%theta_minerl,oxy,temp)
      docr_mineralisation = docr*dom_mineralisation
      donr_mineralisation = donr*dom_mineralisation
      dopr_mineralisation = dopr*dom_mineralisation

      !docr_mineralisation = fdoc_miner(data%use_oxy,data%Rdocr_miner,data%Kdom_minerl,data%theta_minerl,oxy,temp)
      !donr_mineralisation = fdoc_miner(data%use_oxy,data%Rdonr_miner,data%Kdom_minerl,data%theta_minerl,oxy,temp)
      !dopr_mineralisation = fdoc_miner(data%use_oxy,data%Rdopr_miner,data%Kdom_minerl,data%theta_minerl,oxy,temp)
      IF ( data%simPhotolysis ) THEN
         photolysis = photo(vis,cdom,1) + photo(uva,cdom,2) + photo(uvb,cdom,3)
         !# Limit photolysis to 90% of doc pool within 1 hour
         IF(photolysis > 0.9*docr/3.6e3) photolysis = 0.9*docr/3.6e3
      ENDIF
   ENDIF


   !---------------------------------------------------------------------------+
   !# Flux updates for the ODE
   !---------------------------------------------------------------------------+

   !---------------------------------------------------------------------------+
   ! Set temporal derivatives : POM breakdown
   _FLUX_VAR_(data%id_poc)  = _FLUX_VAR_(data%id_poc) - poc_hydrolysis
   _FLUX_VAR_(data%id_pon)  = _FLUX_VAR_(data%id_pon) - pon_hydrolysis
   _FLUX_VAR_(data%id_pop)  = _FLUX_VAR_(data%id_pop) - pop_hydrolysis

   !---------------------------------------------------------------------------+
   ! Set temporal derivatives : DOM metabolism
   _FLUX_VAR_(data%id_doc)  = _FLUX_VAR_(data%id_doc) + poc_hydrolysis - doc_mineralisation
   _FLUX_VAR_(data%id_don)  = _FLUX_VAR_(data%id_don) + pon_hydrolysis - don_mineralisation
   _FLUX_VAR_(data%id_dop)  = _FLUX_VAR_(data%id_dop) + pop_hydrolysis - dop_mineralisation

   !---------------------------------------------------------------------------+
   ! Set temporal derivatives : Refractory Pools
   IF (data%simRPools) THEN
      docr = MAX(1.e-3,docr)
      _FLUX_VAR_(data%id_docr) = _FLUX_VAR_(data%id_docr) - docr_mineralisation - photolysis
      _FLUX_VAR_(data%id_doc)  = _FLUX_VAR_(data%id_doc)  +  &
                                    docr_mineralisation + photolysis*photo_fmin
      _FLUX_VAR_(data%id_donr) = _FLUX_VAR_(data%id_donr) -  &
                                    donr_mineralisation - photolysis*(donr/docr)
      _FLUX_VAR_(data%id_don)  = _FLUX_VAR_(data%id_don)  +  &
                                    donr_mineralisation + photolysis*photo_fmin*(donr/docr)
      _FLUX_VAR_(data%id_dopr) = _FLUX_VAR_(data%id_dopr) -  &
                                    dopr_mineralisation - photolysis*(dopr/docr)
      _FLUX_VAR_(data%id_dop)  = _FLUX_VAR_(data%id_dop)  +  &
                                    dopr_mineralisation + photolysis*photo_fmin*(dopr/docr)
      _FLUX_VAR_(data%id_cpom) = _FLUX_VAR_(data%id_cpom) - cpom_breakdown
      _FLUX_VAR_(data%id_poc)  = _FLUX_VAR_(data%id_poc) + (cpom_breakdown)
      _FLUX_VAR_(data%id_pon)  = _FLUX_VAR_(data%id_pon) + (cpom_breakdown)*data%X_cpom_n
      _FLUX_VAR_(data%id_pop)  = _FLUX_VAR_(data%id_pop) + (cpom_breakdown)*data%X_cpom_p
   ENDIF

   !---------------------------------------------------------------------------+
   !# Set temporal derivatives : OXIDANTS
   doc_min_aerobic   = doc_mineralisation * (oxy/(data%Kdom_minerl+oxy))
   doc_min_anaerobic = doc_mineralisation - doc_min_aerobic

   fnit = one_; ffe = one_; fso4 = one_

   IF (data%use_oxy) THEN
     !# Aerobic mineralisation : O2 consumption
     _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - doc_min_aerobic
   ENDIF
   IF( data%simDenitrification==1 ) THEN
     !# De-nitrification : simple denitrification, using NO3
     denitrification = doc_min_anaerobic * nit/(data%K_nit+nit)
     _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) - denitrification
     fnit = data%K_nit/(data%K_nit+nit)

   ELSEIF( data%simDenitrification==2 ) THEN
     !# De-nitrification : complex model involving NO2, NO3, N2O
     !   We 1st partition metabolism across the pathways to N2:
     !   Denitratation (1st step in denitrification: NO3->NO2)
     !   Denitritation (2nd step in denitrification: NO2->N2O&NH4)
     !   Denitrousation (Last step in denitrification: N2O->N2)
     denitratation    = (doc_mineralisation) * nit/(data%K_nit+nit) * data%Kin_denitrat/(data%Kin_denitrat+oxy) * &
                        data%Klim_denitrous/(data%Klim_denitrous +n2o) * data%Klim_denitrit/(data%Klim_denitrit +no2)
     denitritation    = (doc_mineralisation) * no2/(data%Klim_denitrit + no2) * data%Kin_denitrit/(data%Kin_denitrit+oxy) * &
                        data%Klim_denitrous/(data%Klim_denitrous + n2o)
     denitrousation   = (doc_mineralisation) * n2o/(data%Klim_denitrous+n2o) * data%Kin_denitrous/(data%Kin_denitrous+oxy)

     !# Check the denitrified DOM is not more than the total anaerobic allowed
     denitrification  = denitratation+denitritation+denitrousation
     IF ( denitrification > doc_min_anaerobic ) THEN
       denitratation  = denitratation * doc_min_anaerobic/denitrification
       denitritation  = denitritation * doc_min_anaerobic/denitrification
       denitrousation = denitrousation * doc_min_anaerobic/denitrification
     ENDIF
     dnra                  =  denitritation * data%Kpart_denitrit/(data%Kpart_denitrit+no2)
     nitrous_denitritation =  denitritation * no2/(data%Kpart_denitrit+no2) * 0.5
     ammonium_release      = (denitratation + denitritation + denitrousation) * MIN(don/MAX(doc,1e-2),one_) ! data%X_nc

     _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) - denitratation
     _FLUX_VAR_(data%id_no2) = _FLUX_VAR_(data%id_no2) + denitratation - denitritation
     _FLUX_VAR_(data%id_n2o) = _FLUX_VAR_(data%id_n2o) + nitrous_denitritation - denitrousation
     _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + dnra + ammonium_release

     fnit = data%K_nit/(data%K_nit+(nit+no2))
   ENDIF
   doc_min_anaerobic = doc_min_anaerobic - denitrification

   IF( data%simFeReduction>0 ) THEN
     !# Iron reduction : anaerobic breakdown, after denitrification etc
     fereduction = doc_min_anaerobic * feiii/(data%K_fe+feiii) * fnit
     ffe  = data%K_fe/(data%K_fe+feiii)
     _FLUX_VAR_(data%id_fe) = _FLUX_VAR_(data%id_fe) - fereduction ! stoich
   ENDIF
   doc_min_anaerobic = doc_min_anaerobic - fereduction

   IF( data%simSO4Reduction>0 ) THEN
     !# Sulfate reduction : anaerobic breakdown, after denitrification etc
     so4reduction = doc_min_anaerobic * so4/(data%K_so4+so4) * ffe * fnit
     fso4 = so4/(data%K_so4+so4)
     _FLUX_VAR_(data%id_so4) = _FLUX_VAR_(data%id_so4) - so4reduction ! stoich
   ENDIF
   doc_min_anaerobic = doc_min_anaerobic - so4reduction

   IF( data%simMethanogenesis>0 ) THEN
     !# Methanogenesis : anaerobic breakdown, after all other pathways
     methanogenesis = doc_min_anaerobic * fso4 * ffe * fnit
     _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) - methanogenesis  ! stoich
   ENDIF
   doc_min_anaerobic = doc_min_anaerobic - methanogenesis

   !---------------------------------------------------------------------------+
   ! Set temporal derivatives : PRODUCTS of mineralisation
   if (data%use_dic) THEN
      _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (doc_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_dic) =  &
                                    _FLUX_VAR_(data%id_dic) + photolysis*(1.-photo_fmin)
   ENDIF
   IF (data%use_amm) THEN
      IF( data%simDenitrification==1 ) &
        _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (don_mineralisation)
      !IF( data%simDenitrification==2 ) & !MH needs balacing with Denit fraction
      !  _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (don_mineralisation)
      IF ( data%simRPools )  &
         _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + &
                                       photolysis*(1.-photo_fmin)*MIN(donr/MAX(docr,1e-2),one_)
   ENDIF
   IF (data%use_frp) THEN
      _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + (dop_mineralisation)
      IF ( data%simRPools )  &
         _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + photolysis*(1.-photo_fmin)*(dopr/docr)
   ENDIF


   !---------------------------------------------------------------------------+
   !# Export diagnostic variables
   !---------------------------------------------------------------------------+

   !-- CDOM computed as a function of DOC amount, as empirically defined by
   !   Kostoglidis et al 2005 for Swan-Canning
   IF ( diag_level>0 ) THEN
     IF ( data%simRPools ) THEN
        cdoc = MIN(doc+docr,1.e4)
     ELSE
        cdoc = MIN(doc,1.e4)
     ENDIF
     _DIAG_VAR_(data%id_cdom) = 0.35*exp(0.1922*(cdoc)*(12./1e3))
   ENDIF

   !-- Extra diagnostics
   IF ( diag_level>1 ) THEN
     ! General breakdown rates
     _DIAG_VAR_(data%id_poc_miner) = poc_hydrolysis*secs_per_day
     _DIAG_VAR_(data%id_pon_miner) = pon_hydrolysis*secs_per_day
     _DIAG_VAR_(data%id_pop_miner) = pop_hydrolysis*secs_per_day
     _DIAG_VAR_(data%id_doc_miner) = doc_mineralisation*secs_per_day
     _DIAG_VAR_(data%id_don_miner) = don_mineralisation*secs_per_day
     _DIAG_VAR_(data%id_dop_miner) = dop_mineralisation*secs_per_day

     IF ( data%simRPools ) THEN
       _DIAG_VAR_(data%id_docr_miner) = docr_mineralisation*secs_per_day
       _DIAG_VAR_(data%id_donr_miner) = donr_mineralisation*secs_per_day
       _DIAG_VAR_(data%id_dopr_miner) = dopr_mineralisation*secs_per_day
     ENDIF

     ! BOD5 is computed as the amount of oxygen consumed over a 5-day period (mgO2/L)
     _DIAG_VAR_(data%id_bod) = 32./12.*(doc_mineralisation*secs_per_day*12./1e3)*5.0

     IF ( data%simPhotolysis ) &
        _DIAG_VAR_(data%id_photolysis) = photolysis*secs_per_day

     ! Anaerobic OM mineralisation
     _DIAG_VAR_(data%id_denit)     = denitrification*secs_per_day
     _DIAG_VAR_(data%id_anaerobic) = doc_min_anaerobic*secs_per_day
     IF ( data%simDenitrification>1 ) &
          _DIAG_VAR_(data%id_denit) = denitratation + denitritation + denitrousation
   ENDIF

END SUBROUTINE aed_calculate_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_organic_matter(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink / source terms of AED OGM.
! Everything here is in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, oxy

   ! State
   AED_REAL :: poc, pon, pop
   AED_REAL :: doc, don, dop

   ! Temporary variables
   AED_REAL :: Psed_poc = 0., Psed_pon = 0., Psed_pop = 0.
   AED_REAL :: poc_flux,doc_flux
   AED_REAL :: pon_flux,don_flux
   AED_REAL :: pop_flux,dop_flux
   AED_REAL :: Fsed_poc,Fsed_doc
   AED_REAL :: Fsed_pon,Fsed_don
   AED_REAL :: Fsed_pop,Fsed_dop
   AED_REAL :: fT, fDO, fDOM

!-------------------------------------------------------------------------------
!BEGIN
   Fsed_doc = 0.0 ; Fsed_poc = 0.0 ; Fsed_dop = 0.0 ; Fsed_pop = 0.0 ; Fsed_don = 0.0 ; Fsed_pon = 0.0

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp)    ! local temperature
   oxy  = 300. ; IF (data%use_oxy) oxy = _STATE_VAR_(data%id_oxy)

   ! Retrieve current (local) state variable values.
   poc = _STATE_VAR_(data%id_poc)      ! particulate organic carbon
   doc = _STATE_VAR_(data%id_doc)      ! dissolved organic carbon
   pon = _STATE_VAR_(data%id_pon)      ! particulate organic nitrogen
   don = _STATE_VAR_(data%id_don)      ! dissolved organic nitrogen
   pop = _STATE_VAR_(data%id_pop)      ! particulate organic phosphorus
   dop = _STATE_VAR_(data%id_dop)      ! dissolved organic phosphorus

   ! Set flux rate of dissolved organic matter pools
   fT = data%theta_sed_dom**(temp-20.0)
   fDO = data%Ksed_dom/(data%Ksed_dom+oxy)
   fDOM = 1.

   IF (data%use_Fsed_link_doc) THEN !Get the flux value from the linked variable
     Fsed_doc = _STATE_VAR_S_(data%id_Fsed_doc)
   ELSE                             !Compute directly
     Fsed_doc = data%Fsed_doc * fDO * fT * fDOM
   ENDIF
   IF (data%use_Fsed_link_don) THEN
     Fsed_don = _STATE_VAR_S_(data%id_Fsed_don)
   ELSE
     Fsed_don = data%Fsed_don * fDO * fT * fDOM
   ENDIF
   IF (data%use_Fsed_link_dop) THEN
     Fsed_dop = _STATE_VAR_S_(data%id_Fsed_dop)
   ELSE
     Fsed_dop = data%Fsed_dop * fDO * fT * fDOM
   ENDIF


   ! Set flux rate of particulate organic matter pools
   IF (data%use_Fsed_link_poc) THEN
     Fsed_poc = _STATE_VAR_S_(data%id_Fsed_poc)
   ELSE ! Compute directly
     Fsed_poc = zero_ !data%Fsed_pom * sedimentOMfrac * data%Xsc *(bottom_stress - data%tau_0) / data%tau_r
     IF( data%resuspension>0 ) Fsed_poc = _DIAG_VAR_S_(data%id_l_resus) * data%sedimentOMfrac * data%Xsc
   ENDIF
   IF (data%use_Fsed_link_pon) THEN
     Fsed_pon = _STATE_VAR_S_(data%id_Fsed_pon)
   ELSE ! Compute directly
     Fsed_pon = zero_ !data%Fsed_pom * sedimentOMfrac * data%Xsn *(bottom_stress - data%tau_0) / data%tau_r
     IF( data%resuspension>0 ) Fsed_pon = _DIAG_VAR_S_(data%id_l_resus) * data%sedimentOMfrac * data%Xsn
   ENDIF
   IF (data%use_Fsed_link_pop) THEN
     Fsed_pop = _STATE_VAR_S_(data%id_Fsed_pop)
   ELSE ! Compute directly
     Fsed_pop = zero_ !data%Fsed_pom * sedimentOMfrac * data%Xsp *(bottom_stress - data%tau_0) / data%tau_r
     IF( data%resuspension>0 ) Fsed_pop = _DIAG_VAR_S_(data%id_l_resus) * data%sedimentOMfrac * data%Xsp
   ENDIF

   ! Set bottom fluxes for the pelagic variables, note this doesnt consdier
   ! sedimentation of particles (change per surface area per second)
   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + Fsed_doc
   _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + Fsed_don
   _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) + Fsed_dop
   _FLUX_VAR_(data%id_poc) = _FLUX_VAR_(data%id_poc) + Fsed_poc
   _FLUX_VAR_(data%id_pon) = _FLUX_VAR_(data%id_pon) + Fsed_pon
   _FLUX_VAR_(data%id_pop) = _FLUX_VAR_(data%id_pop) + Fsed_pop

   ! Get sedimentation flux (mmmol/m2/s) loss into the sediment, crossing the SWI
   IF (data%use_Psed_link_poc) THEN
     !-- linked variable requires populating
     Psed_poc = data%w_pom * max(zero_,poc)
   ELSE
     Psed_poc = _DIAG_VAR_(data%id_Psed_poc)
   ENDIF
   IF (data%use_Psed_link_pon) THEN
     !-- linked variable requires populating
     Psed_pon = data%w_pom * max(zero_,pon)
   ELSE
     Psed_pon = _DIAG_VAR_(data%id_Psed_pon)
   ENDIF
   IF (data%use_Psed_link_pop) THEN
     !-- linked variable requires populating
     Psed_pop = data%w_pom * max(zero_,pop)
   ELSE
     !-- diagnostic was set in mobility
     Psed_pop = _DIAG_VAR_(data%id_Psed_pop)
   ENDIF

   ! Store the net sediment fluxes across the SWI as diag variables (mmol/m2/day)
   IF ( diag_level>1 ) THEN
     _DIAG_VAR_S_(data%id_swi_poc) = secs_per_day*(Fsed_poc + Psed_poc) ! resus & settling
     _DIAG_VAR_S_(data%id_swi_doc) = secs_per_day*(Fsed_doc)            ! dissolved flux
     _DIAG_VAR_S_(data%id_swi_pon) = secs_per_day*(Fsed_pon + Psed_pon) ! resus & settling
     _DIAG_VAR_S_(data%id_swi_don) = secs_per_day*(Fsed_don)            ! dissolved flux
     _DIAG_VAR_S_(data%id_swi_pop) = secs_per_day*(Fsed_pop + Psed_pop) ! resus & settling
     _DIAG_VAR_S_(data%id_swi_dop) = secs_per_day*(Fsed_dop)            ! dissolved flux
   ENDIF

   ! Set source & sink terms for the sediment pools (change per surface area per second)
   ! Note that this includes fluxes to and from the pelagic.
   IF( data%simSedimentOM>0 ) THEN
    _FLUX_VAR_B_(data%id_sed_poc) = _FLUX_VAR_B_(data%id_sed_poc) -(Fsed_poc + Psed_poc)
    _FLUX_VAR_B_(data%id_sed_pon) = _FLUX_VAR_B_(data%id_sed_pon) -(Fsed_pon + Psed_pon)
    _FLUX_VAR_B_(data%id_sed_pop) = _FLUX_VAR_B_(data%id_sed_pop) -(Fsed_pop + Psed_pop)
    _FLUX_VAR_B_(data%id_sed_doc) = _FLUX_VAR_B_(data%id_sed_doc) - Fsed_doc
    _FLUX_VAR_B_(data%id_sed_don) = _FLUX_VAR_B_(data%id_sed_don) - Fsed_don
    _FLUX_VAR_B_(data%id_sed_dop) = _FLUX_VAR_B_(data%id_sed_dop) - Fsed_dop
    IF ( diag_level>0 ) THEN
      _DIAG_VAR_S_(data%id_sed_toc) = _STATE_VAR_S_(data%id_sed_poc)+_STATE_VAR_S_(data%id_sed_doc)
      _DIAG_VAR_S_(data%id_sed_ton) = _STATE_VAR_S_(data%id_sed_pon)+_STATE_VAR_S_(data%id_sed_don)
      _DIAG_VAR_S_(data%id_sed_top) = _STATE_VAR_S_(data%id_sed_pop)+_STATE_VAR_S_(data%id_sed_dop)
    ENDIF
  ENDIF

END SUBROUTINE aed_calculate_benthic_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_organic_matter(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: doc, poc, cpom, cdom
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
END SUBROUTINE aed_light_extinction_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_organic_matter(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: vvel, vvel_cpom
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: rho_pom, temp
!
!-------------------------------------------------------------------------------
!BEGIN
     ! settling = 0 : no settling
     ! settling = 1 : constant settling @ w_pom
     ! settling = 2 : constant settling @ w_pom, corrected for variable density
     ! settling = 3 : settling based on Stoke's Law (calculated below)
     SELECT CASE (data%settling)

        CASE ( _MOB_OFF_ )
          ! disable settling by settign vertical velocity to 0
          vvel = zero_
          vvel_cpom = zero_

        CASE ( _MOB_CONST_ )
          ! constant settling velocity using user provided value
          vvel = data%w_pom
          vvel_cpom = data%w_cpom

        CASE ( _MOB_TEMP_ )
          ! constant settling velocity @20C corrected for density changes
          pw = _STATE_VAR_(data%id_rho)
          temp = _STATE_VAR_(data%id_temp)
          mu = water_viscosity(temp)
          mu20 = 0.001002  ! N s/m2
          pw20 = 998.2000  ! kg/m3 (assuming freshwater)
          vvel = data%w_pom*mu20*pw / ( mu*pw20 )
          vvel_cpom = data%w_cpom*mu20*pw / ( mu*pw20 )

        CASE ( _MOB_STOKES_ )
          ! settling velocity based on Stokes Law calculation and cell density

          pw = _STATE_VAR_(data%id_rho)              ! water density
          temp = _STATE_VAR_(data%id_temp)
          mu = water_viscosity(temp)                 ! water dynamic viscosity
          rho_pom = data%rho_pom
          vvel = -9.807*(data%d_pom**2.)*( rho_pom-pw ) / ( 18.*mu )
          IF(data%simRPools) &
          vvel_cpom = -9.807*(data%d_cpom**2.)*( data%rho_cpom-pw ) / ( 18.*mu )
        CASE DEFAULT
          ! unknown settling/migration option selection
          vvel = data%w_pom
          vvel_cpom = data%w_cpom

      END SELECT

      ! set global mobility array (m/s), later used to compute settling
      mobility(data%id_poc) = vvel
      mobility(data%id_pon) = vvel
      mobility(data%id_pop) = vvel
      IF(data%simRPools) mobility(data%id_cpom) = vvel_cpom
      ! set sedimentation flux (mmmol/m2) for later use/reporting
      _DIAG_VAR_(data%id_Psed_poc) = mobility(data%id_poc)*_STATE_VAR_(data%id_poc)
      _DIAG_VAR_(data%id_Psed_pon) = mobility(data%id_pon)*_STATE_VAR_(data%id_pon)
      _DIAG_VAR_(data%id_Psed_pop) = mobility(data%id_pop)*_STATE_VAR_(data%id_pop)
      IF(data%simRPools) _DIAG_VAR_(data%id_Psed_cpom) = mobility(data%id_cpom)*_STATE_VAR_(data%id_cpom)

      IF ( diag_level>1 ) THEN
        _DIAG_VAR_(data%id_pom_vvel) = vvel*secs_per_day
        IF(data%simRPools) _DIAG_VAR_(data%id_cpom_vvel) = vvel_cpom*secs_per_day
      ENDIF


END SUBROUTINE aed_mobility_organic_matter
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


!###############################################################################
PURE AED_REAL FUNCTION R_dom(Rdom_minerl,Ko2_an,f_an,theta,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: Rdom_minerl,Ko2_an,f_an,theta,oxy,temp
   AED_REAL :: f_oxy
!
!-----------------------------------------------------------------------
!BEGIN
   f_oxy = oxy/(Ko2_an+oxy) + f_an*(Ko2_an/(Ko2_an+oxy))
   R_dom = Rdom_minerl * f_oxy * (theta**(temp-20.0))

END FUNCTION R_dom
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION photo(Q, cdom, band)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: Q          ! [W/m2] incoming light intensity
   AED_REAL,INTENT(in) :: cdom       ! [/m] current cdom absorbance, based on DOC
   INTEGER, INTENT(in) :: band       ! [-] bandwidth identifier (VIS, UVA, UVB)
!LOCALS
!  AED_REAL,PARAMETER :: c = 7.52    ! in NML as photo_c and module global c
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

END FUNCTION photo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_organic_matter
