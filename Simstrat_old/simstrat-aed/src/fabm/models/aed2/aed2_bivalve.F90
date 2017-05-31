!###############################################################################
!#                                                                             #
!# aed2_bivalve.F90                                                            #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# In collaboration with :                                                     #
!#     Cornell University, Biology Department                                  #
!#                                                                             #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created January 2015                                                        #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#define _PHYLEN_ 3
#define _PHYMOD_ 'PHY'
#define _OGMPOC_ 'OGM'

MODULE aed2_bivalve
!-------------------------------------------------------------------------------
!  aed2_bivalve --- multi-group bivalve biogeochemical model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,ONLY : find_free_lun,aed2_bio_temp_function,fTemp_function,qsort
!  USE aed2_zoop_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_bivalve_data_t
!
   TYPE type_bivalve_prey
      !State variable name for bivalvelankton prey
      CHARACTER(64) :: bivalve_prey
      !Preference factors for bivalvelankton predators grazing on prey
      AED_REAL      :: Pbiv_prey
   END TYPE type_bivalve_prey


   TYPE type_bivalve_params
      ! General Attributes
      CHARACTER(64) :: name
      AED_REAL :: initial_conc, min
      INTEGER  :: Length

      ! Nutrient parameters
       AED_REAL :: INC, IPC

      ! Growth rate parameters
      AED_REAL :: Rgrz
      AED_REAL :: Ing, WaI, WbI
      AED_REAL :: fassim
      ! Minumum prey concentration parameters
      AED_REAL :: Cmin_grz, Kgrz
      AED_REAL :: minT, Tmin, Tmax, maxT, Dmax, maxD, SSmax, maxSS

      ! Respiration, mortaility and excretion parameters
      AED_REAL :: Rexcr, Regst, gegst, Rresp

      ! Salinity parameters
      INTEGER  :: saltfunc
      AED_REAL :: minS, Smin, Smax, maxS

      AED_REAL :: fR20, War, Wbr, fR

      AED_REAL :: theta_resp
      AED_REAL :: TmaxR, maxTR, Qresp

      AED_REAL :: SDA, Rmort, Rpred
      AED_REAL :: fDO, K_BDO, KDO

      ! The prey
      INTEGER  :: num_prey
      TYPE(type_bivalve_prey) :: prey(MAX_ZOOP_PREY)
   END TYPE

   TYPE,extends(type_bivalve_params) :: type_bivalve_data
      INTEGER  :: id_prey(MAX_ZOOP_PREY)
      INTEGER  :: id_phyIN(MAX_ZOOP_PREY), id_phyIP(MAX_ZOOP_PREY)
   END TYPE


   TYPE,extends(aed2_model_data_t) :: aed2_bivalve_data_t
      !# Variable identifiers
      INTEGER  :: id_biv(MAX_ZOOP_TYPES)
      INTEGER  :: id_Nexctarget,id_Nmorttarget
      INTEGER  :: id_Pexctarget,id_Pmorttarget
      INTEGER  :: id_Cexctarget,id_Cmorttarget
      INTEGER  :: id_DOupttarget
      INTEGER  :: id_SSupttarget
      INTEGER  :: id_tem, id_sal, id_extc
      INTEGER  :: id_grz,id_resp,id_mort
      INTEGER  :: id_excr, id_egst, id_tbiv, id_pbiv, id_3d_grz


      !# Model parameters
      INTEGER  :: num_biv
      TYPE(type_bivalve_data),DIMENSION(:),ALLOCATABLE :: bivalves
      LOGICAL  :: simDNexcr, simDPexcr, simDCexcr
      LOGICAL  :: simPNexcr, simPPexcr, simPCexcr
      LOGICAL  :: simSSupt, preyFeedback

     CONTAINS
         PROCEDURE :: define            => aed2_define_bivalve
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_bivalve
!        PROCEDURE :: mobility          => aed2_mobility_bivalve
!        PROCEDURE :: light_extinction  => aed2_light_extinction_bivalve
!        PROCEDURE :: delete            => aed2_delete_bivalve

   END TYPE

   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed2_bivalve_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count   !Number of bivalve groups
   INTEGER,INTENT(in)          :: list(*) !List of bivalve groups to simulate
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   AED_REAL :: Pbiv_prey(MAX_ZOOP_PREY)

   TYPE(type_bivalve_params)  :: bivalve_param(MAX_ZOOP_TYPES)
   NAMELIST /bivalve_params/ bivalve_param
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD',iostat=status)
    IF (status /= 0) STOP 'Error opening bivalves_params namelist file'
    read(tfil,nml=bivalve_params,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist bivalves_params'

    data%num_biv = count
    allocate(data%bivalves(count))
    DO i=1,count
       ! General
       data%bivalves(i)%name          = bivalve_param(list(i))%name
       data%bivalves(i)%initial_conc  = bivalve_param(list(i))%initial_conc
       data%bivalves(i)%min           = bivalve_param(list(i))%min
       data%bivalves(i)%Length        = bivalve_param(list(i))%Length
       data%bivalves(i)%INC           = bivalve_param(list(i))%INC
       data%bivalves(i)%IPC           = bivalve_param(list(i))%IPC
       ! Filtration & Ingestion
       data%bivalves(i)%Rgrz          = bivalve_param(list(i))%Rgrz/secs_per_day
       data%bivalves(i)%Ing           = bivalve_param(list(i))%Ing
       data%bivalves(i)%WaI           = bivalve_param(list(i))%WaI
       data%bivalves(i)%WbI           = bivalve_param(list(i))%WbI
       data%bivalves(i)%fassim        = bivalve_param(list(i))%fassim
       data%bivalves(i)%Cmin_grz      = bivalve_param(list(i))%Cmin_grz
       data%bivalves(i)%Kgrz          = bivalve_param(list(i))%Kgrz
       data%bivalves(i)%minT          = bivalve_param(list(i))%minT
       data%bivalves(i)%Tmin          = bivalve_param(list(i))%Tmin
       data%bivalves(i)%Tmax          = bivalve_param(list(i))%Tmax
       data%bivalves(i)%maxT          = bivalve_param(list(i))%maxT
       data%bivalves(i)%Dmax          = bivalve_param(list(i))%Dmax
       data%bivalves(i)%maxD          = bivalve_param(list(i))%maxD
       data%bivalves(i)%SSmax         = bivalve_param(list(i))%SSmax
       data%bivalves(i)%maxSS         = bivalve_param(list(i))%maxSS
       ! Excretion & Egestion
       data%bivalves(i)%Rexcr         = bivalve_param(list(i))%Rexcr/secs_per_day
       data%bivalves(i)%Regst         = bivalve_param(list(i))%Regst/secs_per_day
       data%bivalves(i)%gegst         = bivalve_param(list(i))%gegst
       ! Respiration
       data%bivalves(i)%Rresp         = bivalve_param(list(i))%Rresp/secs_per_day
       data%bivalves(i)%saltfunc      = bivalve_param(list(i))%saltfunc
       data%bivalves(i)%minS          = bivalve_param(list(i))%minS
       data%bivalves(i)%Smin          = bivalve_param(list(i))%Smin
       data%bivalves(i)%Smax          = bivalve_param(list(i))%Smax
       data%bivalves(i)%maxS          = bivalve_param(list(i))%maxS
       data%bivalves(i)%fR20          = bivalve_param(list(i))%fR20
       data%bivalves(i)%War           = bivalve_param(list(i))%War
       data%bivalves(i)%Wbr           = bivalve_param(list(i))%Wbr
       data%bivalves(i)%fR            = bivalve_param(list(i))%fR
       data%bivalves(i)%theta_resp    = bivalve_param(list(i))%theta_resp
       data%bivalves(i)%TmaxR         = bivalve_param(list(i))%TmaxR
       data%bivalves(i)%maxTR         = bivalve_param(list(i))%maxTR
       data%bivalves(i)%Qresp         = bivalve_param(list(i))%Qresp
       data%bivalves(i)%SDA           = bivalve_param(list(i))%SDA
       ! Mortality
       data%bivalves(i)%Rmort         = bivalve_param(list(i))%Rmort/secs_per_day
       data%bivalves(i)%Rpred         = bivalve_param(list(i))%Rpred/secs_per_day
       data%bivalves(i)%fDO           = bivalve_param(list(i))%fDO
       data%bivalves(i)%K_BDO         = bivalve_param(list(i))%K_BDO
       data%bivalves(i)%KDO           = bivalve_param(list(i))%KDO

       data%bivalves(i)%num_prey      = bivalve_param(list(i))%num_prey

       !Loop through prey variables assigning a target variable and preference factor
       !First sort in descending order of food preferences
       DO j=1,data%bivalves(i)%num_prey
          sort_i(j) = j
          Pbiv_prey(j) = bivalve_param(list(i))%prey(j)%Pbiv_prey
       ENDDO
       CALL qsort(Pbiv_prey,sort_i,1,data%bivalves(i)%num_prey)
       DO j=1,data%bivalves(i)%num_prey
          data%bivalves(i)%prey(j)%bivalve_prey = &
                   bivalve_param(list(i))%prey(sort_i(data%bivalves(i)%num_prey-j+1))%bivalve_prey
          data%bivalves(i)%prey(j)%Pbiv_prey = &
                   bivalve_param(list(i))%prey(sort_i(data%bivalves(i)%num_prey-j+1))%Pbiv_prey
       ENDDO

       ! Register group as a state variable
       data%id_biv(i) = aed2_define_sheet_variable(                 &
                              bivalve_param(list(i))%name,          &
                              'mmolC/m**2', 'bivalve',              &
                              bivalve_param(list(i))%initial_conc,  &
                              minimum=bivalve_param(list(i))%min)
    ENDDO
!
END SUBROUTINE aed2_bivalve_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_bivalve(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the bivalve biogeochemical model
!
!  Here, the aed2_bivalve namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_biv
   INTEGER  :: the_biv(MAX_ZOOP_TYPES)

   CHARACTER(len=64)  :: dn_target_variable='' !dissolved nitrogen target variable
   CHARACTER(len=64)  :: pn_target_variable='' !particulate nitrogen target variable
   CHARACTER(len=64)  :: dp_target_variable='' !dissolved phosphorus target variable
   CHARACTER(len=64)  :: pp_target_variable='' !particulate phosphorus target variable
   CHARACTER(len=64)  :: dc_target_variable='' !dissolved carbon target variable
   CHARACTER(len=64)  :: pc_target_variable='' !particulate carbon target variable
   CHARACTER(len=64)  :: do_uptake_variable='' !oxy uptake variable
   CHARACTER(len=64)  :: ss_uptake_variable='' !sus. solids uptake variable
   CHARACTER(len=128) :: dbase='aed2_bivalve_pars.nml'

   INTEGER  :: biv_i, prey_i, phy_i

   NAMELIST /aed2_bivalve/ num_biv, the_biv, &
                    dn_target_variable, pn_target_variable, dp_target_variable, &
                    pp_target_variable, dc_target_variable, pc_target_variable, &
                    do_uptake_variable, ss_uptake_variable, dbase
!-----------------------------------------------------------------------
!BEGIN
!print *,'**** Reading /aed2_bivalve/ namelist'
   ! Read the namelist
   read(namlst,nml=aed2_bivalve,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_bivalve'

    data%preyFeedback = .true.
    data%num_biv = 0
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed2_bivalve_load_params to values per second.
   CALL aed2_bivalve_load_params(data, dbase, num_biv, the_biv)


   ! Not required if we use the Spillman quadratic fT
   !CALL aed2_bio_temp_function(data%num_biv,                 &
   !                           data%bivalves%theta_resp_zoo, &
   !                           data%bivalves%Tstd_zoo,       &
   !                           data%bivalves%Topt_zoo,       &
   !                           data%bivalves%Tmax_zoo,       &
   !                           data%bivalves%aTn,            &
   !                           data%bivalves%bTn,            &
   !                           data%bivalves%kTn,            &
   !                           data%bivalves%name)


   !Register link to prey state variables
   DO biv_i = 1,num_biv
      phy_i = 0
      DO prey_i = 1,data%bivalves(biv_i)%num_prey
          data%bivalves(biv_i)%id_prey(prey_i) = aed2_locate_variable( &
                                       data%bivalves(biv_i)%prey(prey_i)%bivalve_prey)
          !If the prey is phytoplankton then also register state dependency on
          !internal nitrogen and phosphorus
          IF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
              phy_i = phy_i + 1
              data%bivalves(biv_i)%id_phyIN(phy_i) = aed2_locate_variable( &
                                       TRIM(data%bivalves(biv_i)%prey(prey_i)%bivalve_prey)//'_IN')
              data%bivalves(biv_i)%id_phyIP(phy_i) = aed2_locate_variable( &
                                       TRIM(data%bivalves(biv_i)%prey(prey_i)%bivalve_prey)//'_IP')

          ENDIF
      ENDDO
   ENDDO

   ! Register link to nutrient pools, if variable names are provided in namelist.
   data%simDNexcr = dn_target_variable .NE. ''
   IF (data%simDNexcr) THEN
     data%id_Nexctarget = aed2_locate_variable(dn_target_variable)
   ENDIF
   data%simDPexcr = dp_target_variable .NE. ''
   IF (data%simDPexcr) THEN
     data%id_Pexctarget = aed2_locate_variable(dp_target_variable)
   ENDIF
   data%simDCexcr = dc_target_variable .NE. ''
   IF (data%simDCexcr) THEN
     data%id_Cexctarget = aed2_locate_variable(dc_target_variable)
   ENDIF

   data%simPNexcr = pn_target_variable .NE. ''
   IF (data%simPNexcr) THEN
     data%id_Nmorttarget = aed2_locate_variable(pn_target_variable)
   ENDIF
   data%simPPexcr = pp_target_variable .NE. ''
   IF (data%simPPexcr) THEN
     data%id_Pmorttarget = aed2_locate_variable(pp_target_variable)
   ENDIF
   data%simPCexcr = pc_target_variable .NE. ''
   IF (data%simPCexcr) THEN
     data%id_Cmorttarget = aed2_locate_variable(pc_target_variable)
   ENDIF

   if (do_uptake_variable .EQ. '') STOP 'bivalve needs DO uptake variable'
   data%id_DOupttarget = aed2_locate_variable(do_uptake_variable)

   data%simSSupt = ss_uptake_variable .NE. ''
   IF (data%simSSupt) THEN
     data%id_SSupttarget = aed2_locate_variable(ss_uptake_variable)
   ENDIF

   ! Register diagnostic variables
   data%id_grz = aed2_define_sheet_diag_variable('grz','mmolC/m**3',  'net bivalve grazing')
   data%id_resp = aed2_define_sheet_diag_variable('resp','mmolC/m**3',  'net bivalve respiration')
   data%id_mort = aed2_define_sheet_diag_variable('mort','mmolC/m**3/d','net bivalve mortality')

   data%id_excr = aed2_define_sheet_diag_variable('excr','mmolC/m**3/d','net bivalve excretion')
   data%id_egst = aed2_define_sheet_diag_variable('egst','mmolC/m**3/d','net bivalve egestion')

   data%id_tbiv = aed2_define_sheet_diag_variable('tbiv','mmolC/m**2','total bivalve mass')
   data%id_pbiv = aed2_define_sheet_diag_variable('pbiv','mmolC/m**2','previous total bivalve mass')
   data%id_3d_grz = aed2_define_diag_variable('tgrz','mmolC/m**3','total bivalve graze')

   ! Register environmental dependencies
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_extc = aed2_locate_global('extc_coef')

END SUBROUTINE aed2_define_bivalve
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_bivalve(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of bivalve biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: biv,temp,salinity,oxy,ss !State variables
   AED_REAL :: prey(MAX_ZOOP_PREY), grazing_prey(MAX_ZOOP_PREY) !Prey state variables
   AED_REAL :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for phytoplankton
   AED_REAL :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL :: pon, pop, poc !Mortaility and fecal pellet state variables
   AED_REAL :: FGrazing_Limitation, f_Temp, f_Salinity, f_SS, I_max
   AED_REAL :: pref_factor, Ctotal_prey !total concentration of available prey
   AED_REAL :: food, grazing, respiration, mortality !Growth & decay functions
   AED_REAL :: grazing_n, grazing_p !Grazing on nutrients
   AED_REAL :: pon_excr, pop_excr, poc_excr !POM excretion rates
   AED_REAL :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   INTEGER  :: biv_i,prey_i,prey_j,phy_i

   !CAB added
   AED_REAL :: f_dens, W, Imax, psuedofaeces, ingestion, excretion, egestion, iteg, R20
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)       ! local temperature
   salinity = _STATE_VAR_(data%id_sal)   ! local salinity
   oxy = _STATE_VAR_(data%id_DOupttarget)! local oxygen
   IF (data%simSSupt) THEN
      ss = _STATE_VAR_(data%id_SSupttarget) ! local suspended solids (inorganic)
   ELSE
      ss = 0.
   ENDIF

   ! Retrieve current (local) state variable values.
   IF (data%simDNexcr)  dn_excr = _STATE_VAR_(data%id_Nexctarget)
   IF (data%simDPexcr)  dp_excr = _STATE_VAR_(data%id_Pexctarget)
   IF (data%simDCexcr)  dc_excr = _STATE_VAR_(data%id_Cexctarget)

   IF (data%simPNexcr)  pon = _STATE_VAR_(data%id_Nmorttarget)
   IF (data%simPPexcr)  pop = _STATE_VAR_(data%id_Pmorttarget)
   IF (data%simPCexcr)  poc = _STATE_VAR_(data%id_Cmorttarget)

   _DIAG_VAR_S_(data%id_tbiv ) = 0
   _DIAG_VAR_S_(data%id_pbiv ) = 0

   DO biv_i=1,data%num_biv

      ! Retrieve this bivalve group
      biv = _STATE_VAR_S_(data%id_biv(biv_i))
      _DIAG_VAR_S_(data%id_pbiv) = _DIAG_VAR_S_(data%id_pbiv) + biv

      grazing       = zero_
      respiration   = zero_
      mortality     = zero_

      !Retrieve prey groups
      Ctotal_prey   = zero_
      DO prey_i=1,data%bivalves(biv_i)%num_prey
         prey(prey_i) = _STATE_VAR_(data%bivalves(biv_i)%id_prey(prey_i))
         Ctotal_prey = Ctotal_prey + prey(prey_i)
      ENDDO

      ! Get the grazing limitation function
      !fGrazing_Limitation = fPrey_Limitation(data%bivalves,biv_i,Ctotal_prey)
      fGrazing_Limitation = min(Ctotal_prey/data%bivalves(biv_i)%Kgrz,one_)


      ! Get the temperature function
      f_Temp = fTemp_function_biv(data,biv_i, temp)
      !f_T = fTemp_function(1, data%bivalves(biv_i)%Tmax,       &
      !                        data%bivalves(biv_i)%Tstd,       &
      !                        data%bivalves(biv_i)%theta_resp, &
      !                        data%bivalves(biv_i)%aTn,        &
      !                        data%bivalves(biv_i)%bTn,        &
      !                        data%bivalves(biv_i)%kTn,        &
      !                        temp)


      ! Get the suspended solids function
      f_SS = fSS_function(data,biv_i,ss)

      ! Get the density limitation function
      f_dens = fD_function(data,biv_i,biv)

      ! Get the final ingestion rate (/ s)
      ! amount grazed in units of mass consumed/mass bivalve/unit time
      IF(data%bivalves(biv_i)%Ing==1) THEN
        W = (0.071/1000.) * data%bivalves(biv_i)%Length**2.8
        Imax = data%bivalves(biv_i)%WaI * W** data%bivalves(biv_i)%WbI
      ELSE
        Imax = data%bivalves(biv_i)%Rgrz
      END IF
      grazing = Imax * fGrazing_Limitation * f_Temp * f_dens * f_SS

      ! Now determine available prey and limit grazing amount to availability of prey
      ! food is total amount of food in units of mass/unit volume/unit time
      food = grazing * biv
      IF (Ctotal_prey < data%bivalves(biv_i)%num_prey * data%bivalves(biv_i)%Cmin_grz ) THEN
         food = zero_
         grazing = zero_
      ELSEIF (food > Ctotal_prey - data%bivalves(biv_i)%num_prey * data%bivalves(biv_i)%Cmin_grz ) THEN
         food = Ctotal_prey - data%bivalves(biv_i)%num_prey * data%bivalves(biv_i)%Cmin_grz
         grazing = food / biv
      ENDIF


      ! Now determine prey composition based on preference factors and availability of prey

      ! Prey has been ordered in grazing preference
      ! So take food in order of preference up to availability minus value of minimum residual
      ! grazing_prey is in units of mass consumed/unit volumne/unit time

      DO prey_i = 1,data%bivalves(biv_i)%num_prey
          !Add up preferences for remaining prey
          pref_factor = zero_
          DO prey_j = prey_i,data%bivalves(biv_i)%num_prey
             pref_factor = pref_factor + data%bivalves(biv_i)%prey(biv_i)%Pbiv_prey
          ENDDO
          IF (food * data%bivalves(biv_i)%prey(prey_i)%Pbiv_prey / pref_factor <= &
                        prey(prey_i) - data%bivalves(biv_i)%Cmin_grz) THEN
             !Take fraction of left over food based on preference factor
             grazing_prey(prey_i) = food * data%bivalves(biv_i)%prey(prey_i)%Pbiv_prey / pref_factor
          ELSEIF (prey(prey_i) > data%bivalves(biv_i)%Cmin_grz) THEN
             grazing_prey(prey_i) = prey(prey_i) - data%bivalves(biv_i)%Cmin_grz
          ELSE
             grazing_prey(prey_i) = zero_
          ENDIF
          !Food remaining after grazing from current prey
          food = food - grazing_prey(prey_i)
      ENDDO

      ! Now determine nutrient composition of food based on prey type
      ! At this stage only the AED model state variables have multiple
      ! nutrients (C,N&P) so assume all others have a single nutrient
      ! and thus not need to calculate nutrient excretion as is taken
      ! care of in the respiration term.  22/12/2011
      ! grazing_n is in units of mass N consumed/unit volume/unit time
      ! grazing_p is in units of mass P consumed/unit volume/unit time

      grazing_n = zero_
      grazing_p = zero_
      phy_i = 0
      DO prey_i = 1,data%bivalves(biv_i)%num_prey
         IF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey .EQ. _OGMPOC_) THEN
            IF (poc > zero_) THEN
                grazing_n = grazing_n + grazing_prey(prey_i) * pon/poc
                grazing_p = grazing_p + grazing_prey(prey_i) * pop/poc
            ELSE
                grazing_n = zero_
                grazing_p = zero_
            ENDIF
         ELSEIF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
            phy_i = phy_i + 1
            phy_INcon(phy_i) = _STATE_VAR_(data%bivalves(biv_i)%id_phyIN(phy_i))
            phy_IPcon(phy_i) = _STATE_VAR_(data%bivalves(biv_i)%id_phyIP(phy_i))
            grazing_n = grazing_n + grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i) /14.0
            grazing_p = grazing_p + grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i) /31.0
         ELSEIF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:15).EQ.'aed2_bivalve') THEN
            grazing_n = grazing_n + grazing_prey(prey_i) * data%bivalves(biv_i)%INC
            grazing_p = grazing_p + grazing_prey(prey_i) * data%bivalves(biv_i)%IPC
         ENDIF
      ENDDO



      psuedofaeces = (one_ - data%bivalves(biv_i)%fassim) * grazing
      ingestion = data%bivalves(biv_i)%fassim * grazing

      IF (biv <= data%bivalves(biv_i)%min) THEN
        ! Don't excrete or die if we are at the min biomass otherwise we have a
        ! mass conservation leak in the C mass balance

        respiration = zero_
        mortality   = zero_
        excretion   = zero_
        egestion    = zero_

      ELSE

        egestion = data%bivalves(biv_i)%Regst * exp(data%bivalves(biv_i)%gegst + &
                   min(Ctotal_prey/data%bivalves(biv_i)%Kgrz,one_)) * ingestion

        ! Get the respiration rate (/ s)
        iteg = ingestion - egestion
        respiration = aed2_bivalve_respiration(data,biv_i,iteg,temp,salinity)

        ! Get the excretion rate (of carbon!) (/s)
        excretion =  data%bivalves(biv_i)%Rexcr * iteg

        ! Get the mortality rate (/ s)
        mortality = data%bivalves(biv_i)%Rmort * f_DO(data,biv_i,oxy)

        ! Add the predation losses to mortality
        mortality = mortality + data%bivalves(biv_i)%Rpred


      ENDIF


      ! Calculate losses into the particulate organic matter pools - Units mmol/s
      poc_excr = (psuedofaeces + egestion + mortality)*biv

      pon_excr = (psuedofaeces * grazing_n / grazing  &
               +  egestion*data%bivalves(biv_i)%INC + mortality*data%bivalves(biv_i)%INC)*biv

      pop_excr = (psuedofaeces * grazing_p / grazing  &
               +  egestion*data%bivalves(biv_i)%IPC + mortality*data%bivalves(biv_i)%IPC)*biv

      ! Now we know the rates of carbon consumption and excretion, calculate rates
      ! of n & p excretion to maintain internal nutrient stores

      ! First, compute rate of change so far of bivalve carbon biomass (mmolC / m2 /s)
      delta_C = (ingestion - respiration - egestion - excretion - mortality) * biv


      ! Then calc nutrient excretion require to balance internal nutrient store
      ! Note pon_excr includes loss due to messy feeding so no need to include assimilation fraction on grazing_n & grazing_p
      don_excr = grazing_n - pon_excr - delta_C * data%bivalves(biv_i)%INC
      dop_excr = grazing_p - pop_excr - delta_C * data%bivalves(biv_i)%IPC

      ! If nutrients are limiting then must excrete doc to maintain balance
      IF ((don_excr < zero_) .AND. (dop_excr < zero_)) THEN
         !Determine which nutrient is more limiting
         IF ((data%bivalves(biv_i)%INC * (grazing_n - pon_excr) - delta_C) .GT. &
            (data%bivalves(biv_i)%IPC * (grazing_p - pop_excr) - delta_C)) THEN
             don_excr = zero_
             doc_excr =  (grazing_n - pon_excr) / data%bivalves(biv_i)%INC - delta_C
             delta_C = delta_C - doc_excr
             dop_excr = grazing_p - pop_excr - delta_C*data%bivalves(biv_i)%IPC
         ELSE
             dop_excr = zero_
             doc_excr = (grazing_p - pop_excr) / data%bivalves(biv_i)%IPC - delta_C
             delta_C = delta_C - doc_excr
             don_excr = grazing_n - pon_excr - delta_C*data%bivalves(biv_i)%INC
         ENDIF
      ELSEIF (don_excr < zero_) THEN !nitrogen limited
         don_excr = zero_
         doc_excr = (grazing_n - pon_excr) / data%bivalves(biv_i)%INC - delta_C
         delta_C = delta_C - doc_excr
         dop_excr = grazing_p - pop_excr - delta_C*data%bivalves(biv_i)%IPC
      ELSEIF (dop_excr < zero_) THEN !phosphorus limited
         dop_excr = zero_
         doc_excr = (grazing_p - pop_excr) / data%bivalves(biv_i)%IPC - delta_C
         delta_C = delta_C - doc_excr
         don_excr = grazing_n - pon_excr - delta_C*data%bivalves(biv_i)%INC
      ELSE !just excrete nutrients no need to balance c
          doc_excr = zero_
      ENDIF


      !write(*,"(4X,'limitations (f_T,f_Salinity): ',2F8.2)")f_T,f_Salinity
      !write(*,"(4X,'sources/sinks (grazing,respiration,mortaility): ',3F8.2)")grazing,excretion,mortality


      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Production / losses in mmolC/s

!print*,'POOP',ingestion,respiration,excretion,egestion,mortality
      _FLUX_VAR_B_(data%id_biv(biv_i)) = _FLUX_VAR_B_(data%id_biv(biv_i)) + ((ingestion - respiration - excretion - egestion - mortality)*biv )


      IF (data%preyFeedback) THEN
         ! Now take food grazed by zooplankton from food pools in mmolC/s
         phy_i = 0
         DO prey_i = 1,data%bivalves(biv_i)%num_prey
            _FLUX_VAR_(data%bivalves(biv_i)%id_prey(prey_i)) = _FLUX_VAR_(data%bivalves(biv_i)%id_prey(prey_i)) + ( -1.0 * grazing_prey(prey_i))
            IF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey .EQ. _OGMPOC_) THEN
               IF (poc > zero_) THEN
                  _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + ( -1.0 * grazing_prey(prey_i) * pon/poc)
                  _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + ( -1.0 * grazing_prey(prey_i) * pop/poc)
               ENDIF
            ELSEIF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
               phy_i = phy_i + 1
               _FLUX_VAR_(data%bivalves(biv_i)%id_phyIN(phy_i)) = _FLUX_VAR_(data%bivalves(biv_i)%id_phyIN(phy_i)) + ( -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i))
               _FLUX_VAR_(data%bivalves(biv_i)%id_phyIP(phy_i)) = _FLUX_VAR_(data%bivalves(biv_i)%id_phyIP(phy_i)) + ( -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i))
            ENDIF
            _DIAG_VAR_(data%id_3d_grz) = _FLUX_VAR_(data%bivalves(biv_i)%id_prey(prey_i))
         ENDDO
      ENDIF


      ! Now manage excretion contributions to DOM
      IF (data%simDCexcr) THEN
         _FLUX_VAR_(data%id_Cexctarget) = _FLUX_VAR_(data%id_Cexctarget) + excretion + doc_excr
      ENDIF
      IF (data%simDNexcr) THEN
         _FLUX_VAR_(data%id_Nexctarget) = _FLUX_VAR_(data%id_Nexctarget) + (don_excr)
      ENDIF
      IF (data%simDPexcr) THEN
         _FLUX_VAR_(data%id_Pexctarget) = _FLUX_VAR_(data%id_Pexctarget) + (dop_excr)
      ENDIF

      ! Now manage psuedofaeces, egestion and mortality contributions to POM
      IF (data%simPCexcr) THEN
         _FLUX_VAR_(data%id_Cmorttarget) = _FLUX_VAR_(data%id_Cmorttarget) + ( poc_excr * mortality )
      ENDIF
      IF (data%simPNexcr) THEN
         _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + ( pon_excr)
      ENDIF
      IF (data%simPPexcr) THEN
         _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + ( pop_excr)
      ENDIF

      ! Export diagnostic variables
      _DIAG_VAR_S_(data%id_grz)  = grazing*secs_per_day
      _DIAG_VAR_S_(data%id_resp) = respiration*secs_per_day
      _DIAG_VAR_S_(data%id_mort) = mortality*secs_per_day

      _DIAG_VAR_S_(data%id_excr) = excretion*secs_per_day
      _DIAG_VAR_S_(data%id_egst) = egestion*secs_per_day

      _DIAG_VAR_S_(data%id_tbiv) = _DIAG_VAR_S_(data%id_tbiv) + biv
   ENDDO
END SUBROUTINE aed2_calculate_benthic_bivalve
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_bivalve_respiration(data,biv_i,iteg,temp,sal) RESULT(resp)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER  :: biv_i ! Invertebrate group
   AED_REAL, INTENT(IN)                :: temp  ! Temp value being used
   AED_REAL, INTENT(IN)                :: sal   ! Salinity value being used
   AED_REAL, INTENT(IN)                :: iteg  ! Ingestion-Egestion
!
!LOCALS
   AED_REAL :: W, TmaxR, maxTR, VV,WW,YY,XX,fT
   AED_REAL :: resp, Q, R20
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Get R20 value
   IF(data%bivalves(biv_i)%fR20==1) THEN
     ! Compute respiration rate coefficient from length
     W = (0.071/1000.) * data%bivalves(biv_i)%Length**2.8
     R20 = data%bivalves(biv_i)%War * W** data%bivalves(biv_i)%Wbr
   ELSE
     ! Use constant respiration rate coefficient
     R20 = data%bivalves(biv_i)%Rresp
   END IF

   ! Now compute actual respiration
   IF(data%bivalves(biv_i)%fR==1) THEN
     ! Option 1) Spillman et al 2008; Bocioniv et al 2013
     resp = R20 * data%bivalves(biv_i)%theta_resp**(temp-20.0)

   ELSEIF(data%bivalves(biv_i)%fR==2) THEN
     ! Option 2) Modified Schneider 1992; Bierman 2005; Gudimov et al. 2015

     TmaxR = data%bivalves(biv_i)%TmaxR
     maxTR = data%bivalves(biv_i)%maxTR
     Q     = data%bivalves(biv_i)%Qresp

     VV = ((TmaxR - temp)/(TmaxR - maxTR))
     WW = log(Q)*(TmaxR - maxTR)
     YY = log(Q)*(TmaxR - maxTR + 2.)
     XX = (WW * (1. + SQRT(1. + (40. / YY))) / 20.)**2
     fT = VV**XX * exp(XX*(1.-VV))

     resp = R20 * fT + data%bivalves(biv_i)%SDA * iteg

   ELSE
     ! Unknown option
     resp = data%bivalves(biv_i)%Rresp
   END IF

!   ! Get the salinity limitation.
!   resp = resp * fSalinity_Limitation(data%bivalves,biv_i,sal)

END FUNCTION aed2_bivalve_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION f_DO(data,biv,oxy)
!-------------------------------------------------------------------------------
! Dissolved oxygen effect on bivalve mortality
!
! Add f(DO) to base mortality, ie. mort is only ever enhanced by low DO
! Two options:
! 1) Abrupt increase in mortality once below a minimum DO threshold is exceeded
! 2) Steady exponential increase as DO decreases, enhancing mortality
!    by K_BDO so that when DO = 0, f(DO) = 1 + K_BDO
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(in) :: data
   INTEGER,INTENT(in)   :: biv              ! Invertebrate group
   AED_REAL, INTENT(IN) :: oxy              ! Dissolved oxygen
!
!LOCALS
   AED_REAL :: bot     ! Work Array
!
!-------------------------------------------------------------------------------
!BEGIN

   IF(data%bivalves(biv)%fDO==0) THEN
     ! Option 0) Abrupt increase below DO threshold
      IF (oxy < data%bivalves(biv)%KDO) THEN
        f_DO = 10.0
      ELSE
        f_DO = 1.0
      END IF

   ELSEIF(data%bivalves(biv)%fDO==1) THEN
     ! Option 1) Exponential increase: steepness determined by KDO and K_BDO
     ! i.e. when DO = 0, final f(DO) = 1 + K_BDO
     bot  = data%bivalves(biv)%KDO + oxy
     f_DO = 1.0 + data%bivalves(biv)%K_BDO * data%bivalves(biv)%KDO/bot

   ELSE
     ! Unknown option
     f_DO = 1.0
   END IF

END FUNCTION f_DO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fTemp_function_biv(data,biv,temp)
!-------------------------------------------------------------------------------
! Temperature growth multiplier for bivalves
!
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                 :: biv   ! Invertebrate group
   AED_REAL, INTENT(IN)                :: temp  ! Temp value being used
!
!LOCALS
   AED_REAL , PARAMETER :: a = 1.00
   AED_REAL             :: MINt,Tmin,Tmax,MAXt

!
!-------------------------------------------------------------------------------
!BEGIN

   MINt = data%bivalves(biv)%MINt
   Tmin = data%bivalves(biv)%Tmin
   Tmax = data%bivalves(biv)%Tmax
   MAXt = data%bivalves(biv)%MAXt

   ! If temp below extreme, temp fn = 0, ie filtration ceases
   IF(temp <= MINt) &
      fTemp_function_biv = zero_

   ! If below min temp for optimal range, limited as below
   IF (temp>MINt .and. temp<Tmin) THEN
      fTemp_function_biv = a*(1.0/(-((Tmin-MINt)*(Tmin-MINt)/ &
            (Tmin*Tmin))+(2.0*(Tmin-MINt)/ &
            Tmin)))*(-((temp-MINt)*(temp-MINt)/&
            (Tmin*Tmin))+2.0*((temp-MINt) / Tmin))
   END IF

   ! If between Tmin and Tmax, ie in the opt temp range, then not limited
   IF (temp >=Tmin .and. temp<=Tmax) THEN
      fTemp_function_biv = a
   END IF

   ! If above max temp for optimal range, limited as below
   IF (temp>Tmax .and. temp <MAXt) THEN
      fTemp_function_biv = a*((-(temp*temp)+2.0*(Tmax*temp)- &
            2.0*(Tmax*MAXt)+(MAXt*MAXt))/ &
            ((Tmax-MAXt)*(Tmax-MAXt)))
   END IF

   ! If temp above extreme, temp fn = 0, ie filtration ceases
   IF (temp>=MAXt) &
       fTemp_function_biv = zero_

 END FUNCTION fTemp_function_biv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 PURE AED_REAL FUNCTION fSS_function(data,biv,SS)
!-------------------------------------------------------------------------------
! Suspended solids function for clams
!
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                 :: biv   ! Invertebrate group
   AED_REAL, INTENT(IN)                :: SS    ! Susp solids value being used
!
!LOCALS
   AED_REAL , PARAMETER :: a = 1.00
   AED_REAL             :: pseudo, maxSS

!
!-------------------------------------------------------------------------------
!BEGIN

   fSS_function = one_

   pseudo =  data%bivalves(biv)%SSmax    ! SS conc where ingestion decreases
   maxSS  =  data%bivalves(biv)%maxSS    ! SS conc where mussel is buggered


   IF(SS <= pseudo) THEN
        fSS_function = one_
   ELSEIF(SS > pseudo .AND. SS <= maxSS) THEN
        ! From Spillman et al 2008
        fSS_function = (-SS*SS + 2.0*pseudo*SS - &
               2.0*pseudo*maxSS + maxSS*maxSS)/  &
               ((pseudo-maxSS) * (pseudo-maxSS))
   ELSE
        fSS_function = zero_
   END IF

   ! Ensure suspended solids function is not negative
   IF (fSS_function <= zero_) &
        fSS_function = zero_

 END FUNCTION fSS_function
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 PURE AED_REAL FUNCTION fD_function(data,biv,dens)
!-------------------------------------------------------------------------------
! Suspended solids function for clams
!
!ARGUMENTS
   CLASS (aed2_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                 :: biv   ! Invertebrate group
   AED_REAL, INTENT(in)                :: dens  ! Density value being used
!
!LOCALS
   AED_REAL             :: maxD, Dmax

!
!-------------------------------------------------------------------------------
!BEGIN

   fD_function = one_

   Dmax =  data%bivalves(biv)%Dmax    ! density where ingestion decreases
   maxD =  data%bivalves(biv)%maxD    ! density where mussel is buggered


   IF(dens <= Dmax) THEN
        fD_function = one_
   ELSEIF(dens > Dmax .AND. dens <= maxD) THEN
        fD_function = (-dens*dens + 2.0*Dmax*dens - &
                       2.0*Dmax*maxD + maxD*maxD)/  &
                      ((Dmax-maxD) * (Dmax-maxD))
   ELSE
        fD_function = zero_
   END IF

   ! Ensure fD is not negative
   IF (fD_function <= zero_) &
        fD_function = zero_

 END FUNCTION fD_function
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if 0
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
! Salinity tolerance of clams/mussels                                   !
!-----------------------------------------------------------------------!
! Directly cut and pasted from Maya's clam model                        !
 FUNCTION clamSalinity(Sbot) RESULT(ret)                                !
   !-- Incoming                                                         !
   REAL (r_wq) :: Sbot(:)            ! Salinity in the bottom layer     !
   !-- Returns the salinity function                                    !
   REAL (r_wq), DIMENSION(SIZE(Sbot)) :: ret
   REAL (r_wq), DIMENSION(SIZE(Sbot)) :: tmp1,tmp2

   ! BSmin minimum sal and BSmax is max sal set in WQcons               !
   ! Don't seem to need Bep, Aep or Sop from consts file                !
   ! Salinity is non-limiting for sals > BSmin and <BSmax               !
   ! clamLowCut is salinity when shells close up and feeding ceases     !
   ! CS 020609 Need to make an exception if FW species and BSmin and    !
   ! clamLowCut are both equal to zero.

   ! Salinity is within the tolerance; no limitation.                   !
   ! If sal in bott cell is > min sal & < max sal set in WQcons, pf=1   !
   WHERE(Sbot >= BSmin .and. Sbot <= BSmax)
     ret = wq_one
   END WHERE

   ! Salinity is greater than the upper bound                           !
   ! maxS is set in caedym_globals at 45psu                             !
   WHERE(Sbot > BSmax)
     ret = (-Sbot*Sbot+2.0*BSmax*Sbot-  &
            2.0*BSmax*maxS+maxS*maxS)/((BSmax-maxS)*(BSmax-maxS))
   END WHERE

   ! Salinity is less than the lower bound but greater than low cut     !
   ! If sal is < min set in WQcons but > clamLowCut set at clamCons.dat)!
   tmp1 = wq_zero
   tmp2 = wq_zero
   WHERE(Sbot < BSmin .AND. Sbot > clamLowCut)
     tmp1 = Sbot-clamLowCut
     tmp2 = BSmin-clamLowCut
     ret =  (2*tmp1/BSmin-(tmp1*tmp1/(BSmin*BSmin)))/ &
            (2*tmp2/BSmin-(tmp2*tmp2/(BSmin*BSmin)))
   END WHERE

   ! Salinity is less than the clamLowCut                               !
   ! If sal < lowest sal (hardwired at start of fn), shells close       !
   WHERE(Sbot <= clamLowCut)
     ret = wq_zero
   END WHERE

   ! If lower bound and low cut are both zero i.e. Freshwater species   !
   ! then need to set f(S) to one                                       !
   IF(BSmin==wq_zero) THEN
     WHERE(Sbot <= BSmin)
         ret =  wq_one
     END WHERE
   ENDIF

   ! Ensure temp function is not negative                               !
   WHERE(ret <= wq_zero)
     ret = wq_zero
   END WHERE

 END FUNCTION clamSalinity
!-----------------------------------------------------------------------!


#endif


END MODULE aed2_bivalve
