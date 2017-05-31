!###############################################################################
!#                                                                             #
!# aed2_vegetation.F90                                                         #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
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

MODULE aed2_vegetation
!-------------------------------------------------------------------------------
!  aed2_vegetation --- multi-group vegetation biogeochemical model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,ONLY : find_free_lun,aed2_bio_temp_function,fTemp_function,qsort

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_vegetation_data_t



   TYPE type_vegetation_params
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

   END TYPE

   TYPE,extends(aed2_model_data_t) :: aed2_vegetation_data_t
      !# Variable identifiers
      INTEGER  :: id_veg(MAX_VEG_TYPES)
      INTEGER  :: id_Nexctarget,id_Nmorttarget
      INTEGER  :: id_Pexctarget,id_Pmorttarget
      INTEGER  :: id_Cexctarget,id_Cmorttarget
      INTEGER  :: id_DOupttarget
      INTEGER  :: id_SSupttarget
      INTEGER  :: id_tem, id_sal, id_extc
      INTEGER  :: id_grz,id_resp,id_mort
      INTEGER  :: id_excr, id_egst

      !# Model parameters
      INTEGER  :: num_veg
      TYPE(type_vegetation_params),DIMENSION(:),ALLOCATABLE :: vegdata
      LOGICAL  :: simDNexcr, simDPexcr, simDCexcr
      LOGICAL  :: simPNexcr, simPPexcr, simPCexcr
      LOGICAL  :: simSSupt, preyFeedback

     CONTAINS
         PROCEDURE :: define             => aed2_define_vegetation
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_vegetation
!        PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_vegetation
!        PROCEDURE :: mobility           => aed2_mobility_vegetation
!        PROCEDURE :: light_extinction   => aed2_light_extinction_vegetation
!        PROCEDURE :: delete             => aed2_delete_vegetation

   END TYPE

   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed2_vegetation_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count   !Number of vegetation groups
   INTEGER,INTENT(in)          :: list(*) !List of vegetation groups to simulate
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   AED_REAL :: Pveg_prey(MAX_ZOOP_PREY)

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   TYPE(type_vegetation_params)  :: vegetation_param(MAX_ZOOP_TYPES)
   NAMELIST /vegetation_params/ vegetation_param
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD',iostat=status)
    IF (status /= 0) STOP 'Error opening vegdata_params namelist file'
    read(tfil,nml=vegetation_params,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist vegdata_params'

    data%num_veg = count
    allocate(data%vegdata(count))
    DO i=1,count
       ! General
       data%vegdata(i)%name          = vegetation_param(list(i))%name
       data%vegdata(i)%initial_conc  = vegetation_param(list(i))%initial_conc
       data%vegdata(i)%min           = vegetation_param(list(i))%min
       data%vegdata(i)%Length        = vegetation_param(list(i))%Length
       data%vegdata(i)%INC           = vegetation_param(list(i))%INC
       data%vegdata(i)%IPC           = vegetation_param(list(i))%IPC
       ! Filtration & Ingestion
       data%vegdata(i)%Rgrz          = vegetation_param(list(i))%Rgrz/secs_pr_day
       data%vegdata(i)%Ing           = vegetation_param(list(i))%Ing
       data%vegdata(i)%WaI           = vegetation_param(list(i))%WaI
       data%vegdata(i)%WbI           = vegetation_param(list(i))%WbI
       data%vegdata(i)%fassim        = vegetation_param(list(i))%fassim
       data%vegdata(i)%Cmin_grz      = vegetation_param(list(i))%Cmin_grz
       data%vegdata(i)%Kgrz          = vegetation_param(list(i))%Kgrz
       data%vegdata(i)%minT          = vegetation_param(list(i))%minT
       data%vegdata(i)%Tmin          = vegetation_param(list(i))%Tmin
       data%vegdata(i)%Tmax          = vegetation_param(list(i))%Tmax
       data%vegdata(i)%maxT          = vegetation_param(list(i))%maxT
       data%vegdata(i)%Dmax          = vegetation_param(list(i))%Dmax
       data%vegdata(i)%maxD          = vegetation_param(list(i))%maxD
       data%vegdata(i)%SSmax         = vegetation_param(list(i))%SSmax
       data%vegdata(i)%maxSS         = vegetation_param(list(i))%maxSS
       ! Excretion & Egestion
       data%vegdata(i)%Rexcr         = vegetation_param(list(i))%Rexcr/secs_pr_day
       data%vegdata(i)%Regst         = vegetation_param(list(i))%Regst/secs_pr_day
       data%vegdata(i)%gegst         = vegetation_param(list(i))%gegst
       ! Respiration
       data%vegdata(i)%Rresp         = vegetation_param(list(i))%Rresp/secs_pr_day
       data%vegdata(i)%saltfunc      = vegetation_param(list(i))%saltfunc
       data%vegdata(i)%minS          = vegetation_param(list(i))%minS
       data%vegdata(i)%Smin          = vegetation_param(list(i))%Smin
       data%vegdata(i)%Smax          = vegetation_param(list(i))%Smax
       data%vegdata(i)%maxS          = vegetation_param(list(i))%maxS
       data%vegdata(i)%fR20          = vegetation_param(list(i))%fR20
       data%vegdata(i)%War           = vegetation_param(list(i))%War
       data%vegdata(i)%Wbr           = vegetation_param(list(i))%Wbr
       data%vegdata(i)%fR            = vegetation_param(list(i))%fR
       data%vegdata(i)%theta_resp    = vegetation_param(list(i))%theta_resp
       data%vegdata(i)%TmaxR         = vegetation_param(list(i))%TmaxR
       data%vegdata(i)%maxTR         = vegetation_param(list(i))%maxTR
       data%vegdata(i)%Qresp         = vegetation_param(list(i))%Qresp
       data%vegdata(i)%SDA           = vegetation_param(list(i))%SDA
       ! Mortality
       data%vegdata(i)%Rmort         = vegetation_param(list(i))%Rmort/secs_pr_day
       data%vegdata(i)%Rpred         = vegetation_param(list(i))%Rpred/secs_pr_day
       data%vegdata(i)%fDO           = vegetation_param(list(i))%fDO
       data%vegdata(i)%K_BDO         = vegetation_param(list(i))%K_BDO
       data%vegdata(i)%KDO           = vegetation_param(list(i))%KDO

!      data%vegdata(i)%num_prey      = vegetation_param(list(i))%num_prey

       !Loop through prey variables assigning a target variable and preference factor
       !First sort in descending order of food preferences
!      DO j=1,data%vegdata(i)%num_prey
!         sort_i(j) = j
!         Pveg_prey(j) = vegetation_param(list(i))%prey(j)%Pveg_prey
!      ENDDO
!      CALL qsort(Pveg_prey,sort_i,1,data%vegdata(i)%num_prey)
!      DO j=1,data%vegdata(i)%num_prey
!         data%vegdata(i)%prey(j)%vegetation_prey = &
!                  vegetation_param(list(i))%prey(sort_i(data%vegdata(i)%num_prey-j+1))%vegetation_prey
!         data%vegdata(i)%prey(j)%Pveg_prey = &
!                  vegetation_param(list(i))%prey(sort_i(data%vegdata(i)%num_prey-j+1))%Pveg_prey
!      ENDDO

       ! Register group as a state variable
       data%id_veg(i) = aed2_define_sheet_variable(                 &
                              vegetation_param(list(i))%name,          &
                              'mmolC/m**2', 'vegetation',              &
                              vegetation_param(list(i))%initial_conc,  &
                              minimum=vegetation_param(list(i))%min)
    ENDDO
!
END SUBROUTINE aed2_vegetation_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_define_vegetation(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the vegetation biogeochemical model
!
!  Here, the aed2_vegetation namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_veg
   INTEGER  :: the_veg(MAX_ZOOP_TYPES)

   CHARACTER(len=64)  :: dn_target_variable='' !dissolved nitrogen target variable
   CHARACTER(len=64)  :: pn_target_variable='' !particulate nitrogen target variable
   CHARACTER(len=64)  :: dp_target_variable='' !dissolved phosphorus target variable
   CHARACTER(len=64)  :: pp_target_variable='' !particulate phosphorus target variable
   CHARACTER(len=64)  :: dc_target_variable='' !dissolved carbon target variable
   CHARACTER(len=64)  :: pc_target_variable='' !particulate carbon target variable
   CHARACTER(len=64)  :: do_uptake_variable='' !oxy uptake variable
   CHARACTER(len=64)  :: ss_uptake_variable='' !sus. solids uptake variable
   CHARACTER(len=128) :: dbase='aed2_vegetation_pars.nml'

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   INTEGER  :: veg_i, prey_i, phy_i

   NAMELIST /aed2_vegetation/ num_veg, the_veg, &
                    dn_target_variable, pn_target_variable, dp_target_variable, &
                    pp_target_variable, dc_target_variable, pc_target_variable, &
                    do_uptake_variable, ss_uptake_variable, dbase
!-----------------------------------------------------------------------
!BEGIN
!print *,'**** Reading /aed2_vegetation/ namelist'
   ! Read the namelist
   read(namlst,nml=aed2_vegetation,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_vegetation'

    data%preyFeedback = .true.
    data%num_veg = 0
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed2_vegetation_load_params to values per second.
   CALL aed2_vegetation_load_params(data, dbase, num_veg, the_veg)


   ! Not required if we use the Spillman quadratic fT
   !CALL aed2_bio_temp_function(data%num_veg,                 &
   !                           data%vegdata%theta_resp_zoo, &
   !                           data%vegdata%Tstd_zoo,       &
   !                           data%vegdata%Topt_zoo,       &
   !                           data%vegdata%Tmax_zoo,       &
   !                           data%vegdata%aTn,            &
   !                           data%vegdata%bTn,            &
   !                           data%vegdata%kTn,            &
   !                           data%vegdata%name)


   !Register link to prey state variables
!  DO veg_i = 1,num_veg
!     phy_i = 0
!     DO prey_i = 1,data%vegdata(veg_i)%num_prey
!         data%vegdata(veg_i)%id_prey(prey_i) = aed2_locate_variable( &
!                                      data%vegdata(veg_i)%prey(prey_i)%vegetation_prey)
!         !If the prey is phytoplankton then also register state dependency on
!         !internal nitrogen and phosphorus
!         IF (data%vegdata(veg_i)%prey(prey_i)%vegetation_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
!             phy_i = phy_i + 1
!             data%vegdata(veg_i)%id_phyIN(phy_i) = aed2_locate_variable( &
!                                      TRIM(data%vegdata(veg_i)%prey(prey_i)%vegetation_prey)//'_IN')
!             data%vegdata(veg_i)%id_phyIP(phy_i) = aed2_locate_variable( &
!                                      TRIM(data%vegdata(veg_i)%prey(prey_i)%vegetation_prey)//'_IP')

!         ENDIF
!     ENDDO
!  ENDDO

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

   if (do_uptake_variable .EQ. '') STOP 'vegetation needs DO uptake variable'
   data%id_DOupttarget = aed2_locate_variable(do_uptake_variable)

   data%simSSupt = ss_uptake_variable .NE. ''
   IF (data%simSSupt) THEN
     data%id_SSupttarget = aed2_locate_variable(ss_uptake_variable)
   ENDIF

   ! Register diagnostic variables
   data%id_grz = aed2_define_sheet_diag_variable('grz','mmolC/m**3',  'net vegetation grazing')
   data%id_resp = aed2_define_sheet_diag_variable('resp','mmolC/m**3',  'net vegetation respiration')
   data%id_mort = aed2_define_sheet_diag_variable('mort','mmolC/m**3/d','net vegetation mortality')

   data%id_excr = aed2_define_sheet_diag_variable('excr','mmolC/m**3/d','net vegetation excretion')
   data%id_egst = aed2_define_sheet_diag_variable('egst','mmolC/m**3/d','net vegetation egestion')

   ! Register environmental dependencies
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_extc = aed2_locate_global('extc_coef')

END SUBROUTINE aed2_define_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian_vegetation(data,column,layer_idx, pc_wet)
!-------------------------------------------------------------------------------
! Right hand sides of vegetation biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in)    :: data
   TYPE (aed2_column_t),          INTENT(inout) :: column(:)
   INTEGER,                       INTENT(in)    :: layer_idx
   AED_REAL,                      INTENT(in)    :: pc_wet
!
!LOCALS
   AED_REAL           :: veg,temp,salinity !State variables
   AED_REAL           :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for veg groups
   AED_REAL           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL           :: pon, pop, poc             !Mortaility and literfall state variables
   AED_REAL           :: FGrowth_Limitation, f_Temp, f_Salinity
   AED_REAL           :: photsynthesis, respiration, mortality !Growth & decay functions
   AED_REAL           :: pon_excr, pop_excr, poc_excr !POM excretion rates
   AED_REAL           :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   INTEGER            :: veg_i

   !CAB added
   AED_REAL :: f_dens, W, Imax, psuedofaeces, ingestion, excretion, egestion, iteg, R20, oxy


   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)       ! local temperature
   salinity = _STATE_VAR_(data%id_sal)   ! local salinity
   oxy = _STATE_VAR_(data%id_DOupttarget)! local oxygen
!  IF (data%simSSupt) THEN
!     ss = _STATE_VAR_(data%id_SSupttarget) ! local suspended solids (inorganic)
!  ELSE
!     ss = 0.
!  ENDIF

   ! Retrieve current (local) state variable values.
   IF (data%simDNexcr)  dn_excr = _STATE_VAR_(data%id_Nexctarget)
   IF (data%simDPexcr)  dp_excr = _STATE_VAR_(data%id_Pexctarget)
   IF (data%simDCexcr)  dc_excr = _STATE_VAR_(data%id_Cexctarget)

   IF (data%simPNexcr)  pon = _STATE_VAR_(data%id_Nmorttarget)
   IF (data%simPPexcr)  pop = _STATE_VAR_(data%id_Pmorttarget)
   IF (data%simPCexcr)  poc = _STATE_VAR_(data%id_Cmorttarget)

   DO veg_i=1,data%num_veg

      ! Retrieve this vegetation group
      veg = _STATE_VAR_S_(data%id_veg(veg_i))

!     grazing       = zero_
      respiration   = zero_
      mortality     = zero_

      !Retrieve prey groups
!     Ctotal_prey   = zero_
!     DO prey_i=1,data%vegdata(veg_i)%num_prey
!        prey(prey_i) = _STATE_VAR_(data%vegdata(veg_i)%id_prey(prey_i))
!        Ctotal_prey = Ctotal_prey + prey(prey_i)
!     ENDDO

      ! Get the grazing limitation function
       !fGrazing_Limitation = fPrey_Limitation(data%vegdata,veg_i,Ctotal_prey)
!      fGrazing_Limitation = min(Ctotal_prey/data%vegdata(veg_i)%Kgrz,one_)


      ! Get the temperature function
       f_Temp = fTemp_function_veg(data,veg_i, temp)
       !f_T = fTemp_function(1, data%vegdata(veg_i)%Tmax,       &
       !                        data%vegdata(veg_i)%Tstd,       &
       !                        data%vegdata(veg_i)%theta_resp, &
       !                        data%vegdata(veg_i)%aTn,        &
       !                        data%vegdata(veg_i)%bTn,        &
       !                        data%vegdata(veg_i)%kTn,        &
       !                        temp)


       ! Get the suspended solids function
!      f_SS = fSS_function(data,veg_i,ss)

       ! Get the density limitation function
       f_dens = fD_function(data,veg_i,veg)

      ! Get the final ingestion rate (/ s)
      ! amount grazed in units of mass consumed/mass vegetation/unit time
      IF(data%vegdata(veg_i)%Ing==1) THEN
        W = (0.071/1000.) * data%vegdata(veg_i)%Length**2.8
        Imax = data%vegdata(veg_i)%WaI * W** data%vegdata(veg_i)%WbI
      ELSE
        Imax = data%vegdata(veg_i)%Rgrz
      END IF
!     grazing = Imax * fGrazing_Limitation * f_Temp * f_dens * f_SS

      ! Now determine available prey and limit grazing amount to availability of prey
      ! food is total amount of food in units of mass/unit volume/unit time
!     food = grazing * veg
!     IF (Ctotal_prey < data%vegdata(veg_i)%num_prey * data%vegdata(veg_i)%Cmin_grz ) THEN
!         food = zero_
!         grazing = zero_
!     ELSEIF (food > Ctotal_prey - data%vegdata(veg_i)%num_prey * data%vegdata(veg_i)%Cmin_grz ) THEN
!         food = Ctotal_prey - data%vegdata(veg_i)%num_prey * data%vegdata(veg_i)%Cmin_grz
!         grazing = food / veg
!     ENDIF


      ! Now determine prey composition based on preference factors and availability of prey

      ! Prey has been ordered in grazing preference
      ! So take food in order of preference up to availability minus value of minimum residual
      ! grazing_prey is in units of mass consumed/unit volumne/unit time

!     DO prey_i = 1,data%vegdata(veg_i)%num_prey
!         !Add up preferences for remaining prey
!         pref_factor = zero_
!         DO prey_j = prey_i,data%vegdata(veg_i)%num_prey
!            pref_factor = pref_factor + data%vegdata(veg_i)%prey(veg_i)%Pveg_prey
!         ENDDO
!         IF (food * data%vegdata(veg_i)%prey(prey_i)%Pveg_prey / pref_factor <= &
!                       prey(prey_i) - data%vegdata(veg_i)%Cmin_grz) THEN
!            !Take fraction of left over food based on preference factor
!            grazing_prey(prey_i) = food * data%vegdata(veg_i)%prey(prey_i)%Pveg_prey / pref_factor
!         ELSEIF (prey(prey_i) > data%vegdata(veg_i)%Cmin_grz) THEN
!            grazing_prey(prey_i) = prey(prey_i) - data%vegdata(veg_i)%Cmin_grz
!         ELSE
!            grazing_prey(prey_i) = zero_
!         ENDIF
!         !Food remaining after grazing from current prey
!         food = food - grazing_prey(prey_i)
!     ENDDO

!     ! Now determine nutrient composition of food based on prey type
!     ! At this stage only the AED model state variables have multiple
!     ! nutrients (C,N&P) so assume all others have a single nutrient
!     ! and thus not need to calculate nutrient excretion as is taken
!     ! care of in the respiration term.  22/12/2011
!     ! grazing_n is in units of mass N consumed/unit volume/unit time
!     ! grazing_p is in units of mass P consumed/unit volume/unit time

!     grazing_n = zero_
!     grazing_p = zero_
!     phy_i = 0
!     DO prey_i = 1,data%vegdata(veg_i)%num_prey
!        IF (data%vegdata(veg_i)%prey(prey_i)%vegetation_prey .EQ. _OGMPOC_) THEN
!           IF (poc > zero_) THEN
!               grazing_n = grazing_n + grazing_prey(prey_i) * pon/poc
!               grazing_p = grazing_p + grazing_prey(prey_i) * pop/poc
!           ELSE
!               grazing_n = zero_
!               grazing_p = zero_
!           ENDIF
!        ELSEIF (data%vegdata(veg_i)%prey(prey_i)%vegetation_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
!           phy_i = phy_i + 1
!           phy_INcon(phy_i) = _STATE_VAR_(data%vegdata(veg_i)%id_phyIN(phy_i))
!           phy_IPcon(phy_i) = _STATE_VAR_(data%vegdata(veg_i)%id_phyIP(phy_i))
!           grazing_n = grazing_n + grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i) /14.0
!           grazing_p = grazing_p + grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i) /31.0
!        ELSEIF (data%vegdata(veg_i)%prey(prey_i)%vegetation_prey(1:15).EQ.'aed2_vegetation') THEN
!           grazing_n = grazing_n + grazing_prey(prey_i) * data%vegdata(veg_i)%INC
!           grazing_p = grazing_p + grazing_prey(prey_i) * data%vegdata(veg_i)%IPC
!        ENDIF
!     ENDDO



!     psuedofaeces = (one_ - data%vegdata(veg_i)%fassim) * grazing
!     ingestion = data%vegdata(veg_i)%fassim * grazing

      IF (veg <= data%vegdata(veg_i)%min) THEN
        ! Don't excrete or die if we are at the min biomass otherwise we have a
        ! mass conservation leak in the C mass balance

        respiration = zero_
        mortality   = zero_
        excretion   = zero_
        egestion    = zero_

      ELSE

!       egestion = data%vegdata(veg_i)%Regst * exp(data%vegdata(veg_i)%gegst + &
!                  min(Ctotal_prey/data%vegdata(veg_i)%Kgrz,one_)) * ingestion

        ! Get the respiration rate (/ s)
        iteg = ingestion - egestion
        respiration = aed2_vegetation_respiration(data,veg_i,iteg,temp,salinity)

        ! Get the excretion rate (of carbon!) (/s)
        excretion =  data%vegdata(veg_i)%Rexcr * iteg

        ! Get the mortality rate (/ s)
!       mortality = data%vegdata(veg_i)%Rmort * f_DO(data,veg_i,oxy)

        ! Add the predation losses to mortality
        mortality = mortality + data%vegdata(veg_i)%Rpred


      ENDIF


      ! Calculate losses into the particulate organic matter pools - Units mmol/s
      poc_excr = (psuedofaeces + egestion + mortality)*veg

!     pon_excr = (psuedofaeces * grazing_n / grazing  &
!              +  egestion*data%vegdata(veg_i)%INC + mortality*data%vegdata(veg_i)%INC)*veg

!     pop_excr = (psuedofaeces * grazing_p / grazing  &
!              +  egestion*data%vegdata(veg_i)%IPC + mortality*data%vegdata(veg_i)%IPC)*veg

      ! Now we know the rates of carbon consumption and excretion, calculate rates
      ! of n & p excretion to maintain internal nutrient stores

      ! First, compute rate of change so far of vegetation carbon biomass (mmolC / m2 /s)
      delta_C = (ingestion - respiration - egestion - excretion - mortality) * veg


      ! Then calc nutrient excretion require to balance internal nutrient store
      ! Note pon_excr includes loss due to messy feeding so no need to include assimilation fraction on grazing_n & grazing_p
!     don_excr = grazing_n - pon_excr - delta_C * data%vegdata(veg_i)%INC
!     dop_excr = grazing_p - pop_excr - delta_C * data%vegdata(veg_i)%IPC

      ! If nutrients are limiting then must excrete doc to maintain balance
!     IF ((don_excr < zero_) .AND. (dop_excr < zero_)) THEN
!        !Determine which nutrient is more limiting
!        IF ((data%vegdata(veg_i)%INC * (grazing_n - pon_excr) - delta_C) .GT. &
!           (data%vegdata(veg_i)%IPC * (grazing_p - pop_excr) - delta_C)) THEN
!            don_excr = zero_
!            doc_excr =  (grazing_n - pon_excr) / data%vegdata(veg_i)%INC - delta_C
!            delta_C = delta_C - doc_excr
!            dop_excr = grazing_p - pop_excr - delta_C*data%vegdata(veg_i)%IPC
!        ELSE
!            dop_excr = zero_
!            doc_excr = (grazing_p - pop_excr) / data%vegdata(veg_i)%IPC - delta_C
!            delta_C = delta_C - doc_excr
!            don_excr = grazing_n - pon_excr - delta_C*data%vegdata(veg_i)%INC
!        ENDIF
!     ELSEIF (don_excr < zero_) THEN !nitrogen limited
!        don_excr = zero_
!        doc_excr = (grazing_n - pon_excr) / data%vegdata(veg_i)%INC - delta_C
!        delta_C = delta_C - doc_excr
!        dop_excr = grazing_p - pop_excr - delta_C*data%vegdata(veg_i)%IPC
!     ELSEIF (dop_excr < zero_) THEN !phosphorus limited
!        dop_excr = zero_
!        doc_excr = (grazing_p - pop_excr) / data%vegdata(veg_i)%IPC - delta_C
!        delta_C = delta_C - doc_excr
!        don_excr = grazing_n - pon_excr - delta_C*data%vegdata(veg_i)%INC
!     ELSE !just excrete nutrients no need to balance c
!         doc_excr = zero_
!     ENDIF


      !write(*,"(4X,'limitations (f_T,f_Salinity): ',2F8.2)")f_T,f_Salinity
      !write(*,"(4X,'sources/sinks (grazing,respiration,mortaility): ',3F8.2)")grazing,excretion,mortality


      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Production / losses in mmolC/s

!print*,'POOP',ingestion,respiration,excretion,egestion,mortality
      _FLUX_VAR_B_(data%id_veg(veg_i)) = _FLUX_VAR_B_(data%id_veg(veg_i)) + ((ingestion - respiration - excretion - egestion - mortality)*veg )


!     IF (data%preyFeedback) THEN
!     ! Now take food grazed by zooplankton from food pools in mmolC/s
!     phy_i = 0
!     DO prey_i = 1,data%vegdata(veg_i)%num_prey
!        _FLUX_VAR_(data%vegdata(veg_i)%id_prey(prey_i)) = _FLUX_VAR_(data%vegdata(veg_i)%id_prey(prey_i)) + ( -1.0 * grazing_prey(prey_i))
!         IF (data%vegdata(veg_i)%prey(prey_i)%vegetation_prey .EQ. _OGMPOC_) THEN
!             IF (poc > zero_) THEN
!                _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + ( -1.0 * grazing_prey(prey_i) * pon/poc)
!                _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + ( -1.0 * grazing_prey(prey_i) * pop/poc)
!             ENDIF
!         ELSEIF (data%vegdata(veg_i)%prey(prey_i)%vegetation_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
!           phy_i = phy_i + 1
!           _FLUX_VAR_(data%vegdata(veg_i)%id_phyIN(phy_i)) = _FLUX_VAR_(data%vegdata(veg_i)%id_phyIN(phy_i)) + ( -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i))
!           _FLUX_VAR_(data%vegdata(veg_i)%id_phyIP(phy_i)) = _FLUX_VAR_(data%vegdata(veg_i)%id_phyIP(phy_i)) + ( -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i))
!        ENDIF
!     ENDDO
!    ENDIF


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
!     _DIAG_VAR_S_(data%id_grz  ) = grazing*secs_pr_day
      _DIAG_VAR_S_(data%id_resp ) = respiration*secs_pr_day
      _DIAG_VAR_S_(data%id_mort ) = mortality*secs_pr_day

      _DIAG_VAR_S_(data%id_excr ) = excretion*secs_pr_day
      _DIAG_VAR_S_(data%id_egst ) = egestion*secs_pr_day
   ENDDO

END SUBROUTINE aed2_calculate_riparian_vegetation
!END SUBROUTINE aed2_calculate_benthic_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
FUNCTION aed2_vegetation_respiration(data,veg_i,iteg,temp,sal) RESULT(resp)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER  :: veg_i ! Invertebrate group
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

     !Make this better
     resp = data%vegdata(veg_i)%Rresp

!   ! Get the salinity limitation.
!   resp = resp * fSalinity_Limitation(data%vegdata,veg_i,sal)


END FUNCTION aed2_vegetation_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





!###############################################################################
PURE AED_REAL FUNCTION fTemp_function_veg(data,veg,temp)
!-------------------------------------------------------------------------------
! Temperature growth multiplier for vegdata
!
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                 :: veg   ! Veg functional group
   AED_REAL, INTENT(IN)                :: temp  ! Temp value being used
!
!LOCALS
   AED_REAL , PARAMETER :: a = 1.00
   AED_REAL             :: MINt,Tmin,Tmax,MAXt

!
!-------------------------------------------------------------------------------
!BEGIN

   MINt = data%vegdata(veg)%MINt
   Tmin = data%vegdata(veg)%Tmin
   Tmax = data%vegdata(veg)%Tmax
   MAXt = data%vegdata(veg)%MAXt

   ! If temp below extreme, temp fn = 0, ie filtration ceases
   IF(temp <= MINt) &
      fTemp_function_veg = zero_

   ! If below min temp for optimal range, limited as below
   IF (temp>MINt .and. temp<Tmin) THEN
      fTemp_function_veg = a*(1.0/(-((Tmin-MINt)*(Tmin-MINt)/ &
            (Tmin*Tmin))+(2.0*(Tmin-MINt)/ &
            Tmin)))*(-((temp-MINt)*(temp-MINt)/&
            (Tmin*Tmin))+2.0*((temp-MINt) / Tmin))
   END IF

   ! If between Tmin and Tmax, ie in the opt temp range, then not limited
   IF (temp >=Tmin .and. temp<=Tmax) THEN
      fTemp_function_veg = a
   END IF

   ! If above max temp for optimal range, limited as below
   IF (temp>Tmax .and. temp <MAXt) THEN
      fTemp_function_veg = a*((-(temp*temp)+2.0*(Tmax*temp)- &
            2.0*(Tmax*MAXt)+(MAXt*MAXt))/ &
            ((Tmax-MAXt)*(Tmax-MAXt)))
   END IF

   ! If temp above extreme, temp fn = 0, ie filtration ceases
   IF (temp>=MAXt) &
       fTemp_function_veg = zero_

 END FUNCTION fTemp_function_veg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





!###############################################################################
 PURE AED_REAL FUNCTION fD_function(data,veg,dens)
!-------------------------------------------------------------------------------
! Density dependence function, for growth limtation
!
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data   ! Module data, with params
   INTEGER, INTENT(in)                 :: veg   ! Veg functional group
   AED_REAL, INTENT(in)                :: dens  ! Density value being used
!
!LOCALS
   AED_REAL             :: maxD, Dmax

!
!-------------------------------------------------------------------------------
!BEGIN

   fD_function = one_

   Dmax =  data%vegdata(veg)%Dmax    ! density where growth decreases
   maxD =  data%vegdata(veg)%maxD    ! density where growth stops

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






END MODULE aed2_vegetation
