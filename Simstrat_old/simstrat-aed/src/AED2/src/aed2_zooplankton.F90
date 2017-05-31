!###############################################################################
!#                                                                             #
!# aed2_zooplankton.F90                                                        #
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
!# Created October 2011                                                        #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#define _PHYLEN_ 3
#define _PHYMOD_ 'PHY'
#define _OGMPOC_ 'OGM'

MODULE aed2_zooplankton
!-------------------------------------------------------------------------------
!  aed2_zooplankton --- multi zooplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,ONLY : find_free_lun,aed2_bio_temp_function, fTemp_function,qsort
   USE aed2_zoop_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_zooplankton_data_t
!

   TYPE,extends(aed2_model_data_t) :: aed2_zooplankton_data_t
      !# Variable identifiers
      INTEGER  :: id_zoo(MAX_ZOOP_TYPES)
      INTEGER  :: id_Nexctarget,id_Nmorttarget
      INTEGER  :: id_Pexctarget,id_Pmorttarget
      INTEGER  :: id_Cexctarget,id_Cmorttarget
      INTEGER  :: id_DOupttarget
      INTEGER  :: id_tem, id_sal, id_extc
      INTEGER  :: id_grz,id_resp,id_mort


      !# Model parameters
      INTEGER  :: num_zoops
      TYPE(type_zoop_data),DIMENSION(:),ALLOCATABLE :: zoops
      LOGICAL  :: simDNexcr, simDPexcr, simDCexcr
      LOGICAL  :: simPNexcr, simPPexcr, simPCexcr

     CONTAINS
         PROCEDURE :: define            => aed2_define_zooplankton
         PROCEDURE :: calculate         => aed2_calculate_zooplankton
!        PROCEDURE :: mobility          => aed2_mobility_zooplankton
!        PROCEDURE :: light_extinction  => aed2_light_extinction_zooplankton
!        PROCEDURE :: delete            => aed2_delete_zooplankton

   END TYPE

   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed2_zooplankton_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_zooplankton_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count !Number of zooplankton groups
   INTEGER,INTENT(in)          :: list(*) !List of zooplankton groups to simulate
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   AED_REAL :: Pzoo_prey(MAX_ZOOP_PREY)

   TYPE(type_zoop_params)  :: zoop_param(MAX_ZOOP_TYPES)
   NAMELIST /zoop_params/ zoop_param
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD',iostat=status)
    IF (status /= 0) STOP 'Error opening zoop_params namelist file'
    read(tfil,nml=zoop_params,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist zoop_params'

    data%num_zoops = count
    allocate(data%zoops(count))
    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%zoops(i)%zoop_name         = zoop_param(list(i))%zoop_name
       data%zoops(i)%zoop_initial      = zoop_param(list(i))%zoop_initial
       data%zoops(i)%min_zoo           = zoop_param(list(i))%min_zoo
       data%zoops(i)%Rgrz_zoo          = zoop_param(list(i))%Rgrz_zoo/secs_per_day
       data%zoops(i)%fassim_zoo        = zoop_param(list(i))%fassim_zoo
       data%zoops(i)%Kgrz_zoo          = zoop_param(list(i))%Kgrz_zoo
       data%zoops(i)%theta_grz_zoo     = zoop_param(list(i))%theta_grz_zoo
       data%zoops(i)%Rresp_zoo         = zoop_param(list(i))%Rresp_zoo/secs_per_day
       data%zoops(i)%Rmort_zoo         = zoop_param(list(i))%Rmort_zoo/secs_per_day
       data%zoops(i)%ffecal_zoo        = zoop_param(list(i))%ffecal_zoo
       data%zoops(i)%fexcr_zoo         = zoop_param(list(i))%fexcr_zoo
       data%zoops(i)%ffecal_sed        = zoop_param(list(i))%ffecal_sed
       data%zoops(i)%theta_resp_zoo    = zoop_param(list(i))%theta_resp_zoo
       data%zoops(i)%Tstd_zoo          = zoop_param(list(i))%Tstd_zoo
       data%zoops(i)%Topt_zoo          = zoop_param(list(i))%Topt_zoo
       data%zoops(i)%Tmax_zoo          = zoop_param(list(i))%Tmax_zoo
       data%zoops(i)%saltfunc_zoo      = zoop_param(list(i))%saltfunc_zoo
       data%zoops(i)%Smin_zoo          = zoop_param(list(i))%Smin_zoo
       data%zoops(i)%Smax_zoo          = zoop_param(list(i))%Smax_zoo
       data%zoops(i)%Sint_zoo          = zoop_param(list(i))%Sint_zoo
       data%zoops(i)%INC_zoo           = zoop_param(list(i))%INC_zoo
       data%zoops(i)%IPC_zoo           = zoop_param(list(i))%IPC_zoo
       data%zoops(i)%simDOlim          = zoop_param(list(i))%simDOlim
       data%zoops(i)%DOmin_zoo         = zoop_param(list(i))%DOmin_zoo
       data%zoops(i)%Cmin_grz_zoo      = zoop_param(list(i))%Cmin_grz_zoo
       data%zoops(i)%num_prey          = zoop_param(list(i))%num_prey
       !Loop through prey variables assigning a target variable and preference factor
       !First sort in decending order of food preferences
       DO j=1,data%zoops(i)%num_prey
          sort_i(j) = j
          Pzoo_prey(j) = zoop_param(list(i))%prey(j)%Pzoo_prey
       ENDDO
       CALL qsort(Pzoo_prey,sort_i,1,data%zoops(i)%num_prey)
       DO j=1,data%zoops(i)%num_prey
          data%zoops(i)%prey(j)%zoop_prey = zoop_param(list(i))%prey(sort_i(data%zoops(i)%num_prey-j+1))%zoop_prey
          data%zoops(i)%prey(j)%Pzoo_prey = zoop_param(list(i))%prey(sort_i(data%zoops(i)%num_prey-j+1))%Pzoo_prey
       ENDDO

       ! Register group as a state variable
       data%id_zoo(i) = aed2_define_variable(           &
                              zoop_param(list(i))%zoop_name,       &
                              'mmolC/m**3', 'zooplankton',         &
                              zoop_param(list(i))%zoop_initial,    &
                              minimum=zoop_param(list(i))%min_zoo)
    ENDDO
!
END SUBROUTINE aed2_zooplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_zooplankton(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the zooplankton biogeochemical model
!
!  Here, the aed2_zooplankton namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_zooplankton_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_zoops
   INTEGER  :: the_zoops(MAX_ZOOP_TYPES)

   CHARACTER(len=64)  :: dn_target_variable='' !dissolved nitrogen target variable
   CHARACTER(len=64)  :: pn_target_variable='' !particulate nitrogen target variable
   CHARACTER(len=64)  :: dp_target_variable='' !dissolved phosphorus target variable
   CHARACTER(len=64)  :: pp_target_variable='' !particulate phosphorus target variable
   CHARACTER(len=64)  :: dc_target_variable='' !dissolved carbon target variable
   CHARACTER(len=64)  :: pc_target_variable='' !particulate carbon target variable
   CHARACTER(len=128) :: dbase='aed2_zoop_pars.nml'

   INTEGER  :: zoop_i, prey_i, phy_i

   NAMELIST /aed2_zooplankton/ num_zoops, the_zoops, &
                    dn_target_variable, pn_target_variable, dp_target_variable, &
                    pp_target_variable, dc_target_variable, pc_target_variable, &
                    dbase
!-----------------------------------------------------------------------
!BEGIN
!print *,'**** Reading /aed2_zooplankton/ namelist'
   ! Read the namelist
   read(namlst,nml=aed2_zooplankton,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_zooplankton'

    data%num_zoops = 0
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed2_zooplankton_load_params to values per second.
   CALL aed2_zooplankton_load_params(data, dbase, num_zoops, the_zoops)


   CALL aed2_bio_temp_function(data%num_zoops,            &
                              data%zoops%theta_resp_zoo, &
                              data%zoops%Tstd_zoo,       &
                              data%zoops%Topt_zoo,       &
                              data%zoops%Tmax_zoo,       &
                              data%zoops%aTn,            &
                              data%zoops%bTn,            &
                              data%zoops%kTn,            &
                              data%zoops%zoop_name)


   !Register link to prey state variables
   DO zoop_i = 1,num_zoops
      phy_i = 0
      DO prey_i = 1,data%zoops(zoop_i)%num_prey
          data%zoops(zoop_i)%id_prey(prey_i) = aed2_locate_variable( &
                                       data%zoops(zoop_i)%prey(prey_i)%zoop_prey)
          !If the zooplankton prey is phytoplankton then also register state dependency on
          !internal nitrogen and phosphorus
          IF (data%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
              phy_i = phy_i + 1
              data%zoops(zoop_i)%id_phyIN(phy_i) = aed2_locate_variable( &
                                       TRIM(data%zoops(zoop_i)%prey(prey_i)%zoop_prey)//'_IN')
              data%zoops(zoop_i)%id_phyIP(phy_i) = aed2_locate_variable( &
                                       TRIM(data%zoops(zoop_i)%prey(prey_i)%zoop_prey)//'_IP')

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


   ! Register diagnostic variables
   data%id_grz = aed2_define_diag_variable('grz','mmolC/m**3/d',  'net zooplankton grazing')
   data%id_resp = aed2_define_diag_variable('resp','mmolC/m**3/d',  'net zooplankton respiration')
   data%id_mort = aed2_define_diag_variable('mort','mmolC/m**3/d','net zooplankton mortality')

   ! Register environmental dependencies
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_extc = aed2_locate_global('extc_coef')
END SUBROUTINE aed2_define_zooplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_zooplankton(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_zooplankton_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: zoo,temp,salinity !State variables
   AED_REAL           :: prey(MAX_ZOOP_PREY), grazing_prey(MAX_ZOOP_PREY) !Prey state variables
   AED_REAL           :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for phytoplankton
   AED_REAL           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL           :: pon, pop, poc !Mortaility and fecal pellet state variables
   AED_REAL           :: FGrazing_Limitation, f_T, f_Salinity
   AED_REAL           :: pref_factor, Ctotal_prey !total concentration of available prey
   AED_REAL           :: food, grazing, respiration, mortality !Growth & decay functions
   AED_REAL           :: grazing_n, grazing_p !Grazing on nutrients
   AED_REAL           :: pon_excr, pop_excr, poc_excr !POM excretion rates
   AED_REAL           :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   INTEGER  :: zoop_i,prey_i,prey_j,phy_i
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)    ! local temperature
   salinity = _STATE_VAR_(data%id_sal)! local salinity

   ! Retrieve current (local) state variable values.
   IF (data%simDNexcr)  dn_excr = _STATE_VAR_(data%id_Nexctarget)
   IF (data%simDPexcr)  dp_excr = _STATE_VAR_(data%id_Pexctarget)
   IF (data%simDCexcr)  dc_excr = _STATE_VAR_(data%id_Cexctarget)

   IF (data%simPNexcr)  pon = _STATE_VAR_(data%id_Nmorttarget)
   IF (data%simPPexcr)  pop = _STATE_VAR_(data%id_Pmorttarget)
   IF (data%simPCexcr)  poc = _STATE_VAR_(data%id_Cmorttarget)

   DO zoop_i=1,data%num_zoops

      ! Retrieve this zooplankton group
      zoo = _STATE_VAR_(data%id_zoo(zoop_i))
      !Retrieve prey groups
      Ctotal_prey   = zero_
      DO prey_i=1,data%zoops(zoop_i)%num_prey
         prey(prey_i) = _STATE_VAR_(data%zoops(zoop_i)%id_prey(prey_i))
         Ctotal_prey = Ctotal_prey + prey(prey_i)
      ENDDO

      grazing       = zero_
      respiration   = zero_
      mortality     = zero_

      ! Get the grazing limitation function
       fGrazing_Limitation = fPrey_Limitation(data%zoops,zoop_i,Ctotal_prey)

      ! Get the temperature function
       f_T = fTemp_function(1, data%zoops(zoop_i)%Tmax_zoo,       &
                               data%zoops(zoop_i)%Tstd_zoo,       &
                               data%zoops(zoop_i)%theta_resp_zoo, &
                               data%zoops(zoop_i)%aTn,            &
                               data%zoops(zoop_i)%bTn,            &
                               data%zoops(zoop_i)%kTn, temp)

      ! Get the salinity limitation.
       f_Salinity = fSalinity_Limitation(data%zoops,zoop_i,salinity)

      ! Get the growth rate (/ s)
      ! grazing is in units of mass consumed/mass zoops/unit time
      grazing = data%zoops(zoop_i)%Rgrz_zoo * fGrazing_Limitation * f_T

      ! Now determine available prey and limit grazing amount to
      ! availability of prey
      ! food is total amount of food in units of mass/unit volume/unit time
      food = grazing * zoo
      IF (Ctotal_prey < data%zoops(zoop_i)%num_prey * data%zoops(zoop_i)%Cmin_grz_zoo ) THEN
          food = zero_
          grazing = food / zoo
      ELSEIF (food > Ctotal_prey - data%zoops(zoop_i)%num_prey * data%zoops(zoop_i)%Cmin_grz_zoo ) THEN
          food = Ctotal_prey - data%zoops(zoop_i)%num_prey * data%zoops(zoop_i)%Cmin_grz_zoo
          grazing = food / zoo
      ENDIF


      ! Now determine prey composition based on preference factors and
      ! availability of prey

      ! Prey has been ordered in grazing preference
      ! So take food in order of preference up to availability minus
      !value of minimum residual
      ! grazing_prey is in units of mass consumed/unit volumne/unit time

      DO prey_i = 1,data%zoops(zoop_i)%num_prey
          !Add up preferences for remaining prey
          pref_factor = zero_
          DO prey_j = prey_i,data%zoops(zoop_i)%num_prey
             pref_factor = pref_factor + data%zoops(zoop_i)%prey(prey_j)%Pzoo_prey
          ENDDO
          IF (food * data%zoops(zoop_i)%prey(prey_i)%Pzoo_prey / pref_factor <= &
                        prey(prey_i) - data%zoops(zoop_i)%Cmin_grz_zoo) THEN
             !Take fraction of left over food based on preference factor
             grazing_prey(prey_i) = food * data%zoops(zoop_i)%prey(prey_i)%Pzoo_prey / pref_factor
          ELSEIF (prey(prey_i) > data%zoops(zoop_i)%Cmin_grz_zoo) THEN
             grazing_prey(prey_i) = prey(prey_i) - data%zoops(zoop_i)%Cmin_grz_zoo
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
      DO prey_i = 1,data%zoops(zoop_i)%num_prey
         IF (data%zoops(zoop_i)%prey(prey_i)%zoop_prey .EQ. _OGMPOC_) THEN
            IF (poc > zero_) THEN
                grazing_n = grazing_n + grazing_prey(prey_i) * pon/poc
                grazing_p = grazing_p + grazing_prey(prey_i) * pop/poc
            ELSE
                grazing_n = zero_
                grazing_p = zero_
            ENDIF
         ELSEIF (data%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
            phy_i = phy_i + 1
            phy_INcon(phy_i) = _STATE_VAR_(data%zoops(zoop_i)%id_phyIN(phy_i))
            phy_IPcon(phy_i) = _STATE_VAR_(data%zoops(zoop_i)%id_phyIP(phy_i))
            grazing_n = grazing_n + grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i) /14.0
            grazing_p = grazing_p + grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i) /31.0
         ELSEIF (data%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:15).EQ.'aed2_zooplankton') THEN
            grazing_n = grazing_n + grazing_prey(prey_i) * data%zoops(zoop_i)%INC_zoo
            grazing_p = grazing_p + grazing_prey(prey_i) * data%zoops(zoop_i)%IPC_zoo
         ENDIF
      ENDDO


      ! Get the respiration rate (/ s)
      respiration = data%zoops(zoop_i)%Rresp_zoo * f_Salinity

      ! Get the mortality rate (/ s)
      mortality = data%zoops(zoop_i)%Rmort_zoo * f_T

      ! Don't excrete or die if we are at the min biomass otherwise we have a
      ! mass conservation leak in the C mass balance
      IF (zoo <= data%zoops(zoop_i)%min_zoo) THEN
        respiration = zero_
        mortality = zero_
      ENDIF

      ! Now we know the rates of carbon consumption and excretion,
      ! calculate rates of n & p excretion to maintain internal
      ! nutrient stores
      ! Calculate excretion of particulate organic matter - Units mmol/s
      poc_excr = ((1 - data%zoops(zoop_i)%fassim_zoo) * grazing + &
               (1 - data%zoops(zoop_i)%ffecal_sed) * data%zoops(zoop_i)%ffecal_zoo * respiration +  &
                                                        mortality) * zoo
      pon_excr = (1 - data%zoops(zoop_i)%fassim_zoo) * grazing_n  + &
               ((1 - data%zoops(zoop_i)%ffecal_sed) * data%zoops(zoop_i)%ffecal_zoo * respiration + &
                              mortality) * data%zoops(zoop_i)%INC_zoo * zoo
      pop_excr = (1 - data%zoops(zoop_i)%fassim_zoo) * grazing_p + &
               ((1 - data%zoops(zoop_i)%ffecal_sed) * data%zoops(zoop_i)%ffecal_zoo * respiration + &
                              mortality) * data%zoops(zoop_i)%IPC_zoo * zoo


      ! Calculate rate of change of zooplankton carbon (mmolC/s)          !
      delta_C = (data%zoops(zoop_i)%fassim_zoo * grazing - respiration - mortality) * zoo
      ! Calculate nutrient excretion require to balance internal nutrient store
      ! Note pon_excr includes loss due to messy feeding so no need to include assimilation fraction on grazing_n & grazing_p
      don_excr = grazing_n - pon_excr - delta_C * data%zoops(zoop_i)%INC_zoo
      dop_excr = grazing_p - pop_excr - delta_C * data%zoops(zoop_i)%IPC_zoo
      !If nutrients are limiting then must excrete doc to maintain balance
      IF ((don_excr < zero_) .AND. (dop_excr < zero_)) THEN
         !Determine which nutrient is more limiting
         IF ((data%zoops(zoop_i)%INC_zoo * (grazing_n - pon_excr) - delta_C) .GT. &
            (data%zoops(zoop_i)%IPC_zoo * (grazing_p - pop_excr) - delta_C)) THEN
             don_excr = zero_
             doc_excr =  (grazing_n - pon_excr) / data%zoops(zoop_i)%INC_zoo - delta_C
             delta_C = delta_C - doc_excr
             dop_excr = grazing_p - pop_excr - delta_C*data%zoops(zoop_i)%IPC_zoo
         ELSE
             dop_excr = zero_
             doc_excr = (grazing_p - pop_excr) / data%zoops(zoop_i)%IPC_zoo - delta_C
             delta_C = delta_C - doc_excr
             don_excr = grazing_n - pon_excr - delta_C*data%zoops(zoop_i)%INC_zoo
         ENDIF
      ELSEIF (don_excr < zero_) THEN !nitrogen limited
         don_excr = zero_
         doc_excr = (grazing_n - pon_excr) / data%zoops(zoop_i)%INC_zoo - delta_C
         delta_C = delta_C - doc_excr
         dop_excr = grazing_p - pop_excr - delta_C*data%zoops(zoop_i)%IPC_zoo
      ELSEIF (dop_excr < zero_) THEN !phosphorus limited
         dop_excr = zero_
         doc_excr = (grazing_p - pop_excr) / data%zoops(zoop_i)%IPC_zoo - delta_C
         delta_C = delta_C - doc_excr
         don_excr = grazing_n - pon_excr - delta_C*data%zoops(zoop_i)%INC_zoo
      ELSE !just excrete nutrients no need to balance c
          doc_excr = zero_
      ENDIF


      !write(*,"(4X,'limitations (f_T,f_Salinity): ',2F8.2)")f_T,f_Salinity
      !write(*,"(4X,'sources/sinks (grazing,respiration,mortaility): ',3F8.2)")grazing,excretion,mortality


      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Zooplankton production / losses in mmolC/s

      _FLUX_VAR_(data%id_zoo(zoop_i)) = _FLUX_VAR_(data%id_zoo(zoop_i)) + ( (data%zoops(zoop_i)%fassim_zoo * grazing - respiration - mortality)*zoo )


      ! Now take food grazed by zooplankton from food pools in mmolC/s
      phy_i = 0
      DO prey_i = 1,data%zoops(zoop_i)%num_prey
         _FLUX_VAR_(data%zoops(zoop_i)%id_prey(prey_i)) = _FLUX_VAR_(data%zoops(zoop_i)%id_prey(prey_i)) + ( -1.0 * grazing_prey(prey_i))
          IF (data%zoops(zoop_i)%prey(prey_i)%zoop_prey .EQ. _OGMPOC_) THEN
              IF (poc > zero_) THEN
                 _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + ( -1.0 * grazing_prey(prey_i) * pon/poc)
                 _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + ( -1.0 * grazing_prey(prey_i) * pop/poc)
              ENDIF
          ELSEIF (data%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
            phy_i = phy_i + 1
            _FLUX_VAR_(data%zoops(zoop_i)%id_phyIN(phy_i)) = _FLUX_VAR_(data%zoops(zoop_i)%id_phyIN(phy_i)) + ( -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i))
            _FLUX_VAR_(data%zoops(zoop_i)%id_phyIP(phy_i)) = _FLUX_VAR_(data%zoops(zoop_i)%id_phyIP(phy_i)) + ( -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i))
         ENDIF
      ENDDO


      ! Now manage excretion contributions to DOM
      IF (data%simDCexcr) THEN
         _FLUX_VAR_(data%id_Cexctarget) = _FLUX_VAR_(data%id_Cexctarget) + (data%zoops(zoop_i)%fexcr_zoo * respiration * zoo + doc_excr)
      ENDIF
      IF (data%simDNexcr) THEN
         _FLUX_VAR_(data%id_Nexctarget) = _FLUX_VAR_(data%id_Nexctarget) + (don_excr)
      ENDIF
      IF (data%simDPexcr) THEN
         _FLUX_VAR_(data%id_Pexctarget) = _FLUX_VAR_(data%id_Pexctarget) + (dop_excr)
      ENDIF

      ! Now manage messy feeding, fecal pellets and mortality contributions to POM
      IF (data%simPCexcr) THEN
         _FLUX_VAR_(data%id_Cmorttarget) = _FLUX_VAR_(data%id_Cmorttarget) + ( poc_excr)
      ENDIF
      IF (data%simPNexcr) THEN
         _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + ( pon_excr)
      ENDIF
      IF (data%simPPexcr) THEN
         _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + ( pop_excr)
      ENDIF

      ! Export diagnostic variables
      _DIAG_VAR_(data%id_grz ) = zoo*grazing*secs_per_day
      _DIAG_VAR_(data%id_resp ) = zoo*respiration*secs_per_day
      _DIAG_VAR_(data%id_mort ) = zoo*mortality*secs_per_day

   ENDDO

END SUBROUTINE aed2_calculate_zooplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_zooplankton
