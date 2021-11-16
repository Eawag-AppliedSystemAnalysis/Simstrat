!###############################################################################
!#                                                                             #
!# aed_zoop_utils.F90                                                          #
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
!# Created October 2011                                                        #
!#                                                                             #
!###############################################################################


#include "aed.h"

MODULE aed_zoop_utils
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

!  PRIVATE   ! By default make everything private
!
!  %% NAMELIST   %%  zoop_prey_t
   TYPE zoop_prey_t
      !State variable name for zooplankton prey
      CHARACTER(64) :: zoop_prey
      !Preference factors for zooplankton predators grazing on prey
      AED_REAL      :: Pzoo_prey
   END TYPE zoop_prey_t
!  %% END NAMELIST   %%  zoop_prey_t


!  %% NAMELIST   %%  zoop_param_t
   TYPE zoop_param_t
      ! General Attributes
      CHARACTER(64) :: zoop_name
      AED_REAL :: zoop_initial, min_zoo
      ! Growth rate parameters
      AED_REAL :: Rgrz_zoo, fassim_zoo, Kgrz_zoo, theta_grz_zoo
      ! Respiration, mortaility and excretion parameters
      AED_REAL :: Rresp_zoo, Rmort_zoo, ffecal_zoo, fexcr_zoo, ffecal_sed
      ! Temperature limitation on zooplankton loss terms
      AED_REAL :: theta_resp_zoo, Tstd_zoo, Topt_zoo, Tmax_zoo
      ! Salinity parameters
      INTEGER  :: saltfunc_zoo
      AED_REAL :: Smin_zoo, Smax_zoo, Sint_zoo
      ! Nutrient parameters
       AED_REAL :: INC_zoo, IPC_zoo
      ! Dissolved oxygen parameters
      AED_REAL :: DOmin_zoo
      ! Minumum prey concentration parameters
      AED_REAL :: Cmin_grz_zoo
      ! Prey information
      INTEGER  :: num_prey
      TYPE(zoop_prey_t)      :: prey(MAX_ZOOP_PREY)
      INTEGER  :: simDOlim
      ! Temperature limitation derived terms
      AED_REAL :: kTn, aTn, bTn
   END TYPE zoop_param_t
!  %% END NAMELIST   %%  zoop_param_t

   TYPE,extends(zoop_param_t) :: zoop_data_t
      INTEGER  :: id_prey(MAX_ZOOP_PREY)
      INTEGER  :: id_phyIN(MAX_ZOOP_PREY),id_phyIP(MAX_ZOOP_PREY)
   END TYPE

CONTAINS
!===============================================================================



!###############################################################################
FUNCTION fPrey_Limitation(zoops,group,C) RESULT(fPlim)
!----------------------------------------------------------------------------!
! Michaelis-Menten formulation for zooplankton grazing on available
! prey is applied.
!----------------------------------------------------------------------------!
!ARGUMENTS
   TYPE(zoop_data_t),DIMENSION(:),INTENT(in) :: zoops
   INTEGER  :: group
   AED_REAL,INTENT(in)                          :: C !total concentration of available prey
!
!LOCALS
   ! Returns the M-M limitation function
   AED_REAL                                 :: fPlim
!
!-------------------------------------------------------------------------------
!BEGIN
   fPlim = 1.0

   fPlim = C/(zoops(group)%Kgrz_zoo+C)

   IF( fPlim<zero_ ) fPlim=zero_

 END FUNCTION fPrey_Limitation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION fSalinity_Limitation(zoops,group,S) RESULT(fSal)
!----------------------------------------------------------------------------!
! Salinity tolerance of zooplankton                                          !
!----------------------------------------------------------------------------!
!ARGUMENTS
   TYPE(zoop_data_t),DIMENSION(:),INTENT(in) :: zoops
   INTEGER  :: group
   AED_REAL,INTENT(in)                          :: S
!
!LOCALS
   AED_REAL  :: fSal ! Returns the salinity function
   AED_REAL  :: Smin,Smax,Sint
!
!-------------------------------------------------------------------------------
!BEGIN
   Smin = zoops(group)%Smin_zoo
   Smax = zoops(group)%Smax_zoo
   Sint = zoops(group)%Sint_zoo

   !Salinity factor represents natural mortality in response to salinity stress.
   ! f(S) = 1 at Smin<=S<=Smax, f(S) = Bep at S=0 & S=2*Smax.
   IF (S < Smin) THEN
      fSal = (Sint-1.0)/(Smin**2.0)*(S**2.0) - 2*(Sint-1.0)/Smin*S + Sint
   ELSEIF(S > Smax) THEN
      fSal = (Sint-1.0)/(Smax**2.0)*(S**2.0) - 2*(Sint-1.0)/Smax*S + Sint
   ELSE
      fSal = 1.0
   ENDIF

   IF( fSal<zero_ ) fSal=zero_
END FUNCTION fSalinity_Limitation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_zoop_utils
