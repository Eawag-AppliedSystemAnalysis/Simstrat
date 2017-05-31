!###############################################################################
!#                                                                             #
!# aed2_chlorophylla.F90                                                       #
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
!# Created August 2011                                                         #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_chlorophylla
!-------------------------------------------------------------------------------
! aed2_chlorphylla --- simple lumped chl-a model
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_chla_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_chla_data_t
      !# Variable identifiers
      INTEGER  :: id_p
      INTEGER  :: id_exctarget,id_morttarget,id_upttarget
      INTEGER  :: id_par
      INTEGER  :: id_I_0
      INTEGER  :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR

      !# Model parameters
      AED_REAL :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      AED_REAL :: dic_per_n
      LOGICAL  :: do_exc,do_mort,do_upt

     CONTAINS
         PROCEDURE :: define            => aed2_define_chla
         PROCEDURE :: calculate         => aed2_calculate_chla
!        PROCEDURE :: mobility          => aed2_mobility_chla
         PROCEDURE :: light_extinction  => aed2_light_extinction_chla
!        PROCEDURE :: delete            => aed2_delete_chla

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_chla(data,  namlst)
!-------------------------------------------------------------------------------
! Initialise the chlorophyl model
!
!  Here, the chla namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_chla_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

   AED_REAL           :: p_initial=0.
   AED_REAL           :: p0=0.0225
   AED_REAL           :: w_p=-1.157407e-05
   AED_REAL           :: i_min=25.
   AED_REAL           :: rmax=1.157407e-05
   AED_REAL           :: alpha=0.3
   AED_REAL           :: rpn=1.157407e-07
   AED_REAL           :: rpdu=2.314814e-07
   AED_REAL           :: rpdl=1.157407e-06
   CHARACTER(len=64)  :: excretion_target_variable=''
   CHARACTER(len=64)  :: mortality_target_variable=''
   CHARACTER(len=64)  :: uptake_target_variable=''

   NAMELIST /aed2_chla/ p_initial,p0,w_p,i_min,rmax,alpha,rpn,rpdu,rpdl, &
                    excretion_target_variable,mortality_target_variable,uptake_target_variable

!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed2_chla,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_chla'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%p0    = p0
   data%i_min = i_min
   data%rmax  = rmax/secs_per_day
   data%alpha = alpha
   data%rpn  = rpn /secs_per_day
   data%rpdu = rpdu/secs_per_day
   data%rpdl = rpdl/secs_per_day

   ! Register state variables
   data%id_p = aed2_define_variable('phy','mmol/m**3','phytoplankton', &
                                    p_initial,minimum=zero_,mobility=w_p/secs_per_day)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   data%do_exc = excretion_target_variable .NE. ''
   IF (data%do_exc) data%id_exctarget = aed2_locate_variable(excretion_target_variable)
   data%do_mort = mortality_target_variable .NE. ''
   IF (data%do_mort) data%id_morttarget = aed2_locate_variable(mortality_target_variable)
   data%do_upt = uptake_target_variable .NE. ''
   IF (data%do_upt) data%id_upttarget = aed2_locate_variable(uptake_target_variable)


   ! Register diagnostic variables
   data%id_GPP = aed2_define_diag_variable('GPP','mmol/m**3/d',  'gross primary production')
   data%id_NCP = aed2_define_diag_variable('NCP','mmol/m**3/d',  'net community production')
   data%id_PPR = aed2_define_diag_variable('PPR','-','phytoplankton p/r ratio (gross)')
   data%id_NPR = aed2_define_diag_variable('NPR','-','phytoplankton p/r ratio (net)')
   data%id_dPAR = aed2_define_diag_variable('PAR','W/m**2',    'photosynthetically active radiation')

   ! Register environmental dependencies
   data%id_par = aed2_locate_global('par')
   data%id_I_0 = aed2_locate_global_sheet('par_sf')


   PRINT *,'AED_CHLA : Note this module has not been completed.'
END SUBROUTINE aed2_define_chla
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_chla(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of chlorophylla model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_chla_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: n,p,par,I_0
   AED_REAL           :: iopt,rpd,primprod

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
   p = _STATE_VAR_(data%id_p)! phytoplankton
   n = _STATE_VAR_(data%id_upttarget)! nutrients

   ! Retrieve current environmental conditions.
   par = _STATE_VAR_(data%id_par) ! local photosynthetically active radiation
   I_0 = _STATE_VAR_S_(data%id_I_0) ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,data%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   IF (par .ge. data%I_min) THEN
      rpd = data%rpdu
   ELSE
      rpd = data%rpdl
   ENDIF

   ! Define some intermediate quantities that will be reused multiple times.
   primprod = fnp(data,n,p,par,iopt)

   ! Set temporal derivatives
   _FLUX_VAR_(data%id_p) = _FLUX_VAR_(data%id_p) + (primprod - data%rpn*p - rpd*p)

   ! If an externally maintained ...
   IF (data%do_upt) THEN
      _FLUX_VAR_(data%id_upttarget) = _FLUX_VAR_(data%id_upttarget) + (-primprod)
   ENDIF
   IF (data%do_mort) THEN
      _FLUX_VAR_(data%id_morttarget) = _FLUX_VAR_(data%id_morttarget) + (rpd*p)
   ENDIF
   IF (data%do_exc) THEN
      _FLUX_VAR_(data%id_exctarget) = _FLUX_VAR_(data%id_exctarget) + (data%rpn*p)
   ENDIF

   ! Export diagnostic variables
   _DIAG_VAR_(data%id_dPAR) = par
   _DIAG_VAR_(data%id_GPP ) = primprod*secs_per_day
   _DIAG_VAR_(data%id_NCP ) = (primprod - data%rpn*p)*secs_per_day
   _DIAG_VAR_(data%id_PPR ) = primprod/(data%rpn*p - rpd*p)
   _DIAG_VAR_(data%id_NPR ) = (primprod - data%rpn*p)/(data%rpn*p - rpd*p)


END SUBROUTINE aed2_calculate_chla
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_chla(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_chla_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: p
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
   p = _STATE_VAR_(data%id_p)! phytoplankton

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   extinction = extinction + (0.0*p)


END SUBROUTINE aed2_light_extinction_chla
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fnp(data,n,p,par,iopt)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nutrient uptake
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_chla_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in)                :: n,p,par,iopt
!
!-------------------------------------------------------------------------------
!BEGIN
   fnp = data%rmax*par/iopt*exp(one_-par/iopt)*n/(data%alpha+n)*(p+data%p0)

END FUNCTION fnp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed2_chlorophylla
