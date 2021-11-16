!###############################################################################
!#                                                                             #
!# aed_geochemistry.F90                                                        #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2012 - 2021 -  The University of Western Australia               #
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
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_geochemistry
!-------------------------------------------------------------------------------
! aed_geochemistry --- geochemistry model
!
! The AED module geochemistry contains equations that describe ....
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_gclib,ONLY : AED_GC_Input,printTables
   USE aed_gcsolver,ONLY : ConfigEquilibriumSolver, &
                            GetListOfGeochemDiagnostics, &
                            InitialiseGCProperties, &
                            UpdateEquilibration, &
                            returnGCDerivedVector, &
                            simManganRedox, &
                            simArsenicRedox, &
                            simSulfurRedox, &
                            simCarbonRedox, &
                            simIronRedox


   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_geochemistry_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_geochemistry_data_t
      !# Variable identifiers
      INTEGER  :: id_comp(MAX_GC_COMPONENTS), id_mins(MAX_GC_MINERALS)
      INTEGER  :: id_cdep(MAX_GC_COMPONENTS), id_mdep(MAX_GC_MINERALS)
      INTEGER  :: id_ubalchg ,id_pH, id_c_pco2, id_o_oxy
      INTEGER  :: id_temp,id_sal
      INTEGER  :: id_sed_dic
      INTEGER  :: id_gcdiag(MAX_GC_COMPONENTS), id_noncon
      INTEGER  :: id_totC
      INTEGER  :: id_feii, id_feiii
      INTEGER  :: id_h2s, id_so4

      !# Model parameters
      INTEGER  :: num_comp, num_mins
      LOGICAL  :: component_linked(MAX_GC_COMPONENTS),mineral_linked(MAX_GC_MINERALS)
      LOGICAL  :: simEq
      AED_REAL :: Riron_red, theta_iron_red, Kiron_red
      AED_REAL :: Riron_aox, Riron_box, theta_iron_ox
      AED_REAL :: Rsulf_red, theta_sulf_red, Ksulf_red
      AED_REAL :: Rsulf_ox, theta_sulf_ox, Ksulf_ox
      AED_REAL :: speciation_dt
      AED_REAL :: Fsed_gch(MAX_GC_COMPONENTS)
      AED_REAL :: Ksed_gch_o2(MAX_GC_COMPONENTS)
      AED_REAL :: Ksed_gch_pH(MAX_GC_COMPONENTS)
      AED_REAL :: w_gch(MAX_GC_MINERALS)
      AED_REAL,          DIMENSION(:), ALLOCATABLE :: DissComp,PartComp
      CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: listDissTransVars,listPartTransVars

     CONTAINS
         PROCEDURE :: define            => aed_define_geochemistry
         PROCEDURE :: initialize        => aed_initialize_geochemistry
         PROCEDURE :: calculate         => aed_calculate_geochemistry
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_geochemistry
         PROCEDURE :: equilibrate       => aed_equilibrate_geochemistry
!        PROCEDURE :: mobility          => aed_mobility_geochemistry
!        PROCEDURE :: light_extinction  => aed_light_extinction_geochemistry
!        PROCEDURE :: delete            => aed_delete_geochemistry

   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS







!###############################################################################
SUBROUTINE aed_define_geochemistry(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER           :: status, i

!  %% NAMELIST   %%  /aed_geochemistry/
!  %% Last Checked 20/08/2021
   INTEGER           :: speciation_dt
   INTEGER           :: num_components,num_minerals
   INTEGER           :: nDissTransportables, nPartTransportables
   LOGICAL           :: simEq = .TRUE.
   AED_REAL          :: min
   AED_REAL          :: dis_initial(MAX_GC_COMPONENTS) = 0.0
   AED_REAL          :: Fsed_gch(MAX_GC_COMPONENTS)
   AED_REAL          :: Ksed_gch_o2(MAX_GC_COMPONENTS)
   AED_REAL          :: Ksed_gch_pH(MAX_GC_COMPONENTS)
   AED_REAL          :: min_initial(MAX_GC_MINERALS) = 0.0
   AED_REAL          :: w_gch(MAX_GC_MINERALS)
   AED_REAL          :: pH_initial = 7.5
   AED_REAL          :: pe_initial = 8.0
   AED_REAL          :: Riron_red, theta_iron_red, Kiron_red
   AED_REAL          :: Riron_aox, Riron_box, theta_iron_ox
   AED_REAL          :: Rsulf_red, theta_sulf_red, Ksulf_red
   AED_REAL          :: Rsulf_ox, theta_sulf_ox, Ksulf_ox
   CHARACTER(len=64) :: geochem_file = ''
   CHARACTER(len=64) :: dis_components(MAX_GC_COMPONENTS) = ''
   CHARACTER(len=64) :: component_link(MAX_GC_COMPONENTS) = ''
   CHARACTER(len=64) :: speciesOutput(10) = ''
   CHARACTER(len=64) :: the_minerals(MAX_GC_MINERALS) = ''
   CHARACTER(len=64) :: mineral_link(MAX_GC_MINERALS) = ''
   CHARACTER(len=64) :: ph_link = 'CAR_pH'
   CHARACTER(len=64) :: pco2_link = 'CAR_pCO2'
   CHARACTER(len=64) :: oxy_link = 'OXY_oxy'
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_geochemistry/

   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: diagnosticList

   NAMELIST /aed_geochemistry/ speciation_dt, geochem_file,                   &
                    num_components, dis_components, component_link, Fsed_gch,  &
                    dis_initial, num_minerals, the_minerals, mineral_link,     &
                    w_gch, min_initial, pH_initial, speciesOutput, simEq,      &
                    Riron_red, theta_iron_red, Kiron_red,                      &
                    Riron_aox, Riron_box, theta_iron_ox, &
                    Rsulf_red, theta_sulf_red, Ksulf_red, &
                    Rsulf_ox, theta_sulf_ox, Ksulf_ox, &
                    ph_link, pco2_link, diag_level
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_geochemistry initialization"
   print *,"         WARNING! aed_geochemistry model is under development"

   ! MH:JOBS
   ! remove pe
   ! species outputs

   !----------------------------------------------------------------------------
   ! Initialise variables
   component_link(:) = ''
   data%component_linked(:) = .FALSE.
   data%mineral_linked(:) = .FALSE.

! Initialisation now done in declaration
!  dis_initial = 0.0  ! default, overwritten by namelist
!  min_initial = 0.0  ! default, overwritten by namelist
!  pH_initial  = 7.5  ! default, overwritten by namelist
!  pe_initial  = 8.0  ! default, not used

!  ph_link     = 'CAR_pH'
!  pco2_link   = 'CAR_pCO2'

   !----------------------------------------------------------------------------
   ! Read the namelist
   read(namlst,nml=aed_geochemistry,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_geochemistry'

   data%simEq = simEq
   data%speciation_dt = speciation_dt  ! Note this is now managed in FV_AED2 or GLM_AED2

   data%Riron_red = Riron_red             ; data%Kiron_red= Kiron_red
   data%theta_iron_red= theta_iron_red    ; data%theta_iron_ox=theta_iron_ox
   data%Riron_aox=Riron_aox/ secs_per_day ; data%Riron_box=Riron_box/ secs_per_day

   data%Rsulf_red = Rsulf_red             ; data%Ksulf_red= Ksulf_red
   data%theta_sulf_red= theta_sulf_red    ; data%theta_sulf_ox=theta_sulf_ox
   data%Rsulf_ox=Rsulf_ox / secs_per_day

   speciesOutput = ''
   speciesOutput(1) = 'NONCON'
   !speciesOutput(1) = 'HCO3-'

   !----------------------------------------------------------------------------
   ! Now load the geochem database

   !CALL aed_geochem_load_params(data, geochem_file, modelinfo)
   CALL AED_GC_Input(geochem_file)
   !CALL printTables()

   !----------------------------------------------------------------------------
   ! Configure the geochemical solver
   CALL ConfigEquilibriumSolver( num_components,  num_minerals,                &
                                 dis_components(1:num_components),             &
                                 the_minerals(1:num_minerals),                 &
                                 nDissTransportables, nPartTransportables,     &
                                 data%listDissTransVars,data%listPartTransVars)

   data%num_comp = nDissTransportables
   data%num_mins = nPartTransportables

   !----------------------------------------------------------------------------
   ! Initialise the module level geochemical values ready for the registration

   ALLOCATE(data%DissComp(nDissTransportables))
   ALLOCATE(data%PartComp(nPartTransportables))

   data%DissComp = zero_
   DO i=1,nDissTransportables
     data%DissComp(i) = dis_initial(i)
     data%Fsed_gch(i) = Fsed_gch(i) / secs_per_day
   END DO
   component_link(num_components+1) = ph_link  ! Special pH var
   data%DissComp(num_components+1) = pH_initial
   data%DissComp(num_components+2) = pe_initial

   data%PartComp = zero_
   DO i=1,num_minerals
     data%PartComp(i) = min_initial(i)
     data%w_gch(i) = w_gch(i)
   END DO

   CALL InitialiseGCProperties(data%DissComp, data%PartComp, 2)

   !print *,'data%DissComp',data%DissComp
   !print *,'data%PartComp',data%PartComp

   !----------------------------------------------------------------------------
   ! Process dis components adding as state vars or dependancies as appropriate
   DO i=1,nDissTransportables
      !print *,'i,',i,component_link(i)
      IF ( component_link(i) .EQ. '' ) THEN
         min = zero_
         IF ( TRIM(data%listDissTransVars(i)) == 'ubalchg' ) min=nan_
         ! Register state variables
         data%id_comp(i) = aed_define_variable(                               &
          !                          TRIM(dis_components(i)),                  &
                                    TRIM(data%listDissTransVars(i)),           &
                                    'mmol/m**3','geochemistry',                &
                                    data%DissComp(i),                          &
                                    minimum=min)
      ELSE
         ! Register external state variable dependencies
         data%id_cdep(i) = aed_locate_variable( TRIM(component_link(i)) )
         data%component_linked(i) = .true.
      ENDIF
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'FeII' ) data%id_feii = data%id_comp(i)
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'FeIII') data%id_feiii= data%id_comp(i)
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'H2S' ) data%id_h2s   = data%id_comp(i)
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'SO4') data%id_so4    = data%id_comp(i)
   ENDDO

   !----------------------------------------------------------------------------
   ! Process minerals adding as state vars or dependancies as appropriate
   DO i=1,num_minerals
      IF ( mineral_link(i) .EQ. '' ) THEN
         ! Register state variables
         data%id_mins(i) = aed_define_variable(                               &
                                    TRIM(the_minerals(i)),                     &
                                    'mmol/m**3','geochemistry',                &
                                    data%PartComp(i),                          &
                                    minimum=zero_,                             &
                                    mobility=data%w_gch(i))
      ELSE
         ! Register external state variable dependencies
         data%id_mdep(i) = aed_locate_variable( mineral_link(i))
         data%mineral_linked(i) = .true.
      ENDIF
   ENDDO


   !----------------------------------------------------------------------------
   ! Register links to other modules

   data%id_o_oxy = aed_locate_variable(oxy_link)

   IF ( ph_link .EQ. '' ) THEN
     ! Register as a diagnostic variable
     data%id_pH = aed_define_diag_variable( 'pH', 'pH', 'pH')
   ELSE
     ! Link to module
     data%id_pH = aed_locate_variable(ph_link)
   ENDIF

   ! solution to get aed_carbon's pCO2 updated for atm exchange ...
   data%id_c_pco2 = aed_locate_variable(pco2_link)

   !----------------------------------------------------------------------------
   ! Register diagnostic variables

   CALL GetListOfGeochemDiagnostics(speciesOutput,diagnosticList)

   DO i=1,SIZE(diagnosticList)
     data%id_gcdiag(i) = aed_define_diag_variable( diagnosticList(i), &
                         '?mmol/m**3?', 'Geochemistry Diagnostic')
   END DO

   data%id_noncon = aed_define_diag_variable( 'noncon_mh', &
                         'niter', 'non-convergence status')

   !----------------------------------------------------------------------------

   ! Register environmental dependencies
   data%id_temp = aed_locate_global( 'temperature' )
   data%id_sal = aed_locate_global( 'salinity' )

   !----------------------------------------------------------------------------

END SUBROUTINE aed_define_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed_initialize_geochemistry(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: temp, tss

   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
!-------------------------------------------------------------------------------
!BEGIN

   !-- Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   !-- Retrieve current (local) state variable values into array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = _STATE_VAR_(data%id_comp(i))
      ELSE
          dissConcs(i) = _STATE_VAR_(data%id_cdep(i))
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = _STATE_VAR_(data%id_mins(i))
      ELSE
          partConcs(i) = _STATE_VAR_(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Redo geochemical equilibration, now spatial initialisation is done
   CALL InitialiseGCProperties(dissConcs, partConcs, 2, inTemp=REAL(temp))

   !-- Copy back into main AED2 arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         _STATE_VAR_(data%id_comp(i)) =  dissConcs(i)
      ELSE
         _STATE_VAR_(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         _STATE_VAR_(data%id_mins(i)) =  partConcs(i)
      ELSE
         _STATE_VAR_(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO


END SUBROUTINE aed_initialize_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_geochemistry model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: reduction,oxidation
   AED_REAL           :: feii,feiii,h2s,so4,oxy,temp
!-------------------------------------------------------------------------------
!BEGIN

   reduction = zero_
   oxidation = zero_

   !-- 1. Iron  ---------------------------------------------------------------!
   IF (simIronRedox) THEN

       oxy   = _STATE_VAR_(data%id_o_oxy) ! oxygen
       temp  = _STATE_VAR_(data%id_temp)
       feii  = _STATE_VAR_(data%id_feii)
       feiii = _STATE_VAR_(data%id_feiii)

       !-- Reduction
       reduction = calcIronReduction(data,feiii,oxy,temp)

       !-- Oxidation
       oxidation = calcIronOxidation(data,feii,oxy,temp)

       !-- Update FeII, O2, FeIII
       _FLUX_VAR_(data%id_feii) = _FLUX_VAR_(data%id_feii) + reduction - oxidation

       _FLUX_VAR_(data%id_o_oxy) = _FLUX_VAR_(data%id_o_oxy) - oxidation *0.25

       _FLUX_VAR_(data%id_feiii) = _FLUX_VAR_(data%id_feiii) + oxidation - reduction

   END IF

   reduction = zero_
   oxidation = zero_
   !-- 2. Sulfur --------------------------------------------------------------!
   IF (simSulfurRedox) THEN

      oxy   = _STATE_VAR_(data%id_o_oxy) ! oxygen
      temp  = _STATE_VAR_(data%id_temp)
      h2s   = _STATE_VAR_(data%id_h2s)
      so4   = _STATE_VAR_(data%id_so4)

      !-- Reduction
      reduction = calcSulfurReduction(data,so4,oxy,temp)

      !-- Oxidation
      oxidation = calcSulfurOxidation(data,h2s,oxy,temp)

      !-- Update
      _FLUX_VAR_(data%id_h2s) = _FLUX_VAR_(data%id_h2s) + reduction - oxidation

      _FLUX_VAR_(data%id_o_oxy) = _FLUX_VAR_(data%id_o_oxy) - oxidation * 0.25

      _FLUX_VAR_(data%id_so4) = _FLUX_VAR_(data%id_so4) + oxidation - reduction

   END IF

END SUBROUTINE aed_calculate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED
! geochemistry. Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
     AED_REAL :: temp,ph
   ! State
     AED_REAL :: dic,oxy
   ! Parameters
     AED_REAL, PARAMETER :: KpHS = 2.0   ! about half of maximum flux at pH5
     AED_REAL, PARAMETER :: KDOs = 2.0*1e3/16.
   ! Temporary variables
     AED_REAL :: gch_flux,oxyEffect,pHEffect
     INTEGER :: i

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp)   ! local temperature

   ! Retrieve current (local) state variable values.
   ph = _STATE_VAR_(data%id_ph)       ! local pH
   oxy = _STATE_VAR_(data%id_o_oxy)

   DO i=1,data%num_comp

     oxyEffect = one_
     pHEffect = one_

     ! Sediment flux dependent on oxygen
      IF( data%Ksed_gch_o2(i)<-1e-8 ) THEN
        oxyEffect = ( data%Ksed_gch_o2(i) / (data%Ksed_gch_o2(i) + oxy) )
      ELSEIF( data%Ksed_gch_o2(i)>1e-8 ) THEN
        oxyEffect = ( oxy / (data%Ksed_gch_o2(i) + oxy) )
      ENDIF
     ! Sediment flux dependent on pH
      IF( data%Ksed_gch_pH(i)<-1e-8 ) THEN
        pHEffect = ( data%Ksed_gch_pH(i) / (data%Ksed_gch_pH(i) + abs(ph-7.0)) )
      ELSEIF( data%Ksed_gch_pH(i)>1e-8 ) THEN
        pHEffect = ( abs(ph-7.0)  / (data%Ksed_gch_pH(i) + abs(ph-7.0)) )
      ENDIF

      gch_flux = data%Fsed_gch(i) * 1.05**(temp-20.0) * oxyEffect * pHEffect

      IF (.NOT.data%component_linked(i)) THEN
        ! geochem module variables, fluxed here as mmol/m2/day
         _FLUX_VAR_(data%id_comp(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ELSE
        ! other module variables, can be fluxed by them
        !_FLUX_VAR_(data%id_cdep(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ENDIF
   ENDDO

END SUBROUTINE aed_calculate_benthic_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, sal
   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
   REAL :: pco2,nc

!-------------------------------------------------------------------------------
!BEGIN

   !-- Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   sal  = _STATE_VAR_(data%id_sal)  ! local salinity

   !-- Retrieve current gch state variable values into work array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = _STATE_VAR_(data%id_comp(i))
      ELSE
          dissConcs(i) = _STATE_VAR_(data%id_cdep(i))
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = _STATE_VAR_(data%id_mins(i))
      ELSE
          partConcs(i) = _STATE_VAR_(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Do geochemical equilibration
   IF (data%simEq) &
      CALL UpdateEquilibration(dissConcs, partConcs, concMode=2, &
                               inTemp=REAL(temp), inSalt=REAL(sal), &
                               stoEq=.true., upDerv=.true.)


   !-- Copy back into main AED2 arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         _STATE_VAR_(data%id_comp(i)) =  dissConcs(i)
      ELSE
         _STATE_VAR_(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         _STATE_VAR_(data%id_mins(i)) =  partConcs(i)
      ELSE
         _STATE_VAR_(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO

   !-- Update diagnostic arrays
   IF( returnGCDerivedVector("pCO2",pco2) > 0 .AND. data%simEq ) THEN
     !print *,'pco2: ',pco2
     _DIAG_VAR_(data%id_c_pco2) = pco2
     !_DIAG_VAR_(data%id_gcdiag(6)) = pco2
   ENDIF
  ! IF( returnGCDerivedVector("NONCON",nc) > 0) THEN
  !   _DIAG_VAR_(data%id_noncon) = nc
  !   _DIAG_VAR_(data%id_gcdiag(5)) = nc
  ! ENDIF


END SUBROUTINE aed_equilibrate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
 FUNCTION calcIronOxidation(data, ironII, oxygen, temp) RESULT(OxidRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: ironII
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   AED_REAL      :: OxidRate
!LOCALS
   AED_REAL, PARAMETER :: ko = 2.14e-5  !@25C
   AED_REAL, PARAMETER :: k1 = 6.78e1   !@25C
   AED_REAL, PARAMETER :: k2 = 2.14e7   !@25C
   AED_REAL, PARAMETER :: FEII_MolWeight = 55.79
   AED_REAL  :: Fe_2, FeOH, FeOH2
   AED_REAL  :: kb, kf
   INTEGER   :: status
!-------------------------------------------------------------------------------
!BEGIN

   kb = data%Riron_box
   kf = data%Riron_aox

   status = returnGCDerivedVector("Fe+2      ", Fe_2)
   IF(status/=1)Fe_2 = zero_
   status = returnGCDerivedVector("FeOH+     ", FeOH)
   IF(status/=1)FeOH = zero_
   status = returnGCDerivedVector("Fe(OH)2   ", FeOH2)
   IF(status/=1)FeOH2 = zero_

   !   OxidRate = oxygen *  (                               &
   !              kf*ko*Fe_2  *FEII_MolWeight*1e3  +        &
   !              kf*k1*FeOH  *FEII_MolWeight*1e3  +        &
   !              kf*k2*FeOH2 *FEII_MolWeight*1e3  +        &
   !              kb*ironII ) * (data%theta_iron_ox**(temp-20.0))

   OxidRate = kb*ironII * (data%theta_iron_ox**(temp-20.0)) * oxygen/(oxygen+100.)

 END FUNCTION calcIronOxidation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
 FUNCTION calcIronReduction(data, ironIII, oxygen, temp) RESULT(RednRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: ironIII
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   !-- Outgoing
   AED_REAL :: RednRate
   !-- Local
   AED_REAL :: light
   AED_REAL :: FeOH2
   AED_REAL :: ko
   INTEGER  :: status
!-------------------------------------------------------------------------------
!BEGIN
   RednRate = zero_

   !-- Biotic rate: iron reducers
   RednRate  = (data%Riron_red  * (data%theta_iron_red**(temp-20.0)) *         &
                data%Kiron_red / (data%Kiron_red + oxygen)) * ironIII

   !-- Photo-reduction
!   IF() THEN
!     light = zero_
!     light = _STATE_VAR_(data%id)
!
!     !-- Photo-reduction rate (/day) is a linear function of PAR
!     ko = data%kFeRpr * light/2.5e3
!
!     FeOH2 = zero_
!~status = returnGCDerivedVector("Fe(OH)2   ", FeOH2)
!~IF(status/=1)FeOH2 = zero_
!     !#MH: Causing ELCD crash 20090817. Temporailiy Disabled
!     !RednRate = RednRate +                                                   &
!     !           ko(vdo) * FeOH2(vdo) * VarDetails(DICHM(FEIII))%MolWeight*1e3
!     RednRate = zero_
!   END IF

 END FUNCTION calcIronReduction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
 FUNCTION calcSulfurReduction(data,sulfate,oxygen,temp) RESULT(RednRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: sulfate
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   !-- Outgoing
   AED_REAL              :: RednRate
!-------------------------------------------------------------------------------
!BEGIN

   RednRate  = zero_
   RednRate  = data%Rsulf_red * (data%theta_sulf_red**(temp-20.0))             &
                * sulfate * (data%Ksulf_red / (data%Ksulf_red + oxygen))

 END FUNCTION calcSulfurReduction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
 FUNCTION calcSulfurOxidation(data,sulfide,oxygen,temp) RESULT(OxdnRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: sulfide
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   !-- Outgoing
   AED_REAL              :: OxdnRate
!-------------------------------------------------------------------------------
!BEGIN

   OxdnRate  = data%Rsulf_ox * (data%theta_sulf_ox**(temp-20.0))               &
                  * sulfide * (oxygen / (data%Ksulf_ox + oxygen))

 END FUNCTION calcSulfurOxidation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_geochemistry
