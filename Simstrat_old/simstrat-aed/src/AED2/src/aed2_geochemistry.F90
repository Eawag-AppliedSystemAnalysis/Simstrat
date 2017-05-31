!###############################################################################
!#                                                                             #
!# aed2_geochemistry.F90                                                       #
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
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"


MODULE aed2_geochemistry
!-------------------------------------------------------------------------------
! aed2_geochemistry --- geochemistry model
!
! The AED module geochemistry contains equations that describe exchange of
! soluable reactive geochemistry across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_gclib,ONLY : AED_GC_Input,printTables
   USE aed2_gcsolver,ONLY : ConfigEquilibriumSolver, &
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
   PUBLIC aed2_geochemistry_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_geochemistry_data_t
      !# Variable identifiers
      INTEGER  :: id_comp(MAX_GC_COMPONENTS), id_mins(MAX_GC_MINERALS)
      INTEGER  :: id_cdep(MAX_GC_COMPONENTS), id_mdep(MAX_GC_MINERALS)
      INTEGER  :: id_ubalchg ,id_pH, id_c_pco2, id_o_oxy
      INTEGER  :: id_temp,id_sal
      INTEGER  :: id_sed_dic
      INTEGER  :: id_gcdiag(MAX_GC_COMPONENTS), id_noncon
      INTEGER  :: id_totC

      !# Model parameters
      INTEGER  :: num_comp, num_mins

      LOGICAL :: component_linked(MAX_GC_COMPONENTS),mineral_linked(MAX_GC_MINERALS)
      LOGICAL :: simEq
      CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: listDissTransVars,listPartTransVars
      AED_REAL,          DIMENSION(:), ALLOCATABLE :: DissComp,PartComp
      AED_REAL :: speciation_dt
      AED_REAL :: Fsed_gch(MAX_GC_COMPONENTS)
      AED_REAL :: w_gch(MAX_GC_MINERALS)

     CONTAINS
         PROCEDURE :: define            => aed2_define_geochemistry
         PROCEDURE :: initialize        => aed2_initialize_geochemistry
         PROCEDURE :: calculate         => aed2_calculate_geochemistry
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_geochemistry
         PROCEDURE :: equilibrate       => aed2_equilibrate_geochemistry
!        PROCEDURE :: mobility          => aed2_mobility_geochemistry
!        PROCEDURE :: light_extinction  => aed2_light_extinction_geochemistry
!        PROCEDURE :: delete            => aed2_delete_geochemistry

   END TYPE

!===============================================================================
CONTAINS







!###############################################################################
SUBROUTINE aed2_define_geochemistry(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status, i
   INTEGER  :: speciation_dt
   INTEGER  :: num_components,num_minerals
   INTEGER  :: nDissTransportables, nPartTransportables
   LOGICAL  :: simEq = .TRUE.
   AED_REAL          :: min
   AED_REAL          :: dis_initial(MAX_GC_COMPONENTS)
   AED_REAL          :: Fsed_gch(MAX_GC_COMPONENTS)
   AED_REAL          :: min_initial(MAX_GC_MINERALS)
   AED_REAL          :: w_gch(MAX_GC_MINERALS)
   AED_REAL          :: pH_initial, pe_initial
   CHARACTER(len=64) :: geochem_file
   CHARACTER(len=64) :: dis_components(MAX_GC_COMPONENTS)
   CHARACTER(len=64) :: component_link(MAX_GC_COMPONENTS)
   CHARACTER(len=64) :: speciesOutput(10)
   CHARACTER(len=64) :: the_minerals(MAX_GC_MINERALS)
   CHARACTER(len=64) :: mineral_link(MAX_GC_MINERALS)
   CHARACTER(len=64) :: carbon_ph_link = 'CAR_pH'
   CHARACTER(len=64) :: carbon_pco2_link = 'CAR_pCO2'
   CHARACTER(len=64) :: oxy_link = 'OXY_oxy'
   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: diagnosticList

   NAMELIST /aed2_geochemistry/ speciation_dt, geochem_file,                             &
                    num_components, dis_components, component_link, Fsed_gch,dis_initial,&
                    num_minerals, the_minerals, mineral_link, w_gch, min_initial,        &
                    carbon_ph_link, pH_initial, speciesOutput, simEq

!-------------------------------------------------------------------------------
!&aed2_geochemistry
!  speciation_dt  = 10
!  geochem_file   = 'geochem_data.dat'
!  num_components = 5,
!  dis_components = 'DIC','Ca','PO4','FeIII','FeII'
!  component_link = 'aed2_carbon_dic','','aed2_phosphorus_frp','',''
!  Fsed_gch       = 0.,0.,0.,0.,0.
!  num_minerals   = 2
!  the_minerals   = 'Calcite', 'FeOH3A'
!  mineral_link   = '',''
!  w_gch          = 0.,0.
!/
!
!-------------------------------------------------------------------------------
!BEGIN

   print *,"WARNING! aed2_geochemistry model is currently under development"
   ! MH:JOBS
   ! remove pe
   ! species outputs
   ! sediment flux
   ! redox

   !----------------------------------------------------------------------------
   ! Initialise variables
   component_link(:) = ''
   data%component_linked(:) = .FALSE.
   data%mineral_linked(:) = .FALSE.
   data%simEq = simEq

   dis_initial = 1.0  ! default, overwritten by namelist
   min_initial = 0.1  ! default, overwritten by namelist
   pH_initial  = 7.5  ! default, overwritten by namelist
   pe_initial  = 8.0

   !----------------------------------------------------------------------------
   ! Read the namelist
   read(namlst,nml=aed2_geochemistry,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_geochemistry'


   data%speciation_dt = speciation_dt  ! Note this is now managed in FV_AED2 or GLM_AED2
   speciesOutput = ''
   speciesOutput(1) = 'NONCON'
   !speciesOutput(1) = 'HCO3-'

   !----------------------------------------------------------------------------
   ! Now load the geochem database

   !CALL aed2_geochem_load_params(data, geochem_file, modelinfo)
   CALL AED_GC_Input(geochem_file)
   !CALL printTables()

   !----------------------------------------------------------------------------
   ! Configure the geochemical solver
!print*," geo 1"
   CALL ConfigEquilibriumSolver( num_components,  num_minerals,                &
                                 dis_components(1:num_components),             &
                                 the_minerals(1:num_minerals),                 &
                                 nDissTransportables, nPartTransportables,     &
                                 data%listDissTransVars,data%listPartTransVars)

!print*," geo 2"
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
   component_link(num_components+1) = carbon_ph_link  ! Special pH var
   data%DissComp(num_components+1) = pH_initial
   data%DissComp(num_components+2) = pe_initial

   data%PartComp = zero_
   DO i=1,num_minerals
     data%PartComp(i) = min_initial(i)
     data%w_gch(i) = w_gch(i)
   END DO

!print*," geo 3"
   CALL InitialiseGCProperties(data%DissComp, data%PartComp, 2)
!print*," geo 4"

   print *,'data%DissComp',data%DissComp
   print *,'data%PartComp',data%PartComp

   !----------------------------------------------------------------------------
   ! Process dis components adding as state vars or dependancies as appropriate
   DO i=1,nDissTransportables
      !print *,'i,',i,component_link(i)
      IF ( component_link(i) .EQ. '' ) THEN
         min = zero_
         IF ( TRIM(data%listDissTransVars(i)) == 'ubalchg' ) min=nan_
         ! Register state variables
         data%id_comp(i) = aed2_define_variable(                               &
          !                          TRIM(dis_components(i)),                  &
                                    TRIM(data%listDissTransVars(i)),           &
                                    'mmol/m**3','geochemistry',                &
                                    data%DissComp(i),                          &
                                    minimum=min)
      ELSE
         ! Register external state variable dependencies
         data%id_cdep(i) = aed2_locate_variable( TRIM(component_link(i)) )
         data%component_linked(i) = .true.
      ENDIF
   ENDDO

   !----------------------------------------------------------------------------
   ! Process minerals adding as state vars or dependancies as appropriate
   DO i=1,num_minerals
      IF ( mineral_link(i) .EQ. '' ) THEN
         ! Register state variables
         data%id_mins(i) = aed2_define_variable(                               &
                                    TRIM(the_minerals(i)),                     &
                                    'mmol/m**3','geochemistry',                &
                                    data%PartComp(i),                          &
                                    minimum=zero_,                             &
                                    mobility=data%w_gch(i))
      ELSE
         ! Register external state variable dependencies
         data%id_mdep(i) = aed2_locate_variable( mineral_link(i))
         data%mineral_linked(i) = .true.
      ENDIF
   ENDDO


   !----------------------------------------------------------------------------
   ! Register diagnostic variables

!print*," geo 5"
   CALL GetListOfGeochemDiagnostics(speciesOutput,diagnosticList)
!print*," geo 6"

   DO i=1,SIZE(diagnosticList)
     data%id_gcdiag(i) = aed2_define_diag_variable( diagnosticList(i), &
                         '?mmol/m**3?', 'Geochemistry Diagnostic')
   END DO

   !MH solution to get aed_carbon's pCO2 updated for atm exchange ...
   data%id_c_pco2 = aed2_locate_global(carbon_pco2_link)

   data%id_noncon = aed2_define_diag_variable( 'noncon_mh', &
                         'niter', 'non-convergence status')

   data%id_o_oxy = aed2_locate_global(oxy_link)

   data%id_pH = aed2_locate_global(carbon_ph_link)
   !----------------------------------------------------------------------------

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global( 'temperature' )
   data%id_sal = aed2_locate_global( 'salinity' )


   !----------------------------------------------------------------------------

END SUBROUTINE aed2_define_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed2_initialize_geochemistry(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
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
   !print *,'re_initializing pH'

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


END SUBROUTINE aed2_initialize_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_geochemistry model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: reduction,oxidation

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
!  dic = _STATE_VAR_(data%id_dic)! geochemistry

   reduction = zero_
   oxidation = zero_

   !-- 1. Iron  ---------------------------------------------------------------!
   IF (simIronRedox) THEN
!
!     !-- Reduction
!     reduction = calcIronReduction(WQ%a3d(vdo3D,DICHM(FEIII)),                 &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     IF(reduction*DDT > WQ3F%a3d(:,DICHM(FEIII)))THEN
!       reduction = WQ3F%a3d(:,DICHM(FEIII))/DDT
!     ENDIF
!
!     WQ3F%a3d(:,DOC(LABILE)) = WQ3F%a3d(:,DOC(LABILE)) - &
!                               WQ3F%a3d(:,DOC(LABILE)) * reduction * 0.21505 * DDT
!
!     !-- Oxidation
!     oxidation = calcIronOxidation(WQ%a3d(vdo3D,DICHM(FEII)),                  &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     IF(oxidation*DDT > WQ3F%a3d(:,DICHM(FEII)))THEN
!       oxidation = WQ3F%a3d(:,DICHM(FEII))/DDT
!     ENDIF
!
!     WQ3F%a3d(:,DO) = WQ3F%a3d(:,DO) &
!                    - MIN(WQ3F%a3d(:,DO), WQ3F%a3d(:,DO) * oxidation * 0.57278 * DDT)
!
!     !-- Update
!     WQ3F%a3d(:,DICHM(FEII))  = WQ3F%a3d(:,DICHM(FEII))  +                     &
!                                                     (reduction - oxidation)*DDT
!
!     WQ3F%a3d(:,DICHM(FEIII)) = WQ3F%a3d(:,DICHM(FEIII)) +                     &
!                                                     (oxidation - reduction)*DDT
!
   END IF

!  _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (diff_dic)


END SUBROUTINE aed2_calculate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_benthic_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED geochemistry.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
     AED_REAL :: temp,ph
   ! State
     AED_REAL :: dic,oxy
   ! Parameters
     AED_REAL, PARAMETER :: KpHS = 2.0   ! about half of maximum flux at pH5.
     AED_REAL, PARAMETER :: KDOs = 2.0*1e3/16.
   ! Temporary variables
     AED_REAL :: gch_flux,oxyEffect,pHEffect
     INTEGER :: i

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   ! Retrieve current (local) state variable values.
   ph = _STATE_VAR_(data%id_ph)     ! local pH
   oxy = _STATE_VAR_(data%id_o_oxy)
   !IF (data%use_oxy) oxy = _STATE_VAR_(data%id_oxy)

   oxyEffect = one_
   pHEffect = one_

   DO i=1,data%num_comp

      ! Sediment flux dependent on oxygen
      ! Special hard-coded oxy dependent release FeII
      IF( i==2 ) THEN
        oxyEffect = 1e-8
        oxyEffect = ( KdoS / (KdoS + oxy) )
      ENDIF

      ! Special hard-coded pH dependent release for FeIII and Al
      IF( i==3 .OR. i==4 ) THEN
        pHEffect = 1e-8
        IF( pH<6.0 ) pHEffect = ( abs(ph-7.0) / (KpHS + abs(ph-7.0)) )
      ENDIF

      gch_flux = data%Fsed_gch(i) * 1.05**(temp-20.0) * oxyEffect * pHEffect

      !gch_flux = gch_flux * data%Fsed_gch(i) / (data%Fsed_gch(i) + oxy)

      IF (.NOT.data%component_linked(i)) THEN
        ! geochem module variables mmol/m2/day
         _FLUX_VAR_(data%id_comp(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ELSE
        ! other module variables
        !_FLUX_VAR_(data%id_cdep(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ENDIF
   ENDDO


END SUBROUTINE aed2_calculate_benthic_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_equilibrate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
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
   sal = _STATE_VAR_(data%id_sal) ! local salinity

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
   IF( returnGCDerivedVector("pCO2",pco2) > 0) THEN
     !print *,'pco2: ',pco2
     _DIAG_VAR_(data%id_c_pco2) = pco2
     _DIAG_VAR_(data%id_gcdiag(6)) = pco2
   ENDIF
   IF( returnGCDerivedVector("NONCON",nc) > 0) THEN
     _DIAG_VAR_(data%id_noncon) = nc
     _DIAG_VAR_(data%id_gcdiag(5)) = nc
   ENDIF


END SUBROUTINE aed2_equilibrate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!
!
!!------------------------------------------------------------------------------!
! SUBROUTINE performRedoxTransformations(WQ3F, WQ, vdo3D)                       !
!   !-- Incoming                                                                !
!   TYPE (WQ_3D), INTENT(INOUT)       :: WQ    ! Old water quality data         !
!   TYPE (WQ_3F), INTENT(INOUT)       :: WQ3F  ! New water quality data         !
!   INTEGER, DIMENSION(:), INTENT(IN) :: vdo3D ! Process map                    !
!   !-- Local                                                                   !
!   REAL (r_wq), DIMENSION(SIZE(vdo3D)) :: oxidation, reduction                 !
!
!   reduction = zero_
!   oxidation = zero_
!
!   !-- 1. Iron  ---------------------------------------------------------------!
!   IF (simIronRedox) THEN
!
!     !-- Reduction
!     reduction = calcIronReduction(WQ%a3d(vdo3D,DICHM(FEIII)),                 &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     WHERE(reduction*DDT > WQ3F%a3d(:,DICHM(FEIII)))
!       reduction = WQ3F%a3d(:,DICHM(FEIII))/DDT
!     ENDWHERE
!
!     WQ3F%a3d(:,DOC(LABILE)) = WQ3F%a3d(:,DOC(LABILE)) - &
!                               WQ3F%a3d(:,DOC(LABILE)) * reduction * 0.21505 * DDT
!
!     !-- Oxidation
!     oxidation = calcIronOxidation(WQ%a3d(vdo3D,DICHM(FEII)),                  &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     WHERE(oxidation*DDT > WQ3F%a3d(:,DICHM(FEII)))
!       oxidation = WQ3F%a3d(:,DICHM(FEII))/DDT
!     ENDWHERE
!
!     WQ3F%a3d(:,DO) = WQ3F%a3d(:,DO) &
!                    - MIN(WQ3F%a3d(:,DO), WQ3F%a3d(:,DO) * oxidation * 0.57278 * DDT)
!
!
!     !-- Update
!     WQ3F%a3d(:,DICHM(FEII))  = WQ3F%a3d(:,DICHM(FEII))  +                     &
!                                                     (reduction - oxidation)*DDT
!
!     WQ3F%a3d(:,DICHM(FEIII)) = WQ3F%a3d(:,DICHM(FEIII)) +                     &
!                                                     (oxidation - reduction)*DDT
!
!   END IF
!
!
!   reduction = zero_
!   oxidation = zero_
!
!   !-- 2. Manganese -----------------------------------------------------------!
!   IF (simManganRedox) THEN
!
!     !-- Reduction
!     reduction = (c%CHM%Manganese%kMnR * &
!                 (c%CHM%Manganese%vMnR**(GetTemp(vdo3D)-20.0)) * &
!                  c%CHM%Manganese%K_MnR / (c%CHM%Manganese%K_MnR+WQ%a3d(vdo3D,DO))) &
!                  *WQ%a3d(vdo3D,DICHM(MNIV))
!
!     !-- Oxidation
!     oxidation = (c%CHM%Manganese%kMnO * (c%CHM%Manganese%vMnO**(HYS%TEM(vdo3D)-20.0)) * WQ%a3d(vdo3D,DO)      &
!               / (c%CHM%Manganese%K_MnO + WQ%a3d(vdo3D,DO))) * WQ%a3d(vdo3D,DICHM(MNII))
!
!     !-- Update
!     WQ3F%a3d(:,DICHM(MNII))  = WQ3F%a3d(:,DICHM(MNII))  +                     &
!                                                     (reduction - oxidation)*DDT
!     WQ3F%a3d(:,DICHM(MNIV)) = WQ3F%a3d(:,DICHM(MNIV)) +                     &
!                                                     (oxidation - reduction)*DDT
!
!   END IF
!
!   !-- 3. Sulfur --------------------------------------------------------------!
!   IF (simSulfurRedox) THEN
!
!     !-- Reduction
!     reduction = calcSulfurReduction(WQ%a3d(vdo3D,DICHM(SO4)),                 &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     WHERE(reduction*DDT > WQ3F%a3d(:,DICHM(SO4)))
!       reduction = WQ3F%a3d(:,DICHM(SO4))/DDT
!     ENDWHERE
!
!     !-- Oxidation
!     oxidation = calcSulfurOxidation(WQ%a3d(vdo3D,DICHM(H2S)),                 &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     WHERE(oxidation*DDT > WQ3F%a3d(:,DICHM(H2S)))
!       oxidation = WQ3F%a3d(:,DICHM(H2S))/DDT
!     ENDWHERE
!
!     !-- Update
!     WQ3F%a3d(:,DICHM(H2S)) = WQ3F%a3d(:,DICHM(H2S))  +                        &
!                                                     (reduction - oxidation)*DDT
!
!     WQ3F%a3d(:,DICHM(SO4)) = WQ3F%a3d(:,DICHM(SO4))  +                        &
!                                                     (oxidation - reduction)*DDT
!
!   END IF
!
!
!   !-- 4. Arsenic -------------------------------------------------------------!
!   IF (simArsenicRedox) THEN
!
!     !-- Reduction
!     reduction = (c%CHM%Arsenic%kAsR * &
!                 (c%CHM%Arsenic%vAsR**(GetTemp(vdo3D)-20.0)) * &
!                  c%CHM%Arsenic%K_AsR / (c%CHM%Arsenic%K_AsR+WQ%a3d(vdo3D,DO))) &
!                  *WQ%a3d(vdo3D,DICHM(ASV))
!
!     !-- Oxidation
!     oxidation = (c%CHM%Arsenic%kAsO * (c%CHM%Arsenic%vAsO**(HYS%TEM(vdo3D)-20.0)) * WQ%a3d(vdo3D,DO)      &
!               / (c%CHM%Arsenic%K_AsO + WQ%a3d(vdo3D,DO))) * WQ%a3d(vdo3D,DICHM(ASIII))
!
!     !-- Update
!     WQ3F%a3d(:,DICHM(ASIII)) = WQ3F%a3d(:,DICHM(ASIII))  +                    &
!                                                     (reduction - oxidation)*DDT
!     WQ3F%a3d(:,DICHM(ASV))   = WQ3F%a3d(:,DICHM(ASV)) +                       &
!                                                     (oxidation - reduction)*DDT
!
!
!   END IF
!
!
!
!
!   !-- 4. Carbon
!   IF (simCarbonRedox) THEN
!
!     !-- No Reduction included
!
!     !-- CH4 Oxidation
!     oxidation = calcMethaneOxidation(WQ%a3d(vdo3D,DICHM(CH4)),                 &
!                                   WQ%a3d(vdo3D,DO), GetTemp(vdo3D), vdo3D)
!     WHERE(oxidation*DDT > WQ3F%a3d(:,DICHM(CH4)))
!       oxidation = WQ3F%a3d(:,DICHM(CH4))/DDT
!     ENDWHERE
!
!     !-- Update
!     WQ3F%a3d(:,DICHM(CH4)) = WQ3F%a3d(:,DICHM(CH4))  - oxidation*DDT
!
!     WQ3F%a3d(:,DICHM(CO3)) = WQ3F%a3d(:,DICHM(CO3))  + oxidation*DDT
!
!   END IF
!
!   reduction = zero_
!   oxidation = zero_
!
!
!
! END SUBROUTINE performRedoxTransformations
!!------------------------------------------------------------------------------!
!
!
!
!
!!------------------------------------------------------------------------------!
! FUNCTION calcIronOxidation(ironII, oxygen, temp, vdo) RESULT(OxidRate)        !
!   !-- Incoming                                                                !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: ironII                            !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: oxygen                            !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: temp                              !
!   INTEGER,     DIMENSION(:), INTENT(IN)  :: vdo                               !
!   !-- Outgoing                                                                !
!   REAL (r_wq), DIMENSION(SIZE(vdo))      :: OxidRate                          !
!   !-- Local                                                                   !
!   REAL (r_wq), PARAMETER :: ko = 2.14e-5  !@25C
!   REAL (r_wq), PARAMETER :: k1 = 6.78e1   !@25C
!   REAL (r_wq), PARAMETER :: k2 = 2.14e7   !@25C
!   REAL (r_wq), DIMENSION(ngrids) :: Fe_2, FeOH, FeOH2
!   REAL (r_wq)                    :: kb, kf
!   INTEGER                        :: status
!
!
!   kb = c%CHM%Iron%kFeOb
!   kf = c%CHM%Iron%kFeOf
!
!   status = returnGCDerivedVector("Fe+2      ", Fe_2)
!   IF(status/=1)Fe_2 = zero_
!   status = returnGCDerivedVector("FeOH+     ", FeOH)
!   IF(status/=1)FeOH = zero_
!   status = returnGCDerivedVector("Fe(OH)2   ", FeOH2)
!   IF(status/=1)FeOH2 = zero_
!
!   OxidRate = oxygen *  (                                                      &
!              kf*ko*Fe_2(vdo) *VarDetails(DICHM(FEII))%MolWeight*1e3  +        &
!              kf*k1*FeOH(vdo) *VarDetails(DICHM(FEII))%MolWeight*1e3  +        &
!              kf*k2*FeOH2(vdo)*VarDetails(DICHM(FEII))%MolWeight*1e3  +        &
!              kb*ironII ) * (c%CHM%Iron%vFeO**(temp-20.0))
!
!
! END FUNCTION calcIronOxidation
!!------------------------------------------------------------------------------!
!
!
!
!
!!------------------------------------------------------------------------------!
! FUNCTION calcIronReduction(ironIII, oxygen, temp, vdo) RESULT(RednRate)       !
!   !-- Incoming                                                                !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: ironIII                           !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: oxygen                            !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: temp                              !
!   INTEGER,     DIMENSION(:), INTENT(IN)  :: vdo                               !
!   !-- Outgoing                                                                !
!   REAL (r_wq), DIMENSION(SIZE(vdo))      :: RednRate                          !
!   !-- Local                                                                   !
!   REAL (r_wq), DIMENSION(ngrids)         :: light                             !
!   REAL (r_wq), DIMENSION(ngrids)         :: FeOH2                             !
!   REAL (r_wq), DIMENSION(ngrids)         :: ko                                !
!   INTEGER     :: status
!
!   RednRate = zero_
!
!
!   !-- Biotic rate: iron reducers
!   RednRate  = (c%CHM%Iron%kFeR * (c%CHM%Iron%vFeR**(temp-20.0)) *             &
!                c%CHM%Iron%K_FeR / (c%CHM%Iron%K_FeR + oxygen))  * ironIII
!
!   !-- Photo-reduction
!   IF(simLight) THEN
!
!     light = zero_
!     light(vdo) = GetLight(vdo,bwidth="TOTAL",units="W/m2")
!
!     !-- Photo-reduction rate (/day) is a linear function of PAR
!     ko(vdo) = c%CHM%Iron%kFeRpr * light/2500
!
!     status = returnGCDerivedVector("FeOH+2     ", FeOH2)
!     IF(status/=1)FeOH2 = zero_
!
!     !#MH: Causing ELCD crash 20090817. Temporailiy Disabled
!     !RednRate = RednRate +                                                     &
!     !           ko(vdo) * FeOH2(vdo) * VarDetails(DICHM(FEIII))%MolWeight*1e3
!     RednRate = zero_
!
!   END IF
!
! END FUNCTION calcIronReduction
!!------------------------------------------------------------------------------!
!
!
!
!!------------------------------------------------------------------------------!
! FUNCTION calcSulfurReduction(sulfate, oxygen, temp, vdo) RESULT(RednRate)     !
!   !-- Incoming                                                                !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: sulfate                           !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: oxygen                            !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: temp                              !
!   INTEGER,     DIMENSION(:), INTENT(IN)  :: vdo                               !
!   !-- Outgoing                                                                !
!   REAL (r_wq), DIMENSION(SIZE(vdo))      :: RednRate                          !
!   !-- Local                                                                   !
!
!   RednRate  = zero_
!   RednRate  = c%CHM%Sulfur%kSuR * (c%CHM%Sulfur%vSuR**(temp-20.0))            &
!                * sulfate * (c%CHM%Sulfur%K_SuR / (c%CHM%Sulfur%K_SuR + oxygen))
!
! END FUNCTION calcSulfurReduction
!!------------------------------------------------------------------------------!
!
!
!
!!------------------------------------------------------------------------------!
! FUNCTION calcSulfurOxidation(sulfide, oxygen, temp, vdo) RESULT(OxdnRate)     !
!   !-- Incoming                                                                !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: sulfide                           !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: oxygen                            !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: temp                              !
!   INTEGER,     DIMENSION(:), INTENT(IN)  :: vdo                               !
!   !-- Outgoing                                                                !
!   REAL (r_wq), DIMENSION(SIZE(vdo))      :: OxdnRate                          !
!   !-- Local                                                                   !
!
!   OxdnRate  = c%CHM%Sulfur%kSuO * (c%CHM%Sulfur%vSuO**(temp-20.0))            &
!                  * sulfide * (oxygen / (c%CHM%Sulfur%K_SuO + oxygen))
!
! END FUNCTION calcSulfurOxidation
!!------------------------------------------------------------------------------!
!
!
!
!!------------------------------------------------------------------------------!
! FUNCTION calcMethaneOxidation(methane, oxygen, temp, vdo) RESULT(OxdnRate)    !
!   !-- Incoming                                                                !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: methane                           !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: oxygen                            !
!   REAL (r_wq), DIMENSION(:), INTENT(IN)  :: temp                              !
!   INTEGER,     DIMENSION(:), INTENT(IN)  :: vdo                               !
!   !-- Outgoing                                                                !
!   REAL (r_wq), DIMENSION(SIZE(vdo))      :: OxdnRate                          !
!   !-- Local                                                                   !
!
!   OxdnRate  = c%CHM%Methane%kCH4O * (c%CHM%Methane%vCH4O**(temp-20.0))        &
!                  * methane * (oxygen / (c%CHM%Methane%K_CH4O + oxygen))
!
! END FUNCTION calcMethaneOxidation
!!------------------------------------------------------------------------------!
!

END MODULE aed2_geochemistry
