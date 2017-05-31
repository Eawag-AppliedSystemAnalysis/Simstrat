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
                           UpdateEquilibration

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_geochemistry_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_geochemistry_data_t
      !# Variable identifiers
      INTEGER  :: id_comp(MAX_GC_COMPONENTS), id_mins(MAX_GC_MINERALS)
      INTEGER  :: id_cdep(MAX_GC_COMPONENTS), id_mdep(MAX_GC_MINERALS)
      INTEGER  :: id_ubalchg ,id_pH
      INTEGER  :: id_temp
      INTEGER  :: id_sed_dic
      INTEGER  :: id_gcdiag(MAX_GC_COMPONENTS)
      INTEGER  :: id_totC

      !# Model parameters
      INTEGER  :: num_comp, num_mins

      LOGICAL :: component_linked(MAX_GC_COMPONENTS),mineral_linked(MAX_GC_MINERALS)
      CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: listDissTransVars,listPartTransVars
      AED_REAL,          DIMENSION(:), ALLOCATABLE :: DissComp,PartComp
      AED_REAL :: speciation_dt
      AED_REAL :: FSed_cgh(MAX_GC_MINERALS)
      AED_REAL :: w_gch(MAX_GC_MINERALS)

     CONTAINS
         PROCEDURE :: define            => aed2_define_geochemistry
         PROCEDURE :: calculate         => aed2_calculate_geochemistry
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_geochemistry
!        PROCEDURE :: mobility          => aed2_mobility_geochemistry
!        PROCEDURE :: light_extinction  => aed2_light_extinction_geochemistry
!        PROCEDURE :: delete            => aed2_delete_geochemistry

   END TYPE

   INTEGER  :: updateStep

!===============================================================================
CONTAINS


!
!-------------------------------------------------------------------------------
!
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
! Basically, we need to allocate the 'dis_components' and 'the_minerals' as
! state variables in this module. However, if there is a non-blank entry in the
! component_link or mineral_link respectively, then we don't allocate it as
! a geochemistry module state variable and instead just get it as a dependency
! and register the variable id we should use for it.
!
! In addition to these  vars, we need one compulsory one called 'chgbal'
!
! Therefore, in the above example, the variables to allocate as state variables
! in the aed2_geochemistry module would be:
!
!     aed2_geochemistry_Ca
!     aed2_geochemistry_FeIII
!     aed2_geochemistry_FeII
!     aed2_geochemistry_Calcite
!     aed2_geochemistry_FeOH3A
!     aed2_geochemistry_chgbal
!
! (note that DIC and PO4 would need to be dependencies, not state variables in
!       this case)
!
! Once this is sorted, we then need to CALL the I/O function(s) that parse the
! "geochemistry database" file. I think you have already got these compiling in
! the existing aed2_geochemistry.f90, which we set up a while back. This should
! go at the end of the aed2_geochemistry_init(). We then CALL ConfigureGeochem(...)
!
! At this stage aed2_calculate_geochmistry() will be blank
!-------------------------------------------------------------------------------



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
   INTEGER  :: status

   INTEGER  :: speciation_dt
   CHARACTER(len=64) :: geochem_file
   INTEGER  :: num_components,i
   CHARACTER(len=64) :: dis_components(MAX_GC_COMPONENTS)
   AED_REAL          :: dis_initial(MAX_GC_COMPONENTS)
   CHARACTER(len=64) :: component_link(MAX_GC_COMPONENTS)
   AED_REAL          :: Fsed_gch(MAX_GC_COMPONENTS)
   CHARACTER(len=64) :: speciesOutput(10)
   INTEGER  :: num_minerals
   CHARACTER(len=64) :: the_minerals(MAX_GC_MINERALS)
   AED_REAL          :: min_initial(MAX_GC_MINERALS)
   CHARACTER(len=64) :: mineral_link(MAX_GC_MINERALS)
   AED_REAL          :: w_gch(MAX_GC_MINERALS)
   AED_REAL          :: pH_initial, pe_initial

   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: diagnosticList

   INTEGER  :: nDissTransportables, nPartTransportables

   AED_REAL,PARAMETER :: secs_pr_day = 86400.

   NAMELIST /aed2_geochemistry/ speciation_dt, geochem_file,                    &
                     num_components, dis_components, component_link, Fsed_gch, &
                     num_minerals, the_minerals, mineral_link, w_gch, speciesOutput

!
!-------------------------------------------------------------------------------
!BEGIN
   !JOBS
   ! dis_initial and min_initial
   ! remove pe
   ! sub timestep
   ! settling vel
   ! species outputs

   component_link(:) = ''
   data%component_linked(:) = .FALSE.
   data%mineral_linked(:) = .FALSE.


   ! Read the namelist
   read(namlst,nml=aed2_geochemistry,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_geochemistry'

   print *,"WARNING! aed2_geochemistry model is currently under development"

   data%speciation_dt = speciation_dt

   updateStep = 0

   dis_initial = 1.0
   min_initial = 0.1
   pH_initial  = 7.5
   pe_initial  = 8.0

   speciesOutput = ''
   !speciesOutput(1) = 'HCO3-'

   !----------------------------------------------------------------------------
   ! Now load the geochem database

   !CALL aed2_geochem_load_params(data, geochem_file, modelinfo)
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
   data%PartComp = zero_
   DO i=1,num_components
     data%DissComp = dis_initial(i)
   END DO
   component_link(num_components+1) = 'aed2_carbon_pH'

   DO i=1,num_minerals
     data%PartComp = min_initial(i)
   END DO
   data%DissComp(num_components+1) = pH_initial
   data%DissComp(num_components+2) = pe_initial

   CALL InitialiseGCProperties(data%DissComp, data%PartComp, 2)

   print *,'data%DissComp',data%DissComp
   print *,'data%PartComp',data%PartComp

   !----------------------------------------------------------------------------
   ! Process dis components adding as state vars or dependancies as appropriate
   DO i=1,nDissTransportables
   print *,'i,',i,component_link(i)
      IF ( component_link(i) .EQ. '' ) THEN
         ! Register state variables
         data%id_comp(i) = aed2_define_variable(                    &
          !                          TRIM(dis_components(i)),                  &
                                    TRIM(data%listDissTransVars(i)),           &
                                    'mmol/m**3','geochemistry',                &
                                    data%DissComp(i),                          &
                                    minimum=zero_)
      ELSE
         ! Register external state variable dependencies
         data%id_cdep(i) = aed2_locate_variable( TRIM(component_link(i)))
         data%component_linked(i) = .true.
      ENDIF
   ENDDO

   !----------------------------------------------------------------------------
   ! Process minerals adding as state vars or dependancies as appropriate
   DO i=1,num_minerals
      IF ( mineral_link(i) .EQ. '' ) THEN
         ! Register state variables
         data%id_mins(i) = aed2_define_variable(                    &
                                    TRIM(the_minerals(i)),                     &
                                    'mmol/m**3','geochemistry',                &
                                    data%PartComp(i),                          &
                                    minimum=zero_)
      ELSE
         ! Register external state variable dependencies
         data%id_mdep(i) = aed2_locate_variable( mineral_link(i))
         data%mineral_linked(i) = .true.
      ENDIF
   ENDDO


   !----------------------------------------------------------------------------
   ! Register diagnostic variables

   CALL GetListOfGeochemDiagnostics(speciesOutput,diagnosticList)

   DO i=1,SIZE(diagnosticList)
     data%id_gcdiag(i) = aed2_define_diag_variable( diagnosticList(i), &
                     'mmol/m**3', 'Geochemistry Diagnostic')
   END DO

   data%id_sed_dic = aed2_define_sheet_diag_variable(          &
                     'sed_dic','mmol/m**2/d', 'Sediment DIC (CO2) flux')

   !----------------------------------------------------------------------------

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global( 'temperature')


   !----------------------------------------------------------------------------

END SUBROUTINE aed2_define_geochemistry
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
!  AED_REAL           :: dic,diff_dic
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
!  dic = _STATE_VAR_(data%id_dic)! geochemistry


   CALL aed2_geochemistry_update_state(data,column,layer_idx)

   ! Set temporal derivatives
!  diff_dic = 0.

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
!  AED_REAL :: temp

   ! State
!  AED_REAL :: dic,oxy

   ! Temporary variables
!  AED_REAL :: dic_flux,nit_flux

   ! Parameters
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
!  temp = _STATE_VAR_(data%id_temp) ! local temperature

!   ! Retrieve current (local) state variable values.
!  dic = _STATE_VAR_(data%id_dic)! geochemistry

!  IF (data%use_oxy) THEN
!     ! Sediment flux dependent on oxygen and temperature
!     oxy = _STATE_VAR_(data%id_oxy)
!     dic_flux = data%Fsed_dic * data%Ksed_dic/(data%Ksed_dic+oxy) * (data%theta_sed_dic**(temp-20.0))
!  ELSE
!     ! Sediment flux dependent on temperature only.
!     dic_flux = data%Fsed_dic * (data%theta_sed_dic**(temp-20.0))
!  ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   !_SET_BOTTOM_FLUX_(data%id_dic,dic_flux/secs_pr_day)
   !_SET_SED_FLUX_(data%id_dic,dic_flux)
!  _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (dic_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_dic) = _FLUX_VAR_B_(data%id_ben_dic) + (-dic_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
!  _DIAG_VAR_S_(data%id_sed_dic) = dic_flux


END SUBROUTINE aed2_calculate_benthic_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_geochemistry_update_state(data,column,layer_idx)
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
   AED_REAL :: temp, tss

   ! State
!  DOUBLETYPE,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
!  DOUBLETYPE,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs

   ! Temporary variables
   INTEGER  :: i


!-------------------------------------------------------------------------------
!BEGIN

   !Check if the interval is right (sub-timestepping)
!   IF(updateStep<data%speciation_dt) THEN
!     updateStep = updateStep+1
!     RETURN
!   END IF
   updateStep = 0


   ! Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values into array for the gcsolver
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

   ! Do geochemical equilibration
   CALL UpdateEquilibration(dissConcs, partConcs, concMode=2, inTemp=REAL(temp), stoEq=.true.)


   ! Copy back into main arrays
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





END SUBROUTINE aed2_geochemistry_update_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed2_geochemistry
