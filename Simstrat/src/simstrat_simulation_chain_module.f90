module simstrat_simulation_chain_module
  use simstrat_kinds
  use simstrat_model_module
  use simstrat_finite_volume_implementation, only: StaggeredFiniteVolumeDiscretization
  use simstrat_lake_physics_module
  use utilities
  implicit none

  private

  type, public :: SimstratLakePhysics
    class(SimstratSimulationChain), pointer :: simulation_chain

  contains
    procedure, pass(self) :: assembleSimulationChain => assembleSimstratSimulationChain
    procedure, nopass :: temperature_terms
    procedure, nopass :: salinity_terms
    procedure, nopass :: u_terms
    procedure, nopass :: v_terms
    procedure, nopass :: k_terms
    procedure, nopass :: eps_terms
  end type

  type, extends(SimstratLakePhysics) :: SimstratAEDLakePhysics
  end type

contains
  !###################################################################
  !Simulation Chain
  !###################################################################

  subroutine assembleSimstratSimulationChain(self, T, S, u, v, k, eps)
    implicit none
    class(SimstratLakePhysics), intent(inout) :: self
    class(SimstratStateVariable), pointer, intent(inout) :: T, S, u, v, k, eps

    !additional steps to PDE state variables
    class(SimstratSimulationChain), pointer :: pre_steps
    !additional step to calculate shear, buoyancy and seiche production
    class(SimstratSimulationChain), pointer :: turbulence

    !allocate presteps
    allocate(SimstratSimulationChain :: pre_steps)
    pre_steps%step_procedure => do_pre_steps

    !allocate turbulence
    allocate(SimstratSimulationChain :: turbulence)
    turbulence%step_procedure => do_turbulence

    !assign procedures
    T%calculate_source_terms => temperature_terms
    S%calculate_source_terms => salinity_terms
    u%calculate_source_terms => u_terms
    v%calculate_source_terms => v_terms
    k%calculate_source_terms => k_terms
    eps%calculate_source_terms => eps_terms

    !assemble simulation chain
    !not yet sure if we get a performance drawback from the chain structure
    self%simulation_chain => pre_steps
    pre_steps%next        => u
    u%next                => v
    v%next                => turbulence
    turbulence%next       => T
    T%next                => S
    S%next                => k
    k%next                => eps

  end subroutine

  !###################################################################
  ! pre steps
  !###################################################################

  subroutine do_pre_steps(datum, model)
    implicit none
    class(SimstratModel), intent(inout) :: model
    real(RK), intent(in) :: datum

    select type(discretization => model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)
        call model%forcing%process_forcing_at_date(datum, model) !ok
        call water_column_stability(model%Stab, model%salctr, model%delsal, discretization%h_centres, model%k, model%eps, model%T, model%S, model%cde, model%rho, model%NN, model%cmue1, model%cmue2, discretization%nz_mfq, discretization%nz_tq) !needs work
        call coriolis(model%u, model%v, model%Cori, model%dt, discretization%nz_mfq) !ok
      class default
        call error('Simulation chain not implemented for this discretization')
        stop
    end select
  end subroutine

  subroutine do_turbulence(datum, model)
    implicit none
    class(SimstratModel), intent(inout) :: model
    real(RK), intent(in) :: datum

    select type(discretization => model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)
        call shear_buoyancy_production(model%u, model%v, model%NN, model%num, model%nuh, discretization%h_faces, model%P, model%B, discretization%nz_mfq)
        call seiche_production(model%E_Seiche,model%P_Seiche,model%a_seiche,model%u10,model%v10,model%CD,model%C10,model%UseWindFilt,model%Wf,model%NN,model%gamma,model%q_NN,model%ModSNorm,discretization%h_faces,discretization%a_faces,discretization%dAdz_faces,model%dt,discretization%nz_mfq)
      class default
        call error('Simulation chain not implemented for this discretization')
        stop
    end select
  end subroutine
end module
