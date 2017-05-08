module strat_statevar
  use strat_kinds
  use strat_simdata
  use strat_grid
  use strat_discretization
  use strat_solver
  implicit none
  private

  type, abstract, public :: ModelVariable
    class(LinSysSolver), pointer :: solver
    class(Discretization), pointer :: disc
    class(StaggeredGrid), pointer :: grid
    class(ModelConfig), pointer :: cfg
    real(RK), dimension(:), pointer :: nu
    real(RK), dimension(:), pointer :: var

  contains
    procedure,  pass(self), public :: init => generic_var_init
    procedure(generic_var_calc_terms), deferred, pass(self), public :: calc_terms
    procedure, pass(self), public :: update  => generic_var_update
  end type


  contains
    subroutine generic_var_init(self, cfg, grid, solver, disc, nu, var)
      class(ModelVariable), intent(inout) :: self
      class(LinSysSolver), target :: solver
      class(StaggeredGrid), target :: grid
      class(Discretization), target :: disc
      class(ModelConfig), target :: cfg

      real(RK), dimension(:), target :: nu, var

      ! Assign pointers
      self%cfg => cfg
      self%grid => grid
      self%solver => solver
      self%disc => disc
      self%nu => nu
      self%var => var
    end subroutine

    subroutine generic_var_calc_terms(self, state, param, sources, boundaries)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(:) ::  sources, boundaries
    end subroutine

    subroutine generic_var_update(self, state, param)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(size(self%var)) :: main_diag, rhs, sources, lower_diag, upper_diag, boundaries

      call self%calc_terms(state, param, sources, boundaries)
      call self%disc%create_LES(self%var, self%nu, sources, boundaries, lower_diag, main_diag, upper_diag , rhs, state%dt)
      call self%solver%solve(lower_diag, main_diag, upper_diag, rhs, self%var)

    end subroutine




end module strat_statevar
