module strat_transport
use strat_kinds
use strat_consts
use strat_simdata
use strat_statevar
use strat_grid
use strat_solver
implicit none
private

type, extends(ModelVariable), public :: TransportModVar
    real(RK), dimension(:), pointer :: dVar
contains
  procedure, pass(self), public :: assign_external_source => transport_assign_external_source
  procedure, pass(self), public :: calc_terms => transport_var_calc_terms
end type

contains
  subroutine transport_assign_external_source(self, dVar)
    class(TransportModVar), intent(inout) :: self
    real(RK), dimension(:), target :: dVar
    self%dVar => dVar
  end subroutine

  subroutine transport_var_calc_terms(self, state, param, sources, boundaries)
    class(TransportModVar), intent(inout) :: self
    class(ModelState), intent(inout) :: state
    class(ModelParam), intent(inout) :: param
    real(RK),  dimension(:) ::  sources, boundaries
    integer :: i
    associate(grid => self%grid, &
              ubnd_fce => self%grid%ubnd_fce, &
              ubnd_vol => self%grid%ubnd_vol)

    !!!!!!!! Define sources !!!!!!!!
    sources = self%dVar  ! We only have a generic, external source!

    ! no explicit boundary conditions
    boundaries(1:ubnd_vol) = 0

  end associate
  end subroutine


end module strat_transport
