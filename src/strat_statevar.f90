! ---------------------------------------------------------------------------------
!     Simstrat a physical 1D model for lakes and reservoirs
!
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     Copyright (C) 2020, Eawag
!
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>. 
! ---------------------------------------------------------------------------------
!<    +---------------------------------------------------------------+
!     | Abstract base class for simulation variables (K, EPS, T, S etc)
!     | The calc_terms and pos_solve methods should be overwritten according
!<    +---------------------------------------------------------------+


module strat_statevar
   use strat_kinds
   use strat_simdata
   use strat_grid
   use strat_discretization
   use strat_solver
   implicit none
   private

   type, abstract, public :: ModelVariable
      ! Things a model variable needs:
      class(LinSysSolver), pointer :: solver      ! Solver
      class(Discretization), pointer :: disc      ! Discretization scheme
      class(StaggeredGrid), pointer :: grid       ! Grid it lives on
      class(ModelConfig), pointer :: cfg          ! Configuration of the model
      real(RK), dimension(:), pointer :: nu       ! nu concerning this variable
      real(RK), dimension(:), pointer :: var      ! pointer to variable data (usually points to some state var in modelstate)
      integer, pointer :: ubnd                    ! Current upper bound of variable
   contains
      procedure, pass(self), public :: init => generic_var_init
      procedure(generic_var_calc_terms), deferred, pass(self), public :: calc_terms
      procedure, pass(self), public :: update => generic_var_update
      procedure, pass(self), public :: post_solve => generic_var_post_solve
   end type

contains
   subroutine generic_var_init(self, cfg, grid, solver, disc, nu, var, ubnd)
      class(ModelVariable), intent(inout) :: self
      class(LinSysSolver), target :: solver
      class(StaggeredGrid), target :: grid
      class(Discretization), target :: disc
      class(ModelConfig), target :: cfg

      real(RK), dimension(:), target :: nu, var
      integer, target :: ubnd

      ! Assign pointers
      self%cfg => cfg
      self%grid => grid
      self%solver => solver
      self%disc => disc
      self%nu => nu
      self%var => var
      self%ubnd => ubnd
   end subroutine

   ! Generic base method for update of a variable
   subroutine generic_var_update(self, state, param)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(size(self%var)) :: main_diag, rhs, sources, lower_diag, upper_diag, boundaries

      ! Calculate source terms
      call self%calc_terms(state, param, sources, boundaries)

      ! Create linear system of equations
      call self%disc%create_LES(self%var, self%nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, state%dt)

      ! Solve LES
      call self%solver%solve(lower_diag, main_diag, upper_diag, rhs, self%var, self%ubnd)

      ! Do post processing (e.g. set boundary values)
      call self%post_solve(state)

   end subroutine

  ! Generic interface for method that calculates source terms
  subroutine generic_var_calc_terms(self, state, param, sources, boundaries)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(:) ::  sources, boundaries
   end subroutine

   ! Generic interface for method that does post solve processing
   subroutine generic_var_post_solve(self, state)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
   end subroutine

end module strat_statevar
