!<    +---------------------------------------------------------------+
!     |  Type definitions for some global interfaces
!<    +---------------------------------------------------------------+

module strat_kinds
   implicit none
   private

   ! Common constants
   integer, parameter, public :: RK = kind(0.d0) !Real kind
   character(len=3), parameter, public :: version = '2.2'
   integer, parameter, public :: n_simstrat = 4

   ! Common Types
   type, abstract, public :: LinSysSolver
   contains
      procedure(generic_linsyssolver_solve), deferred, nopass :: solve
   end type

   type, abstract, public :: SimModule
   contains
      procedure(generic_simmodule_update), deferred, nopass :: update
   end type

   type, abstract, extends(SimModule), public :: StateVariable
contains
   procedure(generic_statevariable_solve), deferred, nopass :: solve
end type

contains

subroutine generic_linsyssolver_solve(ld, md, ud, rhs, x, ubnd)
   implicit none

   ! Arguments
   real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
   real(RK), intent(in) :: rhs(:) ! right-hand side (b)
   real(RK), intent(out) :: x(:) ! solution (x)
   integer, intent(in) :: ubnd

end subroutine

subroutine generic_statevariable_solve()
   implicit none
end subroutine

subroutine generic_simmodule_update()
   implicit none
end subroutine
end module strat_kinds
