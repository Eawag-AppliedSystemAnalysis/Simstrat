module strat_grid
  implicit none
  private

type, public :: StaggeredGrid

end type

contains
  subroutine generic_linsyssolver_solve(ld, md, ud, rhs, x)
        implicit none

        ! Arguments
        real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
        real(RK), intent(in) :: rhs(:)              ! right-hand side (b)
        real(RK), intent(out) :: x(:)               ! solution (x)

  end subroutine

  subroutine generic_statevariable_solve()
        implicit none
  end subroutine

  subroutine generic_simmodule_update()
        implicit none
  end subroutine
end module strat_kinds
