module strat_solver
  use strat_kinds
  implicit none

  private

  ! Type that represents a solver based on the Thomas algorithm
  ! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  type, extends(LinSysSolver), public :: ThomasAlgSolver
     contains
       procedure, nopass :: solve => solve_tridiag_thomas
   end type
contains

  ! Implementation of Thomas algorithm
  subroutine solve_tridiag_thomas(ld, md, ud, rhs, solution, N)
      implicit none

      ! Arguments
      real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
      real(RK), intent(in) :: rhs(:)              ! right-hand side (b)
      real(RK), intent(out) :: x(:)               ! solution (x)

      ! Local variables
      integer :: N
      real(RK), dimension(size(md)) :: ru, qu
      integer :: i

      ru(N) = ld(N-1)/md(N)
      qu(N) = rhs(N)/md(N)

      do i=N-1,2,-1
          ru(i) = ld(i-1) / (md(i)-ud(i)*ru(i+1))
          qu(i) = (rhs(i)-ud(i)*qu(i+1)) / (md(i)-ud(i)*ru(i+1))
      end do

      qu(1) = (rhs(1)-ud(1)*qu(2)) / (md(1)-ud(1)*ru(2))

      solution(1) = qu(1)
      do i=2,N
          solution(i)=qu(i)-ru(i)*solution(i-1)
      end do

      return
  end subroutine

end module strat_solver
