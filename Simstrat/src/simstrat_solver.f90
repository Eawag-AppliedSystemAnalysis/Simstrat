module simstrat_solver
  use simstrat_kinds
  implicit none

  private
  public solve_tridiag_thomas!, solve_tridiag_lapack


contains
  ! Tridiagonal matrix algorithm (solver)
  subroutine solve_tridiag_thomas(ld, md, ud, rhs, solution, N)
      implicit none

      ! Arguments
      integer, intent(in) :: N !dimension
      real(RK), dimension(N), intent(in) :: ld, md, ud, rhs !Lower diagonal, main diagonal, upper diagonal, right-hand side
      real(RK), dimension(N), intent(out) :: solution !solution


      ! Local variables
      real(RK), dimension(N) :: ru, qu
      integer :: i
      
      ru(N) = ld(N)/md(N)
      qu(N) = rhs(N)/md(N)

      do i=N-1,2,-1
          ru(i) = ld(i) / (md(i)-ud(i)*ru(i+1))
          qu(i) = (rhs(i)-ud(i)*qu(i+1)) / (md(i)-ud(i)*ru(i+1))
      end do

      qu(1) = (rhs(1)-ud(1)*qu(2)) / (md(1)-ud(1)*ru(2))

      solution(1) = qu(1)
      do i=2,N
          solution(i)=qu(i)-ru(i)*solution(i-1)
      end do

      return
  end subroutine


  !subroutine solve_tridiag_lapack(ld, md, ud, rhs, solution, N, NRHS)
  !  implicit none

  !  integer, intent(in) :: N, NRHS
  !  real(RK), dimension(N), intent(in) :: ud, md, ld
  !  real(RK), dimension(N, NRHS), intent(in) :: rhs
  !  real(RK), dimension(N, NRHS), intent(out) :: solution
  !  integer :: info

  !  external :: dgtsv

  !  solution(:,:) = rhs(:,:)

  !  call dgtsv(N, NRHS, ld(2:N), md, ud(1:N-1), solution, N, info)

  !  return
  !end subroutine solve_tridiag_lapack

end module simstrat_solver
