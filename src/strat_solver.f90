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
!     |  Solver module
!     | Contains implementation of a tridiagonal matrix solver
!<    +---------------------------------------------------------------+


module strat_solver
   use strat_kinds
   implicit none

   private
   public solve_tridiag_thomas

   ! Type that represents a solver based on the Thomas algorithm
   ! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
   type, extends(LinSysSolver), public :: ThomasAlgSolver
   contains
      procedure, nopass :: solve => solve_tridiag_thomas
   end type
contains

   ! Implementation of Thomas algorithm
   subroutine solve_tridiag_thomas(ld, md, ud, rhs, x, ubnd)
      implicit none

      ! Arguments
      real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
      real(RK), intent(in) :: rhs(:) ! right-hand side (b)
      real(RK), intent(out) :: x(:) ! solution (x)
      integer, intent(in) :: ubnd ! upper bound of x

      ! Local variables
      real(RK), dimension(size(md)) :: ru, qu
      integer :: i

      ru(ubnd) = ud(ubnd)/md(ubnd)
      qu(ubnd) = rhs(ubnd)/md(ubnd)

      do i = ubnd - 1, 2, -1
         ru(i) = ud(i)/(md(i) - ld(i)*ru(i + 1))
         qu(i) = (rhs(i) - ld(i)*qu(i + 1))/(md(i) - ld(i)*ru(i + 1))
      end do

      qu(1) = (rhs(1) - ld(1)*qu(2))/(md(1) - ld(1)*ru(2))

      x(1) = qu(1)
      do i = 2, ubnd
         x(i) = qu(i) - ru(i)*x(i - 1)
      end do

      return
   end subroutine

end module strat_solver
