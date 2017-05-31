module utilities
  use simstrat_type
  use SimstratModel
  implicit none

  private
  public Interp, Interp_nan, Tridiagonal, Integrate

contains
! Linear interpolation of y_grid on z_grid (based on the given y_in on z_in)
!####################################################################
subroutine Interp(z_in,y_in,li_in,z_grid,y_grid,li_grid)
!####################################################################
      implicit none
      real(dp), intent(in) :: z_in(0:), y_in(0:), z_grid(0:)
      real(dp), intent(inout) :: y_grid(0:)
      integer, intent(in) :: li_in, li_grid   ! last index in, last index on grid

      integer posk1,posk2,posi,i

      !Assign closest value if out of given grid
      posk1 = 0
      do while (z_grid(posk1)<=z_in(0))
         y_grid(posk1) = y_in(0)
         posk1 = posk1+1
      end do
      posk2 = li_grid
      do while (z_grid(posk2)>=z_in(li_in))
         y_grid(posk2) = y_in(li_in)
         posk2 = posk2-1
      end do

      !Linear interpolation
      posi = 0
      do i=posk1,posk2
         do while (z_grid(i)>z_in(posi+1))
            posi = posi+1
         end do
         y_grid(i) = y_in(posi) + (z_grid(i)-z_in(posi))*(y_in(posi+1)-y_in(posi))/(z_in(posi+1)-z_in(posi))
      end do

      return
end


!Linear interpolation of y_grid on z_grid (based on the given y_in on z_in)
!####################################################################
subroutine Interp_nan(z_in,y_in,li_in,z_grid,y_grid,li_grid)
!####################################################################
      implicit none
      real(dp), intent(in) :: z_in(0:), y_in(0:), z_grid(0:)
      real(dp), intent(inout) :: y_grid(0:)
      integer, intent(in) :: li_in, li_grid    ! last index in, last index on grid

      integer posk1, posk2, posi, i

      !Assign NaN if out of given grid
      posk1 = 0
      do while (z_grid(posk1)<z_in(0))
         y_grid(posk1) = 0.0_dp
         y_grid(posk1) = 0.0_dp/y_grid(posk1)   ! NaN
         posk1 = posk1+1
      end do
      posk2 = li_grid
      do while (z_grid(posk2)>z_in(li_in))
         y_grid(posk1) = 0.0_dp
         y_grid(posk2) = 0.0_dp/y_grid(posk1)   ! NaN
         posk2 = posk2-1
      end do

      !Linear interpolation
      posi = 0
      do i=posk1,posk2
         do while (z_grid(i)>z_in(posi+1))
            posi = posi+1
         end do
         y_grid(i)=y_in(posi)+(z_grid(i)-z_in(posi))*(y_in(posi+1)-y_in(posi))/(z_in(posi+1)-z_in(posi))
      end do

      return
end subroutine Interp_nan

! Tridiagonal matrix algorithm (solver)
! FB 2016: deleted nz as an argument
!####################################################################
subroutine Tridiagonal(fi,li,au,bu,cu,du,value)
!####################################################################
    !use SimstratModel
    implicit none

    ! Global variables
    integer, intent(in) :: fi,li !First index, last index
    ! Upper diagonal, diagonal, lower diagonal, right-hand side, solution
    real(dp), intent(in) :: au(0:),bu(0:),cu(0:),du(0:)
    real(dp), intent(inout) :: value(0:)

    ! Local variables
    real(dp) :: ru(0:nz),qu(0:nz)
    integer i


    ru(li)=au(li)/bu(li)
    qu(li)=du(li)/bu(li)

    do i=li-1,fi+1,-1
        ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
        qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
    end do

    qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

    value(fi)=qu(fi)
    do i=fi+1,li
        value(i)=qu(i)-ru(i)*value(i-1)
    end do

    return
end subroutine

!Integrate discrete function y[x] using the trapezoidal rule
!####################################################################
subroutine Integrate(x,y,inty,nval)
!####################################################################

      implicit none

      ! Global variables
      real(dp), intent(in) :: x(0:), y(0:)
      real(dp), intent(inout) :: inty(0:)
      integer :: nval

      ! Local variables
      integer :: i


      inty(0) = 0
      do i=1,nval-1
         inty(i) = inty(i-1) + 0.5_dp*(x(i)-x(i-1))*(y(i)+y(i-1))
      end do

      return
end subroutine

end module utilities
