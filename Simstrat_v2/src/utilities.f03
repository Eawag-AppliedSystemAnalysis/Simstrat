!     +---------------------------------------------------------------+
!     | Generic utilities that are used throughout the code
!     +---------------------------------------------------------------+


module utilities
   use strat_kinds
   implicit none

   type, public :: string
      character(len=:), allocatable :: str
   end type

contains

   !Interpolation of yi on grid zi (based on the given y on grid z)
   !####################################################################
   subroutine Interp(z, y, num_z, zi, yi, num_zi)
      !####################################################################
      implicit none
      real(RK), dimension(:), intent(in) :: z, y, zi
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi
      integer :: posk1, posk2, posi, i

      if (num_z == 1) then ! only one value
         yi(1:num_zi) = y(1)
         return
      end if

      !Assign closest value if out of given grid
      posk1 = 1
      do while (zi(posk1) <= z(1))
         yi(posk1) = y(1)
         posk1 = posk1 + 1
      end do
      posk2 = num_zi
      do while (zi(posk2) >= z(num_z))
         yi(posk2) = y(num_z)
         posk2 = posk2 - 1
      end do

      !Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + (zi(i) - z(posi))*(y(posi + 1) - y(posi))/(z(posi + 1) - z(posi))
      end do

      return
   end

   !Interpolation of yi on grid zi (based on the given y on grid z)
   !####################################################################
   pure subroutine Interp_nan(z, y, num_z, zi, yi, num_zi)
      !####################################################################
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(RK), dimension(:), intent(in) :: z, y, zi
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi

      integer posk1, posk2, posi, i

      !Assign NaN if out of given grid
      posk1 = 1
      do while (zi(posk1) < z(1))
         yi(posk1) = 0.0_RK
         yi(posk1) = ieee_value(yi(posk1), ieee_quiet_nan) ! NaN
         posk1 = posk1 + 1
      end do
      posk2 = num_zi
      do while (zi(posk2) > z(num_z))
         yi(posk1) = 0.0_RK
         yi(posk2) = ieee_value(yi(posk2), ieee_quiet_nan) ! NaN
         yi(posk2) = z(num_z)
         posk2 = posk2 - 1
      end do

      !Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + (zi(i) - z(posi))*(y(posi + 1) - y(posi))/(z(posi + 1) - z(posi))
      end do

      return
   end subroutine Interp_nan

   !!Integrate discrete function y[x] using the trapezoidal rule
   !!####################################################################
   subroutine Integrate(x, y, inty, num)
      !!####################################################################
      implicit none

      integer :: num, i
      real(RK), intent(in) :: x(1:num), y(1:num)
      real(RK), intent(inout) :: inty(1:num)

      inty(1) = 0
      do i = 2, num
         inty(i) = inty(i - 1) + 0.5_RK*(x(i) - x(i - 1))*(y(i) + y(i - 1))
      end do
      return
   end

   pure function linspace(x0, xend, n, endpoint) result(x)
      implicit none

      integer, intent(in) :: n
      real(RK), intent(in) :: x0, xend
      logical, optional, intent(in) :: endpoint
      real(RK), dimension(n) :: x

      real(RK) :: dx
      real(RK) :: denom
      integer :: i

      denom = real(n - 1, RK)
      if (present(endpoint) .and. (.not. endpoint)) then
         denom = real(n, RK)
      end if

      dx = (xend - x0)/denom

      x = [(real(i, RK)*dx + x0, i=0, n - 1)]

      return
   end function linspace

   pure subroutine diff(d, a, N)

      implicit none

      integer, intent(in) :: N
      real(RK), dimension(N - 1), intent(out) :: d
      real(RK), dimension(N), intent(in) :: a

      d = a(2:N) - a(1:N - 1)

      return
   end subroutine diff

   subroutine check_file_exists(fname)
      implicit none
      character(len=*), intent(in) :: fname

      logical :: file_exists
      if (fname == '') then
         call error('Filename is empty')
         !write(*,*) 'Filename is empty'
         stop
      else
         inquire (file=fname, exist=file_exists)
         if (.not. file_exists) then
            call error('File '//fname//' does not exist')
            !write(*,*) 'File '//fname//' does not exist'
            stop
         end if
      end if
   end subroutine check_file_exists

   subroutine ok(message)
      implicit none
      character(len=*), intent(in) :: message
      write (*, *) '[OK] '//message
   end subroutine ok

   subroutine error(message)
      implicit none
      character(len=*), intent(in) :: message
      write (*, *) '[ERROR] '//message
   end subroutine error

   subroutine warn(message)
      implicit none
      character(len=*), intent(in) :: message
      write (*, *) '[WARNING] '//message
   end subroutine warn

   pure function find_index_ordered(array, target_value) result(idx)
      implicit none
      real(RK), dimension(:), intent(in) :: array
      real(RK), intent(in) :: target_value

      integer :: idx

      do idx = 1, size(array)
         if (array(idx) > target_value) exit
      end do
   end function

   pure function linear_interpolate(t_start, t_end, v_start, v_end, t) result(v)
      implicit none
      real(RK), intent(in) :: t_start, t_end, v_start, v_end, t
      real(RK) :: v

      v = v_start + t*(v_end - v_start)/(t_end - t_start)
   end function

   pure function convert2height_above_sed(z, z_zero) result(h)
      implicit none
      real(RK), dimension(:), intent(in) :: z
      real(RK), intent(in) :: z_zero

      real(RK), dimension(size(z)) :: h
      integer :: n

      n = size(z)
      h = -z_zero + z(n:1:-1)
   end function


   ! Reverse an array without allocating a second array
   subroutine reverse_in_place(in_arr)
      implicit none
      real(RK), intent(inout) :: in_arr(:)
      real(RK) :: temp

      integer :: first, last, i, len

      first = lbound(in_arr, dim=1)
      last = ubound(in_arr, dim=1)
      len = size(in_arr)

      ! Works for even and odd sized arrays
      !(as len/2 is always integer and not rounded, but cutoff)
      do i = last, first + int(len/2), -1
         temp = in_arr(i)
         in_arr(i) = in_arr(len + 1 - i)
         in_arr(len + 1 - i) = temp
      end do

   end subroutine

end module utilities
