!<    +---------------------------------------------------------------+
!     | Generic utilities that are used throughout the code
!<    +---------------------------------------------------------------+


module utilities
   use strat_kinds
   implicit none

   type, public :: string
      character(len=:), allocatable :: str
   end type

   interface toStr
      module procedure str_int, str_real
   end interface

contains

   !> Interpolation of yi on grid zi (based on the given y on grid z)
   subroutine Interp(z, y, num_z, zi, yi, num_zi)
      implicit none
      real(RK), dimension(:), intent(in) :: z, y, zi
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi
      integer :: posk1, posk2, posi, i

      if (num_z == 1) then ! only one value
         yi(1:num_zi) = y(1)
         return
      end if

      ! Assign closest value if out of given grid
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

      ! Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + (zi(i) - z(posi))*(y(posi + 1) - y(posi))/(z(posi + 1) - z(posi))
      end do

      return
   end

   ! Interpolation of yi on grid zi (based on the given y on grid z)
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

      ! Assign NaN if out of given grid
      posk1 = 1
      do while (zi(posk1) < z(1))
         yi(posk1) = 0.0_RK
         yi(posk1) = ieee_value(yi(posk1), ieee_quiet_nan) ! NaN
         posk1 = posk1 + 1
      end do
      posk2 = num_zi
      do while (zi(posk2) > z(num_z))
         yi(posk2) = 0.0_RK
         yi(posk2) = ieee_value(yi(posk2), ieee_quiet_nan) ! NaN
         posk2 = posk2 - 1
      end do

      ! Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + (zi(i) - z(posi))*(y(posi + 1) - y(posi))/(z(posi + 1) - z(posi))
      end do

      return
   end subroutine Interp_nan

   !! Integrate discrete function y[x] using the trapezoidal rule
   !!####################################################################
   subroutine Integrate(x, y, inty, num)
      !!####################################################################
      implicit none

      integer :: num, i
      real(RK), dimension(:), intent(in) :: x, y
      real(RK), dimension(:), intent(inout) :: inty

      inty(1) = 0
      do i = 2, num
         inty(i) = inty(i - 1) + 0.5_RK*(x(i) - x(i - 1))*(y(i) + y(i - 1))
      end do
      return
   end

      ! Assign nan to values out of current grid
   !####################################################################
   pure subroutine Assign_nan(y, ubnd, ubnd_grid)
      !####################################################################
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(RK), dimension(:), intent(out) :: y
      integer, intent(in) :: ubnd, ubnd_grid

      integer :: i

      ! Assign NaN if out of given grid
      i = ubnd_grid
      do while (i > ubnd)
         y(i) = 0.0_RK
         y(i) = ieee_value(y(i), ieee_quiet_nan) ! NaN
         i = i - 1
      end do

      return
   end subroutine Assign_nan


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
      else
         inquire (file=fname, exist=file_exists)
         if (.not. file_exists) then
            call error('File '//fname//' does not exist')
         end if
      end if
   end subroutine check_file_exists

   subroutine ok(message)
      implicit none
      character(len=*), intent(in) :: message
      write(6, *) '[OK] '//message
   end subroutine ok

   subroutine error(message)
      implicit none
      character(len=*), intent(in) :: message
      write(6, *) '[ERROR] '//message
      stop
   end subroutine error

   subroutine warn(message)
      implicit none
      character(len=*), intent(in) :: message
      write(6, *) '[WARNING] '//message
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

   character(len=20) function str_int(k)
      implicit none
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str_int, '(a)') k
      str_int = adjustl(str_int)
   end function str_int

   character(len=20) function str_real(k)
      implicit none
      ! "Convert an integer to string."
      real(RK), intent(in) :: k
      write (str_real, '(a)') k
      str_real = adjustl(str_real)
   end function str_real

   character(len=20) function real_to_str(k, fmt)
      implicit none
      real(RK), intent(in) :: k
      character(len=*), intent(in) :: fmt
      write (real_to_str, fmt) k
      real_to_str = adjustl(real_to_str)
   end function

   pure logical function is_leap_year(year)
      integer, intent(in) :: year
      is_leap_year = mod(year, 4) == 0 .and. (.not. mod(year, 100) == 0 .or. mod(year, 400) == 0)
   end function

   pure function calc_days_per_month(year) result(days_per_month)
      integer, intent(in) :: year
      integer, dimension(12) :: days_per_month

      if (is_leap_year(year)) then
         days_per_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      else
         days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      end if
   end function

   ! Initialize the calender (day, month and year) based on the given starting year and and start datum
   subroutine init_calendar(start_year, datum, current_year, current_month, current_day)
      implicit none
      integer, intent(in) :: start_year
      real(RK), intent(in) :: datum
      integer, intent(out) :: current_year, current_month
      real(RK), intent(out) :: current_day

      ! Local variables
      real(RK) :: elapsed_days, days_left, days_per_year
      integer :: i
      integer, dimension(12) :: days_per_month
      
      ! Determine current year
      current_year = start_year
      elapsed_days = 0
      do
         if (is_leap_year(current_year)) then
            days_per_year = 366
         else
            days_per_year = 365
         end if

         if ((elapsed_days + days_per_year) < datum) then
            elapsed_days = elapsed_days + days_per_year
            current_year = current_year + 1
         else
            exit
         end if
      end do
      days_left = datum - elapsed_days

      ! Determine current month and day
      days_per_month = calc_days_per_month(current_year)

      do i=1,12
         if (days_left > days_per_month(i)) then
            days_left = days_left - days_per_month(i)
         else
            current_month = i
            current_day = days_left
            exit
         end if
      end do
   end subroutine


   ! Update calendar month and day (used for albedo assignment)
   subroutine update_calendar(current_year, current_month, current_day, dt)
      implicit none
      integer, intent(inout) :: current_year, current_month
      real(RK), intent(inout) :: current_day
      real(RK), intent(in) :: dt

      ! Local variables
      integer, dimension(12) :: days_per_month
      real(RK) :: current_day_new


      ! Prepare day per month arrays
      days_per_month = calc_days_per_month(current_year)
      
      ! Update current day
      current_day_new = current_day + dt/24/60/60

      ! If new month is reached
      if (ceiling(current_day_new) > days_per_month(current_month)) then
         current_month = current_month + 1
         current_day = current_day_new - floor(current_day_new)

         ! If new year is reached
         if (current_month > 12) then
            current_month = 1
            current_year = current_year + 1
         end if
      else
         ! If not a new month, just go on counting
         current_day = current_day_new
      end if

   end subroutine

   subroutine save_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine read_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine


   subroutine save_matrix(output_unit, matrix)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:, :), allocatable, intent(in) :: matrix

      write(output_unit) lbound(matrix, 1), ubound(matrix, 1), lbound(matrix, 2), ubound(matrix, 2)
      write(output_unit) matrix
   end subroutine

   subroutine read_matrix(input_unit, matrix)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:, :), allocatable, intent(inout) :: matrix
      integer :: matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2

      read(input_unit) matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2
      if (.not. allocated(matrix)) then
         allocate (matrix(matrix_lbound_1:matrix_ubound_1, matrix_lbound_2:matrix_ubound_2))
      end if
      read(input_unit) matrix(matrix_lbound_1:matrix_ubound_1, matrix_lbound_2:matrix_ubound_2)
   end subroutine

   pure real(RK) function datum(start_datum, simulation_time)
      implicit none
      real(RK), intent(in) :: start_datum
      integer(8), intent(in) :: simulation_time

      datum = start_datum + real(simulation_time, RK) / SECONDS_PER_DAY
   end function
end module utilities
