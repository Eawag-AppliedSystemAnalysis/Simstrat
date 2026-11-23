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
!     | Generic utilities that are used throughout the code
!<    +---------------------------------------------------------------+


module utilities
   use strat_kinds
   use strat_consts
   implicit none

   !> Save/read an array or matrix (allocatable or pointer) to/from a snapshot file, with its bounds
   interface save_array
      module procedure save_real_array, save_integer_array, save_logical_array, save_array_pointer
   end interface
   interface read_array
      module procedure read_real_array, read_integer_array, read_logical_array, read_array_pointer
   end interface
   interface save_matrix
      module procedure save_real_matrix, save_matrix_pointer
   end interface
   interface read_matrix
      module procedure read_real_matrix, read_matrix_pointer
   end interface

   !> Number as string for messages: toStr(integer), toStr(real) or toStr(real, format)
   interface toStr
      module procedure str_int, str_real
   end interface

contains

   !> Linear interpolation of yi on grid zi (based on the given y on the increasing grid z).
   !> Outside of z, the closest value is used, or NaN if nan_outside is true.
   pure subroutine Interp(z, y, num_z, zi, yi, num_zi, nan_outside)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      implicit none
      real(RK), dimension(:), intent(in) :: z, y, zi
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi
      logical, intent(in), optional :: nan_outside
      integer :: posk1, posk2, posi, i
      real(RK) :: below, above
      logical :: nan

      nan = .false.
      if (present(nan_outside)) nan = nan_outside
      below = y(1); above = y(num_z)
      if (nan) then
         below = ieee_value(below, ieee_quiet_nan); above = below
      end if

      if (num_z == 1) then ! only one value
         yi(1:num_zi) = y(1)
         return
      end if

      ! Points outside of the given grid (points on the boundary get the boundary value, also if z starts or ends
      ! with a repeated depth)
      posk1 = 1
      do while (zi(posk1) < z(1) .or. (.not. nan .and. zi(posk1) == z(1)))
         yi(posk1) = below
         posk1 = posk1 + 1
         if (posk1 > num_zi) exit ! All points below given grid
      end do
      posk2 = num_zi
      do while (zi(posk2) > z(num_z) .or. (.not. nan .and. zi(posk2) == z(num_z)))
         yi(posk2) = above
         posk2 = posk2 - 1
         if (posk2 < 1) exit ! All points above given grid
      end do

      ! Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + ((zi(i) - z(posi)) / (z(posi + 1) - z(posi))) * (y(posi + 1) - y(posi))
      end do
   end subroutine Interp

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
      write(*, *) '[OK] '//message
   end subroutine ok

   subroutine error(message)
      implicit none
      character(len=*), intent(in) :: message
      write(6, *) '[ERROR] '//message
      error stop
   end subroutine error

   subroutine warn(message)
      implicit none
      character(len=*), intent(in) :: message
      write(*, *) '[WARNING] '//message
   end subroutine warn



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

   !> Integer as string without blanks, e.g. for messages
   function str_int(k) result(s)
      implicit none
      integer, intent(in) :: k
      character(len=:), allocatable :: s
      character(len=20) :: buf
      write (buf, '(I0)') k
      s = trim(buf)
   end function

   !> Real as string without blanks, e.g. for messages: with format fmt if given, otherwise with 3 decimals
   !> (or in exponential notation for very large/small numbers)
   function str_real(x, fmt) result(s)
      implicit none
      real(RK), intent(in) :: x
      character(len=*), intent(in), optional :: fmt
      character(len=:), allocatable :: s
      character(len=64) :: buf
      if (present(fmt)) then
         write (buf, fmt) x
      else if (abs(x) >= 1e5_RK .or. (abs(x) < 1e-2_RK .and. x /= 0)) then
         write (buf, '(ES10.3)') x
      else
         write (buf, '(F0.3)') x
      end if
      s = trim(adjustl(buf))
      ! Leading zero (F0.d writes e.g. ".5")
      if (s(1:1) == '.') s = '0'//s
      if (len(s) > 1) then
         if (s(1:2) == '-.') s = '-0'//s(2:)
      end if
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

   ! Density of freshwater as a function of T and S (and optionally P)
   ! according to Chen and Millero (1986), Limnol. Oceanogr. 31(3), 657-662
   ! T [°C] in-situ temperature, S [g/kg] salinity, P [bar] pressure above atmospheric pressure (optional)
   !
   ! Without P (default): density at atmospheric pressure. This is used everywhere except for N2
   !    (state%rho, plunging inflows in strat_lateral).
   ! With P: in-situ density at pressure P. It is ONLY used in update_NN to compare two water parcels at the
   !    same pressure (local potential density). The absolute in-situ density must not be used for N2 directly,
   !    because compression alone would then add ~g^2/c^2 = 4e-5 s-2 to N2.
   subroutine calc_density(rho, T, S, P)
      implicit none

      ! Arguments
      real(RK), intent(in) :: T, S
      real(RK), intent(in), optional :: P
      real(RK), intent(out) :: rho

      ! Local variables
      real(RK) :: rho0t, rho0st, K

      ! Density at atmospheric pressure
      rho0t = 0.9998395_RK + T*(6.7914e-5_RK + T*(-9.0894e-6_RK + T*(1.0171e-7_RK + T*(-1.2846e-9_RK + T*(1.1592e-11_RK + T*(-5.0125e-14_RK))))))
      rho0st = (8.181e-4_RK + T*(-3.85e-6_RK + T*(4.96e-8_RK)))*S

      rho = rho_0*(rho0t + rho0st)

      ! Pressure dependence with the secant bulk modulus K [bar]: rho(S,T,P) = rho(S,T,0)/(1 - P/K(S,T,P))
      ! Because K depends on T and S, the thermal expansivity increases with pressure (by ~30 % at 300 m for 6 °C)
      ! and the temperature of maximum density decreases (by ~0.2 K per 100 m)
      if (present(P)) then
         K = 19652.17_RK + T*(148.113_RK + T*(-2.293_RK + T*(1.256e-2_RK + T*(-4.18e-5_RK)))) &
             + (3.2726_RK + T*(-2.147e-4_RK + T*1.128e-4_RK))*P + (53.238_RK - 0.313_RK*T + 5.728e-3_RK*P)*S
         rho = rho/(1 - P/K)
      end if
   end subroutine


   subroutine save_real_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine save_integer_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      integer, dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine save_logical_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      logical, dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine save_array_pointer(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:), pointer :: array ! No intent(in): otherwise not distinguishable from save_real_array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine read_real_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine

   subroutine read_integer_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      integer, dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine

   subroutine read_logical_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      logical, dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine

   subroutine read_array_pointer(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:), pointer, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine


   subroutine save_real_matrix(output_unit, matrix)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:, :), allocatable, intent(in) :: matrix

      write(output_unit) lbound(matrix, 1), ubound(matrix, 1), lbound(matrix, 2), ubound(matrix, 2)
      write(output_unit) matrix
   end subroutine

   subroutine save_matrix_pointer(output_unit, matrix)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:, :), pointer :: matrix ! No intent(in): otherwise not distinguishable from save_real_matrix

      write(output_unit) lbound(matrix, 1), ubound(matrix, 1), lbound(matrix, 2), ubound(matrix, 2)
      write(output_unit) matrix
   end subroutine

   subroutine read_real_matrix(input_unit, matrix)
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

   subroutine read_matrix_pointer(input_unit, matrix)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:, :), pointer, intent(inout) :: matrix
      integer :: matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2

      read(input_unit) matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2
      read(input_unit) matrix(matrix_lbound_1:matrix_ubound_1, matrix_lbound_2:matrix_ubound_2)
   end subroutine

   pure real(RK) function datum(start_datum, simulation_time)
      implicit none
      real(RK), intent(in) :: start_datum
      integer(8), dimension(2), intent(in) :: simulation_time

      ! Simulation time is a tuple of 2 integers (days, seconds)
      datum = start_datum + real(simulation_time(1), RK) + real(simulation_time(2), RK) / SECONDS_PER_DAY
   end function

end module utilities
