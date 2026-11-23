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
!     |  Checks of the inflow/outflow files (issue #73)
!     |     Format, consistency between the files, plausibility of the
!     |     implied temperatures/concentrations and water balance.
!     |     Called once before the simulation starts; the files are read
!     |     as text, independently of the reading in strat_lateral.
!<    +---------------------------------------------------------------+


module strat_inflow_check
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_lateral, only: GenericLateralModule
   use utilities
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private

   public :: check_inflow_files

   ! Content of one inflow/outflow file, read for the input checks (issue #73)
   type :: InflowFileData
      character(len=:), allocatable :: name
      integer :: nval_deep = 0, nval_surface = 0, nval = 0, n_lines = 0
      real(RK), allocatable :: z(:)       ! Depths of the columns (as in the file)
      real(RK), allocatable :: t(:)       ! Times [d]
      real(RK), allocatable :: v(:,:)     ! Values (column, line)
      integer, allocatable :: line_nr(:)  ! Line number in the file of each data line
   end type

contains

   ! ------------------------------------------------------------------------------------------------
   ! Check all inflow/outflow files of the lateral module self. Called once before the simulation starts.
   ! The files are read completely as text, independently of the reading routines used during the
   ! simulation (which are not changed). Clear format errors stop the simulation; suspicious but
   ! possible content gives a warning.
   !
   ! File format: line 1 header; line 2 number of deep and surface columns; line 3 a dummy value and
   ! the depths of all columns [m] (deep: relative to the top of the morphology, surface: relative to
   ! the lake surface); further lines: time [d] and one value per column.
   ! Normal inflow mode (1): all columns are given per metre of depth [unit*m2/s] and integrated over
   ! depth, so each deep or surface group needs at least two depths.
   ! Density-driven inflow mode (2): the deep columns of Qin, Tin, Sin and AED2 files are single
   ! plunging inflows [m3/s], [°C], [‰], [conc.]; surface columns and outflows are as in mode 1.
   ! ------------------------------------------------------------------------------------------------
   subroutine check_inflow_files(self, state, sim_cfg)
      implicit none
      class(GenericLateralModule) :: self
      class(ModelState) :: state
      class(SimConfig) :: sim_cfg

      type(InflowFileData), allocatable :: f(:)
      character(len=:), allocatable :: fname
      integer :: i
      logical :: rho_mode

      rho_mode = (self%cfg%inflow_mode == 2)
      allocate(f(self%n_vars))

      do i = 1, self%n_vars
         if (i > n_simstrat) then
            fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
         else
            fname = trim(self%simstrat_path(i))
         end if
         call read_inflow_file(fname, f(i), self%max_n_inflows)
         call check_single_file(self, f(i), i, rho_mode, sim_cfg)
      end do

      ! Qin, Tin, Sin and AED2 inflow files must have the same columns
      do i = 3, self%n_vars
         call check_same_columns(f(1), f(i), rho_mode)
      end do

      ! Temperatures/concentrations implied by the inflow files
      do i = 3, self%n_vars
         if (i > n_simstrat) then
            if (i - n_simstrat == state%n_pH) then
               call check_implied_values(f(1), f(i), rho_mode, 0.0_RK, 14.0_RK, 'pH')
            else
               call check_implied_values(f(1), f(i), rho_mode, 0.0_RK, huge(1.0_RK), 'concentration')
            end if
         else if (i == 3) then
            call check_implied_values(f(1), f(i), rho_mode, -1.0_RK, 40.0_RK, 'temperature [°C]')
         else
            call check_implied_values(f(1), f(i), rho_mode, 0.0_RK, huge(1.0_RK), 'salinity [‰]')
         end if
      end do

      ! Water balance over the simulation period
      call check_water_balance(self, f(1), f(2), rho_mode, sim_cfg)
   end subroutine

   ! Read a complete inflow/outflow file and check its format (stops with an error for format problems)
   subroutine read_inflow_file(fname, d, max_columns)
      implicit none
      character(len=*), intent(in) :: fname
      type(InflowFileData), intent(out) :: d
      integer, intent(in) :: max_columns

      character(len=:), allocatable :: line, bad
      real(RK), allocatable :: vals(:), t_tmp(:), v_tmp(:,:)
      integer, allocatable :: l_tmp(:)
      integer :: unit, status, line_nr, n, cap
      logical :: eof

      d%name = trim(fname)
      open (newunit=unit, file=fname, action='read', status='old', iostat=status)
      if (status /= 0) call error('File '//trim(fname)//' not found.')
      allocate(vals(16))
      line_nr = 0

      ! Line 1: header
      call read_line(unit, line, eof); line_nr = line_nr + 1
      if (eof) call error(trim(fname)//': file is empty.')

      ! Line 2: number of deep and surface columns. A file with only a header line is allowed (e.g. no outflow: the lake
      ! overflows); the simulation then warns that no data were found and sets the values to zero.
      allocate(d%z(0), d%t(0), d%v(0, 0), d%line_nr(0))
      call next_data_line(unit, line, eof, line_nr)
      if (eof) then
         close (unit)
         return
      end if
      deallocate(d%z, d%t, d%v, d%line_nr)
      call parse_numbers(line, vals, n, bad)
      if (len(bad) > 0) call error(trim(fname)//', line '//toStr(line_nr)//': "'//bad//'" is not a number.')
      if (n < 2) call error('Define the number of columns for deep and surface inflows of '//trim(fname)// &
                            ' (line '//toStr(line_nr)//' must contain two numbers).')
      if (n > 2) call error(trim(fname)//', line '//toStr(line_nr)//': must contain exactly two numbers '// &
                            '(number of deep and surface columns), found '//toStr(n)//'.')
      if (any(vals(1:2) /= real(nint(vals(1:2)), RK)) .or. any(vals(1:2) < 0)) &
         call error(trim(fname)//', line '//toStr(line_nr)//': the numbers of deep and surface columns must be integers >= 0.')
      d%nval_deep = nint(vals(1)); d%nval_surface = nint(vals(2)); d%nval = d%nval_deep + d%nval_surface
      if (d%nval > max_columns) call error(trim(fname)//': more columns ('//toStr(d%nval)// &
                                           ') than ModelConfig.MaxLengthInputData ('//toStr(max_columns)//').')

      ! Line 3: dummy value and depths
      call next_data_line(unit, line, eof, line_nr)
      if (eof) call error(trim(fname)//': line 3 (depths) is missing.')
      call parse_numbers(line, vals, n, bad)
      if (len(bad) > 0) call error(trim(fname)//', line '//toStr(line_nr)//': "'//bad//'" is not a number.')
      if (n /= d%nval + 1) call error(trim(fname)//', line '//toStr(line_nr)//': expected '//toStr(d%nval + 1)// &
         ' numbers (a dummy value and '//toStr(d%nval)//' depths as defined in line 2), found '//toStr(n)//'.')
      allocate(d%z(d%nval)); d%z = vals(2:n)
      if (any(ieee_is_nan(d%z)) .or. .not. all(ieee_is_finite(d%z))) &
         call error(trim(fname)//', line '//toStr(line_nr)//': depths must be finite numbers.')

      ! Data lines
      cap = 1024
      allocate(d%t(cap), d%v(d%nval, cap), d%line_nr(cap))
      d%n_lines = 0
      do
         call next_data_line(unit, line, eof, line_nr)
         if (eof) exit
         call parse_numbers(line, vals, n, bad)
         if (len(bad) > 0) call error(trim(fname)//', line '//toStr(line_nr)//': "'//bad//'" is not a number.')
         if (n /= d%nval + 1) call error(trim(fname)//', line '//toStr(line_nr)//': expected '//toStr(d%nval + 1)// &
            ' numbers (time and '//toStr(d%nval)//' values as defined in line 2), found '//toStr(n)//'.')
         if (any(ieee_is_nan(vals(1:n))) .or. .not. all(ieee_is_finite(vals(1:n)))) &
            call error(trim(fname)//', line '//toStr(line_nr)//': NaN or infinite value.')
         if (d%n_lines == cap) then ! Grow arrays
            cap = 2*cap
            allocate(t_tmp(cap), v_tmp(d%nval, cap), l_tmp(cap))
            t_tmp(1:d%n_lines) = d%t(1:d%n_lines); v_tmp(:, 1:d%n_lines) = d%v(:, 1:d%n_lines)
            l_tmp(1:d%n_lines) = d%line_nr(1:d%n_lines)
            call move_alloc(t_tmp, d%t); call move_alloc(v_tmp, d%v); call move_alloc(l_tmp, d%line_nr)
         end if
         d%n_lines = d%n_lines + 1
         d%t(d%n_lines) = vals(1)
         d%v(:, d%n_lines) = vals(2:n)
         d%line_nr(d%n_lines) = line_nr
      end do
      close (unit)
   end subroutine

   ! Checks of one file that do not need other files
   subroutine check_single_file(self, d, i, rho_mode, sim_cfg)
      implicit none
      class(GenericLateralModule) :: self
      type(InflowFileData), intent(in) :: d
      integer, intent(in) :: i
      logical, intent(in) :: rho_mode
      class(SimConfig) :: sim_cfg

      integer :: k, nd
      logical :: plunging

      nd = d%nval_deep
      plunging = rho_mode .and. i /= 2 ! Deep columns are single plunging inflows

      ! Depths of integrated columns must be increasing (repeated depths are allowed and give a step in the profile).
      ! Plunging inflows are independent of each other and can be given in any order and at the same depth.
      if (nd > 1 .and. .not. plunging) then
         if (any(d%z(2:nd) < d%z(1:nd - 1))) call error(d%name//', line 3: the depths of the deep columns must be '// &
            'increasing (e.g. -20 -10 -5).')
      end if
      if (d%nval_surface > 1) then
         if (any(d%z(nd + 2:d%nval) < d%z(nd + 1:d%nval - 1))) call error(d%name//', line 3: the depths of the '// &
            'surface columns must be increasing (e.g. -2 0).')
      end if
      if (nd > 0) then
         if (any(self%grid%z_zero + d%z(1:nd) < -1e-6_RK)) call error(d%name//', line 3: deep depth below the lake '// &
            'bottom (deepest point of the morphology: '//toStr(-self%grid%z_zero)//' m).')
         if (any(d%z(1:nd) > 1e-6_RK)) call warn(d%name//', line 3: deep depth above the top of the morphology (> 0 m).')
      end if
      if (d%nval_surface > 0) then
         if (any(d%z(nd + 1:d%nval) > 1e-6_RK)) call error(d%name//', line 3: surface depths are relative to the lake '// &
            'surface and must be <= 0.')
      end if

      ! A group with a single depth gives zero flow (values are integrated over depth), except for plunging inflows
      if (nd == 1 .and. .not. plunging) call warn(d%name//': only one deep depth. Values are given per metre of '// &
         'depth and integrated, so this column has no effect. Use at least two depths (e.g. -21 -19).')
      if (d%nval_surface == 1) call warn(d%name//': only one surface depth. Values are given per metre of depth and '// &
         'integrated, so this column has no effect. Use at least two depths (e.g. -1 0).')

      if (d%n_lines == 0) return ! The simulation warns that no data were found and sets the values to zero

      ! Times must not decrease (a repeated time is allowed and gives a step change)
      do k = 2, d%n_lines
         if (d%t(k) < d%t(k - 1)) call error(d%name//', line '//toStr(d%line_nr(k))//': time '//toStr(d%t(k))// &
            ' is smaller than the time in the previous line. Times must be increasing.')
      end do

      ! Coverage of the simulation period (outside the data, the first/last values are used)
      if (d%t(1) > sim_cfg%start_datum + 1e-6_RK) call warn(d%name//': data start at day '//toStr(d%t(1))// &
         ', after the start of the simulation (day '//toStr(sim_cfg%start_datum)//'). The first values are used before.')
      if (d%t(d%n_lines) < sim_cfg%end_datum - 1e-6_RK) call warn(d%name//': data end at day '//toStr(d%t(d%n_lines))// &
         ', before the end of the simulation (day '//toStr(sim_cfg%end_datum)//'). The last values are used after.')

      ! Signs of inflows and outflows
      if (d%nval > 0) then
         if (i == 1 .and. any(d%v(:, 1:d%n_lines) < 0)) call warn(d%name//': negative inflow values (e.g. line '// &
            toStr(first_line(d, -1))//'). Outflows belong into the outflow file.')
         if (i == 2 .and. any(d%v(:, 1:d%n_lines) > 0)) call warn(d%name//': positive outflow values (e.g. line '// &
            toStr(first_line(d, 1))//'). Outflows must be negative.')
      end if
   end subroutine

   ! Qin, Tin, Sin and AED2 inflow files must have the same columns (same numbers of columns and same depths)
   subroutine check_same_columns(q, d, rho_mode)
      implicit none
      type(InflowFileData), intent(in) :: q, d
      logical, intent(in) :: rho_mode
      character(len=:), allocatable :: msg
      logical :: same

      if (q%nval == 0) return ! No inflow
      same = (d%nval_deep == q%nval_deep .and. d%nval_surface == q%nval_surface)
      if (same) same = all(abs(d%z - q%z) < 1e-6_RK)
      if (.not. same) then
         msg = d%name//' does not have the same columns (number of deep/surface columns and depths) as '//q%name
         if (rho_mode) then
            call error(msg//'. In the density-driven inflow mode, each column must correspond to an inflow column.')
         else
            call warn(msg//'.')
         end if
      end if
   end subroutine

   ! Check the temperatures/concentrations implied by an inflow file d together with the inflow rates q.
   ! Values per metre of depth (normal mode, surface columns) are divided by the inflow rate at the same depth and time.
   ! Plunging inflows (deep columns in the density-driven mode) are given directly as temperature/concentration.
   subroutine check_implied_values(q, d, rho_mode, vmin, vmax, what)
      implicit none
      type(InflowFileData), intent(in) :: q, d
      logical, intent(in) :: rho_mode
      real(RK), intent(in) :: vmin, vmax
      character(len=*), intent(in) :: what

      integer :: k, j, kq, n_bad, line_bad
      real(RK) :: qv, x, w, x_bad
      character(len=:), allocatable :: hint

      if (d%nval /= q%nval .or. d%nval == 0 .or. d%n_lines == 0 .or. q%n_lines == 0) return
      if (d%nval_deep /= q%nval_deep) return
      n_bad = 0; kq = 1; line_bad = 0; x_bad = 0
      do k = 1, d%n_lines
         ! Inflow rates at the time of line k (linear interpolation, constant outside the data)
         do while (kq < q%n_lines - 1)
            if (q%t(kq + 1) >= d%t(k)) exit
            kq = kq + 1
         end do
         do j = 1, d%nval
            if (q%n_lines == 1 .or. d%t(k) <= q%t(1)) then
               qv = q%v(j, 1)
            else if (d%t(k) >= q%t(q%n_lines)) then
               qv = q%v(j, q%n_lines)
            else
               w = 1
               if (q%t(kq + 1) > q%t(kq)) w = (d%t(k) - q%t(kq))/(q%t(kq + 1) - q%t(kq))
               qv = (1 - w)*q%v(j, kq) + w*q%v(j, kq + 1)
            end if
            if (qv <= 1e-12_RK) cycle ! No inflow at this depth and time
            if (rho_mode .and. j <= d%nval_deep) then
               x = d%v(j, k)
            else
               x = d%v(j, k)/qv
            end if
            if (x < vmin .or. x > vmax) then
               n_bad = n_bad + 1
               if (n_bad == 1) then
                  line_bad = d%line_nr(k); x_bad = x
               end if
            end if
         end do
      end do
      if (n_bad > 0) then
         if (rho_mode) then
            hint = 'Plunging (deep) inflows are given directly, surface inflows as value times inflow per metre of depth.'
         else
            hint = 'Values must be given as value times inflow per metre of depth (e.g. °C m2/s), like the inflow.'
         end if
         if (vmax < huge(1.0_RK)) then
            call warn(d%name//': implied '//what//' outside ['//toStr(vmin)//', '//toStr(vmax)//'] in '// &
               toStr(n_bad)//' values (first: line '//toStr(line_bad)//', '//toStr(x_bad)//'). '//hint)
         else
            call warn(d%name//': implied '//what//' below '//toStr(vmin)//' in '//toStr(n_bad)// &
               ' values (first: line '//toStr(line_bad)//', '//toStr(x_bad)//'). '//hint)
         end if
      end if
   end subroutine

   ! Mean inflow and outflow over the simulation period and the implied change of the lake level
   subroutine check_water_balance(self, q_in, q_out, rho_mode, sim_cfg)
      implicit none
      class(GenericLateralModule) :: self
      type(InflowFileData), intent(in) :: q_in, q_out
      logical, intent(in) :: rho_mode
      class(SimConfig) :: sim_cfg

      real(RK) :: mean_in, mean_out, dlevel
      real(RK), parameter :: seconds_per_year = 365.25_RK*86400

      mean_in = mean_total_flow(q_in, rho_mode, sim_cfg%start_datum, sim_cfg%end_datum)
      mean_out = mean_total_flow(q_out, .false., sim_cfg%start_datum, sim_cfg%end_datum)
      dlevel = (mean_in + mean_out)*seconds_per_year/self%grid%Az(self%grid%ubnd_fce)
      if (q_out%nval == 0) then
         call ok('Water balance of the inflow/outflow files: mean inflow '//toStr(mean_in)//' m3/s, no outflow data '// &
                 '(the lake overflows when it is full)')
         return
      end if
      call ok('Water balance of the inflow/outflow files: mean inflow '//toStr(mean_in)//' m3/s, mean outflow '// &
              toStr(mean_out)//' m3/s, net '//toStr(mean_in + mean_out)//' m3/s ('//toStr(dlevel)//' m/year)')
      ! A surplus is not a problem: when the lake is full, the excess water leaves by overflow (e.g. without outflow file).
      ! A deficit lowers the lake level until outflows above the lake surface fall dry.
      if (dlevel < -1.0_RK) call warn('The outflows exceed the inflows: the lake level would sink by '// &
         toStr(-dlevel)//' m per year on average (until outflows above the lake surface fall dry). Check the outflow file.')
   end subroutine

   ! Time-mean of the total flow [m3/s] of a file over [t_start, t_end] (constant outside the data, as in the simulation)
   real(RK) function mean_total_flow(d, plunging, t_start, t_end) result(mean)
      implicit none
      type(InflowFileData), intent(in) :: d
      logical, intent(in) :: plunging
      real(RK), intent(in) :: t_start, t_end

      real(RK), allocatable :: qt(:)
      real(RK) :: a, b, qa, qb, integral
      integer :: k

      mean = 0
      if (d%n_lines == 0 .or. d%nval == 0 .or. t_end <= t_start) return
      allocate(qt(d%n_lines))
      do k = 1, d%n_lines
         qt(k) = total_flow(d, k, plunging)
      end do
      if (d%n_lines == 1) then
         mean = qt(1)
         return
      end if
      integral = 0
      ! Before the first and after the last data line
      if (t_start < d%t(1)) integral = integral + qt(1)*(min(t_end, d%t(1)) - t_start)
      if (t_end > d%t(d%n_lines)) integral = integral + qt(d%n_lines)*(t_end - max(t_start, d%t(d%n_lines)))
      ! Linear segments within the data
      do k = 1, d%n_lines - 1
         a = max(t_start, d%t(k)); b = min(t_end, d%t(k + 1))
         if (b <= a) cycle
         qa = qt(k) + (qt(k + 1) - qt(k))*(a - d%t(k))/(d%t(k + 1) - d%t(k))
         qb = qt(k) + (qt(k + 1) - qt(k))*(b - d%t(k))/(d%t(k + 1) - d%t(k))
         integral = integral + 0.5_RK*(qa + qb)*(b - a)
      end do
      mean = integral/(t_end - t_start)
   end function

   ! Total flow [m3/s] of data line k: plunging inflows are summed, values per metre of depth are integrated over depth
   real(RK) function total_flow(d, k, plunging) result(q)
      implicit none
      type(InflowFileData), intent(in) :: d
      integer, intent(in) :: k
      logical, intent(in) :: plunging
      integer :: nd

      nd = d%nval_deep
      q = 0
      if (plunging) then
         q = sum(d%v(1:nd, k))
      else if (nd > 1) then
         q = sum(0.5_RK*(d%v(2:nd, k) + d%v(1:nd - 1, k))*(d%z(2:nd) - d%z(1:nd - 1)))
      end if
      if (d%nval_surface > 1) then
         q = q + sum(0.5_RK*(d%v(nd + 2:d%nval, k) + d%v(nd + 1:d%nval - 1, k))*(d%z(nd + 2:d%nval) - d%z(nd + 1:d%nval - 1)))
      end if
   end function

   ! Line number of the first value with the given sign (-1: negative, 1: positive)
   integer function first_line(d, sgn) result(line_nr)
      implicit none
      type(InflowFileData), intent(in) :: d
      integer, intent(in) :: sgn
      integer :: k
      line_nr = 0
      do k = 1, d%n_lines
         if ((sgn < 0 .and. any(d%v(:, k) < 0)) .or. (sgn > 0 .and. any(d%v(:, k) > 0))) then
            line_nr = d%line_nr(k)
            return
         end if
      end do
   end function

   ! Read one line of arbitrary length
   subroutine read_line(unit, line, eof)
      implicit none
      integer, intent(in) :: unit
      character(len=:), allocatable, intent(out) :: line
      logical, intent(out) :: eof
      character(len=4096) :: buffer
      integer :: status, n_read

      eof = .false.
      ! Usual case: the line fits into the buffer
      read (unit, '(A)', advance='no', iostat=status, size=n_read) buffer
      line = buffer(1:n_read)
      if (is_iostat_eor(status)) return
      if (status /= 0) then ! End of file (or read error)
         if (n_read == 0) eof = .true.
         return
      end if
      ! Longer lines
      do
         read (unit, '(A)', advance='no', iostat=status, size=n_read) buffer
         line = line//buffer(1:n_read)
         if (is_iostat_eor(status)) exit
         if (status /= 0) then ! End of file (or read error)
            if (len(line) == 0) eof = .true.
            exit
         end if
      end do
   end subroutine

   ! Read the next non-empty line (empty lines are skipped, as by the list-directed reads of the simulation).
   ! Tabs, commas and carriage returns are replaced by blanks.
   subroutine next_data_line(unit, line, eof, line_nr)
      implicit none
      integer, intent(in) :: unit
      character(len=:), allocatable, intent(out) :: line
      logical, intent(out) :: eof
      integer, intent(inout) :: line_nr
      do
         call read_line(unit, line, eof)
         if (eof) return
         line_nr = line_nr + 1
         line = blank_separators(line)
         if (len_trim(line) > 0) return
      end do
   end subroutine

   ! Split a line (separated by blanks, see next_data_line) into numbers. bad is the first token that is not a number.
   subroutine parse_numbers(s, vals, n, bad)
      implicit none
      character(len=*), intent(in) :: s
      real(RK), allocatable, intent(inout) :: vals(:)
      integer, intent(out) :: n
      character(len=:), allocatable, intent(out) :: bad

      real(RK), allocatable :: tmp(:)
      real(RK) :: x
      integer :: p, q, l, status

      l = len(s); bad = ''

      n = 0; p = 1
      do
         do while (p <= l)
            if (s(p:p) /= ' ') exit
            p = p + 1
         end do
         if (p > l) exit
         q = p
         do while (q <= l)
            if (s(q:q) == ' ') exit
            q = q + 1
         end do
         if (scan(s(p:q - 1), '/*') > 0) then ! List-directed separators/repeat counts are not allowed
            bad = s(p:q - 1)
            return
         end if
         ! Plain decimal numbers are parsed directly (much faster than a Fortran read); others (e.g. NaN) with read
         if (.not. fast_real(s(p:q - 1), x)) then
            read (s(p:q - 1), *, iostat=status) x
            if (status /= 0) then
               bad = s(p:q - 1)
               return
            end if
         end if
         n = n + 1
         if (n > size(vals)) then
            allocate(tmp(2*size(vals))); tmp(1:size(vals)) = vals; call move_alloc(tmp, vals)
         end if
         vals(n) = x
         p = q
      end do
   end subroutine

   ! Parse a plain decimal number ([sign]digits[.digits][e|E|d|D[sign]digits]). Returns .false. for anything else.
   logical function fast_real(tok, x) result(success)
      implicit none
      character(len=*), intent(in) :: tok
      real(RK), intent(out) :: x
      integer, parameter :: i8 = selected_int_kind(18)
      integer(i8) :: mantissa
      integer :: k, l, n_digits, n_decimals, expo, exp_sign, c
      logical :: negative, after_point, has_digits, has_exp_digits

      success = .false.; x = 0
      l = len(tok); k = 1
      negative = .false.
      if (tok(1:1) == '-' .or. tok(1:1) == '+') then
         negative = (tok(1:1) == '-'); k = 2
      end if
      mantissa = 0; n_digits = 0; n_decimals = 0; after_point = .false.; has_digits = .false.
      do while (k <= l)
         c = ichar(tok(k:k))
         if (c >= 48 .and. c <= 57) then
            if (n_digits >= 18) return ! Too many digits for the integer mantissa
            mantissa = 10*mantissa + (c - 48)
            has_digits = .true.
            if (mantissa > 0) n_digits = n_digits + 1
            if (after_point) n_decimals = n_decimals + 1
         else if (tok(k:k) == '.' .and. .not. after_point) then
            after_point = .true.
         else
            exit
         end if
         k = k + 1
      end do
      if (.not. has_digits) return
      expo = 0
      if (k <= l) then
         if (scan(tok(k:k), 'eEdD') == 0) return
         k = k + 1; exp_sign = 1
         if (k <= l) then
            if (tok(k:k) == '-' .or. tok(k:k) == '+') then
               if (tok(k:k) == '-') exp_sign = -1
               k = k + 1
            end if
         end if
         has_exp_digits = .false.
         do while (k <= l)
            c = ichar(tok(k:k))
            if (c < 48 .or. c > 57) return
            if (expo > 9999) return
            expo = 10*expo + (c - 48); has_exp_digits = .true.
            k = k + 1
         end do
         if (.not. has_exp_digits) return
         expo = exp_sign*expo
      end if
      expo = expo - n_decimals
      if (abs(expo) > 300) return
      x = real(mantissa, RK)*10.0_RK**expo
      if (negative) x = -x
      success = .true.
   end function

   ! Replace tabs, commas and carriage returns by blanks
   function blank_separators(line) result(s)
      implicit none
      character(len=*), intent(in) :: line
      character(len=len(line)) :: s
      integer :: k
      s = line
      do k = 1, len(s)
         if (s(k:k) == char(9) .or. s(k:k) == ',' .or. s(k:k) == char(13)) s(k:k) = ' '
      end do
   end function

end module strat_inflow_check
