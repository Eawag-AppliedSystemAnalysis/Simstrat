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
!     |  Interface and implementation of output data loggers
!<    +---------------------------------------------------------------+

module strat_outputfile
   use strat_kinds
   use strat_simdata
   use strat_grid
   use simstrat_aed2
   use utilities
   use csv_module
   implicit none

   private

   ! Main interface for loggers
   type, abstract, public :: SimstratOutputLogger
      !###########################################
      class(OutputConfig), public, pointer   :: output_config
      class(SimConfig), public, pointer :: sim_config
      class(ModelConfig), public, pointer :: model_config
      class(AED2Config), public, pointer :: aed2_config
      class(StaggeredGrid), public, pointer ::grid
      type(csv_file), dimension(:), allocatable :: output_files
      integer, public :: n_depths
      integer, public :: n_vars, n_vars_Simstrat, n_vars_AED2_state, n_vars_AED2_diagnostic
      integer, public :: counter = 0
      integer(8), public, dimension(2) :: simulation_time_for_next_output = 0
      real(RK), dimension(:,:), allocatable :: last_iteration_data

   contains
      procedure(generic_log_init), deferred, pass(self), public :: initialize
      procedure(generic_init_files), deferred, pass(self), public :: init_files
      procedure(generic_calculate_simulation_time_for_next_output), deferred, pass(self), public :: calculate_simulation_time_for_next_output
      procedure(generic_log_start), deferred, pass(self), public :: start
      procedure(generic_log), deferred, pass(self), public :: log   ! Method that is called each loop to write data
      procedure(generic_log_close), deferred, pass(self), public :: close
      procedure(generic_log_save), deferred, pass(self), public :: save
      procedure(generic_log_load), deferred, pass(self), public :: load
   end type

   ! Logger that interpolates on a specified output grid at specific times.
   type, extends(SimstratOutputLogger), public :: InterpolatingLogger
      private
   contains
      procedure, pass(self), public :: initialize => log_init_interpolating
      procedure, pass(self), public :: init_files => init_files_interpolating
      procedure, pass(self), public :: calculate_simulation_time_for_next_output => calculate_simulation_time_for_next_output_interpolating
      procedure, pass(self), public :: start => log_start
      procedure, pass(self), public :: log => log_interpolating
      procedure, pass(self), public :: close => log_close
      procedure, pass(self), public :: save => log_save
      procedure, pass(self), public :: load => log_load
   end type

   type, private :: OutputHelper
      real(RK) :: w0, w1, output_datum
      real(RK), dimension(:), allocatable :: interpolated_data
      logical write_to_file
   contains
      procedure, pass :: init => output_helper_init
      procedure, pass :: add_datum => output_helper_add_datum
      procedure, pass :: add_data_array => output_helper_add_data_array
      procedure, pass :: add_data_scalar => output_helper_add_data_scalar
      procedure, pass :: next_row => output_helper_next_row
   end type


contains

   ! Abstract interface definitions

   subroutine generic_log_init(self, state, sim_config, model_config, aed2_config, output_config, grid, snapshot_file_exists)

      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(ModelState), target :: state
      class(SimConfig), target :: sim_config
      class(ModelConfig), target :: model_config
      class(AED2Config), target :: aed2_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid
      logical, intent(in) :: snapshot_file_exists

   end subroutine

   subroutine generic_init_files(self, output_config, grid, snapshot_file_exists)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid
      logical, intent(in) :: snapshot_file_exists

   end subroutine

   subroutine generic_calculate_simulation_time_for_next_output(self, simulation_time)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
      integer(8), dimension(2), intent(in) :: simulation_time
   end subroutine

   subroutine generic_log(self, simdata)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
      class(SimulationData) :: simdata
   end subroutine

   subroutine generic_log_start(self)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
   end subroutine

   subroutine generic_log_close(self)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
   end subroutine

   subroutine generic_log_save(self)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
   end subroutine

   subroutine generic_log_load(self)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
   end subroutine

   !************************* Init logging ****************************

   ! Init logging for interpolating logger

   subroutine log_init_interpolating(self, state, sim_config, model_config, aed2_config, output_config, grid, snapshot_file_exists)

      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(ModelState), target :: state
      class(SimConfig), target :: sim_config
      class(ModelConfig), target :: model_config
      class(AED2Config), target :: aed2_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid
      logical, intent(in) :: snapshot_file_exists

      integer :: i, j, n_output_times

      self%sim_config => sim_config
      self%model_config => model_config
      self%aed2_config => aed2_config
      self%output_config => output_config
      self%grid => grid
      self%n_vars_Simstrat = size(output_config%output_vars)

      if (self%model_config%couple_aed2) then
        ! allocate AED2 output structure for state variables
        allocate (output_config%output_vars_aed2_state) ! We don't know yet how many variables
        output_config%output_vars_aed2_state%names => state%AED2_state_names
        output_config%output_vars_aed2_state%values => state%AED2_state
        self%n_vars_AED2_state = state%n_AED2_state

        ! Allocate AED2 output structure for diagnostic variables if necessary
        if (aed2_config%output_diagnostic_variables) then
          allocate (output_config%output_vars_aed2_diagnostic) ! We don't know yet how many variables
          output_config%output_vars_aed2_diagnostic%names => state%AED2_diagnostic_names
          output_config%output_vars_aed2_diagnostic%values => state%AED2_diagnostic
          self%n_vars_AED2_diagnostic = state%n_AED2_diagnostic
        else
         self%n_vars_AED2_diagnostic = 0
        end if
      else
         self%n_vars_AED2_state = 0
      end if
      self%n_vars = self%n_vars_Simstrat + self%n_vars_AED2_state + self%n_vars_AED2_diagnostic

      ! If output times are given in file
      if (output_config%thinning_interval == 0) then
         ! Number of output times specified in file
         n_output_times = size(output_config%tout)
         ! If Simulation start larger than output times, abort.
         if (sim_config%start_datum > output_config%tout(1)) then
            call error('Simulation start time is larger than first output time.')
         end if
      end if

      if (output_config%thinning_interval>1) write(6,*) 'Interval [days]: ',output_config%thinning_interval*sim_config%timestep/real(SECONDS_PER_DAY, RK)

      call ok('Output times successfully read')


      ! If output depth interval is given
      if (output_config%depth_interval > 0) then

         ! Allocate zout
         allocate(output_config%zout(ceiling(grid%max_depth/output_config%depth_interval + 1e-6)))

         output_config%zout(1) = 0
         i = 2

         if (output_config%output_depth_reference == 'bottom') then
            do while ((grid%max_depth + 1e-6) > (output_config%depth_interval + output_config%zout(i - 1)))
               output_config%zout(i) = output_config%zout(i - 1) + output_config%depth_interval
               i = i + 1
            end do

         else if (output_config%output_depth_reference == 'surface') then
            ! Add output depth every "output_config%depth_interval" meters until max_depth is reached
            do while ((grid%max_depth + 1e-6) > (output_config%depth_interval - output_config%zout(i - 1)))
               output_config%zout(i) = output_config%zout(i - 1) - output_config%depth_interval
               i = i + 1
            end do

            call reverse_in_place(output_config%zout)
         end if

      else if (output_config%depth_interval == 0) then
         allocate(output_config%zout(size(output_config%zout_read)))
         output_config%zout = output_config%zout_read
      end if
      call ok('Output depths successfully read')

      call self%init_files(output_config, grid, snapshot_file_exists)
   end subroutine

   !************************* Init files ****************************

   ! Initialize files for interpolating logger
   subroutine init_files_interpolating(self, output_config, grid, snapshot_file_exists)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid
      logical, intent(in) :: snapshot_file_exists

      logical :: status_ok, exist_output_folder
      integer :: i, exitstat
      character(len=256) :: mkdirCmd

      self%n_depths = size(output_config%zout)
      allocate (self%last_iteration_data(self%n_vars, self%n_depths))
      self%last_iteration_data = 0

      ! Check if output directory exists
      inquire(file=output_config%PathOut,exist=exist_output_folder)

      ! Create Simstrat output folder if it does not exist
      if(.not.exist_output_folder) then
         call warn('Result folder does not exist, create folder according to config file...')
         mkdirCmd = 'mkdir '//trim(output_config%PathOut)
         call execute_command_line(mkdirCmd, exitstat = exitstat)

         ! mkdir does not seem to accept a path to a folder in execute_command_line, thus a default result folder "Results" will be generated in this case.
         if (exitstat==1) then
            call warn('Result path specified in config file could not be generated. Default result folder "Results" was generated instead.')
            call execute_command_line('mkdir Results')
            output_config%PathOut = 'Results'
         end if
      end if

      call open_files(self, output_config, grid, snapshot_file_exists)

   end subroutine

   subroutine open_files(self, output_config, grid, snapshot_file_exists)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid
      logical, intent(in) :: snapshot_file_exists
      integer :: i, exitstat
      character(len=:), allocatable :: file_path
      logical :: status_ok, append, exist_output_folder
      character(len=256) :: mkdirCmd

      if (allocated(self%output_files)) deallocate (self%output_files)
      allocate (self%output_files(1:self%n_vars))

      do i = 1, self%n_vars_Simstrat
         file_path = output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat'
         inquire (file=file_path, exist=append)
         append = append .and. snapshot_file_exists
         if (self%output_config%output_vars(i)%volume_grid) then
            !Variable on volume grid
            call self%output_files(i)%open(file_path, n_cols=self%n_depths+1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('Datetime')
               call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Variable on face grid
            call self%output_files(i)%open(file_path, n_cols=self%n_depths+1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('Datetime')
               call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         else
            !Variable at surface
            call self%output_files(i)%open(file_path, n_cols=1 + 1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('Datetime')
               call self%output_files(i)%add(grid%z_face(grid%ubnd_fce), real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         end if
      end do

      ! AED2 part
      if (self%model_config%couple_aed2) then
         do i = self%n_vars_Simstrat + 1, self%n_vars
            if (i < (self%n_vars_Simstrat + self%n_vars_AED2_state + 1)) then
               file_path = output_config%PathOut//'/'//trim(self%output_config%output_vars_aed2_state%names(i - self%n_vars_Simstrat))//'_out.dat'
            else
               file_path = output_config%PathOut//'/'//trim(self%output_config%output_vars_aed2_diagnostic%names(i - self%n_vars_Simstrat - self%n_vars_AED2_state))//'_out.dat'
            end if
            inquire (file=file_path, exist=append)

            append = append .and. snapshot_file_exists
            call self%output_files(i)%open(file_path, n_cols=self%n_depths+1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('')
               call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         end do
      end if

   end subroutine

   !************************* Calculate next time point ****************************

   ! Calculate time point for next output (if necessar)
   subroutine calculate_simulation_time_for_next_output_interpolating(self, simulation_time)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      integer(8), dimension(2), intent(in) :: simulation_time
      integer :: i

      associate (thinning_interval => self%output_config%thinning_interval, &
                 timestep => self%sim_config%timestep, &
                 simulation_times_for_output => self%output_config%simulation_times_for_output)

         if (thinning_interval == 0) then
            do i = self%counter, size(simulation_times_for_output,2)
               if (simulation_times_for_output(1,i) > simulation_time(1) .or. &
                  simulation_times_for_output(1,i) == simulation_time(1) .and. simulation_times_for_output(2,i) > simulation_time(2)) then
                  self%simulation_time_for_next_output = simulation_times_for_output(:,i)
                  exit
               end if
            end do
         end if
      end associate
   end subroutine

   !************************* Log ****************************

   ! Log current state on interpolated grid at specific times
   subroutine log_interpolating(self, simdata)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(SimulationData) :: simdata
      !workaround gfortran bug => cannot pass allocatable array to csv file
      real(RK), dimension(self%n_depths) :: values_on_zout
      type(OutputHelper) :: output_helper
      integer :: i, j

      call output_helper%init(simdata%model%simulation_time, simdata%model%simulation_time_old, self%sim_config%timestep, self%simulation_time_for_next_output, &
                              self%sim_config%start_datum, self%output_config%thinning_interval, &
                              self%output_config%simulation_times_for_output, self%counter)
      ! Standard display: display when logged: datum, lake surface, T(1), T(surf)
      if (self%sim_config%disp_simulation == 1 .and. output_helper%write_to_file) then
         write(6,'(F12.4,F16.4,F20.4,F20.4)') simdata%model%datum, simdata%grid%lake_level, &
                                              simdata%model%T(simdata%grid%nz_occupied), simdata%model%T(1)
      ! Extra display: display every iteration: datum, lake surface, T(1), T(surf)
      else if (self%sim_config%disp_simulation == 2) then
         write(6,'(F12.4,F20.4,F15.4,F15.4)') simdata%model%datum, simdata%grid%lake_level, &
                                              simdata%model%T(simdata%grid%ubnd_vol), simdata%model%T(1)
      end if
      do i = 1, self%n_vars_Simstrat
         call output_helper%add_datum(self%output_files(i), "(F12.4)")
         ! If on volume or faces grid
         if (self%output_config%output_vars(i)%volume_grid) then
            ! Interpolate state on volume grid
            call self%grid%interpolate_from_vol(self%output_config%output_vars(i)%values, self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            call output_helper%add_data_array(self%output_files(i), i, self%last_iteration_data, values_on_zout, "(ES14.4E3)")
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Interpolate state on face grid
            call self%grid%interpolate_from_face(self%output_config%output_vars(i)%values, self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            call output_helper%add_data_array(self%output_files(i), i, self%last_iteration_data, values_on_zout, "(ES14.4E3)")
         else
            call output_helper%add_data_scalar(self%output_files(i), i, self%last_iteration_data, self%output_config%output_vars(i)%values_surf, "(ES14.4E3)")
         end if
         call output_helper%next_row(self%output_files(i))
      end do

      ! AED2 part
      if (self%model_config%couple_aed2) then
         do i = self%n_vars_Simstrat + 1, self%n_vars
            ! Write datum
            call output_helper%add_datum(self%output_files(i), "(F12.4)")
            ! Interpolate state on volume grid
            if (i < (self%n_vars_Simstrat + self%n_vars_AED2_state + 1)) then
               call self%grid%interpolate_from_vol(self%output_config%output_vars_aed2_state%values(:,i - self%n_vars_Simstrat), self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            else
               call self%grid%interpolate_from_vol(self%output_config%output_vars_aed2_diagnostic%values(:,i - self%n_vars_Simstrat - self%n_vars_AED2_state), self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            end if
            ! Write state
            call output_helper%add_data_array(self%output_files(i), i, self%last_iteration_data, values_on_zout, "(ES14.4E3)")
            ! Advance to next row
            call output_helper%next_row(self%output_files(i))
         end do
      end if
   end subroutine

   subroutine output_helper_init(self, simulation_time, simulation_time_old, timestep, simulation_time_for_next_output, start_datum, &
                                 thinning_interval, simulation_times_for_output, counter)
      implicit none
      class(OutputHelper), intent(inout) :: self
      integer(8), dimension(2), intent(in) :: simulation_time
      integer(8), dimension(2), intent(in) :: simulation_time_old
      integer, intent(in) :: timestep
      integer(8), dimension(2), intent(inout) :: simulation_time_for_next_output
      real(RK), intent(in) :: start_datum
      integer, intent(in) :: thinning_interval
      integer(8), dimension(:,:), allocatable, intent(in) :: simulation_times_for_output
      integer, intent(inout) :: counter

      ! Local variables
      integer :: i
      logical :: write_condition1, write_condition2

      ! Write condition 1: the next output time is larger than the old simulation time
      write_condition1 = (simulation_time_for_next_output(1) > simulation_time_old(1) .or. &
               simulation_time_for_next_output(1) == simulation_time_old(1) .and. simulation_time_for_next_output(2) > simulation_time_old(2))
      
      ! Write condition 2: the next output time is smaller or equal to the current simulation time
      write_condition2 = (simulation_time_for_next_output(1) < simulation_time(1) .or. &
               simulation_time_for_next_output(1) == simulation_time(1) .and. simulation_time_for_next_output(2) <= simulation_time(2))
      
      self%write_to_file = (write_condition1 .and. write_condition2)

      ! If both writing conditions are fulfilled
      if (self%write_to_file) then
         counter = counter + 1

         ! w1 and w0 are weights for the case where the output time needs to be interpolated from surrounding timesteps
         self%w1 = ((simulation_time(1) - simulation_time_for_next_output(1))*SECONDS_PER_DAY + &
          simulation_time(2) - simulation_time_for_next_output(2)) / real(timestep, RK)
         self%w0 = 1 - self%w1
         ! Save current output time
         self%output_datum = datum(start_datum, simulation_time_for_next_output)
         if (thinning_interval == 0) then
            ! Look for next output time
            do i = counter, size(simulation_times_for_output,2)
               if (simulation_times_for_output(1,i) > simulation_time(1) .or. &
                  simulation_times_for_output(1,i) == simulation_time(1) .and. simulation_times_for_output(2,i) > simulation_time(2)) then
                  simulation_time_for_next_output = simulation_times_for_output(:,i)
                  exit
               end if
            end do
         else
            ! Regular output time spacing
            simulation_time_for_next_output(2) = simulation_time_for_next_output(2) + thinning_interval * timestep
            if (simulation_time_for_next_output(2) >= SECONDS_PER_DAY) then
               simulation_time_for_next_output(2) = simulation_time_for_next_output(2) - SECONDS_PER_DAY
               simulation_time_for_next_output(1) = simulation_time_for_next_output(1) + 1
            end if
         end if
      end if
   end subroutine

   subroutine output_helper_add_datum(self, output_file, format)
      implicit none
      class(OutputHelper), intent(inout) :: self
      type(csv_file), intent(inout) :: output_file
      character(len=*), intent(in) :: format

      if (self%write_to_file) then
         call output_file%add(self%output_datum, real_fmt=format)
      end if
   end subroutine

   subroutine output_helper_add_data_array(self, output_file, index, last_iteration_data, data, format)
      implicit none
      class(OutputHelper), intent(inout) :: self
      type(csv_file), intent(inout) :: output_file
      integer, intent(in) :: index
      real(RK), dimension(:,:), intent(inout) :: last_iteration_data
      real(RK), dimension(:), intent(in) :: data
      character(len=*), intent(in) :: format

      if (self%write_to_file) then
         if (self%w1 == 0) then
            call output_file%add(data, real_fmt=format)
         else
            if (.not. allocated(self%interpolated_data)) then
               allocate (self%interpolated_data(size(data)))
            end if
            self%interpolated_data = self%w1 * last_iteration_data(index, :) + self%w0 * data
            call output_file%add(self%interpolated_data, real_fmt=format)
         end if
      end if
      last_iteration_data(index, :) = data
   end subroutine

   subroutine output_helper_add_data_scalar(self, output_file, index, last_iteration_data, data, format)
      implicit none
      class(OutputHelper), intent(inout) :: self
      type(csv_file), intent(inout) :: output_file
      integer, intent(in) :: index
      real(RK), dimension(:,:), intent(inout) :: last_iteration_data
      real(RK), intent(in) :: data
      character(len=*), intent(in) :: format

      if (self%write_to_file) then
         if (self%w1 == 0) then
            call output_file%add(data, real_fmt=format)
         else
            call output_file%add(self%w1 * last_iteration_data(index, 1) + self%w0 * data, real_fmt=format)
         end if
      end if
      last_iteration_data(index, 1) = data
   end subroutine

   subroutine output_helper_next_row(self, output_file)
      implicit none
      class(OutputHelper), intent(inout) :: self
      type(csv_file), intent(inout) :: output_file

      if (self%write_to_file) then
         call output_file%next_row()
      end if
   end subroutine

   !************************* Start ****************************

   ! Start logging
   subroutine log_start(self)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self

      if (self%sim_config%disp_simulation /= 0) then
         write(6,*)
         write(6,*) ' -------------------------- '
         write(6,*) '   SIMULATION IN PROGRESS   '
         write(6,*) ' -------------------------- '
         write(6,*)
         write(6,'(A12, A20, A20, A20)') 'Time [d]','Surface level [m]','T_surf [degC]','T_bottom [degC]'
      end if
   end subroutine

   !************************* Close ****************************

   ! Close all files
   subroutine log_close(self)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      integer :: i
      logical :: status_ok

      do i = 1, self%n_vars
         call self%output_files(i)%close (status_ok)
      end do
      if (self%sim_config%disp_simulation /= 0) then
         write(6,*)
         write(6,*) ' -------------------------- '
         write(6,*) '    SIMULATION COMPLETED    '
         write(6,*) ' -------------------------- '
         write(6,*)
      end if
   end subroutine

   !************************* Save ****************************

   ! Save logger state
   subroutine log_save(self)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      write(80) self%counter, self%simulation_time_for_next_output(1), self%simulation_time_for_next_output(2)
      call save_matrix(80, self%last_iteration_data)
   end subroutine


   !************************* Load ****************************

   ! Load logger state
   subroutine log_load(self)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      read(81) self%counter, self%simulation_time_for_next_output(1), self%simulation_time_for_next_output(2)
      call read_matrix(81, self%last_iteration_data)
   end subroutine

end module
