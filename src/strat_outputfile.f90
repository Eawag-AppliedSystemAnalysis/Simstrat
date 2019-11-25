!<    +---------------------------------------------------------------+
!     |  Interface and implementation of output data loggers
!<    +---------------------------------------------------------------+

module strat_outputfile
   use strat_kinds
   use strat_simdata
   use strat_grid
   use utilities
   use csv_module
   implicit none

   private

   ! Main interface for loggers
   type, abstract, public :: SimstratOutputLogger
      !###########################################
      class(OutputConfig), public, pointer   :: output_config
      class(SimConfig), public, pointer :: sim_config
      class(StaggeredGrid), public, pointer ::grid
      type(csv_file), dimension(:), allocatable :: output_files
      integer, public :: n_depths
      integer, public :: n_vars
      integer, public :: counter = 0
      integer(8), public :: simulation_time_for_next_output = 0
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
   subroutine generic_log_init(self, sim_config, output_config, grid, snapshot_file_exists)
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(SimConfig), target :: sim_config
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
      integer(8), intent(in) :: simulation_time
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
   subroutine log_init_interpolating(self, sim_config, output_config, grid, snapshot_file_exists)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(SimConfig), target :: sim_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid
      logical, intent(in) :: snapshot_file_exists

      integer :: i, j, n_output_times
      real(RK), dimension(size(output_config%tout)) :: tout_test ! Array to test if computed tout with
      ! adjusted timestep is the same as the tout given in file

      self%sim_config => sim_config
      self%output_config => output_config
      self%grid => grid
      self%n_vars = size(output_config%output_vars)

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

      ! Create output folder if it does not exist
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
      integer :: i
      character(len=:), allocatable :: file_path
      logical :: status_ok, append

      if (allocated(self%output_files)) deallocate (self%output_files)
      allocate (self%output_files(1:self%n_vars))

      do i = 1, self%n_vars
         file_path = output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat'
         inquire (file=file_path, exist=append)
         append = append .and. snapshot_file_exists
         if (self%output_config%output_vars(i)%volume_grid) then
            !Variable on volume grid
            call self%output_files(i)%open(file_path, n_cols=self%n_depths+1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('')
               call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Variable on face grid
            call self%output_files(i)%open(file_path, n_cols=self%n_depths+1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('')
               call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         else
            !Variable at surface
            call self%output_files(i)%open(file_path, n_cols=1 + 1, append=append, status_ok=status_ok)
            if (.not. append) then
               call self%output_files(i)%add('')
               call self%output_files(i)%add(grid%z_face(grid%ubnd_fce), real_fmt='(F12.3)')
               call self%output_files(i)%next_row()
            end if
         end if
      end do
   end subroutine

   !************************* Calculate next time point ****************************

   ! Calculate time point for next output (if necessary)
   subroutine calculate_simulation_time_for_next_output_interpolating(self, simulation_time)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      integer(8), intent(in) :: simulation_time
      integer :: i

      associate (thinning_interval => self%output_config%thinning_interval, &
                 timestep => self%sim_config%timestep, &
                 simulation_times_for_output => self%output_config%simulation_times_for_output)
         if (thinning_interval == 0) then
            do i = self%counter, size(simulation_times_for_output)
               if (simulation_times_for_output(i) > simulation_time) then
                  self%simulation_time_for_next_output = simulation_times_for_output(i)
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
      integer :: i

      call output_helper%init(simdata%model%simulation_time, self%sim_config%timestep, self%simulation_time_for_next_output, &
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
      do i = 1, self%n_vars
         call output_helper%add_datum(self%output_files(i), "(F12.4)")
         ! If on volume or faces grid
         if (self%output_config%output_vars(i)%volume_grid) then
            ! Interpolate state on volume grid
            call self%grid%interpolate_from_vol(self%output_config%output_vars(i)%values, self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            call output_helper%add_data_array(self%output_files(i), i, self%last_iteration_data, values_on_zout, "(ES14.4)")
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Interpolate state on face grid
            call self%grid%interpolate_from_face(self%output_config%output_vars(i)%values, self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            call output_helper%add_data_array(self%output_files(i), i, self%last_iteration_data, values_on_zout, "(ES14.4)")
         else
            call output_helper%add_data_scalar(self%output_files(i), i, self%last_iteration_data, self%output_config%output_vars(i)%values_surf, "(ES14.4)")
         end if
         call output_helper%next_row(self%output_files(i))
      end do
   end subroutine

   subroutine output_helper_init(self, simulation_time, timestep, simulation_time_for_next_output, start_datum, &
                                 thinning_interval, simulation_times_for_output, counter)
      implicit none
      class(OutputHelper), intent(inout) :: self
      integer(8), intent(in) :: simulation_time
      integer, intent(in) :: timestep
      integer(8), intent(inout) :: simulation_time_for_next_output
      real(RK), intent(in) :: start_datum
      integer, intent(in) :: thinning_interval
      integer(8), dimension(:), allocatable, intent(in) :: simulation_times_for_output
      integer, intent(inout) :: counter
      integer :: i

      self%write_to_file = (simulation_time_for_next_output > simulation_time - timestep &
                      .and. simulation_time_for_next_output <= simulation_time)
      if (self%write_to_file) then
         counter = counter + 1
         self%w1 = (simulation_time - simulation_time_for_next_output) / real(timestep, RK)
         self%w0 = 1 - self%w1
         self%output_datum = datum(start_datum, simulation_time_for_next_output)
         if (thinning_interval == 0) then
            do i = counter, size(simulation_times_for_output)
               if (simulation_times_for_output(i) > simulation_time) then
                  simulation_time_for_next_output = simulation_times_for_output(i)
                  exit
               end if
            end do
         else
            simulation_time_for_next_output = simulation_time_for_next_output + thinning_interval * timestep
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
      write(80) self%counter, self%simulation_time_for_next_output
      call save_matrix(80, self%last_iteration_data)
   end subroutine


   !************************* Load ****************************

   ! Load logger state
   subroutine log_load(self)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      read(81) self%counter, self%simulation_time_for_next_output
      call read_matrix(81, self%last_iteration_data)
   end subroutine

end module
