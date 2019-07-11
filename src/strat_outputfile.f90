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
      integer, public :: n_vars
      integer, public :: n_vars_AED2

   contains
      procedure(generic_log_init), deferred, pass(self), public :: initialize
      procedure(generic_init_files), deferred, pass(self), public :: init_files
      procedure(generic_log), deferred, pass(self), public :: log   ! Method that is called each loop to write data
      procedure(generic_log_close), deferred, pass(self), public :: close
   end type

   ! Logger that interpolates on a specified output grid at specific times.
   type, extends(SimstratOutputLogger), public :: InterpolatingLogger
      private
   contains
      procedure, pass(self), public :: initialize => log_init_interpolating
      procedure, pass(self), public :: init_files => init_files_interpolating
      procedure, pass(self), public :: log => log_interpolating
      procedure, pass(self), public :: close => log_close
   end type

contains

   ! Abstract interface definitions

   subroutine generic_log_init(self, state, sim_config, model_config, aed2_config, output_config, grid)
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(ModelState), target :: state
      class(SimConfig), target :: sim_config
      class(ModelConfig), target :: model_config
      class(AED2Config), target :: aed2_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

   end subroutine

   subroutine generic_init_files(self, output_config, grid)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

   end subroutine

   subroutine generic_log(self, datum)
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      real(RK), intent(in) :: datum
   end subroutine

   subroutine generic_log_close(self)
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
   end subroutine

   !************************* Init logging ****************************

   ! Init logging for interpolating logger
   subroutine log_init_interpolating(self, state, sim_config, model_config, aed2_config, output_config, grid)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(ModelState), target :: state
      class(SimConfig), target :: sim_config
      class(ModelConfig), target :: model_config
      class(AED2Config), target :: aed2_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

      integer :: i, j, n_output_times
      real(RK), dimension(size(output_config%tout)) :: tout_test ! Array to test if computed tout with
      ! adjusted timestep is the same as the tout given in file

      self%sim_config => sim_config
      self%model_config => model_config
      self%aed2_config => aed2_config
      self%output_config => output_config
      self%grid => grid
      self%n_vars = size(output_config%output_vars)

      if (self%model_config%couple_aed2) then
        ! allocate AED2 output structure
        allocate (output_config%output_vars_aed2) ! We don't know yet how many variables
        output_config%output_vars_aed2%names => state%AED2_names
        output_config%output_vars_aed2%values => state%AED2_state
        self%n_vars_AED2 = state%n_AED2
      end if

      ! If output times are given in file
      if (output_config%thinning_interval == 0) then
      ! Number of output times specified in file
      n_output_times = size(output_config%tout)
         ! If Simulation start larger than output times, abort.
         if (sim_config%start_datum > output_config%tout(1)) then
            call error('Simulation start time is larger than first output time.')
         end if

         ! Allocate arrays for number of timesteps between output times and adjusted timestep
         allocate(output_config%n_timesteps_between_tout(n_output_times), output_config%adjusted_timestep(n_output_times))

         ! Compute number of timesteps between simulation start and first output time
         output_config%n_timesteps_between_tout(1) = (output_config%tout(1) - sim_config%start_datum)*86400/sim_config%timestep

         ! If number of timesteps = 0
         if (int(output_config%n_timesteps_between_tout(1)) == 0) then
            output_config%adjusted_timestep(1) = (output_config%tout(1) - sim_config%start_datum)*86400
            tout_test(1) = sim_config%start_datum + output_config%adjusted_timestep(1)/86400

            ! Set to 1, as the adjusted timestep has to be used once, otherwise the model does not advance
            output_config%n_timesteps_between_tout(1) = 1
            call warn('First output time is equal to simulation start time')
         else
            ! If number of timesteps > 0
            output_config%adjusted_timestep(1) = ((output_config%tout(1) - sim_config%start_datum)*86400)/int(output_config%n_timesteps_between_tout(1))
            tout_test(1) = sim_config%start_datum
            ! Add up adjusted timestep, the resulting tout_test(1) should be equal to tout(1)
            do j = 1, int(output_config%n_timesteps_between_tout(1))
                tout_test(1) = tout_test(1) + output_config%adjusted_timestep(1)/86400
            end do
         end if

         ! Compute number of timesteps between subsequent output times
         do i=2,n_output_times
            output_config%n_timesteps_between_tout(i) = (output_config%tout(i) - tout_test(i-1))*86400/sim_config%timestep

            ! If number of timesteps = 0
            if (int(output_config%n_timesteps_between_tout(i))==0) then
               output_config%adjusted_timestep(i) = (output_config%tout(i) - tout_test(i-1))*86400
               tout_test(i) = tout_test(i-1) + output_config%adjusted_timestep(i)/86400

               ! Set to 1, as the adjusted timestep has to be used once, otherwise the model does not advance
               output_config%n_timesteps_between_tout(i) = 1
               call warn('At least one time interval for the model output is smaller than the simulation timestep')
            else
               ! If number of timesteps > 0
               output_config%adjusted_timestep(i) = ((output_config%tout(i) - tout_test(i-1))*86400)/int(output_config%n_timesteps_between_tout(i))
               tout_test(i) = tout_test(i-1)
               ! Add up adjusted timesteps. The resulting tout_test(i) should be equal to tout(i)
               do j = 1, int(output_config%n_timesteps_between_tout(i))
                  tout_test(i) = tout_test(i) + output_config%adjusted_timestep(i)/86400
               end do
            end if
         end do
      end if

      if (output_config%thinning_interval>1) write(6,*) 'Interval [days]: ',output_config%thinning_interval*sim_config%timestep/86400.
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

      call self%init_files(output_config, grid)
   end subroutine

   !************************* Init files ****************************

   ! Initialize files for interpolating logger
   subroutine init_files_interpolating(self, output_config, grid)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

      logical :: status_ok, exist_output_folder
      integer :: i, ppos
      character(len=256) :: mkdirCmd, output_folder

      self%n_depths = size(output_config%zout)

      ! Check if output directory exists
      inquire(file=output_config%PathOut,exist=exist_output_folder)

      ! Create Simstrat output folder if it does not exist
      if(.not.exist_output_folder) then
         call warn('Result folder does not exist, create folder...')
         mkdirCmd = 'mkdir '//trim(output_config%PathOut)
         call execute_command_line(mkdirCmd)
      end if

      if (allocated(self%output_files)) deallocate (self%output_files)
      allocate (self%output_files(1:self%n_vars + self%n_vars_AED2))

      do i = 1, self%n_vars
         if (self%output_config%output_vars(i)%volume_grid) then
            !Variable on volume grid
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=self%n_depths+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Variable on face grid
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=self%n_depths+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(self%output_config%zout, real_fmt='(F12.3)')
         else
            !Variable at surface
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=1 + 1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_face(grid%ubnd_fce), real_fmt='(F12.3)')
         end if
         call self%output_files(i)%next_row()
      end do

      ! AED2 part
      if (self%model_config%couple_aed2) then
        ! Check if output directory exists
        inquire(file=self%aed2_config%path_aed2_output//'/',exist=exist_output_folder)
        ! Create AED2 output folder if it does not exist
        if(.not. exist_output_folder) then
          call warn('AED2 result folder does not exist, create folder...')
          ppos = scan(trim(self%aed2_config%path_aed2_output),"/", BACK= .true.)
          if ( ppos > 0 ) output_folder = self%aed2_config%path_aed2_output(1:ppos - 1)
          mkdirCmd = 'mkdir '//trim(output_folder)
          write(6,*) mkdirCmd
          call execute_command_line(mkdirCmd)
        end if

        do i = 1, self%n_vars_AED2
            call self%output_files(i + self%n_vars)%open(self%aed2_config%path_aed2_output//'/'//trim(self%output_config%output_vars_aed2%names(i))//'_out.dat', n_cols=self%n_depths + 1, status_ok=status_ok)
            call self%output_files(i + self%n_vars)%add('')
            call self%output_files(i + self%n_vars)%add(self%output_config%zout, real_fmt='(F12.3)')
            call self%output_files(i + self%n_vars)%next_row()
            end do
        end if

   end subroutine

   !************************* Log ****************************

   ! Log current state on interpolated grid at specific times
   subroutine log_interpolating(self, datum)
      implicit none

      class(InterpolatingLogger), intent(inout) :: self
      real(RK), intent(in) :: datum
      !workaround gfortran bug => cannot pass allocatable array to csv file
      real(RK), dimension(self%n_depths) :: values_on_zout
      integer :: i, j

      do i = 1, self%n_vars
         ! Write datum
         call self%output_files(i)%add(datum, real_fmt='(F12.4)')
         ! If on volume or faces grid
         if (self%output_config%output_vars(i)%volume_grid) then
            ! Interpolate state on volume grid
            call self%grid%interpolate_from_vol(self%output_config%output_vars(i)%values, self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            ! Write state
            call self%output_files(i)%add(values_on_zout, real_fmt='(ES14.4)')
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Interpolate state on face grid
            call self%grid%interpolate_from_face(self%output_config%output_vars(i)%values, self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
            ! Write state
            call self%output_files(i)%add(values_on_zout, real_fmt='(ES14.4)')
         else
            ! If only value at surface
            call self%output_files(i)%add(self%output_config%output_vars(i)%values_surf, real_fmt='(ES14.4)')
         end if
         ! Advance to next row
         call self%output_files(i)%next_row()
      end do

      if (self%model_config%couple_aed2) then
        do i = 1, self%n_vars_AED2
          ! Write datum
          call self%output_files(i + self%n_vars)%add(datum, real_fmt='(F12.4)')
          ! Interpolate state on volume grid
          call self%grid%interpolate_from_vol(self%output_config%output_vars_aed2%values(:,i), self%output_config%zout, values_on_zout, self%n_depths, self%output_config%output_depth_reference)
          ! Write state
          call self%output_files(i + self%n_vars)%add(values_on_zout, real_fmt='(ES14.4)')
          ! Advance to next row
          call self%output_files(i + self%n_vars)%next_row()
        end do
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

   end subroutine

end module
