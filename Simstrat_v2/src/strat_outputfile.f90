!     +---------------------------------------------------------------+
!     |  Interface and implementation of output data loggers
!     +---------------------------------------------------------------+

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
      integer, public :: thinning_interval, current_i
   contains
      procedure(generic_log_init), deferred, pass(self), public :: initialize
      procedure(generic_log), deferred, pass(self), public :: log   ! Method that is called each loop to write data
      procedure(generic_log_close), deferred, pass(self), public :: close
      procedure(generic_init_files), deferred, pass(self), public :: init_files
   end type

   !Simple logger = Logger that just writes out current state of variables,
   ! without any interpolation etc
   type, extends(SimstratOutputLogger), public :: SimpleLogger
      private
   contains
      procedure, pass(self), public :: initialize => log_init_simple
      procedure, pass(self), public :: init_files => init_files_simple
      procedure, pass(self), public :: log => log_simple
      procedure, pass(self), public :: close => log_close
   end type

   ! Logger that interpolates on a specified output grid. Implementation not finished
   type, extends(SimpleLogger), public :: InterpolatingLogger
      private
   contains
      procedure, pass(self), public :: initialize => log_init_interpolating
      procedure, pass(self), public :: init_files => init_files_interpolating
      procedure, pass(self), public :: log => log_interpolating
   end type

contains

  ! Abstract interface definitions
   subroutine generic_log_init(self, sim_config, output_config, grid)
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(SimConfig), target :: sim_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

   end subroutine

   subroutine generic_init_files(self, output_config, grid)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

   end subroutine

   ! Simple logger initialiation
   subroutine init_files_simple(self, output_config, grid)
      implicit none
      class(SimpleLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

      logical :: status_ok
      integer :: n_depths, i
      n_depths = grid%l_fce
      self%n_depths = n_depths

      ! Allocate output files
      if (allocated(self%output_files)) deallocate (self%output_files)
      allocate (self%output_files(self%n_vars))

      ! For each configured variable, create file and write header
      do i = 1, self%n_vars
         if (self%output_config%output_vars(i)%volume_grid) then
            !Variable on volume grid
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=grid%l_vol+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_volume(1:grid%ubnd_vol), real_fmt='(F12.3)')
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Variable on face grid
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=grid%l_fce+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_face(1:grid%ubnd_fce), real_fmt='(F12.3)')
         else  
            !Variable at surface
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=1 + 1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_face(grid%ubnd_fce), real_fmt='(F12.3)')             
         end if
         call self%output_files(i)%next_row()
      end do
      !self%thinning_interval = 1
      self%current_i = 0

   end subroutine

   subroutine init_files_interpolating(self, output_config, grid)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

      logical :: status_ok
      integer :: n_depths, i
      n_depths = size(output_config%zout)
      self%n_depths = n_depths
      if (allocated(self%output_files)) deallocate (self%output_files)
      allocate (self%output_files(1:self%n_vars))

      do i = 1, self%n_vars
         if (self%output_config%output_vars(i)%volume_grid) then
            !Variable on volume grid
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=grid%l_vol+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_volume(1:grid%ubnd_vol), real_fmt='(F12.3)')
         else if (self%output_config%output_vars(i)%face_grid) then
            ! Variable on face grid
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=grid%l_fce+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_face(1:grid%ubnd_fce), real_fmt='(F12.3)')
         else  
            !Variable at surface
            call self%output_files(i)%open(output_config%PathOut//'/'//trim(self%output_config%output_vars(i)%name)//'_out.dat', n_cols=1 + 1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_face(grid%ubnd_fce), real_fmt='(F12.3)')             
         end if
         call self%output_files(i)%next_row()
      end do

   end subroutine

   subroutine log_init_simple(self, sim_config, output_config, grid)
      implicit none
      class(SimpleLogger), intent(inout) :: self
      class(SimConfig), target :: sim_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

      integer :: i, n_vars, n_depths

      self%sim_config => sim_config
      self%output_config => output_config
      self%grid => grid
      n_vars = size(output_config%output_vars)
      self%n_vars = n_vars

      call self%init_files(output_config, grid)
   end subroutine

   subroutine log_init_interpolating(self, sim_config, output_config, grid)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(SimConfig), target :: sim_config
      class(OutputConfig), target :: output_config
      class(StaggeredGrid), target :: grid

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
            write(6,*) 'Error: simulation start time is larger than first output time'
            stop
        end if

        ! Allocate arrays for number of timesteps between output times and adjusted timestep
        allocate(output_config%n_timesteps_between_tout(n_output_times), output_config%adjusted_timestep(n_output_times))
        
        ! Compute number of timesteps between simulation start and first output time
        output_config%n_timesteps_between_tout(1) = (output_config%tout(1) - sim_config%start_datum)*86400/sim_config%timestep

        ! If number of timesteps = 0
        if (int(output_config%n_timesteps_between_tout(1)) == 0) then
            output_config%adjusted_timestep(1) = (output_config%tout(1) - sim_config%start_datum)*86400
            tout_test(1) = sim_config%start_datum + output_config%adjusted_timestep(1)
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
                tout_test(i) = tout_test(i-1) + output_config%adjusted_timestep(i)
                write(6,*) 'Warning: time interval for output is smaller than dt for iteration'
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
      write(6,*) 'Output times successfully read'
      write(6,*)

    call self%init_files(output_config, grid)
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


   subroutine log_interpolating(self, datum)
      implicit none

      class(InterpolatingLogger), intent(inout) :: self
      real(RK), intent(in) :: datum
      !workaround gfortran bug => cannot pass allocatable array to csv file
      real(RK), dimension(:), allocatable :: values, values_on_zout, test
      integer :: i

      allocate (values_on_zout(self%n_depths))
      allocate (test(self%n_depths))

      do i = 1, self%n_depths
         test(i) = self%grid%z_zero + self%output_config%zout(self%n_depths - i + 1)
      end do

      do i = 1, self%n_vars
        call self%output_files(i)%add(datum, real_fmt='(F12.4)')
        if (self%output_config%output_vars(i)%volume_grid) then
          call self%grid%interpolate_from_vol(test, self%output_config%output_vars(i)%values, self%n_depths, values_on_zout)
          call self%output_files(i)%add(self%output_config%output_vars(i)%values, real_fmt='(ES14.4)')
        else if (self%output_config%output_vars(i)%face_grid) then
          call self%grid%interpolate_from_face(test, self%output_config%output_vars(i)%values, self%n_depths, values_on_zout)
          call self%output_files(i)%add(self%output_config%output_vars(i)%values, real_fmt='(ES14.4)')
        else
          call self%output_files(i)%add(self%output_config%output_vars(i)%values_surf, real_fmt='(ES14.4)')
        end if
        call self%output_files(i)%next_row()
      end do

      deallocate (values_on_zout)
   end subroutine


   ! log current state
   subroutine log_simple(self, datum)
      implicit none

      class(SimpleLogger), intent(inout) :: self
      real(RK), intent(in) :: datum
      !workaround gfortran bug => cannot pass allocatable array to csv file
      !real(RK), dimension(:), allocatable :: values, values_on_zout,test
      integer :: i

      if (mod(self%current_i, self%thinning_interval) /= 0) then
         self%current_i = self%current_i + 1
         return ! dont log
      else
         self%current_i = self%current_i + 1
      end if

      ! For each variable, write state
      do i = 1, self%n_vars
        call self%output_files(i)%add(datum, real_fmt='(F12.4)')
        if (self%output_config%output_vars(i)%volume_grid .or. self%output_config%output_vars(i)%face_grid) then
           call self%output_files(i)%add(self%output_config%output_vars(i)%values, real_fmt='(ES14.4)')   
        else
           call self%output_files(i)%add(self%output_config%output_vars(i)%values_surf, real_fmt='(ES14.4)')
        end if
           call self%output_files(i)%next_row()
      end do

   end subroutine

   ! Close all files
   subroutine log_close(self)
      implicit none
      class(SimpleLogger), intent(inout) :: self

      integer :: i
      logical :: status_ok
      do i = 1, self%n_vars
         call self%output_files(i)%close (status_ok)
      end do

   end subroutine

end module
