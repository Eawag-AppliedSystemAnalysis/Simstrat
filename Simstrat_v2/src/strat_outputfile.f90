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
      class(OutputConfig), public, pointer   :: config
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
      procedure, pass(self), public :: initialize => log_init
      procedure, pass(self), public :: init_files => simple_init_files
      procedure, pass(self), public :: log => log_simple
      procedure, pass(self), public :: close => log_close
   end type

   ! Logger that interpolates on a specified output grid. Implementation not finished
   type, extends(SimpleLogger), public :: InterpolatingLogger
      private
   contains
      procedure, pass(self), public :: init_files => interpolating_init_files
      procedure, pass(self), public :: log => log_interpolating
   end type

contains

  ! Abstract interface definitions
   subroutine generic_log_init(self, config, grid)
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(OutputConfig), target :: config
      class(StaggeredGrid), target :: grid

   end subroutine

   subroutine generic_init_files(self, config, grid)
      implicit none
      class(SimstratOutputLogger), intent(inout) :: self
      class(OutputConfig), target :: config
      class(StaggeredGrid), target :: grid

   end subroutine

   ! Simple logger initialiation
   subroutine simple_init_files(self, config, grid)
      implicit none
      class(SimpleLogger), intent(inout) :: self
      class(OutputConfig), target :: config
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
         if (self%config%output_vars(i)%volume_grid) then
            !Variable on volume grid
            call self%output_files(i)%open(config%PathOut//'/'//trim(self%config%output_vars(i)%name)//'_out.dat', n_cols=grid%l_vol+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_volume(1:grid%ubnd_vol), real_fmt='(F12.3)')
         else
            ! Variable on face grid
            call self%output_files(i)%open(config%PathOut//'/'//trim(self%config%output_vars(i)%name)//'_out.dat', n_cols=grid%l_fce+1, status_ok=status_ok)
            call self%output_files(i)%add('')
            call self%output_files(i)%add(grid%z_face(1:grid%ubnd_fce), real_fmt='(F12.3)')
         end if
         call self%output_files(i)%next_row()
      end do
      self%thinning_interval = 1
      self%current_i = 0

   end subroutine

   subroutine interpolating_init_files(self, config, grid)
      implicit none
      class(InterpolatingLogger), intent(inout) :: self
      class(OutputConfig), target :: config
      class(StaggeredGrid), target :: grid

      logical :: status_ok
      integer :: n_depths, i
      n_depths = size(config%zout)
      self%n_depths = n_depths
      if (allocated(self%output_files)) deallocate (self%output_files)
      allocate (self%output_files(1:self%n_vars))

      do i = 1, self%n_vars
      call self%output_files(i)%open(config%PathOut//'/'//trim(self%config%output_vars(i)%name)//'_out.dat', n_cols=self%n_depths+1, status_ok=status_ok)
         call self%output_files(i)%add('')
         call self%output_files(i)%add(self%config%zout(self%n_depths:1:-1), real_fmt='(F12.3)')
         call self%output_files(i)%next_row()
      end do

   end subroutine

   subroutine log_init(self, config, grid)
      implicit none
      class(SimpleLogger), intent(inout) :: self
      class(OutputConfig), target :: config
      class(StaggeredGrid), target :: grid

      integer :: i, n_vars, n_depths

      self%config => config
      self%grid => grid
      n_vars = size(config%output_vars)
      self%n_vars = n_vars

      call self%init_files(config, grid)
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

   !################################################
   ! Specific implementations for OnTheFly Output
   !################################################

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
         test(i) = self%grid%z_zero + self%config%zout(self%n_depths - i + 1)
      end do

      do i = 0, self%n_vars - 1
         if (self%config%output_vars(i)%volume_grid) then
            call self%grid%interpolate_from_vol(test, self%config%output_vars(i)%values, self%n_depths, values_on_zout)
         else
            call self%grid%interpolate_from_face(test, self%config%output_vars(i)%values, self%n_depths, values_on_zout)
         end if

         where (abs(values_on_zout) < 1E-20_RK) values_on_zout = 0.0_RK
         call self%output_files(i)%add(datum, real_fmt='(F12.4)')
         call self%output_files(i)%add(values_on_zout, real_fmt='(ES12.4)')
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
         where (abs(self%config%output_vars(i)%values) < 1E-20_RK) &
            self%config%output_vars(i)%values = 0.0_RK
         call self%output_files(i)%add(datum, real_fmt='(F12.4)')
         call self%output_files(i)%add(self%config%output_vars(i)%values, real_fmt='(ES14.6)')
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
