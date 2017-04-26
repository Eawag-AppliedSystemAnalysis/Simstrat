module simstrat_output_module
  use simstrat_kinds
  use simstrat_model_module
  use utilities
  use csv_module
  implicit none

  private


  type, abstract, extends(SimstratOutputLogger), public :: SimstratThinningLogger
    private
    integer :: n_vars, n_depths, n_timepoints
    integer :: thinning_cycle_idx
    real(RK), dimension(:), allocatable :: out_depths
    real(RK), dimension(:), allocatable :: out_hight_above_sed
    character(len=:), allocatable :: PathOut
    type(string), dimension(:), allocatable :: out_names
    type(csv_file), dimension(:), allocatable :: output_files
    contains
      procedure, pass(self), public :: log => thinning_log
      procedure, pass(self), public :: close => log_close
      procedure(generic_write), deferred, pass(self) :: write
  end type

  type, extends(SimstratThinningLogger), public :: SimstratMemoryBoundThinningLogger
    private
    real(RK), dimension(:,:,:), allocatable :: memory_log
    integer :: curr_log_idx
    contains
      procedure, pass(self), public :: initialize => log_init_memory_bound
      procedure, pass(self) :: write => write_memory_bound
      procedure, pass(self), public :: close => log_close_memory_bound
  end type

  type, extends(SimstratThinningLogger), public :: SimstratOnTheFlyThinningLogger
    private
    contains
      procedure, pass(self), public :: initialize => log_init_on_the_fly
      procedure, pass(self) :: write => write_on_the_fly
  end type

  abstract interface
    subroutine generic_write(self, datum, state_vars)
      import SimstratThinningLogger, RK
      implicit none

      class(SimstratThinningLogger), intent(inout) :: self
      real(RK), intent(in) :: datum
      real(RK), dimension(:,:), intent(in) :: state_vars
    end subroutine
  end interface

contains
  !############################################
  ! General implementations for ThinningOutput
  !############################################

  subroutine load_out_depths(zoutName, depths)
    implicit none

    character(len=:), allocatable, intent(in) :: zoutName
    real(RK), dimension(:), allocatable, intent(out) :: depths

    type(csv_file) :: f
    logical :: status_ok

    call check_file_exists(zoutName)

    call f%read(zoutName, header_row=1, status_ok=status_ok)
    if(.not.status_ok) then
      call error('Unable to read output depths: '//zoutName)
      stop
    end if
    call f%get(1, depths, status_ok)
    call f%destroy()
  end subroutine

  subroutine log_init(self, model)
    implicit none
    class(SimstratThinningLogger), intent(inout) :: self
    class(SimstratModel), intent(in) :: model
    
    integer :: i

    self%n_vars = size(model%output_vars)
    self%PathOut = model%PathOut
    call load_out_depths(model%zoutName, self%out_depths)
    self%n_depths = size(self%out_depths)
    self%out_hight_above_sed = convert2hight_above_sed(self%out_depths, model%z_zero)
    self%n_timepoints = int((model%t_end-model%t_start)/(model%dt*model%thinning_interval/86400))
    self%thinning_cycle_idx = 1

    if(allocated(self%output_files)) deallocate(self%output_files)
    if(allocated(self%out_names)) deallocate(self%out_names)

    allocate(self%output_files(self%n_vars))
    allocate(self%out_names(self%n_vars))

    do i=1,self%n_vars
      self%out_names(i)%str = model%output_vars(i)%name
    end do
  end subroutine

  subroutine thinning_log(self, datum, model)
    implicit none

    class(SimstratThinningLogger), intent(inout) :: self
    class(SimstratModel), intent(in) :: model
    real(RK), intent(in) :: datum

    real(RK), dimension(self%n_vars, self%n_depths) :: state_vars
    integer :: i

    if(self%thinning_cycle_idx == model%thinning_interval) then
      do i=1,self%n_vars
        state_vars(i, 1:self%n_depths) = model%discretization%exportMFQ(self%out_hight_above_sed, model%output_vars(i)%value)
      end do
      call self%write(datum, state_vars)
      self%thinning_cycle_idx = 1
    else
      self%thinning_cycle_idx = self%thinning_cycle_idx +1
    end if

  end subroutine

  subroutine log_open(self)
    implicit none

    class(SimstratThinningLogger), intent(inout) :: self
    integer :: i
    logical :: status_ok

    real(RK), dimension(self%n_depths+1) :: header

    do i=1,self%n_vars
      call self%output_files(i)%open(self%PathOut//'/'//trim(self%out_names(i)%str)//'_out.dat', n_cols=self%n_depths+1, status_ok=status_ok)
      call self%output_files(i)%add('')
      call self%output_files(i)%add(self%out_depths(self%n_depths:1:-1), real_fmt='(F12.3)')
      call self%output_files(i)%next_row()
    end do
  end subroutine
  
  subroutine log_close(self)
    implicit none

    class(SimstratThinningLogger), intent(inout) :: self

    integer :: i
    logical :: status_ok
    do i=1,self%n_vars
      call self%output_files(i)%close(status_ok)
    end do
  end subroutine
  
  !################################################
  ! Specific implementations for MemoryBound Output
  !################################################


  subroutine log_init_memory_bound(self, model)
    implicit none

    class(SimstratMemoryBoundThinningLogger), intent(inout) :: self
    class(SimstratModel), intent(in) :: model
    logical :: status_ok

    call log_init(self, model)
    self%curr_log_idx = 1
    if(allocated(self%memory_log)) deallocate(self%memory_log)
    allocate(self%memory_log(self%n_timepoints, self%n_vars, self%n_depths+1))
  end subroutine

  subroutine write_memory_bound(self, datum, state_vars)
    implicit none

    class(SimstratMemoryBoundThinningLogger), intent(inout) :: self
    real(RK), intent(in) :: datum
    real(RK), dimension(:, :), intent(in) :: state_vars

    integer :: i

    do i=1,self%n_vars
      self%memory_log(self%curr_log_idx, i, 1) = datum
      self%memory_log(self%curr_log_idx, i, 2:self%n_depths+1) = state_vars(i,1:self%n_depths)
    end do

    self%curr_log_idx = self%curr_log_idx+1
  end subroutine

  subroutine log_close_memory_bound(self)
    implicit none

    class(SimstratMemoryBoundThinningLogger), intent(inout) :: self
    integer :: n, i

    write(*,*) ' Writing results...'

    call log_open(self)

    where(abs(self%memory_log) < 1E-20_RK) self%memory_log = 0.0_RK

    do n=1,self%n_timepoints
      do i=1,self%n_vars
        call self%output_files(i)%add(self%memory_log(n, i, 1), real_fmt='(F12.4)')
        call self%output_files(i)%add(self%memory_log(n, i, 2:self%n_depths+1), real_fmt='(ES12.4)')
        call self%output_files(i)%next_row()
      end do
    end do
    !flush all data

    call log_close(self)
  end subroutine

  !################################################
  ! Specific implementations for OnTheFly Output
  !################################################

  subroutine log_init_on_the_fly(self, model)
    implicit none

    class(SimstratOnTheFlyThinningLogger), intent(inout) :: self
    class(SimstratModel), intent(in) :: model

    call log_init(self, model)
    call log_open(self)

  end subroutine

  subroutine write_on_the_fly(self, datum, state_vars)
    implicit none

    class(SimstratOnTheFlyThinningLogger), intent(inout) :: self
    real(RK), intent(in) :: datum
    real(RK), dimension(:, :), intent(in) :: state_vars
    !workaround gfortran bug => cannot pass allocatable array to csv file
    real(RK), dimension(:), allocatable :: dummy


    integer :: i
    allocate(dummy(self%n_depths))
    do i=1,self%n_vars
      dummy(1:self%n_depths) = state_vars(i, 1:self%n_depths)
      where(abs(dummy) < 1E-20_RK) dummy = 0.0_RK
      call self%output_files(i)%add(datum, real_fmt='(F12.4)')
      call self%output_files(i)%add(dummy, real_fmt='(ES12.4)')
      call self%output_files(i)%next_row()
    end do
  end subroutine

end module
