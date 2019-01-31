!<    +---------------------------------------------------------------+
!     |  Absorption module
!     |  - reads and processes absorption input file
!<    +---------------------------------------------------------------+

module strat_absorption
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use utilities
   implicit none
   private

   type, public :: AbsorptionModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param
      character(len=:), allocatable  :: file

      ! Variables that are used in between iteration. These used to be "save" variables
      real(RK) :: tb_start, tb_end !Start and end time
      real(RK), dimension(:), allocatable :: z_absorb !Read depths
      real(RK), dimension(:), allocatable :: absorb_start, absorb_end !Interpolated start and end values
      integer :: eof, nval

   contains
      procedure, pass :: init => absorption_init
      procedure, pass :: update => absorption_update
   end type

contains

   subroutine absorption_init(self, model_config, model_param, absorption_file, grid)
   !! initialize absorption

      implicit none
      class(AbsorptionModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param
      character(len=:), allocatable :: absorption_file

      self%cfg => model_config
      self%param => model_param
      self%grid => grid
      self%file = absorption_file

      ! Allocate arrays (used to be "save" variables)
      allocate (self%z_absorb(grid%max_length_input_data))
      allocate (self%absorb_start(grid%nz_grid))
      allocate (self%absorb_end(grid%nz_grid))

   end subroutine

   subroutine absorption_update(self, state)
   !! update absorption
      implicit none
      class(AbsorptionModule) :: self
      class(ModelState) :: state

      ! Local Variables
      real(RK) :: dummy !Read depths
      real(RK) :: absorb_read_start(self%grid%max_length_input_data), absorb_read_end(self%grid%max_length_input_data) !Read start and end values
      integer :: i

      ! Associations for easier readability / comparability to old code
      associate (tb_start=>self%tb_start, &
                 tb_end=>self%tb_end, &
                 z_absorb=>self%z_absorb, &
                 absorb_start=>self%absorb_start, &
                 absorb_end=>self%absorb_end, &
                 eof=>self%eof, &
                 nval=>self%nval, &
                 nz=>self%grid%nz_occupied)

         if (state%first_timestep) then ! First iteration
            open (30, status='old', file=self%file)
            eof = 0

            !Read depths: first line are names, second line is number of depths available
            read (30, *, end=9)
            read (30, *, end=9) nval
            read (30, *, end=9) dummy, (z_absorb(i), i=1, nval)

            !Make depths positives
            do i = 1, nval
               z_absorb(i) = abs(z_absorb(i))
            end do

            !Read first values
            read (30, *, end=9) tb_start, (absorb_read_start(i), i=1, nval)
            if (state%datum < tb_start) call warn('First light attenuation date after simulation start time.')

            !Interpolate absorb_read_start on z_absorb onto faces of grid
            call self%grid%interpolate_to_face(z_absorb, absorb_read_start, nval, absorb_start)

            read (30, *, end=7) tb_end, (absorb_read_end(i), i=1, nval)

            ! Write to console that file was successfully read
            call ok('Absorption input file successfully read')

            ! Do the same for absorb_read_end
            call self%grid%interpolate_to_face(z_absorb, absorb_read_end, nval, absorb_end)

         end if

         if (state%datum <= tb_start .or. eof == 1) then !If datum before first date or end of file reached
            goto 8
         else
            do while (state%datum > tb_end) !Move to appropriate interval to get correct value
               tb_start = tb_end
               absorb_start(1:nz) = absorb_end(1:nz)
               !Read next value
               read (30, *, end=7) tb_end, (absorb_read_end(i), i=1, nval)
               call self%grid%interpolate_to_face(z_absorb, absorb_read_end, nval, absorb_end)
            end do
            !Linearly interpolate value at correct datum (for all depths)
            state%absorb(1:nz) = absorb_start(1:nz) + (state%datum - tb_start)*(absorb_end(1:nz) &
                                 - absorb_start(1:nz))/(tb_end - tb_start)
         end if
         return

7        eof = 1
         if(state%datum>tb_start) call warn('Last light attenuation date before simulation end time.')

8        state%absorb(1:nz) = absorb_start(1:nz)           !Take first value of current interval
         return

9        call error('Reading light attenuation file (no data found).')

      end associate
   end subroutine

end module
