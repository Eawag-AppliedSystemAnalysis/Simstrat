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
      integer :: unit = 0 ! Unit of the opened absorption file (0 if not open)

      ! Variables that are used in between iteration. These used to be "save" variables
      real(RK) :: tb_start, tb_end !Start and end time
      real(RK), dimension(:), allocatable :: z_absorb !Read depths
      real(RK), dimension(:), allocatable :: absorb_start, absorb_end !Start and end values at the read depths z_absorb
      integer :: number_of_lines_read = 0
      integer :: eof, nval

   contains
      procedure, pass :: init => absorption_init
      procedure, pass :: update => absorption_update
      procedure, pass :: save => absorption_save
      procedure, pass :: load => absorption_load
      procedure, pass :: deallocate => absorption_deallocate
   end type

contains

   subroutine absorption_init(self, model_config, model_param, absorption_file, grid)
   ! Initialize absorption

      implicit none
      class(AbsorptionModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param
      character(len=:), allocatable :: absorption_file

      ! Close files and free memory of a previous initialization (if Simstrat is initialized several times)
      call self%deallocate()
      self%number_of_lines_read = 0

      self%cfg => model_config
      self%param => model_param
      self%grid => grid
      self%file = absorption_file

      ! Allocate arrays (used to be "save" variables)
      allocate (self%z_absorb(grid%max_length_input_data))
      allocate (self%absorb_start(grid%max_length_input_data))
      allocate (self%absorb_end(grid%max_length_input_data))
   end subroutine

   subroutine absorption_save(self)
      implicit none
      class(AbsorptionModule) :: self

      write (80) self%number_of_lines_read
      write (80) self%tb_start, self%tb_end
      write (80) self%eof, self%nval
      call save_array(80, self%z_absorb)
      call save_array(80, self%absorb_start)
      call save_array(80, self%absorb_end)
   end subroutine

   subroutine absorption_load(self)
      implicit none
      class(AbsorptionModule) :: self

      read (81) self%number_of_lines_read
      read (81) self%tb_start, self%tb_end
      read (81) self%eof, self%nval
      call read_array(81, self%z_absorb)
      call read_array(81, self%absorb_start)
      call read_array(81, self%absorb_end)
   end subroutine

   ! Deallocate all data allocated in self and close the absorption file
   subroutine absorption_deallocate(self)
      class(AbsorptionModule) :: self
      logical :: is_open

      if (self%unit /= 0) then
         inquire (unit=self%unit, opened=is_open)
         if (is_open) close (self%unit)
         self%unit = 0
      end if

      if (allocated(self%file)) deallocate(self%file)
      if (allocated(self%z_absorb)) deallocate(self%z_absorb)
      if (allocated(self%absorb_start)) deallocate(self%absorb_start)
      if (allocated(self%absorb_end)) deallocate(self%absorb_end)
   end subroutine

   subroutine absorption_update(self, state)
   ! Update absorption
      implicit none
      class(AbsorptionModule) :: self
      class(ModelState) :: state

      ! Local Variables
      real(RK) :: dummy !Read depths
      real(RK) :: absorb_cur(self%grid%max_length_input_data) !Values at the current time at the read depths
      integer :: i

      ! Associations for easier readability / comparability to old code
      associate (tb_start=>self%tb_start, &
                 tb_end=>self%tb_end, &
                 z_absorb=>self%z_absorb, &
                 absorb_start=>self%absorb_start, &
                 absorb_end=>self%absorb_end, &
                 eof=>self%eof, &
                 nval=>self%nval)

         if (state%first_timestep) then ! First iteration
            open (newunit=self%unit, status='old', file=self%file)
            if (self%number_of_lines_read > 0) then
               do i = 1, self%number_of_lines_read
                  read (self%unit, *, end=9) ! skip over already read lines
               end do
            else
               eof = 0

               !Read depths: first line are names, second line is number of depths available
               read (self%unit, *, end=9)
               call count_read(self)
               read (self%unit, *, end=9) nval
               call count_read(self)
               read (self%unit, *, end=9) dummy, (z_absorb(i), i=1, nval)
               call count_read(self)

               !Make depths positives
               do i = 1, nval
                  z_absorb(i) = abs(z_absorb(i))
               end do

               !Read first values
               read (self%unit, *, end=9) tb_start, (absorb_start(i), i=1, nval)
               call count_read(self)

               if (state%datum < tb_start) call warn('First light attenuation date after simulation start time.')

               read (self%unit, *, end=7) tb_end, (absorb_end(i), i=1, nval)
               call count_read(self)

               ! Write to console that file was successfully read
               call ok('Absorption input file successfully read')
            end if
         end if

         if (state%datum <= tb_start .or. eof == 1) then !If datum before first date or end of file reached
            absorb_cur(1:nval) = absorb_start(1:nval) !Take first value of current interval
         else
            do while (state%datum > tb_end) !Move to appropriate interval to get correct value
               tb_start = tb_end
               absorb_start(1:nval) = absorb_end(1:nval)
               !Read next value
               read (self%unit, *, end=7) tb_end, (absorb_end(i), i=1, nval)
               call count_read(self)
            end do
            !Linearly interpolate value at correct datum (for all depths)
            absorb_cur(1:nval) = absorb_start(1:nval) + (state%datum - tb_start)/(tb_end - tb_start)*(absorb_end(1:nval) - absorb_start(1:nval))
         end if
         goto 8

7        eof = 1
         if(state%datum>tb_start) call warn('Last light attenuation date before simulation end time.')
         absorb_cur(1:nval) = absorb_start(1:nval) !Take first value of current interval

         ! Interpolate from the read depths (below the current surface) onto the occupied faces and volumes (grid order)
8        absorb_cur(1:nval) = self%param%p_absorb*absorb_cur(1:nval)
         call self%grid%interpolate_depth_to_face(z_absorb, absorb_cur, nval, state%absorb)
         call self%grid%interpolate_depth_to_vol(z_absorb, absorb_cur, nval, state%absorb_vol)
         return

9        call error('Reading absorption file (no data found).')

      end associate
   end subroutine

   subroutine count_read(self)
      implicit none
      class(AbsorptionModule) :: self
      self%number_of_lines_read = self%number_of_lines_read + 1
   end subroutine

end module
