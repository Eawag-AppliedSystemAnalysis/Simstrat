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
      class(FABMConfig), pointer :: fabm_cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param
      character(len=:), allocatable  :: file

      ! Variables that are used in between iteration. These used to be "save" variables
      real(RK) :: tb_start, tb_end !Start and end time
      real(RK) :: inp !Fraction of input absorption added to FABM bioshade feedback
      real(RK), dimension(:), allocatable :: z_absorb !Read depths
      real(RK), dimension(:), allocatable :: absorb_start, absorb_end !Interpolated start and end values
      real(RK), dimension(:), allocatable :: bg, bg_vol !Background absorption added to FABM bioshade feedback
      integer :: number_of_lines_read = 0
      integer :: eof, nval
      logical :: fabm_contribution

   contains
      procedure, pass :: init => absorption_init
      procedure, pass :: update => absorption_update
      procedure, pass :: save => absorption_save
      procedure, pass :: load => absorption_load
   end type

contains

   subroutine absorption_init(self, model_config, fabm_config, model_param, absorption_file, grid)
   ! Initialize absorption

      implicit none
      class(AbsorptionModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(FABMConfig), target :: fabm_config
      class(ModelParam), target :: model_param
      character(len=:), allocatable :: absorption_file

      self%cfg => model_config
      self%fabm_cfg => fabm_config
      self%param => model_param
      self%grid => grid
      self%file = absorption_file

      ! Allocate arrays (used to be "save" variables)
      allocate (self%z_absorb(grid%max_length_input_data))
      allocate (self%absorb_start(grid%nz_grid + 1))
      allocate (self%absorb_end(grid%nz_grid + 1))

      ! Check whether there is FABM contribution (bioshade feedback)
      ! No FABM contribution if input contribution is 1
      if (model_config%couple_fabm) then
         if ((fabm_config%input_extinction < 0.0_RK) .or. (fabm_config%input_extinction > 1.0_RK)) then
            call error('FABM Input Extinction Factor must be between 0 and 1')
         else if (fabm_config%input_extinction == 1.0_RK) then
            self%fabm_contribution = .false.
         else
            self%fabm_contribution = .true.
            if ((fabm_config%background_extinction < 0.0_RK) .or. (fabm_config%background_extinction > 1.0_RK)) then
               call error('FABM Background Extinction must be between 0 and 1')
            end if
         end if
      else
         self%fabm_contribution = .false.
      end if

      ! Get input contribution and background extinction in case of FABM contribution
      if (self%fabm_contribution) then
         allocate (self%bg(grid%nz_grid + 1))
         allocate (self%bg_vol(grid%nz_grid))
         self%inp = fabm_config%input_extinction
         self%bg(:) = fabm_config%background_extinction
         self%bg_vol(:) = fabm_config%background_extinction
      end if
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

   subroutine absorption_update(self, state)
   ! Update absorption
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

         ! In case of FABM contribution
         if (self%fabm_contribution) then      
            ! If there is no input contribution, set absorb to the biogeochemical contribution + background extinction
            ! Absorption passed to FABM is just background extinction
            ! No further calculation necessary
            if (self%inp == 0.0_RK) then
               state%absorb(1:nz+1) = state%absorb_from_fabm(1:nz+1) + self%bg(1:nz+1)
               state%absorb_to_fabm(1:nz) = self%bg_vol(1:nz)
               return
            end if
         end if

         if (state%first_timestep) then ! First iteration
            open (30, status='old', file=self%file)
            if (self%number_of_lines_read > 0) then
               do i = 1, self%number_of_lines_read
                  read (30, *, end=9) ! skip over already read lines
               end do
            else
               eof = 0

               !Read depths: first line are names, second line is number of depths available
               read (30, *, end=9)
               call count_read(self)
               read (30, *, end=9) nval
               call count_read(self)
               read (30, *, end=9) dummy, (z_absorb(i), i=1, nval)
               call count_read(self)

               !Make depths positives
               do i = 1, nval
                  z_absorb(i) = abs(z_absorb(i))
               end do

               !Read first values
               read (30, *, end=9) tb_start, (absorb_read_start(i), i=1, nval)
               call count_read(self)

               if (state%datum < tb_start) call warn('First light attenuation date after simulation start time.')

               !Interpolate absorb_read_start on z_absorb onto faces of grid
               call self%grid%interpolate_to_face(z_absorb, absorb_read_start, nval, absorb_start)

               read (30, *, end=7) tb_end, (absorb_read_end(i), i=1, nval)
               call count_read(self)

               ! Write to console that file was successfully read
               call ok('Absorption input file successfully read')

               ! Do the same for absorb_read_end
               call self%grid%interpolate_to_face(z_absorb, absorb_read_end, nval, absorb_end)
            end if
         end if

         if (state%datum <= tb_start .or. eof == 1) then !If datum before first date or end of file reached
            goto 8
         else
            do while (state%datum > tb_end) !Move to appropriate interval to get correct value
               tb_start = tb_end
               absorb_start(1:nz+1) = absorb_end(1:nz+1)
               !Read next value
               read (30, *, end=7) tb_end, (absorb_read_end(i), i=1, nval)
               call count_read(self)
               call self%grid%interpolate_to_face(z_absorb, absorb_read_end, nval, absorb_end)
            end do
            !Linearly interpolate value at correct datum (for all depths)
            state%absorb(1:nz+1) = absorb_start(1:nz+1) + (state%datum - tb_start)/(tb_end - tb_start)*(absorb_end(1:nz+1)  - absorb_start(1:nz+1))
            state%absorb(1:nz+1) = self%param%p_absorb*state%absorb(1:nz+1)
         end if
         ! In case of FABM contribution add the contribution and calculate absorption passed to FABM
         if (self%fabm_contribution) then
            ! Interpolate to volume for absorption passed to FABM array, interpolate_to_vol handles boundaries
            call self%grid%interpolate_to_vol(self%grid%z_face, state%absorb, nz + 1, state%absorb_to_fabm)
            ! Add biogeochemical contribution (absorb_from_fabm) and background extinction to absorption
            ! Linear regression between input absorbtion and biogeochemical absorption with factor inp
            state%absorb(1:nz+1) = self%inp * state%absorb(1:nz+1) + (1 - self%inp) * state%absorb_from_fabm(1:nz+1) + self%bg(1:nz+1)
            ! Absorption passed to FABM such that adding bgc contribution once (done by bgc model) results in absorption above
            ! Assuming that biogeochemical contribution does not change in one timestep
            state%absorb_to_fabm(1:nz) = self%inp * state%absorb_to_fabm(1:nz) - self%inp * state%absorb_from_fabm_vol(1:nz) + self%bg_vol(1:nz)
         end if
         return

7        eof = 1
         if(state%datum>tb_start) call warn('Last light attenuation date before simulation end time.')

8        state%absorb(1:nz+1) = absorb_start(1:nz+1)           !Take first value of current interval
         ! In case of FABM contribution add the contribution and calculate absorption passed to FABM
         if (self%fabm_contribution) then
            call self%grid%interpolate_to_vol(self%grid%z_face, state%absorb, nz + 1, state%absorb_to_fabm)
            state%absorb(1:nz+1) = self%inp * state%absorb(1:nz+1) + (1 - self%inp) * state%absorb_from_fabm(1:nz+1) + self%bg
            state%absorb_to_fabm(1:nz) = self%inp * state%absorb_to_fabm(1:nz) - self%inp * state%absorb_from_fabm_vol(1:nz) + self%bg_vol(1:nz)
         end if
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
