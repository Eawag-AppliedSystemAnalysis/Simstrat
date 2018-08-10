!     +---------------------------------------------------------------+
!     |  Lateral module
!     |     Reads and processes inflows/outflows such that
!     |     Advection can be calculated in the next step
!     |     There are at least two different implementations possible:
!     |         - LateralRhoModule : Inflow plunges according to density
!     |         - LateralModule:     Inflow affects layer as configured in file
!     +---------------------------------------------------------------+


module strat_lateral
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use utilities
   implicit none
   private

   ! Generic base class for both modules (LateralRho and Normal)
   type, abstract, public :: GenericLateralModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param

      ! Variables that where either marked with "save" before, or that have been
      ! global, but only used in the lateral environment:
      real(RK), dimension(:, :), allocatable   :: z_Inp, Q_start, Qs_start, Q_end, Qs_end, Q_read_start, Q_read_end
      real(RK), dimension(:, :), allocatable   :: Inp_read_start, Inp_read_end, Qs_read_start, Qs_read_end
      real(RK) :: tb_start(1:4), tb_end(1:4) ! Input depths, start time, end time
      integer :: eof(1:4)
      integer :: nval(1:4), nval_deep(1:4), nval_surface(1:4) ! Number of values

   contains
      procedure, pass :: init => lateral_generic_init
      procedure(lateral_generic_update), deferred, pass :: update
      !procedure, pass :: surface_flow => lateral_generic_surface_flow
   end type

   ! Subclasses
   type, extends(GenericLateralModule), public :: LateralRhoModule
   contains
      procedure, pass, public :: update => lateral_rho_update
   end type

   type, extends(GenericLateralModule), public:: LateralModule
   contains
      procedure, pass, public :: update => lateral_update
   end type

contains
   subroutine lateral_generic_update(self, state)
      implicit none
      class(GenericLateralModule) :: self
      class(ModelState) :: state
   end subroutine

   subroutine lateral_generic_init(self, model_config, model_param, grid)
      implicit none
      class(GenericLateralModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param

      self%cfg => model_config
      self%param => model_param
      self%grid => grid
   end subroutine

   ! Calculate surface flow
!    subroutine lateral_generic_surface_flow(self, Q_surface, Q_total, depth_surface, i)
!       implicit none
!       class(GenericLateralModule) :: self
!       ! Global declarations
!       real(RK), intent(in) :: Q_surface, depth_surface
!       real(RK), intent(inout) :: Q_total(1:)
!       integer, intent(in) :: i

!       ! Local variables
!       real(RK) :: Q_surf_integral, add_surf
!       integer :: j

!       Q_surf_integral = Q_surface

!       do j = 1, self%grid%ubnd_fce

!          if (((Q_surf_integral > 0) .and. (.not. i == 2)) .or. ((Q_surf_integral < 0) .and. (i == 2))) then
!             Q_total(self%grid%ubnd_fce + 1 - j) = Q_total(self%grid%ubnd_fce + 1 - j) + Q_surf_integral
!          else
!             exit
!          end if

!          add_surf = Q_surface*self%grid%h(self%grid%ubnd_fce + 1 - j)/depth_surface
!          Q_surf_integral = Q_surf_integral - add_surf
!       end do
!       return
!    end subroutine

   ! Implementation for lateral rho
   subroutine lateral_rho_update(self, state)
      implicit none
      class(LateralRhoModule) :: self
      class(ModelState) :: state
      write (*, *) "lateral_rho_update not implemented"
   end subroutine

   ! "Normal" Implementation
   ! Note that this method has been taken from old simstrat and needs refactoring
   subroutine lateral_update(self, state)
      implicit none
      class(LateralModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      !real(RK) :: Q_read_start(1:4, 1:state%nz_input), Q_read_end(1:4, 1:state%nz_input) ! Integrated input
      !real(RK) :: Qs_read_start(1:4, 1:state%nz_input), Qs_read_end(1:4, 1:state%nz_input) ! Integrated surface input
      real(RK) :: dummy

      integer :: i, j
      integer :: fnum(1:4) ! File number
      character(len=20) :: fname(1:4)

      associate (datum=>state%datum, &
                 idx=>state%model_step_counter, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step, Q_vert is the integrated net water input
                 Q_vert=>state%Q_vert, &
                 grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)

         fname = ['inflow           ', 'outflow          ', 'input temperature', 'input salinity   ']
         fnum = [41, 42, 43, 44]

         ! FB 2016: Major revision to include surface inflow
         ! Do this for inflow, outflow, temperature and salinity
         do i = 1, 4
            if (idx == 1) then ! First iteration

               ! Allocate arrays if not already done, this saves memory compared to declaring with max_length_input_data
               if (.not. allocated(self%z_Inp)) allocate (self%z_Inp(1:4, 1:state%nz_input)) ! Input depths
               if (.not. allocated(self%Inp_read_start)) allocate (self%Inp_read_start(1:4, 1:state%nz_input))
               if (.not. allocated(self%Inp_read_end)) allocate (self%Inp_read_end(1:4, 1:state%nz_input))
               if (.not. allocated(self%Q_read_start)) allocate (self%Q_read_start(1:4, 1:state%nz_input))
               if (.not. allocated(self%Q_read_end)) allocate (self%Q_read_end(1:4, 1:state%nz_input))               
               if (.not. allocated(self%Qs_read_start)) allocate (self%Qs_read_start(1:4, 1:state%nz_input))
               if (.not. allocated(self%Qs_read_end)) allocate (self%Qs_read_end(1:4, 1:state%nz_input))
               if (.not. allocated(self%Q_start)) allocate (self%Q_start(1:4, 1:grid%nz_grid + 1)) ! Input interpolated on grid
               if (.not. allocated(self%Q_end)) allocate (self%Q_end(1:4, 1:grid%nz_grid + 1)) ! Input interpolated on grid
               if (.not. allocated(self%Qs_start)) allocate (self%Qs_start(1:4, 1:grid%nz_grid + 1)) ! Surface input interpolated on grid
               if (.not. allocated(self%Qs_end)) allocate (self%Qs_end(1:4, 1:grid%nz_grid + 1)) ! Surface input interpolated on grid

               self%Qs_start(i, :) = 0.0_RK
               self%Qs_end(i, :) = 0.0_RK

               ! Open file and start to read
               self%eof(i) = 0
               read (fnum(i), *, end=9) ! Skip first row: description of columns

               ! Read number of deep and surface inflows
               read (fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)

               ! Total number of values to read
               self%nval(i) = self%nval_deep(i) + self%nval_surface(i)

               ! Read input depths and convert coordinate system
               read (fnum(i), *, end=9) dummy, (self%z_Inp(i, j), j=1, self%nval(i))

               self%z_Inp(i, 1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i, 1:self%nval_deep(i))
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))
               !write(6,*) self%nval_deep(i), self%nval_surface(i), self%nval(i), grid%lake_level
               !write(6,*) self%z_Inp(i,1:self%nval_deep(i)), self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i))
               ! Allocate array for depth of the surface inflow (only done for i=1)
               !if (.not. allocated(self%depth_surfaceFlow)) allocate (self%depth_surfaceFlow(1:4, 1:(self%nval_surface(i) + 2)))

               ! Read first line
               read (fnum(i), *, end=9) self%tb_start(i), (self%Inp_read_start(i, j), j=1, self%nval(i))

               ! Integrate the inflow (direct interpolation of inflow is not correct)
               call Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), self%Q_read_start(i, :), self%nval_deep(i))
               call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))
               ! Very important: once the inflowing quantitiy is integrated, it necessarily has to be
               ! interpolated to the face grid
               if (self%nval_surface(i)>0) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               end if


               ! Read second line and treatment of deep inflow
               read (fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
               call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
               call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))

               if (self%nval_surface(i)>0) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if
               ! Add surface flow for both in- and outflow

!                if (self%nval_surface(i) > 0) then
!                   do j = 1, self%nval_surface(i)
!                      self%depth_surfaceFlow(i, j) = self%z_Inp(i, self%nval_deep(i) - 1 + j)
!                      call self%surface_flow(self%Inp_read_start(i, self%nval_deep(i) - 1 + j), self%Q_start(i, :), self%depth_surfaceFlow(i, j), i)
!                      call self%surface_flow(self%Inp_read_end(i, self%nval_deep(i) - 1 + j), self%Q_end(i, :), self%depth_surfaceFlow(i, j), i)
!                   end do
!                end if
               write(6,*) "--Input file successfully read: ",fname(i)
            end if ! idx==1

            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep them at the surface
            if ((.not. grid%lake_level == grid%lake_level_old) .and. (self%nval_surface(i) > 0)) then

            self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level
               ! Recalculate Q_start from deep inflows
               !call Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), Q_read_start(i, :), self%nval_deep(i))
               !call grid%interpolate_to_face_from_second(self%z_Inp(i, :), Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))
            call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
            call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))               
            !write(6,*) self%Qs_start(i,:)
            !write(6,*) self%Qs_read_start(i,1:2)
            !write(6,*) self%Qs_start(i,grid%nz_grid)
            !write(6,*) self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i))
               ! Recalculate Q_end from deep inflows
               !call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), Q_read_end(i, :), self%nval_deep(i))
               !call grid%interpolate_to_face_from_second(self%z_Inp(i, :), Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))

               ! Add surface flow
!                if (self%nval_surface(i) > 0) then
!                   do j = 1, self%nval_surface(i)
!                      call self%surface_flow(self%Inp_read_start(i, self%nval_deep(i) - 1 + j), self%Q_start(i, :), self%depth_surfaceFlow(i, j), i)
!                      call self%surface_flow(self%Inp_read_end(i, self%nval_deep(i) - 1 + j), self%Q_end(i, :), self%depth_surfaceFlow(i, j), i)
!                   end do
!                end if
            end if ! end if not lake_level...

            ! Temporal treatment of inflow
            if ((datum <= self%tb_start(i)) .or. (self%eof(i) == 1)) then ! if datum before first date or end of file reached
               goto 8
            else
               do while (.not. ((datum >= self%tb_start(i)) .and. (datum <= self%tb_end(i)))) ! Do until datum between dates
                  self%tb_start(i) = self%tb_end(i) ! Move one step in time
                  self%Q_start(i, :) = self%Q_end(i, :)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)

                  read (fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))

                  call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))

                  if (self%nval_surface(i)>0) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
                  !write(6,*) self%Q_start(i,:)
                  !write(6,*) self%Q_end(i,:)
                  ! Treat deep inflow
                  !call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                  !call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))

                  ! Add surface flow
!                   if (self%nval_surface(i) > 0) then
!                      do j = 1, self%nval_surface(i)
!                         call self%surface_flow(self%Inp_read_end(i, self%nval_deep(i) - 1 + j), self%Q_end(i, :), self%depth_surfaceFlow(i, j), i)
!                      end do
!                   end if

               end do ! end do while
            end if

            ! Linearly interpolate value at correct datum
            ! Note that Q_inp is on the face grid
            !self%Q_start(i,1)=0
            !self%Q_end(i,1)=0
            do j = 1, ubnd_fce
               Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i)) &
               * (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))/(self%tb_end(i)-self%tb_start(i))
            end do
            goto 11

 7          self%eof(i) = 1
 8          state%Q_inp(i,1:ubnd_fce) = self%Q_start(i,1:ubnd_fce) + self%Qs_start(i,1:ubnd_fce)     ! Set to closest available value
            goto 11

 9          write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            state%Q_inp(i, 1:ubnd_fce) = 0.0_RK
            self%Q_start(i, 1:ubnd_fce) = 0.0_RK
            self%Qs_start(i, 1:ubnd_fce) = 0.0_RK
            11        continue

         end do ! end do i=1,4
         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1)=0
         Q_vert(2:ubnd_fce) = Q_inp(1, 2:ubnd_fce) + Q_inp(2, 2:ubnd_fce)

         ! Set all Q to the differences (from the integrals)
         ! The new Q_inp is located on the volume grid, element 1 remains unchanged since element 1 is 0
         do i = 1, 4
            Q_inp(i, 1) = 0
            do j = 1, ubnd_vol
               Q_inp(i, j + 1) = Q_inp(i, j + 2) - Q_inp(i, j + 1)
            end do
               Q_inp(i,ubnd_vol + 1) = 0
               !write(6,*) Q_inp(i,ubnd_vol-5:ubnd_vol)
         end do

      end associate
   end subroutine

end module
