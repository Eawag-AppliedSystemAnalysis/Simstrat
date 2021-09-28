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
!     |  Lateral module
!     |     Reads and processes inflows/outflows such that
!     |     Advection can be calculated in the next step
!     |     There are at least two different implementations possible:
!     |         - LateralRhoModule : Inflow plunges according to density
!     |         - LateralModule:     Inflow affects layer as configured in file
!<    +---------------------------------------------------------------+


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
      real(RK), dimension(:), allocatable  :: tb_start, tb_end ! Start time, end time
      integer, dimension(:), allocatable  :: eof, nval, nval_deep, nval_surface, fnum
      integer, dimension(:), allocatable :: number_of_lines_read
      logical, dimension(:), allocatable :: has_surface_input, has_deep_input
      integer :: n_vars, max_n_inflows
      logical :: couple_aed2
      character(len=100) :: simstrat_path(n_simstrat), aed2_path

   contains
      procedure, pass :: init => lateral_generic_init
      procedure, pass :: save => lateral_generic_save
      procedure, pass :: load => lateral_generic_load
      procedure(lateral_generic_update), deferred, pass :: update
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

   subroutine lateral_generic_init(self, state, model_config, input_config, aed2_config, model_param, grid)
      implicit none
      class(GenericLateralModule) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(InputConfig), target :: input_config
      class(AED2Config), target :: aed2_config
      class(ModelParam), target :: model_param

      ! Locals
      integer :: i

      self%cfg => model_config
      self%param => model_param
      self%grid => grid

      self%n_vars = n_simstrat
      self%simstrat_path(1) = input_config%QinpName
      self%simstrat_path(2) = input_config%QoutName
      self%simstrat_path(3) = input_config%TinpName
      self%simstrat_path(4) = input_config%SinpName

      self%couple_aed2 = model_config%couple_aed2
      if (self%couple_aed2) then
         self%n_vars = self%n_vars + state%n_AED2_state
         self%aed2_path = aed2_config%path_aed2_inflow
      end if

      self%max_n_inflows = model_config%max_length_input_data

      allocate(self%eof(self%n_vars))
      allocate(self%nval(self%n_vars))
      allocate(self%nval_deep(self%n_vars))
      allocate(self%nval_surface(self%n_vars))
      allocate(self%tb_start(self%n_vars))
      allocate(self%tb_end(self%n_vars))
      allocate(self%fnum(self%n_vars))
      allocate(self%number_of_lines_read(self%n_vars))
      self%number_of_lines_read = 0

      allocate (self%z_Inp(1:self%n_vars,1:self%max_n_inflows)) ! Input depths read from file
      allocate (self%Inp_read_start(1:self%n_vars,1:self%max_n_inflows))  ! Input read from file
      allocate (self%Inp_read_end(1:self%n_vars,1:self%max_n_inflows))  ! Input read from file
      allocate (self%Q_read_start(1:self%n_vars, 1:self%max_n_inflows)) ! Integrated input
      allocate (self%Q_read_end(1:self%n_vars, 1:self%max_n_inflows)) ! Integrated input
      allocate (self%Qs_read_start(1:self%n_vars,1:self%max_n_inflows))  ! Integrated surface input
      allocate (self%Qs_read_end(1:self%n_vars,1:self%max_n_inflows))  ! Integrated surface input
      allocate (self%Q_start(1:self%n_vars,1:grid%nz_grid + 1)) ! Input interpolated on grid
      allocate (self%Q_end(1:self%n_vars,1:grid%nz_grid + 1)) ! Input interpolated on grid
      allocate (self%Qs_start(1:self%n_vars,1:grid%nz_grid + 1)) ! Surface input interpolated on grid
      allocate (self%Qs_end(1:self%n_vars,1:grid%nz_grid + 1)) ! Surface input interpolated on grid

      allocate(state%Q_inp(1:self%n_vars,1:grid%nz_grid + 1))
      allocate(self%has_surface_input(1:self%n_vars))
      allocate(self%has_deep_input(1:self%n_vars))

      ! Get location of pH in AED2 array
      if (self%couple_aed2) then
         do i = 1, state%n_AED2_state
            select case(trim(state%AED2_state_names(i)))
            case('CAR_pH')
               state%n_pH = i
            end select
         end do
      end if
   end subroutine

   subroutine lateral_generic_save(self)
      implicit none
      class(GenericLateralModule) :: self
      logical :: has_allocated

      has_allocated = allocated(self%z_Inp)

      if (has_allocated) then
         write (80) has_allocated
         call save_integer_array(80, self%number_of_lines_read)
         call save_integer_array(80, self%eof)
         call save_integer_array(80, self%nval)
         call save_integer_array(80, self%nval_deep)
         call save_integer_array(80, self%nval_surface)
         call save_integer_array(80, self%fnum)
         call save_logical_array(80, self%has_surface_input)
         call save_logical_array(80, self%has_deep_input)
         call save_array(80, self%tb_start)
         call save_array(80, self%tb_end)
         call save_matrix(80, self%z_Inp)
         call save_matrix(80, self%Q_start)
         call save_matrix(80, self%Qs_start)
         call save_matrix(80, self%Q_end)
         call save_matrix(80, self%Qs_end)
         call save_matrix(80, self%Q_read_start)
         call save_matrix(80, self%Q_read_end)
         call save_matrix(80, self%Inp_read_start)
         call save_matrix(80, self%Inp_read_end)
         call save_matrix(80, self%Qs_read_start)
         call save_matrix(80, self%Qs_read_end)
      end if
   end subroutine

   subroutine lateral_generic_load(self)
      implicit none
      class(GenericLateralModule) :: self
      logical :: has_allocated

      has_allocated = allocated(self%z_Inp)

      if (has_allocated) then
         read (81) has_allocated
         call read_integer_array(81, self%number_of_lines_read)
         call read_integer_array(81, self%eof)
         call read_integer_array(81, self%nval)
         call read_integer_array(81, self%nval_deep)
         call read_integer_array(81, self%nval_surface)
         call read_integer_array(81, self%fnum)
         call read_logical_array(81, self%has_surface_input)
         call read_logical_array(81, self%has_deep_input)
         call read_array(81, self%tb_start)
         call read_array(81, self%tb_end)
         call read_matrix(81, self%z_Inp)
         call read_matrix(81, self%Q_start)
         call read_matrix(81, self%Qs_start)
         call read_matrix(81, self%Q_end)
         call read_matrix(81, self%Qs_end)
         call read_matrix(81, self%Q_read_start)
         call read_matrix(81, self%Q_read_end)
         call read_matrix(81, self%Inp_read_start)
         call read_matrix(81, self%Inp_read_end)
         call read_matrix(81, self%Qs_read_start)
         call read_matrix(81, self%Qs_read_end)
      end if
   end subroutine
      
      
   ! Implementation for lateral rho
   subroutine lateral_rho_update(self, state)
      implicit none
      class(LateralRhoModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: Inp(1:self%n_vars,1:self%max_n_inflows)
      real(RK) :: dummy
      real(RK) :: Q_in(1:self%grid%ubnd_vol), h_in(1:self%grid%ubnd_vol)
      real(RK) :: T_in, S_in, co2_in, ch4_in, rho_in, CD_in, g_red, slope, Ri, E, Q_inp_inc
      real(RK) :: AED2_in(state%n_AED2_state)
      integer :: i, j, k, i1, i2, l, status
      character(len=100) :: fname


      associate (datum=>state%datum, &
                 idx=>state%first_timestep, &
                 number_of_lines_read=>self%number_of_lines_read, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input at each depth (integrated inflow - outflow)
                 grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)

         do i=1, self%n_vars
            if (idx) then  ! If first timestep
               if (self%number_of_lines_read(i) == 0) then ! If start is not from snapshot
                  ! max_n_inflows was set to 1000 automatically. To reduce the size, it is redetermined here.
                  self%max_n_inflows = 0

                  ! Read inflow files
                  if (i > n_simstrat) then
                     fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
                  else
                     fname = trim(self%simstrat_path(i))
                  end if

                  self%fnum(i) = i + 60  ! Should find a better way to manage unit numbers
                  open(self%fnum(i), action='read', status='old', file=fname)
                  
                  if (status .ne. 0) then
                     call error('File '//fname//' not found.')
                  else
                     write(6,*) 'Reading ', fname
                  end if

                  ! Default values
                  self%Q_start(i,:) = 0.0_RK
                  self%Q_end(i,:) = 0.0_RK
                  self%Qs_start(i,:) = 0.0_RK
                  self%Qs_end(i,:) = 0.0_RK

                  ! End of file is not reached
                  self%eof(i) = 0

                  ! Read input depths
                  read(self%fnum(i),*,end=9)
                  call count_read(self, i)

                  ! Read number of deep and surface columns
                  read(self%fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)
                  call count_read(self, i)

                  ! Check whether there is deep inflow (fixed) or surface inflow (varies with lake level)
                  if (self%nval_deep(i) > 0) then
                     self%has_deep_input(i) = .true.
                  else
                     self%has_deep_input(i) = .false.
                  end if
                  if (self%nval_surface(i) > 0) then
                     self%has_surface_input(i) = .true.
                  else
                     self%has_surface_input(i) = .false.
                  end if

                  ! Total number of values to read
                  self%nval(i) = self%nval_deep(i) + self%nval_surface(i)

                  if (self%nval(i) > self%max_n_inflows) self%max_n_inflows = self%nval(i)

                  ! Read input depths
                  read(self%fnum(i),*,end=9) dummy, (self%z_Inp(i,j),j=1,self%nval(i))
                  call count_read(self, i)

                  ! Convert deep input depths
                  self%z_Inp(i,1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i,1:self%nval_deep(i))
                  ! Convert surface input depths
                  self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))

                  !Read first input values
                  read(self%fnum(i),*,end=9) self%tb_start(i),(self%Inp_read_start(i,j),j=1,self%nval(i))
                  call count_read(self, i)

                  ! If there is deep outflow (i==2)
                  if (i==2 .and. self%has_deep_input(i)) then
                     call Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_start(i,1:self%nval_deep(i)),self%Q_read_start(i,:),self%nval_deep(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_start(i,:),self%nval_deep(i),self%Q_start(i,:))
                  end if
                  ! If there is any surface inflow
                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  end if

                  ! Read next line
                  read(self%fnum(i),*,end=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))
                  call count_read(self, i)

                  ! If there is deep outflow (i==2)
                  if (i==2 .and. self%has_deep_input(i)) then
                     call Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_end(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i),self%Q_end(i,:))
                  end if

                  ! If there is any surface inflow
                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
                  call ok('Input file successfully read: '//fname)
               else ! if start from snapshot
                  ! Open inflow files
                  if (i > n_simstrat) then
                     fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
                  else
                     fname = trim(self%simstrat_path(i))
                  end if
                  open(self%fnum(i), action='read', status='old', file=fname)
                  
                  if (status .ne. 0) then
                     call error('File '//fname//' not found.')
                  else
                     write(6,*) 'Reading ', fname
                  end if
                  do l = 1, self%number_of_lines_read(i)
                     read (self%fnum(i), *, end=9) ! Skip already read and processed lines
                  end do
                  call ok('Input file successfully opened: '//fname)
               end if
            end if ! if        

            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep relative inflow depth constant
            if ((.not. grid%lake_level == grid%lake_level_old) .and. self%has_surface_input(i)) then

               ! Readjust surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

               ! Adjust surface inflow to new lake level
               if (self%has_surface_input(i)) then
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if
            end if ! end if not lake_level...


            if ((datum<=self%tb_start(i)).or.(self%eof(i)==1)) then    ! if datum before first date or end of file reached
               goto 8
            else
               do while (.not.((datum>=self%tb_start(i)).and.(datum<=self%tb_end(i)))) ! do until datum between dates
                  ! Move one step in time
                  self%tb_start(i) = self%tb_end(i)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)

                  ! For outflow, take Q_start; for inflow, temperature and salinity, take Inp for the plunging algorithm
                  if (i==2) then
                     self%Q_start(i, :) = self%Q_end(i, :)
                  else
                     self%Inp_read_start(i,1:self%nval_deep(i)) = self%Inp_read_end(i,1:self%nval_deep(i))
                  end if

                  ! Read next line
                  read(self%fnum(i),*,end=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))
                  call count_read(self, i)

                  ! If there is deep outflow (i==2)
                  if (i==2 .and. self%has_deep_input(i)) then
                     call Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_end(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i),self%Q_end(i,:))
                  end if

                  ! If there is any surface inflow
                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
               end do

               if(self%tb_end(i)<=self%tb_start(i)) then
                  call error('Dates in '//trim(fname)//' file must always be increasing.')
               end if

               ! Linearly interpolate value at correct datum
               if (i/=2) then
                  ! For plunging input, Inp will be needed later
                  Inp(i,1:self%nval_deep(i)) = self%Inp_read_start(i,1:self%nval_deep(i)) + (datum-self%tb_start(i))/(self%tb_end(i) - self%tb_start(i)) &
                  * (self%Inp_read_end(i,1:self%nval_deep(i)) - self%Inp_read_start(i,1:self%nval_deep(i)))


                  ! Surface input is already added to Q_inp; the plunging algorithm will add the deep input further below
                  do j=1,ubnd_fce
                     Q_inp(i,j) = (self%Qs_start(i,j)) + (datum - self%tb_start(i))/(self%tb_end(i) - self%tb_start(i))* &
                     (self%Qs_end(i,j) - self%Qs_start(i,j))
                  end do
               else
                  ! For outflow (i==2), both surface and deep inputs are added to Q_inp
                  do j=1,ubnd_fce
                     Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i))/(self%tb_end(i)-self%tb_start(i))* &
                     (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))
                  end do
               end if
            end if
            goto 11

7           self%eof(i) = 1
8           if(i/=2) Inp(i,1:self%nval_deep(i)) = self%Inp_read_start(i,1:self%nval_deep(i)) ! Set to closest available value
            Q_inp(i,1:ubnd_vol) = self%Q_start(i,1:ubnd_vol) + self%Qs_start(i,1:ubnd_vol) ! Set to closest available value
            goto 11

9           write(6,*) '[WARNING] ','No data found in ',trim(fname),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            if(i/=2) Inp(i,1:self%nval_deep(i)) = 0.0_RK
            if(i/=2) self%Inp_read_start(i,1) = 0.0_RK
            if(i==2) Q_inp(i,1:ubnd_vol) = 0.0_RK
            if(i==2) self%Q_start(i,1:ubnd_fce) = 0.0_RK
            self%Qs_start(i,1:ubnd_fce) = 0.0_RK

11          continue
         end do      ! end do i=1,n_vars

         !Set Q_inp to the differences (from the integrals)
         do i = 1,self%n_vars
            do j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            end do
            Q_inp(i,ubnd_vol + 1) = 0
         end do

         ! Only if biochemistry enabled: Transform pH to [H] for physical mixing processes
         if (self%couple_aed2 .and.state%n_pH > 0) then
            ! current pH profile
            state%AED2_state(:,state%n_pH) = 10.**(-state%AED2_state(:,state%n_pH))
            do i=1,ubnd_vol
               if (Q_inp(n_simstrat + state%n_pH,i) > 0) then
                  ! Surface inflows: pH is given as pH*m2/s, so before transforming to [H], we need to get rid of the m2/s temporarily
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)/Q_inp(1,i)
                  Q_inp(n_simstrat + state%n_pH,i) = 10.**(-Q_inp(n_simstrat + state%n_pH,i))
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)*Q_inp(1,i)
               end if
            end do
         end if

         ! Plunging algorithm
         do j = 1,self%nval_deep(1)  ! nval_deep needs to be the same for all i
            if (Inp(1,j) > 1E-15) then
               k = ubnd_vol
               do while (grid%z_volume(k) > self%z_Inp(1,j)) ! Find the place where the plunging inflow enters the lake (defined in file)
                  k = k - 1
               end do

               ! Get initial Q, T and S for the plunging inflow (before entrainment of ambient water)
               Q_in(k) = Inp(1,j) !Inflow flow rate [m3/s]
               T_in = Inp(3,j) !Inflow temperature [°C]
               S_in = Inp(4,j) !Inflow salinity [‰]

               ! Only if biochemistry enabled
               if (self%couple_aed2) then
                  ! Get AED2 values for the plunging inflow (before entrainment of ambient water)
                  AED2_in = Inp(n_simstrat + 1 : self%n_vars,j)
                  ! Transform pH to [H] for physical mixing processes
                  if (state%n_pH > 0) AED2_in(state%n_pH) = 10.**(-AED2_in(state%n_pH))
               end if
               ! Compute density as a function of T and S
               call calc_density(rho_in, T_in, S_in)
               g_red = g*(rho_in - state%rho(k))/rho_in !Reduced gravity [m/s2]

               slope = pi/72 !Slope of inflow
               !hang = pi/3 !Stream half-angle
               CD_in = self%param%CD*10 !Inflow drag coefficient
               !Ri = CD_in*(1+0.21*CD_in**0.5*sin(hang))/(sin(hang)*tan(slope)) !Richardson number
               !Ri = CD_in/tan(slope)*(1/sin(hang)+0.21*CD_in**0.5) !Richardson number
               Ri = CD_in/tan(slope)*(1.15 + 0.21*CD_in**0.5) !Richardson number (assuming an inflow half-angle of pi/3)
               E = 1.6*CD_in**1.5/Ri !Entrainment coefficient
               h_in(k) = (2*Q_in(k)**2*Ri*tan(slope)**2/abs(g_red))**0.2 !Inflow thickness [m]

               if (g_red > 0) then !Inflow plunges
                  do while ((rho_in > state%rho(k)).and.(k > 1))
                     h_in(k - 1) = 1.2*E*(grid%z_volume(k) - grid%z_volume(k-1))/sin(slope) + h_in(k)
                     Q_in(k - 1) = Q_in(k)*(h_in(k - 1)/h_in(k))**(5./3.)
                     Q_inp(2,k) = Q_inp(2,k) - (Q_in(k-1) - Q_in(k))
                     T_in = (T_in*Q_in(k) + state%T(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     S_in = (S_in*Q_in(k) + state%S(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     if (self%couple_aed2) then
                        AED2_in = (AED2_in*Q_in(k) + state%AED2_state(k,:)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     end if
                     rho_in = (rho_in*Q_in(k) + state%rho(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     k = k - 1
                  end do
                  i2 = k
                  do i1 = k,ubnd_vol !extend upwards
                     if(i1 == ubnd_vol) exit
                     if(grid%z_volume(i1 + 1) > (grid%z_volume(k) + h_in(k))) exit
                  end do
               else if (g_red < 0) then !Inflow rises
                  do while ((rho_in < state%rho(k)) .and. (k < ubnd_vol))
                     h_in(k + 1) = 1.2*E*(grid%z_volume(k + 1) - grid%z_volume(k))/sin(slope) + h_in(k)
                     Q_in(k + 1) = Q_in(k)*(h_in(k + 1)/h_in(k))**(5./3.)
                     Q_inp(2,k) = Q_inp(2,k) - (Q_in(k + 1) - Q_in(k))
                     T_in = (T_in*Q_in(k) + state%T(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     S_in = (S_in*Q_in(k) + state%S(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     if (self%couple_aed2) then
                        AED2_in = (AED2_in*Q_in(k) + state%AED2_state(k,:)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     end if
                     rho_in = (rho_in*Q_in(k) + state%rho(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     k = k + 1
                  end do
                  i1 = k
                  do i2 = k,1,-1 !extend downwards
                     if(i2 == 1) exit
                     if(grid%z_volume(i2 - 1) < (grid%z_volume(k) - h_in(k))) exit
                  end do
               end if

               ! Deep plunging input is added to Q_inp for i=1,3,4 (inflow, temperature, salinity)
               do i = i2,i1
                  Q_inp_inc = Q_in(k)/(grid%z_face(i1 + 1) - grid%z_face(i2))*grid%h(i)
                  Q_inp(1,i) = Q_inp(1,i) + Q_inp_inc
                  Q_inp(3,i) = Q_inp(3,i) + T_in*Q_inp_inc
                  Q_inp(4,i) = Q_inp(4,i) + S_in*Q_inp_inc
                  if (self%couple_aed2) Q_inp(n_simstrat + 1 : self%n_vars,i) = Q_inp(n_simstrat + 1 : self%n_vars,i) + AED2_in*Q_inp_inc
               end do
            end if
         end do

         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1) = 0
         do i = 2,ubnd_fce
            Q_vert(i) = Q_vert(i - 1) + Q_inp(1,i - 1) + Q_inp(2,i - 1)
         end do
      end associate
   end subroutine

   ! "Normal" Implementation
   subroutine lateral_update(self, state)
      implicit none
      class(LateralModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: dummy
      integer :: i, j, l, status
      character(len=100) :: fname


      associate (datum=>state%datum, &
                 idx=>state%first_timestep, &
                 number_of_lines_read=>self%number_of_lines_read, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input
                 grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_Fce=>self%grid%ubnd_fce)

         do i=1, self%n_vars
            if (idx) then  ! If first timestep
               if (self%number_of_lines_read(i) == 0) then  ! If not started from snapshot
                  if (i > n_simstrat) then
                     fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
                  else
                     fname = trim(self%simstrat_path(i))
                  end if

                  ! Read inflow files
                  self%fnum(i) = i + 60  ! Should find a better way to manage unit numbers
                  open(self%fnum(i), action='read', status='old', file=fname)
                  
                  if (status .ne. 0) then
                     call error('File '//fname//' not found.')
                  else
                     write(6,*) 'Reading ', fname
                  end if

                   ! Default values
                  self%Q_start(i,:) = 0.0_RK
                  self%Q_end(i,:) = 0.0_RK
                  self%Qs_start(i, :) = 0.0_RK
                  self%Qs_end(i, :) = 0.0_RK

                  ! Open file and start to read
                  self%eof(i) = 0
                  read (self%fnum(i), *, end=9) ! Skip first row: description of columns
                  call count_read(self, i)

                  ! Number of deep (fixed) and surface inputs to read
                  read (self%fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)
                  call count_read(self, i)

                  ! Check whether there is deep inflow (fixed) or surface inflow (varies with lake level)
                  if (self%nval_deep(i) > 0) then
                     self%has_deep_input(i) = .true.
                  else
                     self%has_deep_input(i) = .false.
                  end if
                  if (self%nval_surface(i) > 0) then
                     self%has_surface_input(i) = .true.
                  else
                     self%has_surface_input(i) = .false.
                  end if

                  ! Total number of values to read
                  self%nval(i) = self%nval_deep(i) + self%nval_surface(i)

                  ! Read input depths
                  read (self%fnum(i), *, end=9) dummy, (self%z_Inp(i, j), j=1, self%nval(i))
                  call count_read(self, i)

                  ! Convert input depths
                  self%z_Inp(i, 1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i, 1:self%nval_deep(i))

                  if (self%has_surface_input(i)) then
                     ! Convert surface input depths
                     self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))
                  end if

                  ! Read first input line
                  read (self%fnum(i), *, end=9) self%tb_start(i), (self%Inp_read_start(i, j), j=1, self%nval(i))
                  call count_read(self, i)

                  if (self%has_deep_input(i)) then
                     ! Cumulative integration of input
                     call Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), self%Q_read_start(i, :), self%nval_deep(i))
                     ! Interpolation on face grid
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))
                  end if

                  ! If there is surface input, integrate and interpolate
                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  end if

                  ! Read second line and treatment of deep inflow
                  read (self%fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
                  call count_read(self, i)

                  if (self%has_deep_input(i)) then
                     call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
                  end if
                  ! If there is surface input, integrate and interpolate
                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if

                  call ok('Input file successfully read: '//fname)
               else ! if start from snapshot
                  ! Open inflow files
                  if (i > n_simstrat) then
                     fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
                  else
                     fname = trim(self%simstrat_path(i))
                  end if
                  open(self%fnum(i), action='read', status='old', file=fname)

                  do l = 1, self%number_of_lines_read(i)
                     read (self%fnum(i), *, end=9) ! Skip already read and processed lines
                  end do
                  call ok('Input file successfully opened: '//fname)
               end if
            end if ! idx = 1

            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep them at the surface
            if ((.not. grid%lake_level == grid%lake_level_old) .and. (self%has_surface_input(i))) then

               ! Readjust surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

               ! Adjust surface inflow to new lake level
               call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))

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

                  read (self%fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
                  call count_read(self, i)

                  if (self%has_deep_input(i)) then
                    call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                    call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
                  end if

                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
               end do ! end do while
            end if

            ! Linearly interpolate value at correct datum (Q_inp is on face grid)
            do j = 1, ubnd_fce
               Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i))/(self%tb_end(i)-self%tb_start(i))* &
               (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))
            end do
            !write(6,*) i, self%z_Inp(i, self%nval(i)), self%Qs_read_end(i, self%nval(i)), self%nval_surface(i), self%Qs_end(i, ubnd_vol)
            goto 11

            ! If end of file reached, set to closest available value
 7          self%eof(i) = 1
 8          Q_inp(i,1:ubnd_fce) = self%Q_start(i,1:ubnd_fce) + self%Qs_start(i,1:ubnd_fce)
            goto 11

            ! If no data available
 9          write(6,*) '[WARNING] ','No data found in ',trim(fname),' file. Check number of depths. Values set to zero.'

            self%eof(i) = 1
            Q_inp(i, 1:ubnd_fce) = 0.0_RK
            self%Q_start(i, 1:ubnd_fce) = 0.0_RK
            self%Qs_start(i, 1:ubnd_fce) = 0.0_RK
11          continue

         end do ! end do i=1,self%n_vars
         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1)=0
         Q_vert(2:ubnd_fce) = Q_inp(1, 2:ubnd_fce) + Q_inp(2, 2:ubnd_fce)

         ! The final Q_inp is located on the volume grid
         do i = 1, self%n_vars
            do j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            end do
            Q_inp(i,ubnd_vol + 1) = 0
         end do

         ! Only if biochemistry enabled: Transform pH to [H] for physical mixing processes
         if (self%couple_aed2 .and. state%n_pH > 0) then
            ! current pH profile
            state%AED2_state(:,state%n_pH) = 10.**(-state%AED2_state(:,state%n_pH))
            do i=1,ubnd_vol
               if (Q_inp(n_simstrat + state%n_pH,i) > 0) then
                  ! Surface inflows: pH is given as pH*m2/s, so before transforming to [H], we need to get rid of the m2/s temporarily
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)/Q_inp(1,i)
                  Q_inp(n_simstrat + state%n_pH,i) = 10.**(-Q_inp(n_simstrat + state%n_pH,i))
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)*Q_inp(1,i)
               end if
            end do
         end if

      end associate
   end subroutine

   subroutine count_read(self, i)
      implicit none
      class(GenericLateralModule) :: self
      integer i
      
      self%number_of_lines_read(i) = self%number_of_lines_read(i) + 1
   end subroutine

end module
