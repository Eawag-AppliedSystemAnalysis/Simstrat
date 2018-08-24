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


   ! Implementation for lateral rho
   subroutine lateral_rho_update(self, state)
      implicit none
      class(LateralRhoModule) :: self
      class(ModelState) :: state
   
      ! Local Declarations
      real(RK) :: Q_read_start(1:4,1:state%nz_input), Q_read_end(1:4,1:state%nz_input) ! Integrated outflow
      real(RK) :: Inp(1:4,1:state%nz_input)
      real(RK) :: dummy
      real(RK) :: Q_in(1:self%grid%ubnd_vol), h_in(1:self%grid%ubnd_vol)
      real(RK) :: T_in, S_in, rho_in, CD_in, g_red, slope, Ri, E, Q_inp_inc
     
      integer :: i, j, k, i1, i2
      integer :: fnum(1:4) ! File number
      character(len=20) :: fname(1:4)
     
      associate (datum=>state%datum, &
                 idx=>state%model_step_counter, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input
                 grid=>self%grid)
        
        fname = ['inflow           ','outflow          ','input temperature','input salinity   ']
        fnum = [41,42,43,44]
        do i=1,4
          if (idx==1) then
            if (i==1) then
               ! Allocate arrays if not already done, this saves memory compared to declaring with max_length_input_data
               allocate (self%z_Inp(1:4,1:state%nz_input)) ! Input depths
               allocate (self%Inp_read_start(1:4,1:state%nz_input))
               allocate (self%Inp_read_end(1:4,1:state%nz_input))
               allocate (self%Q_start(1:4,1:grid%nz_grid+1)) ! Input interpolated on grid
               allocate (self%Q_end(1:4,1:grid%nz_grid+1)) ! Input interpolated on grid
            end if
            
            if(self%cfg%disp_diagnostic==2) write(6,*) 'Starting to read '//trim(fname(i))//' file...'
            self%eof(i) = 0
            !Skip first three rows (except for Qout)
            read(fnum(i),*,end=9)
            read(fnum(i),*,end=9) self%nval(i)
            read(fnum(i),*,end=9) dummy, (self%z_Inp(i,j),j=1,self%nval(i))
            self%z_Inp(i,1:self%nval(i)) = grid%z_zero + self%z_Inp(i,1:self%nval(i))
            if (i==2 .and. self%nval(i)>1) then
               if (any(self%z_Inp(i,2:self%nval(i))<=self%z_Inp(i,1:(self%nval(i)-1)))) then
                  write(6,*) 'Error: outflow depths in ',trim(fname(i)),' file must be strictly increasing.'
                  stop
               end if
            end if
            if (i==3 .or. i==4) then
               if (any(self%z_Inp(i,1:self%nval(i))/=self%z_Inp(1,1:self%nval(i)))) then
                  write(6,*) 'Error: inflow depths in ',trim(fname(i)),' file must match the ones in inflow file.'
                  stop
               end if
            end if
            !Read first values
            read(fnum(i),*,end=9) self%tb_start(i),(self%Inp_read_start(i,j),j=1,self%nval(i))
            if(i==2) call Integrate(self%z_Inp(i,:),self%Inp_read_start(i,:),Q_read_start(i,:),self%nval(i))
            if(i==2) call grid%interpolate_to_face_from_second(self%z_Inp(i,:),Q_read_start(i,:),self%nval(i),self%Q_start(i,:))
            read(fnum(i),*,end=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))
            if(i==2) call Integrate(self%z_Inp(i,:),self%Inp_read_end(i,:),Q_read_end(i,:),self%nval(i))
            if(i==2) call grid%interpolate_to_face_from_second(self%z_Inp(i,:),Q_read_end(i,:),self%nval(i),self%Q_end(i,:))
          end if

          if ((datum<=self%tb_start(i)).or.(self%eof(i)==1)) then    ! if datum before first date or end of file reached
             goto 8
          else
             do while (.not.((datum>=self%tb_start(i)).and.(datum<=self%tb_end(i)))) ! do until datum between dates
                self%tb_start(i) = self%tb_end(i)             ! move one step in time
                if(i/=2) self%Inp_read_start(i,:) = self%Inp_read_end(i,:)
                if(i==2) self%Q_start(i,:) = self%Q_end(i,:)
                read(fnum(i),*,end=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))
                if(i==2) call Integrate(self%z_Inp(i,:),self%Inp_read_end(i,:),Q_read_end(i,:),self%nval(i))
                if(i==2) call grid%interpolate_to_face_from_second(self%z_Inp(i,:),Q_read_end(i,:),self%nval(i),self%Q_end(i,:))
             end do
             !Linearly interpolate value at correct datum
             if (i/=2) then
               Inp(i,1:self%nval(i)) = self%Inp_read_start(i,1:self%nval(i)) +&
                      (datum-self%tb_start(i)) * (self%Inp_read_end(i,1:self%nval(i))-self%Inp_read_start(i,1:self%nval(i)))/(self%tb_end(i)-self%tb_start(i))
             else
                if(self%tb_end(i)<=self%tb_start(i)) then
                  write(6,*) 'Error: dates in ',trim(fname(i)),' file must always be increasing.'
                  stop
                end if
                do j=1,grid%ubnd_vol
                   Q_inp(i,j) = self%Q_start(i,j) + (datum-self%tb_start(i)) * (self%Q_end(i,j)-self%Q_start(i,j))/(self%tb_end(i)-self%tb_start(i))
                end do
             end if
          end if
          goto 11

7         self%eof(i) = 1
8         if(i/=2) Inp(i,:) = self%Inp_read_start(i,:) ! Set to closest available value
          if(i==2) Q_inp(i,1:grid%ubnd_vol) = self%Q_start(i,1:grid%ubnd_vol) ! Set to closest available value
          goto 11

9         write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
          self%eof(i) = 1
          if(i/=2) Inp(i,1:self%nval(i)) = 0.0_RK
          if(i/=2) self%Inp_read_start(i,1) = 0.0_RK
          if(i==2) Q_inp(i,1:grid%ubnd_vol) = 0.0_RK
          if(i==2) self%Q_start(i,1:grid%ubnd_fce) = 0.0_RK

11        continue
        end do      ! end do i=1,4

        !Set Qout to the differences (from the integrals)
        do j=1,grid%ubnd_vol-1
           Q_inp(2,grid%ubnd_vol-j+1) = Q_inp(2,grid%ubnd_vol-j+1)-Q_inp(2,grid%ubnd_vol-j)
        end do
        Q_inp(1,:) = 0
        Q_inp(3:4,:) = 0

        do j=1,self%nval(1)
           if (Inp(1,j)>1E-15) then
              k=grid%ubnd_vol
              do while (grid%z_volume(k)>self%z_Inp(1,j))
               k=k-1
              end do
              Q_in(k) = Inp(1,j) !Inflow flow rate [m3/s]
              T_in = Inp(3,j) !Inflow temperature [°C*m3/s]
              S_in = Inp(4,j) !Inflow salinity [‰*m3/s]
              rho_in = rho_0*(0.9998395+T_in*(6.7914e-5+T_in*(-9.0894e-6+T_in*&
                     (1.0171e-7+T_in*(-1.2846e-9+T_in*(1.1592e-11+T_in*(-5.0125e-14))))))+&
                     (8.181e-4+T_in*(-3.85e-6+T_in*(4.96e-8)))*S_in) !Inflow density [kg/m3]
              g_red = g*(rho_in-state%rho(k))/rho_in !Reduced gravity [m/s2]

              slope = pi/72 !Slope of inflow
              !hang = pi/3 !Stream half-angle
              CD_in = self%param%CD*10 !Inflow drag coefficient
              !Ri = CD_in*(1+0.21*CD_in**0.5*sin(hang))/(sin(hang)*tan(slope)) !Richardson number
              !Ri = CD_in/tan(slope)*(1/sin(hang)+0.21*CD_in**0.5) !Richardson number
              Ri = CD_in/tan(slope)*(1.15+0.21*CD_in**0.5) !Richardson number (assuming an inflow half-angle of pi/3)
              E = 1.6*CD_in**1.5/Ri !Entrainment coefficient
              h_in(k) = (2*Q_in(k)**2*Ri*tan(slope)**2/abs(g_red))**0.2 !Inflow thickness [m]

              if (g_red>0) then !Inflow plunges
                 do while ((rho_in>state%rho(k)).and.(k>1))
                    h_in(k-1) = 1.2*E*(grid%z_volume(k)-grid%z_volume(k-1))/sin(slope) + h_in(k)
                    Q_in(k-1) = Q_in(k)*(h_in(k-1)/h_in(k))**(5./3.)
                    Q_inp(2,k) = Q_inp(2,k) - (Q_in(k-1)-Q_in(k))
                    T_in = (T_in*Q_in(k)+state%T(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                    S_in = (S_in*Q_in(k)+state%S(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                    rho_in = (rho_in*Q_in(k)+state%rho(k)*(Q_in(k-1)-Q_in(k)))/Q_in(k-1)
                    k=k-1
                 end do
                 i2 = k
                 do i1=k,grid%ubnd_vol !extend upwards
                    if((i1==grid%ubnd_vol).or.(grid%z_volume(i1+1)>(grid%z_volume(k)+h_in(k)))) exit
                 end do
              else if (g_red<0) then !Inflow rises
                 do while ((rho_in<state%rho(k)).and.(k<grid%ubnd_vol))
                    h_in(k+1) = 1.2*E*(grid%z_volume(k+1)-grid%z_volume(k))/sin(slope) + h_in(k)
                    Q_in(k+1) = Q_in(k)*(h_in(k+1)/h_in(k))**(5./3.)
                    Q_inp(2,k) = Q_inp(2,k) - (Q_in(k+1)-Q_in(k))
                    T_in = (T_in*Q_in(k)+state%T(k)*(Q_in(k+1)-Q_in(k)))/Q_in(k+1)
                    S_in = (S_in*Q_in(k)+state%S(k)*(Q_in(k+1)-Q_in(k)))/Q_in(k+1)
                    rho_in = (rho_in*Q_in(k)+state%rho(k)*(Q_in(k+1)-Q_in(k)))/Q_in(k+1)
                    k=k+1
                 end do
                 i1 = k
                 do i2=k,1,-1 !extend downwards
                    if((i2==1).or.(grid%z_volume(i2-1)<(grid%z_volume(k)-h_in(k)))) exit
                 end do
              end if

              do i=i2,i1
                 Q_inp_inc = Q_in(k)/(grid%z_face(i1)-grid%z_face(i2)+grid%h(i))*grid%h(i)
                 Q_inp(1,i) = Q_inp(1,i) + Q_inp_inc
                 Q_inp(3,i) = Q_inp(3,i) + T_in*Q_inp_inc
                 Q_inp(4,i) = Q_inp(4,i) + S_in*Q_inp_inc
              end do
           end if
        end do

        Q_vert(1) = Q_inp(1,1) + Q_inp(2,1)
        do i=2,grid%ubnd_fce
           Q_vert(i) = Q_vert(i-1) + Q_inp(1,i) + Q_inp(2,i)
        end do
      end associate
   end subroutine

   ! "Normal" Implementation
   ! Note that this method has been taken from old simstrat and needs refactoring
   subroutine lateral_update(self, state)
      implicit none
      class(LateralModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: dummy

      integer :: i, j
      integer :: fnum(1:4) ! File number
      character(len=20) :: fname(1:4)

      associate (datum=>state%datum, &
                 idx=>state%model_step_counter, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input
                 grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_Fce=>self%grid%ubnd_fce)

         fname = ['inflow           ', 'outflow          ', 'input temperature', 'input salinity   ']
         fnum = [41, 42, 43, 44]

         ! FB 2016: Major revision to include surface inflow
         do i = 1, 4 ! Do this for inflow, outflow, temperature and salinity
            if (idx==1) then ! First iteration
               if (i==1) then

                  ! Allocate arrays for very first iteration
                  allocate (self%z_Inp(1:4, 1:state%nz_input)) ! Input depths
                  allocate (self%Inp_read_start(1:4, 1:state%nz_input)) ! Raw input read
                  allocate (self%Inp_read_end(1:4, 1:state%nz_input)) ! Raw input read
                  allocate (self%Q_read_start(1:4, 1:state%nz_input)) ! Integrated input
                  allocate (self%Q_read_end(1:4, 1:state%nz_input)) ! Integrated input           
                  allocate (self%Qs_read_start(1:4, 1:state%nz_input))  ! Integrated surface input
                  allocate (self%Qs_read_end(1:4, 1:state%nz_input))  ! Integrated surface input
                  allocate (self%Q_start(1:4, 1:grid%nz_grid+1)) ! Input interpolated on grid
                  allocate (self%Q_end(1:4, 1:grid%nz_grid+1)) ! Input interpolated on grid
                  allocate (self%Qs_start(1:4, 1:grid%nz_grid+1)) ! Surface input interpolated on grid
                  allocate (self%Qs_end(1:4, 1:grid%nz_grid+1)) ! Surface input interpolated on grid
               end if

               ! Default value of surface inflow
               self%Qs_start(i, :) = 0.0_RK
               self%Qs_end(i, :) = 0.0_RK

               ! Open file and start to read
               self%eof(i) = 0
               read (fnum(i), *, end=9) ! Skip first row: description of columns

               ! Read number of deep and surface inflows
               read (fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)

               ! Total number of values to read
               self%nval(i) = self%nval_deep(i) + self%nval_surface(i)

               ! Read input depths
               read (fnum(i), *, end=9) dummy, (self%z_Inp(i, j), j=1, self%nval(i))

               ! Convert input depths
               self%z_Inp(i, 1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i, 1:self%nval_deep(i))
               ! Convert surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))

               ! Read first input line
               read (fnum(i), *, end=9) self%tb_start(i), (self%Inp_read_start(i, j), j=1, self%nval(i))

               ! Cumulative integration of input
               call Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), self%Q_read_start(i, :), self%nval_deep(i))
               ! Interpolation on face grid
               call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))

               ! If there is surface input, integrate and interpolate
               if (self%nval_surface(i)>0) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               end if


               ! Read second line and treatment of deep inflow
               read (fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
               call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
               call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))

               ! If there is surface input, integrate and interpolate
               if (self%nval_surface(i)>0) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if

               write(6,*) "--Input file successfully read: ",fname(i)
            end if ! idx==1



            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep them at the surface
            if ((.not. grid%lake_level == grid%lake_level_old) .and. (self%nval_surface(i) > 0)) then

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

                  read (fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))

                  call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))

                  if (self%nval_surface(i)>0) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
               end do ! end do while
            end if

            ! Linearly interpolate value at correct datum (Q_inp is on face grid)
            do j = 1, ubnd_fce
               Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i)) &
               * (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))/(self%tb_end(i)-self%tb_start(i))
            end do
            goto 11

            ! If end of file reached, set to closest available value
 7          self%eof(i) = 1
 8          Q_inp(i,1:ubnd_fce) = self%Q_start(i,1:ubnd_fce) + self%Qs_start(i,1:ubnd_fce)
            goto 11

            ! If no data available
 9          write(6,*) 'No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            Q_inp(i, 1:ubnd_fce) = 0.0_RK
            self%Q_start(i, 1:ubnd_fce) = 0.0_RK
            self%Qs_start(i, 1:ubnd_fce) = 0.0_RK
            11        continue

         end do ! end do i=1,4
         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1)=0
         Q_vert(2:ubnd_fce) = Q_inp(1, 2:ubnd_fce) + Q_inp(2, 2:ubnd_fce)

         ! The final Q_inp is located on the volume grid
         do i = 1, 4
            do j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            end do
               Q_inp(i,ubnd_vol + 1) = 0
         end do

      end associate
   end subroutine

end module
