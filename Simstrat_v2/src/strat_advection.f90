!     +---------------------------------------------------------------+
!     |  Advection module
!     |  - Based on already read in/outflows, calculates Advection
!     |  - Might grow or shrink grid (methods merge/add_box)
!     +---------------------------------------------------------------+

module strat_advection
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use utilities
   implicit none
   private

   type, public :: AdvectionModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param

   contains
      procedure, pass :: init => advection_init
      procedure, pass :: update => advection_update
      procedure, pass :: merge_box => advection_merge_box
      procedure, pass :: add_box => advection_add_box
   end type

contains

   subroutine advection_init(self, model_config, model_param, grid)
      implicit none
      class(AdvectionModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param

      self%cfg => model_config
      self%param => model_param
      self%grid => grid

   end subroutine

   ! A lot of code that is hard to test - might be refactored in the future
   subroutine advection_update(self, state)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state

      real(RK) :: top_z, top_h
      real(RK) :: top
      real(RK) :: dh, dh_i(1:2), h_div_2, h_mult_2 ! depth differences
      real(RK) :: dU(self%grid%nz_grid), dV(self%grid%nz_grid), dTemp(self%grid%nz_grid), dS(self%grid%nz_grid)
      real(RK) :: dt_i(1:2) ! first and second time step
      real(RK) :: AreaFactor_adv(1:self%grid%nz_grid)
      integer :: i, t_i

      associate (grid=>self%grid, &
                 nz_occupied=>self%grid%nz_occupied, &
                 dt=>state%dt, &
                 h=>self%grid%h, &
                 Q_vert=>state%Q_vert, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)

         !Depth difference compared to previous timestep
         top_z = grid%z_face(ubnd_fce)
         top_h = grid%h(ubnd_vol)

         dh = state%Q_vert(ubnd_fce)/grid%Az(ubnd_fce)*state%dt
         h_div_2 = 0.5_RK*h(nz_occupied - 1) ! Take second highest box since the top box might not be at the full height
         h_mult_2 = 2_RK*h(nz_occupied - 1)

         ! Calculate timestep splitting
         !Split timestep depending on situation
         if (dh == 0.) then ! If volume does not change, take one normal time step
            dt_i(1) = dt
         else if ((dh + top_z) >= grid%max_depth) then ! If surface level reached, take a step until surface
            dt_i(1) = (grid%max_depth - top_z)/dh*dt
         else if (((dh + top_h) > h_div_2) .and. & ! If top box>0.5*lower box and <2*lower box, take one time step
                  ((dh + top_h) < h_mult_2)) then
            dt_i(1) = dt
         else if ((dh + top_h) <= h_div_2) then ! If top box<=0.5*lower box, first step until top box=0.5*lower box
            dt_i(1) = abs((top_h - h_div_2)/dh)*dt
         else ! If top box>=2*lower box, first step until top box = 2*lower box
            dt_i(1) = abs((2*h(nz_occupied - 1) - top_h)/dh)*dt
         end if
         dt_i(2) = dt - dt_i(1) ! Rest of timestep

         ! FB 2016: Revision
         do t_i = 1, 2 !First and (if needed) second timestep
            AreaFactor_adv(1:nz_occupied) = dt_i(t_i)/(grid%Az(2:nz_occupied+1)*grid%h(1:nz_occupied)) ! Area factor for dt(t_i)
            dh_i(t_i) = dh*dt_i(t_i)/dt ! Depth difference for dt(t_i)

            ! Calculate changes
            do i = 1, ubnd_vol
               if (i == ubnd_vol .and. Q_vert(i+1) > 0) then
                  top = 0
               else
                  top = 1
               end if
               ! Advective flow out of box i, always negative
               dU(i) = -top*abs(Q_vert(i+1))*state%U(i)
               dV(i) = -top*abs(Q_vert(i+1))*state%V(i)
               dTemp(i) = -top*abs(Q_vert(i+1))*state%T(i)
               dS(i) = -top*abs(Q_vert(i+1))*state%S(i)
               if (i > 2 .and. Q_vert(i) > 0) then ! Advective flow into box i, from above
                  dU(i) = dU(i) + Q_vert(i)*state%U(i - 1)
                  dV(i) = dV(i) + Q_vert(i)*state%V(i - 1)
                  dTemp(i) = dTemp(i) + Q_vert(i)*state%T(i - 1)
                  dS(i) = dS(i) + Q_vert(i)*state%S(i - 1)
               end if
               if (i < ubnd_vol) then
                  if (Q_vert(i + 2) < 0) then ! Advective flow into box i, from below
                     dU(i) = dU(i) - Q_vert(i + 2)*state%U(i + 1)
                     dV(i) = dV(i) - Q_vert(i + 2)*state%V(i + 1)
                     dTemp(i) = dTemp(i) - Q_vert(i + 2)*state%T(i + 1)
                     dS(i) = dS(i) - Q_vert(i + 2)*state%S(i + 1)
                  end if
               end if
            end do

            ! Add change to state variables
            ! dT = dT(vertical advection) + dT(inflow) + dT(negative outflow), units: °C*m^3/s
            dTemp(1:ubnd_vol) = dTemp(1:ubnd_vol) + state%Q_inp(3, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%T(1:ubnd_vol)
            ! dS = dS(vertical advection) + dS(inflow) + dT(negative outflow), units: ‰*m^3/s
            dS(1:ubnd_vol) = dS(1:ubnd_vol) + state%Q_inp(4, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%S(1:ubnd_vol)
            ! Add change to the state variable
            state%U(1:ubnd_vol) = state%U(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dU(1:ubnd_vol)
            state%V(1:ubnd_vol) = state%V(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dV(1:ubnd_vol)
            state%T(1:ubnd_vol) = state%T(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dTemp(1:ubnd_vol)
            state%S(1:ubnd_vol) = state%S(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dS(1:ubnd_vol)

            ! Variation of variables due to change in volume
            state%U(ubnd_vol) = state%U(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))
            state%V(ubnd_vol) = state%V(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))
            state%T(ubnd_vol) = state%T(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))
            state%S(ubnd_vol) = state%S(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))

            ! Adjust boxes (Horrible if/else construction - replace!)
            if (dh == 0) then ! If volume does not change, return
               return
            else if ((dh + top_z) >= grid%max_depth) then ! If surface level reached
               if (.not. (top_z == grid%max_depth)) then
                  call grid%modify_top_box(dh_i(t_i))
                  return
               else
                  return
               end if

            else if (((dh_i(t_i) + top_h) > h_div_2) .and. &
                     ((dh_i(t_i) + top_h) < (h_mult_2))) then ! and top box<2*lower box
               call grid%modify_top_box(dh_i(t_i))
               return
            else if (t_i == 1 .and. (dh + top_h) <= h_div_2) then ! If top box<=0.5*lower box, merge 2 boxes
               call self%merge_box(state, dh_i(t_i))
            else if (t_i == 1 .and. (dh + top_h) >= h_mult_2) then ! If top box>=2*lower box, add one box
               call self%add_box(state, dh_i(t_i))
            end if ! dh==0

         end do !end do t_i=1,2
      end associate
   end subroutine

   ! Merges two boxes
   ! - Takes care of calculating the new state variable for this box
   ! - Calls grid methods to modify grid spacing etc
   subroutine advection_merge_box(self, state, dh)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK) :: dh
      real(RK) :: w_a, w_b
      associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol)

         ! New values of the state variables are weighted averages
         !determine weighting an normalization connstant
         w_a = 0.5_RK*self%grid%Az(ubnd_fce)
         w_b = self%grid%Az(ubnd_fce - 1)

         ! shrink grid by one (this also updates ubnd_fce/vol)
         call self%grid%shrink(dh)

         ! update quantities in new top box (based on former top box and current value)
         state%U(ubnd_vol) = (w_a*state%U(ubnd_vol + 1) + w_b*state%U(ubnd_vol))/(w_a + w_b)
         state%V(ubnd_vol) = (w_a*state%V(ubnd_vol + 1) + w_b*state%V(ubnd_vol))/(w_a + w_b)
         state%T(ubnd_vol) = (w_a*state%T(ubnd_vol + 1) + w_b*state%T(ubnd_vol))/(w_a + w_b)
         state%S(ubnd_vol) = (w_a*state%S(ubnd_vol + 1) + w_b*state%S(ubnd_vol))/(w_a + w_b)

         state%k(ubnd_fce) = (w_a*state%k(ubnd_fce + 1) + w_b*state%k(ubnd_fce))/(w_a + w_b)
         state%eps(ubnd_fce) = (w_a*state%eps(ubnd_fce + 1) + w_b*state%eps(ubnd_fce))/(w_a + w_b)
         state%Q_vert(ubnd_fce) = (w_a*state%Q_vert(ubnd_fce + 1) + w_b*state%Q_vert(ubnd_fce))/(w_a + w_b)

         ! update area factors
         call self%grid%update_area_factors()

      end associate
   end subroutine

   ! Adds a new box
   ! - Takes care of calculating the new state variable for this box
   ! - Calls grid methods to modify grid spacing etc
   subroutine advection_add_box(self, state, dh)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK) :: dh
      associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol)

         ! extend grid by one (also updates ubnd_vol etc)
         call self%grid%grow(dh)
         
         ! Update quantities in new grid element
         state%U(ubnd_vol) = state%U(ubnd_vol - 1)
         state%V(ubnd_vol) = state%V(ubnd_vol - 1)
         state%T(ubnd_vol) = state%T(ubnd_vol - 1)
         state%S(ubnd_vol) = state%S(ubnd_vol - 1)
         state%Q_vert(ubnd_fce) = state%Q_vert(ubnd_fce - 1) ! Vertical discharge of new box

         state%k(ubnd_fce) = state%k(ubnd_fce - 1)
         state%eps(ubnd_fce) = state%eps(ubnd_fce - 1)

         call self%grid%update_area_factors()

      end associate
   end subroutine

end module
