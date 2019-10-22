!<    +---------------------------------------------------------------+
!     |  Advection module
!     |  - Based on already read in/outflows, calculates Advection
!     |  - Might grow or shrink grid (methods merge/add_box)
!<    +---------------------------------------------------------------+

module strat_advection
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use simstrat_aed2
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

   subroutine advection_init(self, state, model_config, model_param, grid)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param

      ! Local variables
      integer :: i

      self%cfg => model_config
      self%param => model_param
      self%grid => grid

      if (self%cfg%couple_aed2) then
         do i = 1, state%n_AED2
            select case(trim(state%AED2_names(i)))
            case('CAR_pH')
               state%n_pH = i
            end select
         end do
      end if

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
                 AED2_state=>state%AED2_state, &
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
         else if (top_z == grid%max_depth) then ! If we are already at the maximum lake level
            dt_i(1) = dt
         else if ((dh + top_z) >= grid%max_depth) then ! If the full timestep would lead to a lake level higher than maximum allowed lake level, split the timestep.
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

         ! FB 2016/2019: Revisions
         do t_i = 1, 2 !First and (if needed) second timestep
            AreaFactor_adv(1:nz_occupied) = dt_i(t_i)/((grid%Az(1:nz_occupied) + grid%Az(2:nz_occupied+1))/2*grid%h(1:nz_occupied)) ! Area factor for dt(t_i)
            dh_i(t_i) = dh*dt_i(t_i)/dt ! Depth difference for dt(t_i)

            ! Update Simstrat variables U, V, T and S
            call do_update_statvars(self, state, AreaFactor_adv(1:nz_occupied), dh_i(t_i))

            ! Update AED2 variables
            if(self%cfg%couple_aed2) call do_update_statvars_AED2(self, state, AreaFactor_adv(1:nz_occupied), dh_i(t_i))

            ! Adjust boxes (Horrible if/else construction - replace!)
            if (t_i == 1) then
               if (dh == 0) then ! If volume does not change, return
                  return
               else if ((dh + top_z) >= grid%max_depth) then ! If surface level reached
                  call grid%modify_top_box(grid%max_depth - top_z)
                  return
               else if (((dh_i(t_i) + top_h) > h_div_2) .and. &
                        ((dh_i(t_i) + top_h) < (h_mult_2))) then ! and top box<2*lower box
                  call grid%modify_top_box(dh_i(t_i))
                  return
               else if (t_i == 1 .and. (dh + top_h) <= h_div_2) then ! If top box<=0.5*lower box, merge 2 boxes
                  call self%merge_box(state, dh_i(t_i))
               else if (t_i == 1 .and. (dh + top_h) >= h_mult_2) then ! If top box>=2*lower box, add one box
                  call self%add_box(state, dh_i(t_i))
               end if ! dh==0
            end if

         end do !end do t_i=1,2
      end associate
   end subroutine

   subroutine do_update_statvars(self, state, AreaFactor_adv, dh)
      ! Arguments
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK), dimension(:) :: AreaFactor_adv
      real(RK) :: dh

      ! Local variables
      integer :: i, top
      real(RK) :: dU(self%grid%nz_grid), dV(self%grid%nz_grid), dTemp(self%grid%nz_grid), dS(self%grid%nz_grid)
      integer :: outflow_above, outflow_below

      associate(ubnd_vol => self%grid%ubnd_vol, &
         Q_vert => state%Q_vert, &
         h => self%grid%h)

            ! Calculate changes
            do i = 1, ubnd_vol
               ! For the top-most cell, if Q_vert at the upper face is positive, there is still no outflow (the cell is simply growing, but this is done elsewhere)
               if ((i == ubnd_vol) .and. Q_vert(i + 1) > 0) then
                  top = 0
               else
                  top = 1
               end if

               ! If Q_vert at the upper face of cell i is positive, then there is outflow to the cell above
               if (Q_vert(i + 1) > 0) then
                     outflow_above = 1
               else 
                     outflow_above = 0
               end if

               ! If Q_vert at the lower face of cell i is negative, then there is outflow to the cell below
               if (Q_vert(i) < 0) then
                     outflow_below = 1
               else
                     outflow_below = 0
               end if

               ! Calculate advective flow out of cell (thus negative sign in the front) i to the cells above and below
               dU(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%U(i)
               dV(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%V(i)
               dTemp(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%T(i)
               dS(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%S(i)

               ! Calculate the advective flow into cell i from below
               if (i > 1 .and. Q_vert(i ) > 0) then
                  dU(i) = dU(i) + Q_vert(i)*state%U(i - 1)
                  dV(i) = dV(i) + Q_vert(i)*state%V(i - 1)
                  dTemp(i) = dTemp(i) + Q_vert(i)*state%T(i - 1)
                  dS(i) = dS(i) + Q_vert(i)*state%S(i - 1)
               end if

               ! Calculate the advective flow into cell i from above (- sign in front because Q_vert is negative if there is inflow)
               if (i < ubnd_vol .and. Q_vert(i + 1) < 0) then
                  dU(i) = dU(i) - Q_vert(i + 1)*state%U(i + 1)
                  dV(i) = dV(i) - Q_vert(i + 1)*state%V(i + 1)
                  dTemp(i) = dTemp(i) - Q_vert(i + 1)*state%T(i + 1)
                  dS(i) = dS(i) - Q_vert(i + 1)*state%S(i + 1)
               end if
            end do

            ! Add change to state variables
            ! dT = dT(vertical advection) + dT(inflow) + dT(outflow), units: °C*m^3/s
            dTemp(1:ubnd_vol) = dTemp(1:ubnd_vol) + state%Q_inp(3, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%T(1:ubnd_vol)
            ! dS = dS(vertical advection) + dS(inflow) + dS(outflow), units: ‰*m^3/s
            dS(1:ubnd_vol) = dS(1:ubnd_vol) + state%Q_inp(4, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%S(1:ubnd_vol)

            ! Add change to the state variable
            state%U(1:ubnd_vol) = state%U(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dU(1:ubnd_vol)
            state%V(1:ubnd_vol) = state%V(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dV(1:ubnd_vol)
            state%T(1:ubnd_vol) = state%T(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dTemp(1:ubnd_vol)
            state%S(1:ubnd_vol) = state%S(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dS(1:ubnd_vol)

            ! Variation of variables due to change in volume
            state%U(ubnd_vol) = state%U(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%V(ubnd_vol) = state%V(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%T(ubnd_vol) = state%T(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%S(ubnd_vol) = state%S(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
      end associate
   end subroutine

   subroutine do_update_statvars_AED2(self, state, AreaFactor_adv, dh)
      ! Arguments
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK), dimension(:) :: AreaFactor_adv(1:self%grid%ubnd_vol)
      real(RK) :: dh

      ! Local variables
      integer :: i, top, outflow_above, outflow_below
      real(RK) :: dAED2(self%grid%nz_grid, state%n_AED2)

      associate(ubnd_vol => self%grid%ubnd_vol, &
         Q_vert => state%Q_vert, &
         AED2_state => state%AED2_state, &
         h => self%grid%h)

         ! Calculate changes
         do i = 1, ubnd_vol
            ! For the top-most cell, if Q_vert at the upper face is positive, there is still no outflow (the cell is simply growing, but this is done elsewhere)
            if ((i == ubnd_vol) .and. Q_vert(i + 1) > 0) then
               top = 0
            else
               top = 1
            end if

            ! If Q_vert at the upper face of cell i is positive, then there is outflow to the cell above
            if (Q_vert(i + 1) > 0) then
                  outflow_above = 1
            else 
                  outflow_above = 0
            end if

            ! If Q_vert at the lower face of cell i is negative, then there is outflow to the cell below
            if (Q_vert(i) < 0) then
                  outflow_below = 1
            else
                  outflow_below = 0
            end if

            ! Calculate advective flow out of cell (thus negative sign in the front) i to the cells above and below
            dAED2(i,:) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*AED2_state(i,:)

            ! Calculate the advective flow into cell i from below
            if (i > 1 .and. Q_vert(i) > 0) then
               dAED2(i,:) = dAED2(i,:) + Q_vert(i)*AED2_state(i - 1,:)
            end if
            ! Calculate the advective flow into cell i from above (- sign in front because Q_vert is negative if there is inflow)
            if (i < ubnd_vol .and. Q_vert(i + 1) < 0) then
               dAED2(i,:) = dAED2(i,:) - Q_vert(i + 1)*AED2_state(i + 1,:)
            end if
         end do

         ! Add change to state variables
         ! dAED2 = dAED2(vertical advection) + dAED2(inflow) + Outflow(negative)*AED2, units: C*m^3/s

         ! Add change to the state variable
         do i=1,state%n_AED2
            dAED2(1:ubnd_vol,i) = dAED2(1:ubnd_vol,i) + state%Q_inp(n_simstrat + i, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*AED2_state(1:ubnd_vol,i)
            AED2_state(1:ubnd_vol,i) = AED2_state(1:ubnd_vol,i) + AreaFactor_adv(1:ubnd_vol)*dAED2(1:ubnd_vol,i)
         end do

         ! Variation of variables due to change in volume
         AED2_state(ubnd_vol,:) = AED2_state(ubnd_vol,:)*h(ubnd_vol)/(h(ubnd_vol) + dh)

         ! Transform [H] back to pH
         if(self%cfg%couple_aed2) then
            AED2_state(:,state%n_pH) = -log10(AED2_state(:,state%n_pH))
         end if

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
      associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol, AED2_state=>state%AED2_state)

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

         ! AED2
         if (self%cfg%couple_AED2) AED2_state(ubnd_vol,:) = (w_a*AED2_state(ubnd_vol + 1,:) + w_b*AED2_state(ubnd_vol,:))/(w_a + w_b)

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

         if (self%cfg%couple_AED2) state%AED2_state(ubnd_vol,:) = state%AED2_state(ubnd_vol - 1,:)

         call self%grid%update_area_factors()

      end associate
   end subroutine

end module
