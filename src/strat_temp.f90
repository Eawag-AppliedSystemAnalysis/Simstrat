!<    +---------------------------------------------------------------+
!     | Implementation of a statevar for a temperature variable
!<    +---------------------------------------------------------------+


module strat_temp
   use strat_kinds
   use strat_consts
   use strat_simdata
   use strat_statevar
   use strat_grid
   use strat_solver
   implicit none
   private

   type, extends(ModelVariable), public :: TempModelVar
   contains
      procedure, pass(self), public :: calc_terms => temp_var_calc_terms
      procedure, pass(self), public :: post_solve => temp_var_post_solve
   end type

contains

  ! Only calc_terms overwritten
   subroutine temp_var_calc_terms(self, state, param, sources, boundaries)
      class(TempModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(:) ::  sources, boundaries
      integer :: i
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         ! Radiation reaching top layer
         state%rad(self%grid%ubnd_fce) = state%rad0/rho_0/cp ![°C*m/s]  , rad is on faces

         ! Radiation reaching a layer is equal to radiation in the layer above minus absorption
         do i = ubnd_fce - 1, 1, -1
            state%rad(i) = state%rad(i + 1)*exp(-grid%h(i)*(state%absorb(ubnd_fce - i)+state%absorb(ubnd_fce + 1 - i))/2) !Attenuated by absorption
         end do

         !!!!!!!! Define sources !!!!!!!!
         ! Add Hsol Term to sources (Eq 1, Goudsmit(2002))
         sources(1:ubnd_vol) = (state%rad(2:ubnd_fce) - state%rad(1:ubnd_fce - 1))/grid%h(1:ubnd_vol)

         ! Set boundary heat flux at surface (Eq 25, Goudsmit(2002))
         sources(ubnd_vol) = sources(ubnd_vol) + state%heat/rho_0/cp/grid%h(ubnd_vol)

         ! No explicit boundary conditions
         boundaries(1:ubnd_vol) = 0

         ! Forcing mode 1 for temp is done in post solve method
         if (self%cfg%forcing_mode==1) then
         !   sources(ubnd_vol) = (state%SST-state%T(ubnd_vol))/state%dt
         !   boundaries(ubnd_vol) = 0
         !    bu(nz) = 1.0_dp
         !    au(nz) = 0.0_dp
         !    cu(nz) = 0.0_dp
         !    du(nz) = SST
         end if

         ! Add geothermal flux to sources (Eq 1, Goudsmit(2002))
         if (param%fgeo /= 0) sources(1:ubnd_vol) = sources(1:ubnd_vol) + state%fgeo_add(1:ubnd_vol)

      end associate
   end subroutine

   subroutine temp_var_post_solve(self, state)
      class(TempModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state

      real(RK) :: Az_vol(self%grid%nz_grid), rho(self%grid%length_vol), depth(size(self%grid%z_volume))
      real(RK) :: rho0t(self%grid%length_vol), rho0st(self%grid%length_vol)
      integer :: i
      real(RK) :: volume, sum_z, zv, sumy

      if (self%cfg%forcing_mode==1) then
         state%T(self%grid%ubnd_vol) = state%SST
      end if


      ! Calculation of heat (per layer), schmidt stability and mixing depht for model calibration
      call self%grid%interpolate_to_vol(self%grid%z_face,self%grid%Az, self%grid%nz_grid+1, Az_vol)

      volume = 0
      sum_z = 0

      do i = 1, self%grid%ubnd_vol
         rho0t(i) = 0.9998395_RK + state%T(i)*(6.7914e-5_RK + state%T(i)*(-9.0894e-6_RK + state%T(i)* &
                                             (1.0171e-7_RK + state%T(i)*(-1.2846e-9_RK + state%T(i)*(1.1592e-11_RK + state%T(i)*(-5.0125e-14_RK))))))
         rho0st(i) = (8.181e-4_RK + state%T(i)*(-3.85e-6_RK + state%T(i)*(4.96e-8_RK)))*state%S(i)
         rho(i) = rho_0*(rho0t(i) + rho0st(i))

         depth(i) = self%grid%z_volume(self%grid%ubnd_vol) - self%grid%z_volume(i)

         volume = volume + Az_vol(i)*self%grid%h(i)
         sum_z = sum_z + depth(i)*Az_vol(i)*self%grid%h(i)

         state%heat_per_layer(i) = Az_vol(i)*state%T(i)*cp*rho(i)/1e12
      end do
      zv = sum_z/volume

      state%schmidt_stability = 9.81/self%grid%Az(self%grid%ubnd_fce)*sum((depth - zv)*rho*Az_vol*self%grid%h)

      sumy = 0
      do i = 1, self%grid%ubnd_vol
         sumy = sumy + (depth(i) - zv)*rho(i)*Az_vol(i)*self%grid%h(i)
      end do
      state%schmidt_stability = 9.81/self%grid%Az(self%grid%ubnd_fce)*sumy

      do i = self%grid%ubnd_vol,1,-1
         if ((state%T(self%grid%ubnd_vol) - state%T(i)) > 1) then
            state%mixing_depth = depth(i)
            return
         end if
      end do

   end subroutine

end module strat_temp
