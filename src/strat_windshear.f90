!<    +---------------------------------------------------------------+
!     | Implementation of a statevar for U/V variable
!<    +---------------------------------------------------------------+


module strat_windshear
   use strat_kinds
   use strat_consts
   use strat_simdata
   use strat_statevar
   use strat_grid
   use strat_solver
   implicit none
   private

   type, extends(ModelVariable), public :: UVModelVar
      real(RK), pointer :: stress_t
   contains
      procedure, pass(self), public :: calc_terms => uv_var_calc_terms
      procedure, pass(self), public :: assign_shear_stress => uv_var_assign_shear_stress
   end type

contains
  ! Method to assign shear stress variable acting on this state variable
   subroutine uv_var_assign_shear_stress(self, stress_t)
      class(UVModelVar), intent(inout) :: self
      real(RK), target ::stress_t
      self%stress_t => stress_t
   end subroutine

   ! Calculate source terms for u or v
   subroutine uv_var_calc_terms(self, state, param, sources, boundaries)
      class(UVModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK) :: integr, length
      real(RK), dimension(:) ::  sources, boundaries
      real(RK), dimension(size(self%var)) ::  uv_norm
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         if (self%cfg%pressure_gradients == 1) then
            integr = sum(self%var(1:ubnd_vol))
            length = sqrt(grid%Az(ubnd_fce))
         else
            integr = 0
            length = 0
         end if

         uv_norm = sqrt(state%U**2 + state%V**2)

         !!!!!!!! Define sources !!!!!!!!
         sources = 0
         if (self%cfg%pressure_gradients == 1) then !Svensson 1978
            sources(2:ubnd_vol - 1) = -pi**2*rho_0*g*integr/ubnd_vol*grid%max_depth/length**2
         elseif (self%cfg%pressure_gradients == 2) then !???
            sources(2:ubnd_vol - 1) = -state%drag*self%var(2:ubnd_vol - 1)* &
                                      uv_norm(2:ubnd_vol - 1)* &
                                      grid%dAz(2:ubnd_vol - 1)/grid%Az(2:ubnd_vol - 1)
         end if

         ! Set surface condition based on shear stress variable
         sources(ubnd_vol) = self%stress_t/grid%h(ubnd_vol)

         !!!!!!!! Explicit boundary conditions !!!!!
         boundaries(1) = state%drag*uv_norm(1)/grid%h(1)
         boundaries(2:ubnd_vol) = 0

      end associate
   end subroutine

end module strat_windshear
