!<    +---------------------------------------------------------------+
!     |  Stability Module
!     |  - Contains methods to update cmue_cn /qe and NN
!<    +---------------------------------------------------------------+

module strat_stability
   use strat_kinds
   use strat_consts
   use strat_grid
   use strat_simdata
   use utilities
   implicit none
   private

   ! Common Types
   type, public :: StabilityModule
      class(StaggeredGrid), pointer :: grid
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param

   contains
      procedure, pass :: init => stability_module_init
      procedure, pass :: update => stability_module_update
      procedure, pass :: update_cmue_cn => stability_module_update_cmue_cn
      procedure, pass :: update_cmue_qe => stability_module_update_cmue_qe
      procedure, pass :: update_NN => stability_module_update_NN

   end type
contains

   subroutine stability_module_init(self, grid, model_cfg, model_param)
      implicit none
      class(StabilityModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param

      self%grid => grid
      self%model_cfg => model_cfg
      self%model_param => model_param
   end subroutine

   ! Update state variables
   subroutine stability_module_update(self, state)
      implicit none
      class(StabilityModule) :: self
      class(ModelState) :: state
      real(RK), dimension(self%grid%ubnd_fce) :: beta

      ! Do buoyancy update (update NN)
      call self%update_NN(state%T, state%S, state%rho, state%NN)

      ! Update cmue depending on selected stabilty function
      if (self%model_cfg%stability_func == 1) then
         call self%update_cmue_cn(state%cmue1, state%cmue2)

      else if (self%model_cfg%stability_func == 2) then
         beta(1:self%grid%ubnd_fce) = state%NN(1:self%grid%ubnd_fce)*(state%k(1:self%grid%ubnd_fce)/state%eps(1:self%grid%ubnd_fce))**2

         beta(1) = 0
         beta(self%grid%ubnd_fce) = 0

         call self%update_cmue_qe(beta, state%cmue1, state%cmue2, state%cde)

      end if

   end subroutine

   ! Compute NN from T and salinity
   subroutine stability_module_update_NN(self, T, S, rho, NN)
      implicit none
      class(StabilityModule) :: self

      ! Global variables
      real(RK), dimension(:), intent(in) :: T, S
      real(RK), dimension(:), intent(inout) :: NN, rho

      ! Local variables
      real(RK) :: buoy(self%grid%length_fce)
      real(RK) :: rho0t(self%grid%length_fce), rho0st(self%grid%length_fce)
      integer :: i

      associate (grd=>self%grid)

         do i = 1, grd%ubnd_vol
            call calc_density(rho(i),T(i),S(i))
            buoy(i) = -g*(rho(i) - rho_0)/rho_0
         end do

         NN(2:grd%ubnd_fce - 1) = grd%meanint(1:grd%ubnd_vol - 1)*(buoy(2:grd%ubnd_fce - 1) - buoy(1:grd%ubnd_fce - 2))
         NN(1) = NN(2)
         NN(grd%ubnd_fce) = NN(grd%ubnd_fce - 1)

      end associate
   end subroutine

   subroutine stability_module_update_cmue_cn(self, cmue1, cmue2)
      implicit none

      ! Global variables
      class(StabilityModule) :: self
      real(RK), dimension(:), intent(inout) :: cmue1, cmue2

      ! Standard version of k-eps model
      cmue1 = cmue
      cmue2 = cmue/Prndtl

      ! Set boundaries
      cmue1(1) = cmue1(2)
      cmue2(1) = cmue2(2)
      cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce - 1)
      cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce - 1)
   end subroutine

   subroutine stability_module_update_cmue_qe(self, beta, cmue1, cmue2, cde)
      implicit none
      class(StabilityModule) :: self
      real(RK), dimension(:), intent(in) :: beta
      real(RK), dimension(:), intent(inout) ::cmue1, cmue2
      real(RK) :: cde
      real(RK) :: gh, sm, sh
      integer :: i

      do i = 2, self%grid%ubnd_fce - 1
         gh = -cde**2*0.5_RK*beta(i)
         if (gh > 0.02) gh = gh - (gh - 0.02_RK)**2/(gh + 0.0233_RK - 2*0.02_RK)
         if (gh < -0.28) gh = -0.28_RK

         sm = 1.0_RK - 3*c1 - 6*a1/b1 - 3*a2*gh*((b2 - 3*a2)*(1.0_RK - 6*a1/b1) - 3*c1*(b2 + 6*a1))
         sm = a1*sm/((1.0_RK - 3*a2*gh*(6*a1 + b2))*(1.0_RK - 9*a1*a2*gh))
         sh = a2*(1.0_RK - 6*a1/b1)/(1.0_RK - 3*a2*gh*(6*a1 + b2))

         cmue1(i) = sqrt(2.0_RK)*cde*sm
         cmue2(i) = sqrt(2.0_RK)*cde*sh
      end do

      ! Set boundaries
      cmue1(1) = cmue1(2)
      cmue2(1) = cmue2(2)
      cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce - 1)
      cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce - 1)
   end subroutine

end module
