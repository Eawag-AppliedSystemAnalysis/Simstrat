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
         ! Only inner faces (boundary values are set to zero below, eps can be zero there)
         beta(2:self%grid%ubnd_fce - 1) = state%NN(2:self%grid%ubnd_fce - 1)*(state%k(2:self%grid%ubnd_fce - 1)/state%eps(2:self%grid%ubnd_fce - 1))**2

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
      real(RK) :: P, rho_above, rho_below, Tf, Sf, alpha
      integer :: i

      associate (grd=>self%grid)

         ! Density at atmospheric pressure (state%rho, used for plunging inflows and AED2; not used for N2)
         do i = 1, grd%ubnd_vol
            call calc_density(rho(i),T(i),S(i))
         end do

         ! Brunt-Vaisala frequency N2 with pressure effects (issue #87), following M. Schmid (Lake Baikal version of k-eps):
         !    N2 = g*alpha(T,S,p)*(dT/dz + Gamma) - g*beta(T,S,p)*dS/dz
         ! 1) Pressure-dependent expansivity alpha(T,S,p) (and contraction coefficient beta):
         !    The two cells adjacent to face i are both evaluated at the pressure of face i (local potential density),
         !    so that the density difference only reflects T and S, with alpha and beta valid at that depth.
         !    Equivalent to alpha from Chen and Millero (1986) with the pressure term (0.3682-1.520e-2*T+1.91e-4*T^2)*p,
         !    but exact for large T differences and including the salinity gradient.
         ! 2) Adiabatic lapse rate Gamma = g*alpha*(T+273.15)/cp: T is in-situ temperature, so a parcel displaced
         !    downwards warms (alpha > 0) or cools (alpha < 0) adiabatically. This adds (g*alpha)^2*(T+273.15)/cp >= 0.
         !    Relevant for deep lakes near the temperature of maximum density, where N2 is very small.
         do i = 2, grd%ubnd_fce - 1
            P = rho_0*g*(grd%z_face(grd%ubnd_fce) - grd%z_face(i))*1e-5_RK ! Pressure above atmospheric pressure at face i [bar]
            ! 1) Density difference of the two adjacent cells at the same pressure P
            call calc_density(rho_above, T(i), S(i), P)
            call calc_density(rho_below, T(i - 1), S(i - 1), P)
            NN(i) = -g*grd%meanint(i - 1)*(rho_above - rho_below)/rho_0
            ! 2) Adiabatic term, with alpha at face i from a central difference of the density (+-0.01 K) at pressure P
            Tf = 0.5_RK*(T(i) + T(i - 1)); Sf = 0.5_RK*(S(i) + S(i - 1))
            call calc_density(rho_above, Tf + 0.01_RK, Sf, P)
            call calc_density(rho_below, Tf - 0.01_RK, Sf, P)
            alpha = -(rho_above - rho_below)/(0.02_RK*rho_0)
            NN(i) = NN(i) + (g*alpha)**2*(Tf + 273.15_RK)/cp
         end do
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
