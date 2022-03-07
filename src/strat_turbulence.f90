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
!     | Turbulence module
!     | Updates state variables
!<    +---------------------------------------------------------------+




module strat_turbulence
   use strat_kinds
   use strat_consts
   use strat_grid
   use strat_simdata
   implicit none
   private

   ! Common Types
   type, public :: TurbulenceModule
      class(StaggeredGrid), pointer :: grid
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param

   contains
      procedure, pass :: init => turbulence_module_init
      procedure, pass :: update => turbulence_module_update

      procedure, pass :: do_production => turbulence_module_do_production
      procedure, pass :: do_seiche => turbulence_module_do_seiche

   end type
contains

   subroutine turbulence_module_init(self, state, grid, model_cfg, model_param)
      implicit none
      class(TurbulenceModule) :: self
      class(ModelState), target :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param

      self%grid => grid
      self%model_cfg => model_cfg
      self%model_param => model_param

      state%E_Seiche = model_param%seiche_ini
   end subroutine

   subroutine turbulence_module_update(self, state, param)
      implicit none
      class(TurbulenceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param

      call self%do_production(state)

      call self%do_seiche(state, param)


   end subroutine

   subroutine turbulence_module_do_production(self, state)
      !####################################################################
      implicit none
      class(TurbulenceModule) :: self
      class(ModelState) :: state

      associate (grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)

         ! Equation 5 (left) of Goudsmit, 2002
         ! P is defined on the inner faces
         state%P = 0
         state%P(2:ubnd_fce - 1) = (state%U(2:ubnd_vol) - state%U(1:ubnd_vol - 1))**2 + (state%V(2:ubnd_vol) - state%V(1:ubnd_vol - 1))**2
         state%P(2:ubnd_fce - 1) = state%P(2:ubnd_fce - 1)*state%num(2:ubnd_fce - 1)*grid%meanint(1:ubnd_vol - 1)**2
         ! Equation 5 (right) of Goudsmit, 2002
         state%B = 0
         state%B(2:ubnd_fce - 1) = -state%nuh(2:ubnd_fce - 1)*state%NN(2:ubnd_fce - 1)

         return
      end associate
   end subroutine

   subroutine turbulence_module_do_seiche(self, state, param)
      !####################################################################
      implicit none
      class(TurbulenceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param

      ! Local variables
      real(RK) :: W10, PS, PW, f_norm, minNN, a_seiche_local
      real(RK) :: distrib(self%grid%ubnd_fce)
      integer :: i

      associate (grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)
         minNN = 0
         distrib = 0

         ! Update distrib on inner faces
         do i = 2, ubnd_fce - 1
            distrib(i) = max(state%NN(i)**param%q_NN, minNN)/grid%Az(i)*grid%dAz(i - 1)
         end do

         ! Determine a_seiche
         if (self%model_cfg%split_a_seiche) then   ! If a_seiche is splitted seasonally

            ! If maximum stratification (N2) is higher than threshold
            if (maxval(state%NN(2:ubnd_fce - 1)) >= param%strat_sumr) then
               a_seiche_local = param%a_seiche

            ! If maximum stratification (N2) is lower than threshold
            else if (maxval(state%NN(2:ubnd_fce - 1)) < param%strat_sumr) then
               a_seiche_local = param%a_seiche_w
            end if
         else  ! If a_seiche is not splitted
            a_seiche_local = param%a_seiche
         end if

         ! Exit function if a_seiche is 0
         if (a_seiche_local == 0) then
            state%P_seiche = 0.0_RK
            return
         end if

         ! Calculate Seiche normalization factor
         f_norm = 0.0_RK
         if (self%model_cfg%seiche_normalization == 1) then ! max NN
            f_norm = maxval(state%NN(2:ubnd_fce - 1))

            f_norm = (f_norm**param%q_NN)*grid%Az(ubnd_fce)*rho_0
         else if (self%model_cfg%seiche_normalization == 2) then ! integral
            do i = 2, ubnd_fce - 1
               f_norm = f_norm + distrib(i)*grid%Az(i)*grid%h(i - 1)

            end do

            f_norm = f_norm*rho_0
         end if

         ! todo: direct float comparison...? OK?
         ! why is this code here?
         if (f_norm == 0.) then
            do i = 2, ubnd_fce - 1
               distrib(i) = 1/grid%h(i - 1)
            end do
            f_norm = grid%Az(ubnd_fce)*rho_0
         end if

         ! Adjust wind params based on configuration
         if (self%model_cfg%use_filtered_wind) then !use filtered wind (AG 2014)
            PW = a_seiche_local*grid%Az(ubnd_fce)*rho_air*state%C10*state%Wf**3
         else ! Use real wind
            W10 = sqrt(state%u10**2 + state%v10**2)
            PW = a_seiche_local*grid%Az(ubnd_fce)*rho_air*state%C10*W10**3
         end if

         ! Update E_Seiche
         PS = state%E_Seiche**(1.5_RK)*state%gamma
         state%E_Seiche = state%E_Seiche + (PW - PS)*state%dt

         ! Limit so that E_Seiche does not become negative
         if (state%E_Seiche < 0.) then
            PS = (PS*state%dt + state%E_Seiche)/state%dt
            state%E_Seiche = 0.0_RK
         end if

         ! Update P_Seiche
         ! Equation 24 in Goudsmit, 2002
         do i = 2, ubnd_fce - 1

            state%P_Seiche(i) = 1.0_RK/f_norm*distrib(i)*PS*(1.0_RK - 10*sqrt(param%CD))
         end do
         state%P_Seiche(1) = 0.0_RK
         state%P_Seiche(ubnd_fce) = 0.0_RK

         return
      end associate
   end subroutine

end module
