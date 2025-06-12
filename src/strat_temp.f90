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
      !real(RK) :: T_soil
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         ! Radiation reaching top layer
         state%rad(self%grid%ubnd_fce) = state%rad0/rho_0/cp ![°C*m/s]  , rad is on faces
         state%swr(self%grid%ubnd_fce) = state%rad0 ![J/s/m2]
         state%par(self%grid%ubnd_fce) = state%rad0*swr_par ![J/s/m2]

         ! Radiation reaching a layer is equal to radiation in the layer above minus absorption
         ! Absorption is different for swr and par so the last calculation is not exactly correct
         do i = ubnd_fce - 1, 1, -1
            state%rad(i) = state%rad(i + 1)*exp(-grid%h(i)*(state%absorb(ubnd_fce - i)+state%absorb(ubnd_fce + 1 - i))/2) !Attenuated by absorption
            state%swr(i) = state%rad(i)*state%rho_0*cp
            state%par(i) = state%rad(i)*swr_par
         end do

         ! Radiation staying in a layer (including the heat which goes into the sediment, therefore the weighting using Az)
         ! Absorption is different for swr and par so the last calculation is not exactly correct
         state%rad_vol(1:ubnd_vol) = (state%rad(2:ubnd_fce)*grid%Az(2:ubnd_fce) - state%rad(1:ubnd_fce - 1)*grid%Az(1:ubnd_fce - 1))/(grid%Az(2:ubnd_fce) + grid%Az(1:ubnd_fce - 1))*2
         state%swr_vol(1:ubnd_vol) = state%rad_vol(1:ubnd_vol)*rho_0*cp
         state%par_vol(1:ubnd_vol) = state%swr_vol(1:ubnd_vol)*rho_0*cp

         !!!!!!!! Define sources !!!!!!!!
         ! Add Hsol Term to sources (Eq 1, Goudsmit(2002))
         sources(1:ubnd_vol) = state%rad_vol(1:ubnd_vol)/grid%h(1:ubnd_vol)

         ! Experimenting with influence of soil temperature on water temperature using a sinusoidal soil temperature         
         !T_soil = 5 + 5*cos(2*pi/365*(state%current_month*30 + state%current_day - 270))

         !sources(1:ubnd_vol) = sources(1:ubnd_vol) + 2.4*(T_soil - state%T(1:ubnd_vol))/0.5*(grid%Az(2:ubnd_fce) - grid%Az(1:ubnd_fce-1))/rho_0/cp/grid%h(1:ubnd_vol)/(grid%Az(2:ubnd_fce) + grid%Az(1:ubnd_fce - 1))*2

         ! Set boundary heat flux at surface (Eq 25, Goudsmit(2002))
         ! FB, 2021: Correction of relevant area (accoring to tests by Ana Ayala)
         sources(ubnd_vol) = sources(ubnd_vol) + state%heat/rho_0/cp/grid%h(ubnd_vol)*grid%Az(ubnd_fce)/(grid%Az(ubnd_fce) + grid%Az(ubnd_fce - 1))*2

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

      if (self%cfg%forcing_mode==1) then
         state%T(self%grid%ubnd_vol) = state%SST
      end if

   end subroutine

end module strat_temp
