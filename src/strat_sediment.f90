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
!     |  Sediment heat model (ModelConfig.SedimentHeat = 1, issue #88)
!     |     Heat storage in the sediment and heat exchange with the water
!<    +---------------------------------------------------------------+
!
! Below each grid cell there is a sediment column (n_sed layers, from dz_sed_top = 0.02 m at the top to ~4 m at the
! bottom, 8.9 m in total; the annual temperature wave penetrates ~2 m). The area of the column is the sediment area of
! the cell: Az(i+1) - Az(i) (sloping lake bed), plus Az(1) (flat bottom) for the lowest cell.
! Each column exchanges heat with
!    - the water of its cell (conductive flux between the water temperature and the top sediment layer),
!    - the shortwave radiation reaching the lake bed of the cell (absorbed in the top layer; state%sed_sw from strat_temp,
!      which removes it from the water),
!    - the geothermal heat flux (param%fgeo, entering at the bottom of the column; not added to the water directly).
! The sediment is updated every dt_sed (1 h, strat_consts) after the temperature equation: the water of the cell is the
! top node of the column (heat capacity of the cell per sediment area) and the column and the water are solved together
! implicitly. This is stable for any time step and conserves energy exactly. The absorbed shortwave radiation is
! accumulated between the updates.
! Columns of cells above the lake surface keep their temperature until the cells are flooded again.
! Thermal properties (lambda_sed, rhoc_sed) of water-saturated fine sediment are defined in strat_consts.

module strat_sediment
   use strat_kinds
   use strat_consts
   use strat_simdata
   use strat_grid
   use utilities
   implicit none
   private

   type, public :: SedimentModule
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param
      real(RK) :: dz(n_sed)                     ! Layer thicknesses [m]
      real(RK) :: g(0:n_sed - 1)                ! Conductances [W/m2/K]: g(0) water - layer 1, g(k) layer k - layer k+1
      real(RK), allocatable :: T_sed(:, :)      ! Sediment temperature [°C] (layer, cell)
      real(RK), allocatable :: flux(:)          ! Heat flux from the sediment into the water of each cell [W/m2]
      real(RK), allocatable :: sw_sum(:)        ! Shortwave radiation absorbed by the sediment since the last update [J]
      integer :: n_steps = 1                    ! Number of model time steps per sediment time step
      integer :: step = 0                       ! Model time steps since the last sediment update
   contains
      procedure, pass :: init => sediment_init
      procedure, pass :: update => sediment_update
      procedure, pass :: save => sediment_save
      procedure, pass :: load => sediment_load
      procedure, pass :: deallocate => sediment_deallocate
   end type

contains

   subroutine sediment_init(self, state, grid, param)
      implicit none
      class(SedimentModule) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelParam), target :: param
      integer :: k, i

      self%grid => grid
      self%param => param

      ! Layer thicknesses and conductances
      self%dz(1) = dz_sed_top
      do k = 2, n_sed
         self%dz(k) = self%dz(k - 1)*dz_sed_factor
      end do
      self%g(0) = lambda_sed/(0.5_RK*self%dz(1))
      do k = 1, n_sed - 1
         self%g(k) = lambda_sed/(0.5_RK*(self%dz(k) + self%dz(k + 1)))
      end do

      ! Initial sediment temperature: water temperature of the cell (top water temperature above the lake surface).
      ! The deep layers adjust within 1-2 years (spin-up). Overwritten when continuing from a snapshot.
      call sediment_deallocate(self)
      allocate (self%T_sed(n_sed, grid%nz_grid), self%flux(grid%nz_grid), self%sw_sum(grid%nz_grid))
      do i = 1, grid%nz_grid
         self%T_sed(:, i) = state%T(min(i, grid%ubnd_vol))
      end do
      self%flux = 0.0_RK
      self%sw_sum = 0.0_RK
      self%n_steps = max(1, nint(dt_sed/state%dt))
      self%step = 0

      call ok('Sediment heat model on (heat storage, absorbed light and geothermal flux in the sediment)')
   end subroutine

   ! Update the sediment temperature and the water temperature of the cells (after the temperature equation, every dt_sed)
   subroutine sediment_update(self, state)
      implicit none
      class(SedimentModule) :: self
      class(ModelState) :: state

      ! Node 0: water of the cell, nodes 1..n_sed: sediment layers
      real(RK) :: a(0:n_sed), b(0:n_sed), c(0:n_sed), r(0:n_sed), cap(0:n_sed), x(0:n_sed), w
      real(RK) :: a_sed, vol, dt
      integer :: i, k

      associate (grid => self%grid, T => state%T, T_sed => self%T_sed)
         ! Accumulate the absorbed shortwave radiation; update only every n_steps time steps
         self%sw_sum(1:grid%ubnd_vol) = self%sw_sum(1:grid%ubnd_vol) + state%sed_sw(1:grid%ubnd_vol)*state%dt
         self%step = self%step + 1
         if (self%step < self%n_steps) return
         dt = self%step*state%dt
         self%step = 0

         cap(1:n_sed) = rhoc_sed*self%dz/dt ! Heat capacity per time step [W/m2/K]

         do i = 1, grid%ubnd_vol
            a_sed = max(grid%Az(i + 1) - grid%Az(i), 0.0_RK)
            if (i == 1) a_sed = a_sed + grid%Az(1) ! Flat bottom
            if (a_sed <= 0) then
               self%flux(i) = 0.0_RK
               cycle
            end if
            vol = 0.5_RK*(grid%Az(i) + grid%Az(i + 1))*grid%h(i)
            cap(0) = rho_0*cp*vol/a_sed/dt ! Water of the cell per sediment area

            ! Implicit heat conduction: cap*(T_new - T_old) = conduction + sources
            do k = 0, n_sed
               a(k) = 0.0_RK; c(k) = 0.0_RK
               b(k) = cap(k)
               if (k == 0) then
                  r(k) = cap(k)*T(i)
               else
                  r(k) = cap(k)*T_sed(k, i)
               end if
               if (k > 0) then
                  a(k) = -self%g(k - 1); b(k) = b(k) + self%g(k - 1)
               end if
               if (k < n_sed) then
                  c(k) = -self%g(k); b(k) = b(k) + self%g(k)
               end if
            end do
            r(1) = r(1) + self%sw_sum(i)/a_sed/dt   ! Absorbed shortwave radiation (mean over dt)
            r(n_sed) = r(n_sed) + self%param%fgeo   ! Geothermal heat flux from below

            ! Thomas algorithm
            do k = 1, n_sed
               w = a(k)/b(k - 1)
               b(k) = b(k) - w*c(k - 1)
               r(k) = r(k) - w*r(k - 1)
            end do
            x(n_sed) = r(n_sed)/b(n_sed)
            do k = n_sed - 1, 0, -1
               x(k) = (r(k) - c(k)*x(k + 1))/b(k)
            end do

            ! Heat flux from the sediment into the water (mean over dt) and new temperatures
            self%flux(i) = cap(0)*(x(0) - T(i))
            T(i) = x(0)
            T_sed(:, i) = x(1:n_sed)
         end do
         self%sw_sum = 0.0_RK
      end associate
   end subroutine

   ! Sediment temperatures are written after the other snapshot data (see save_snapshot in simstrat.f90)
   subroutine sediment_save(self)
      implicit none
      class(SedimentModule) :: self
      write(80) self%T_sed, self%sw_sum, self%step
   end subroutine

   ! Returns .false. if the snapshot does not contain sediment temperatures (the initial values are then kept)
   logical function sediment_load(self) result(success)
      implicit none
      class(SedimentModule) :: self
      real(RK), dimension(size(self%T_sed, 1), size(self%T_sed, 2)) :: tmp
      real(RK), dimension(size(self%sw_sum)) :: tmp_sw
      integer :: status, tmp_step
      read(81, iostat=status) tmp, tmp_sw, tmp_step
      success = (status == 0)
      if (success) then
         self%T_sed = tmp; self%sw_sum = tmp_sw; self%step = tmp_step
      end if
   end function

   subroutine sediment_deallocate(self)
      implicit none
      class(SedimentModule) :: self
      if (allocated(self%T_sed)) deallocate(self%T_sed)
      if (allocated(self%flux)) deallocate(self%flux)
      if (allocated(self%sw_sum)) deallocate(self%sw_sum)
   end subroutine

end module strat_sediment
