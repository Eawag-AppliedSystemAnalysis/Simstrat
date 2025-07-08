! ---------------------------------------------------------------------------------
!     Simstrat a physical 1D model for lakes and reservoirs
!
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     Copyright (C) 2020, Eawag
!     FABM: Copyright (C) 2014, Bolding & Bruggeman ApS
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
!     |  Simstrat - FABM interface
!<    +---------------------------------------------------------------+

module simstrat_fabm
   use fabm ! main FABM module
   use strat_simdata
   use strat_kinds
   use strat_grid
   use strat_solver
   use utilities
   use, intrinsic :: ieee_arithmetic

   implicit none
   private

   ! Type for use of FABM
   type, public :: SimstratFABM
      class(type_fabm_model), pointer :: fabm_model ! holds metadata on all bgc models active within FABM

      ! Arrays to hold the values of all biogeochemical state variables
      ! Where this memory resides and how it is laid out is typically host-specific

      ! Array for interior tracer source terms and vertical velocities
      real(RK), dimension(:,:), allocatable :: sms_int, velocity
      ! Arrays for fluxes and tracer source terms at bottom
      ! -> Consider that there is a bottom at every layer
      real(RK), dimension(:), allocatable :: flux_bt, sms_bt
      ! Arrays for fluxes and tracer source terms at surface
      real(RK), dimension(:), allocatable :: flux_sf, sms_sf

      ! Variables for validity of state
      logical :: repair = .true. ! Whether to clip all state variables to valid range from bgc models when update is called
      logical :: valid_int, valid_sf, valid_bt
   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
   end type SimstratFABM

contains

   ! Initialize: called once in within the initialization of Simstrat
   ! Sets up memory, reads FABM configuration from fabm.yaml, links the external Simstrat variables and sets the initial conditions of FABM variables
   subroutine init(self, state, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState) :: state
      class(StaggeredGrid) :: grid
      !class(type_fabm_model), pointer :: fabm_model

      ! Local variables
      integer :: ivar

      ! Create an alias of the model
      !fabm_model => self%fabm_model

      ! Initialize the model (reads run-time FABM model configuration from fabm.yaml and stores it in fabm_model)
      ! After this the number of biogeochemical variables is fixed
      ! (access variable metadata in fabm_model%interior_state_variables, fabm_model%interior_diagnostic_variables)
      ! The model interacts iwth FABM and describes properties of all bgc variables and parameters
      self%fabm_model => fabm_create_model()
      if (.not. associated(self%fabm_model)) then
         call error("FABM model creation failed")
      end if

      ! Provide extents of the spatial domain (number of layers nz for a 1D column)
      ! Used to allocate memory for FABM-managed spatially explicit fields
      ! Because the entire extent is given some layers might be calculated that are not physical, be aware of that in the output
      call self%fabm_model%set_domain(grid%nz_grid)

      ! At this point (after the call to fabm_create_model), memory should be
      ! allocated to hold the values of all size(fabm_model%*_state_variables) state variables
      ! All state variable values are combined in an array *_state with shape grid%nz_grid, state%n_fabm_*_state
      ! Interior state variables
      state%n_fabm_interior_state = size(self%fabm_model%interior_state_variables)
      allocate(state%fabm_interior_state(grid%nz_grid, state%n_fabm_interior_state))
      ! Bottom state variables
      state%n_fabm_bottom_state = size(self%fabm_model%bottom_state_variables)
      if (state%n_fabm_bottom_state > 0) then
         allocate(state%fabm_bottom_state(state%n_fabm_bottom_state))
      end if
      ! Surface state variables
      state%n_fabm_surface_state = size(self%fabm_model%surface_state_variables)
      if (state%n_fabm_surface_state > 0) then
         allocate(state%fabm_surface_state(state%n_fabm_surface_state))
      end if

      ! Total amount of states
      state%n_fabm_state = state%n_fabm_interior_state + state%n_fabm_bottom_state + state%n_fabm_surface_state

      ! Allocate memory for additional information (name, units, long_name, minimum, maximum, missing_value) on state variables
      allocate(state%fabm_state_names(state%n_fabm_state))

      ! Point FABM to fields that hold state variable data and get additional information

      ! Interior state variables
      if (state%n_fabm_interior_state > 0) then
         do ivar = 1, state%n_fabm_interior_state
            call self%fabm_model%link_interior_state_data(ivar, state%fabm_interior_state(:,ivar))
            state%fabm_state_names(ivar) = self%fabm_model%interior_state_variables(ivar)%name
            print *, "int: ", self%fabm_model%interior_state_variables(ivar)%minimum
         end do
      end if
      ! Bottom state variables
      if (state%n_fabm_bottom_state > 0) then
         do ivar = 1, state%n_fabm_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(ivar))
            state%fabm_state_names(state%n_fabm_interior_state + ivar) = self%fabm_model%bottom_state_variables(ivar)%name
            print *, "bot: ", self%fabm_model%bottom_state_variables(ivar)%minimum
         end do
      end if
      ! Surface state variables
      if (state%n_fabm_surface_state > 0) then
         do ivar = 1, state%n_fabm_surface_state
            call self%fabm_model%link_surface_state_data(ivar, state%fabm_surface_state(ivar))
            state%fabm_state_names(state%n_fabm_interior_state + state%n_fabm_bottom_state + ivar) = self%fabm_model%surface_state_variables(ivar)%name
         end do
      end if

      ! Get id for standard variable <variable>: if memory location of <variable> changes send updated pointers
      !type(type_fabm_interior_variable_id) :: id_var
      !id_var = fabm_model%get_interior_variable_id(fabm_standard_variables%variable)
      !call fabm_model%link_varkind(id_var, var)
      ! Determine whether a particular variable is needed by the bgc models in FABM
      !call fabm_model%varibale_needs_values(fabm_standard_variables%variable)

      ! Point FABM to fields that contain values for environmental data, all variables are assumed to be allocated
      ! Do this for all variables on FABM's standard variable list that the model can provide
      ! For this list, visit https://github.com/fabm-model/fabm/wiki/List-of-standard-variables
      ! Listed below are all non biogeochemical variables from that list
      
      ! Interior data (arrays)
      ! Attenuation coefficient of photosynthetically active radiative (PAR) flux, a fraction of shortwave radiative (SWR) flux [m-1]
      call self%fabm_model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, state%absorb_vol)
      ! Attenuation coefficient of SWR flux [m-1]
      call self%fabm_model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_shortwave_flux, state%absorb_vol)
      ! Note: Passing the same absorption coefficient for SWR and PAR is not exactly correct 
      ! since the absorption coefficient for PAR can be higher than for SWR
      ! see for example: https://www.sciencedirect.com/science/article/pii/S1001074217317734
      ! Thickness of each layer [m]
      call self%fabm_model%link_interior_data(fabm_standard_variables%cell_thickness, grid%h(1:grid%nz_grid))
      ! Density of each layer [kg m-3]
      call self%fabm_model%link_interior_data(fabm_standard_variables%density, state%rho)
      ! PAR flux [W m-2]
      call self%fabm_model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux, state%par_vol)
      ! SWR flux [W m-2]
      call self%fabm_model%link_interior_data(fabm_standard_variables%downwelling_shortwave_flux, state%swr_vol)
      ! Salinity at each layer [1e-3]
      call self%fabm_model%link_interior_data(fabm_standard_variables%practical_salinity, state%S)
      ! Pressure at each layer [dbar]: equal to layer depth, assuming one meter in depth is one dbar in pressure
      call self%fabm_model%link_interior_data(fabm_standard_variables%pressure, grid%layer_depth)
      ! Temperature at each layer [°C]
      call self%fabm_model%link_interior_data(fabm_standard_variables%temperature, state%T)
      ! Net rate of SWR energy absorption at each layer [W m-2], not defined in Simstrat
      !call fabm_model%link_interior_data(fabm_standard_variables%net_rate_of_absorption_of_shortwave_energy_in_layer)
      ! Vertical tracer diffusity [m2 s-1]: defined in GOTM, not defined in Simstrat
      !link_interior_data(type_interior_standard_variable(name='vertical_tracer_diffusivity', units='m2 s-1'))

      ! Scalars (horizontal data in a 1D model)
      ! Absolute depth [m]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%bottom_depth, grid%z_zero)
      ! Bottom stress [Pa]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%bottom_stress, state%u_taub)
      ! Cloud area fraction [-]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%cloud_area_fraction, state%Cloud)
      ! Ice area fraction [-]: 1 as soon as ice height is larger than ice_tolerance, 0 else
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%ice_area_fraction, state%ice_area_fraction)
      ! Surface air pressure [Pa]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_air_pressure, state%p_air)
      ! Surface albedo [-]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_albedo, state%wat_albedo)
      ! PAR flux at surface [W m-2]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux, state%par0)
      ! SWR flux at surface [W m-2]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, state%rad0)
      ! Surface drag coefficient [-]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_drag_coefficient_in_air, state%C10)
      ! Surface specific humidity [-]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_specific_humidity, state%qa)
      ! Surface temperature [°C]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_temperature, state%T_atm)
      ! Total wind speed [m s-1]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%wind_speed, state%uv10)
      ! latitude [degree_north]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%latitude, state%Lat)
      ! Depth relative to surface [m], not foud in FABM alhough in list, bot provided in GOTM
      ! call self%fabm_model%link_scalar(fabm_standard_variables%depth, grid%max_depth)
      ! Secchi depth [m], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%secchi_depth)
      ! Depth below geoid [m], provided in GOTM but only necessary when coupling with geodetic data, not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%bottom_depth_below_geoid)
      ! Bottom roughness [m], provided in GOTM, constant in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%bottom_roughness_length)
      ! PAR flux in air [W m-2], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux_in_air)
      ! SWR flux in air [W m-2], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux_in_air)
      ! Longitude [degree_east], provided in GOTM, not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%longitude)
      ! Number of days since start of the year [days], provided in GOTM, could be calculated in Simstrat but only if albedo is calculated
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%number_of_days_since_start_of_the_year)

      ! Complete initialization and check whether FABM has all dependencies fulfilled
      ! (i.e., whether all required calls to fabm_model%link_*_data have been made and all required data have been provided)
      ! Stop with fatal error if not
      ! Selection of diagnostics that FABM will compute and store becomes frozen
      call self%fabm_model%start()

      ! -> Program received signal SIGSEGV: Segmentation fault - invalid memory reference.
      ! Set FABM-provided initial values for state variables (tracers), typically space-independent.
      ! This sets the values of arrays sent to fabm_model%link_*_state_data,
      ! in this case those contained in *_state
      ! -> If model not initialized with custom or previously stored state (needs to be added as argument)
      call self%fabm_model%initialize_interior_state(1, grid%nz_grid)
      call self%fabm_model%initialize_bottom_state()
      call self%fabm_model%initialize_surface_state()

      ! Do the following sequence once to ensure diagnostics called in the model state valdation have been assigned a value

      ! 1. Prepare all fields FABM needs to compute source terms (e.g., light)
      call self%fabm_model%prepare_inputs()

      ! 2. Retrieve sources and fluxes across whole domain: order of call and of processing grid points is up to host
      
      ! 2a. Allocate and initialize with 0 (and then retrieve) interior tracer source terms (tracer units s-1)
      if (state%n_fabm_interior_state > 0) then
         allocate(self%sms_int(grid%nz_grid, state%n_fabm_interior_state))
         self%sms_int = 0.0
      end if
      !call self%fabm_model%get_interior_sources(1, grid%nz_grid, self%sms_int)

      ! 2b. Allocate and initialize with 0 (and then retrieve) fluxes and tracer source terms at bottom
      if (state%n_fabm_interior_state > 0) then
         allocate(self%flux_bt(state%n_fabm_interior_state))
         self%flux_bt = 0.0
      end if
      if (state%n_fabm_bottom_state > 0) then
         allocate(self%sms_bt(state%n_fabm_bottom_state))
         self%sms_bt = 0.0
      end if
      !call self%fabm_model%get_bottom_sources(self%flux_bt, self%sms_bt)

      ! 2c. Allocate and initialize with 0 (and then retrieve) fluxes and tracer source terms at surface
      if (state%n_fabm_interior_state > 0) then
         allocate(self%flux_sf(state%n_fabm_interior_state))
         self%flux_sf = 0.0
      end if
      if (state%n_fabm_surface_state > 0) then
         allocate(self%sms_sf(state%n_fabm_surface_state))
         self%sms_sf = 0.0
      end if
      !call self%fabm_model%get_surface_sources(self%flux_sf, self%sms_sf)

      ! 3. Allocate and initialize with 0 (and then retrieve) vertical velocities (sinking, floating, active movement) in m s-1
      allocate(self%velocity(grid%nz_grid, state%n_fabm_interior_state))
      self%velocity = 0.0
      !call self%fabm_model%get_vertical_movement(1,grid%nz_grid, self%velocity)

      ! 4. Compute any remaining diagnostics
      call self%fabm_model%finalize_outputs()

      ! At this point, initialization is complete
      ! Assign local alias back to self
      !self%fabm_model => fabm_model
   end subroutine init

   ! The update function is called in the main loop of simstrat (in simstrat.f90) at every time step
   ! Particle atmospheric, pelagic and benthic fluxes and diffusion are computed to update bgc state variable values
   subroutine update(self, state, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(StaggeredGrid) :: grid
      !class(type_fabm_model), pointer :: fabm_model

      ! Local variables
      integer :: ivar

      ! Create an alias of the model
      !fabm_model => self%fabm_model

      ! Follow steps 1 to 5 before output and time integration

      ! 1. Validate the model state, stop if variables are not valid
      ! Only run after first time step has been run
      ! -> maybe do this after the update again
      call self%fabm_model%check_interior_state(1 ,grid%nz_grid, self%repair, self%valid_int)
      call self%fabm_model%check_bottom_state(self%repair, self%valid_bt)
      call self%fabm_model%check_surface_state(self%repair, self%valid_sf)
      if (.not. (self%valid_int .and. self%valid_bt .and. self%valid_sf)) then
         if (self%repair) then
            ! call warn("FABM Variable repaired")
         else
            call error("FABM Variable out of bounds")
         end if
      end if

      ! 2. Prepare all fields (e.g. light attenuation) FABM needs to compute source terms (e.g. light)
      ! Operates on entire active spatial domain
      call self%fabm_model%prepare_inputs()

      ! 3. Retrieve sources [var_unit s-1] and fluxes [var_unit m s-1] across entire domain
      ! Order of call and of processing grid points is up to host
      ! >0: flux into water
      ! <0: flux out of water

      ! 3a. Retrieve interior tracer source terms
      call self%fabm_model%get_interior_sources(1, grid%nz_grid, self%sms_int)
      if (allocated(self%sms_int)) then
         if (any(ieee_is_nan(self%sms_int))) then
            ! call warn("FABM Interior Source contains NaNs, set to 0")
            where (ieee_is_nan(self%sms_int))
               self%sms_int = 0.0
            end where
         end if
      end if

      ! 3b. Retrieve fluxes over pelagic-benthic interface and bottom tracer source terms
      call self%fabm_model%get_bottom_sources(self%flux_bt, self%sms_bt)
      if (allocated(self%flux_bt)) then
         if (any(ieee_is_nan(self%flux_bt))) then
            ! call warn("FABM Bottom Flux contains NaN, set to 0")
            where (ieee_is_nan(self%flux_bt))
               self%flux_bt = 0.0
            end where
         end if
      end if
      if (allocated(self%sms_bt)) then
         if (any(ieee_is_nan(self%sms_bt))) then
            ! call warn("FABM Bottom Source contains NaN, set to 0")
            where (ieee_is_nan(self%sms_bt))
               self%sms_bt = 0.0
            end where
         end if
      end if

      ! 3c. Retrieve fluxes over air-water surface and surface tracer source terms
      call self%fabm_model%get_surface_sources(self%flux_sf, self%sms_sf)
      if (allocated(self%flux_sf)) then
         if (any(ieee_is_nan(self%flux_sf))) then
            ! call warn("FABM Surface Flux contains NaN, set to 0")
            where (ieee_is_nan(self%flux_sf))
               self%flux_sf = 0.0
            end where
         end if
      end if
      if (allocated(self%sms_sf)) then
         if (any(ieee_is_nan(self%sms_sf))) then
            ! call warn("FABM Surface Source contains NaN, set to 0")
            where (ieee_is_nan(self%sms_sf))
               self%sms_sf = 0.0
            end where
         end if
      end if

      ! 4. Retrieve local vertical velocities of pelagic state variables [m s-1]
      ! Movement through water, independent of water flow
      ! >0: upward movement (floating, active movement)
      ! <0: downward ovement (sinking, sedimentation, active movement)
      call self%fabm_model%get_vertical_movement(1, grid%nz_grid, self%velocity)
      if (allocated(self%velocity)) then
         if (any(ieee_is_nan(self%velocity))) then
            ! call warn("FABM Interior Velocity contains NaN, set to 0")
            where (ieee_is_nan(self%velocity))
               self%velocity = 0.0
            end where
         end if
      end if

      ! 5. Compute any remaining diagnostics not computed by preceding routines
      ! Operates on entire active spatial domain
      call self%fabm_model%finalize_outputs()

      ! This would be the ideal moment for the output: 
      ! All variables (enviromental, state, mask, diagnostics, source, fluxes, vertical velocities) are in sync.
      ! -> Simstrat does output after update_fabm is done: the state variables are one step further (does not matter if only those are output)

      ! Time-integrate the advection-diffusion-reaction equations (method up to host)
      ! of all tracers, combining the Simstrat transport terms with the FABM biogeochemical source
      ! terms and fluxes (sms, flux) and vertical velocities (velocity). This results in an updated interior_state.
      do ivar = 1, state%n_fabm_interior_state
         call diffusion_FABM_interior_state(self, state, grid, ivar)
      end do

      ! Direct time integration of source terms to update bottom_state and surface_state
      if (state%n_fabm_bottom_state > 0) state%fabm_bottom_state = state%fabm_bottom_state + state%dt * self%sms_bt
      if (state%n_fabm_surface_state > 0) state%fabm_surface_state = state%fabm_surface_state + state%dt * self%sms_sf

      ! Assign local alias back to self
      !self%fabm_model => fabm_model
   end subroutine update

   ! Diffusion algorithm for interior state variables: Simstrat transport terms and FABM biogeochemical terms integrated simultaneously
   ! Assuming small enough dt such that only fluxes between neighbouring layers are relevant
   subroutine diffusion_fabm_interior_state(self, state, grid, ivar)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(StaggeredGrid) :: grid
      integer, intent(in) :: ivar
      
      ! Total Flux
         ! real(RK), dimension(0:grid%nz_grid+1) :: flux_tot ! Include air and benthic layer
         ! real(RK), dimension(grid%nz_grid) :: boundaries, source, lower_diag, main_diag, upper_diag, rhs
         ! flux_tot(1:grid%nz_grid) = self%velocity(:, ivar) * state%fabm_interior_state(:, ivar) ! Local flux
         ! if (self%flux_bt(ivar) > 0) then
         !    flux_tot(0) = self%flux_bt(ivar) ! flux into water
         ! else
         !    flux_tot(0) = 0 ! -> verify: flux out of water is same as negative local flux at 1
         ! end if
         ! if (self%flux_sf(ivar) > 0) then
         !    flux_tot(grid%nz_grid+1) = -self%flux_sf(ivar) ! flux into water
         ! else
         !    flux_tot(grid%nz_grid+1) = 0 ! -> verify: flux out of water is same as positive local flux at grid%nz_grid
         ! end if
      ! Total Flux

      ! Local variables
      real(RK), dimension(grid%nz_grid) :: velocity_up, velocity_down, sources, lower_diag, main_diag, upper_diag, rhs
      real(RK), dimension(grid%nz_grid-1) :: AreaFactor_ext_1, AreaFactor_ext_2

      ! Create linear system of equations with transport and source terms
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}

      ! Build diagonals of A with Simstrat transport terms (code from discretization module)
      ! With diffusivity for temperature nuh
      ! Downward flux to/from each layer
      upper_diag(1) = 0.0_RK
      upper_diag(2:grid%nz_grid) = state%dt*state%nuh(2:grid%nz_grid)*grid%AreaFactor_1(2:grid%nz_grid)
      ! Upward flux to/from each layer
      lower_diag(1:grid%nz_grid-1) = state%dt*state%nuh(2:grid%nz_grid)*grid%AreaFactor_2(1:grid%nz_grid - 1)
      lower_diag(grid%nz_grid) = 0.0_RK
      ! 1 - downward flux - upward flux
      main_diag(:) = 1.0_RK - upper_diag(:) - lower_diag(:)

      ! Add FABM fluxes to A as residual verical advection terms
      ! Area factors for external fluxes
      AreaFactor_ext_1(1:grid%nz_grid-1) = grid%Az(2:grid%nz_grid)/grid%Az_vol(2:grid%nz_grid)/grid%h(2:grid%nz_grid)
      AreaFactor_ext_2(1:grid%nz_grid-1) = grid%Az(2:grid%nz_grid)/grid%Az_vol(1:grid%nz_grid-1)/grid%h(1:grid%nz_grid-1)
      ! Upward and downward movement rates for each layer (positive)
      velocity_up = self%velocity(:, ivar)
      where (velocity_up < 0.0_RK)
         velocity_up = 0.0_RK
      end where
      velocity_down = -self%velocity(:, ivar)
      where (velocity_down < 0.0_RK)
         velocity_down = 0.0_RK
      end where
      ! Downward flux to each layer
      upper_diag(2:grid%nz_grid) = upper_diag(2:grid%nz_grid) + state%dt*velocity_down(2:grid%nz_grid)*AreaFactor_ext_1(1:grid%nz_grid-1)
      ! Upward flux to each layer
      lower_diag(1:grid%nz_grid-1) = lower_diag(1:grid%nz_grid-1) + state%dt*velocity_up(1:grid%nz_grid-1)*AreaFactor_ext_2(1:grid%nz_grid-1)
      ! Flux away from each layer
      main_diag(:) = main_diag(:) - state%dt*self%velocity(:, ivar)/grid%h(1:grid%nz_grid)

      ! Get source S^{n}
      ! Source at each layer
      sources = self%sms_int(:, ivar)
      ! Add pelagic-benthic and air-water flux [var_unit m s-1] as source [var_unit s-1]
      ! Convert to source by division by height of bottomost / uppermost layer [m]
      sources(1) = self%flux_bt(ivar) / grid%h(1)
      sources(grid%nz_grid) = self%flux_sf(ivar) / grid%h(grid%nz_grid)

      ! Calculate RHS (phi^{n}+dt*S^{n})
      rhs(:) = state%fabm_interior_state(:, ivar) + state%dt*sources(:)

      ! Solve LES to get phi^{n+1}
      call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%fabm_interior_state(:, ivar), grid%nz_grid)
   end subroutine diffusion_fabm_interior_state

   ! Deallocate memory
   subroutine deallocate_fabm(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self

      ! Deallocate model
      if (associated(self%fabm_model)) deallocate(self%fabm_model)
      ! Deallocate internal arrays
      if (allocated(self%sms_int)) deallocate(self%sms_int)
      if (allocated(self%flux_bt)) deallocate(self%flux_bt)
      if (allocated(self%sms_bt)) deallocate(self%sms_bt)
      if (allocated(self%flux_sf)) deallocate(self%flux_sf)
      if (allocated(self%sms_sf)) deallocate(self%sms_sf)
      if (allocated(self%velocity)) deallocate(self%velocity)
   end subroutine deallocate_fabm

   ! Process feedbacks from bgc to physics (absorption, albedo, wind drag changes, ...)

   ! -> Light absorption feedback by FAMB variables
   ! subroutine absorption_update_fabmS(self, state)

      ! ! Arguments
      ! class(SimstratFABM) :: self
      ! class(ModelState) :: state

      ! ! Local variables
      ! integer :: i
      ! real(RK) :: bio_extinction

      ! do i=self%grid%nz_occupied, 1, -1
      !    bio_extinction = 0.0_RK
      !    call fabm_light_extinction(self%column, i, bio_extinction) : change for fabm
      !    state%absorb_vol(i) = self%fabm_cfg%background_extinction + bio_extinction

      ! end do
      ! ! Interpolate to faces to be compatible with Simstrat temperature module
      ! call self%grid%interpolate_to_face(self%grid%z_volume, state%absorb_vol, self%grid%nz_occupied, state%absorb)

   ! end subroutine

   ! -> Calculate photosynthetically active radiation (PAR) and short wave
   ! radiation (SWR) over entire column, using surface short wave radiation,
   ! and background and biotic extinction.
   !     GOTM: Copyright by the GOTM-team under the GNU Public License - www.gnu.org
   !     Original author(s): Jorn Bruggeman (rem.: for every subroutine)
   ! subroutine light(nlev)
      ! !INPUT PARAMETERS:
      ! integer, intent(in) :: nlev
      ! !LOCAL VARIABLES:
      ! integer :: i
      ! real(RK) :: bioext
      ! real(RK) :: localexts(1:nlev)
      ! bioext = 0
      ! call fabm_get_light_extinction(model,1,nlev,localexts(1:nlev))
      ! do i=nlev,1,-1
      !    ! Add the extinction of the first half of the grid box.
      !    bioext = bioext+localexts(i)*curh(i)/2
      !    ! Calculate photosynthetically active radiation (PAR), shortwave radiation, and PAR attenuation.
      !    par(i) = I_0*(1-A)*exp(-z(i)/g2-bioext)
      !    swr(i) = par(i)+I_0*A*exp(-z(i)/g1)
      !    k_par(i) = 1/g2+localexts(i)
      !    ! Add the extinction of the second half of the grid box.
      !    bioext = bioext+localexts(i)*curh(i)/2
      ! end do
   ! end subroutine light

end module simstrat_fabm