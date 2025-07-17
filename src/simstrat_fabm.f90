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
      ! two-dimensional if there is a bottom at every layer
      real(RK), dimension(:,:), allocatable :: flux_bt, sms_bt
      ! Arrays for fluxes and tracer source terms at surface
      real(RK), dimension(:), allocatable :: flux_sf, sms_sf

      ! Variables for validity of state
      logical :: valid_int, valid_sf, valid_bt

      ! Index of current bottom location (pelagic-benthic interface)
      integer, dimension(:), pointer :: bottom_index
      ! Variable to enumerate over all layers if bottom_everywhere is true
      ! 1 if bottom_everywhere is false, nz_grid if bottom_everywhere is true
      ! All associated code copied and adapted from https://gitlab.com/wateritech-public/waterecosystemstool/gotm
      integer :: kmax_bot

      integer :: att_index ! Index of attenuation_coefficient_of_photosynthetic_radiative_flux in FABM diagnostic variables
      integer, dimension(:), allocatable :: diagnostic_index ! Index of FABM diagnostic variables
   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
   end type SimstratFABM

contains

   ! Initialize: called once in within the initialization of Simstrat
   ! Sets up memory, reads FABM configuration from fabm_cfg_file, links the external Simstrat variables and sets the initial conditions of FABM variables
   subroutine init(self, state, fabm_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState) :: state
      class(FABMConfig) :: fabm_cfg
      class(StaggeredGrid) :: grid
      !class(type_fabm_model), pointer :: fabm_model

      ! Local variables
      integer :: ivar, k

      ! Create an alias of the model
      !fabm_model => self%fabm_model

      ! make sure everything is deallocated
      call deallocate_fabm(self)

      ! Initialize the model: Reads run-time FABM model configuration from fabm_cfg%fabm_cfg_file (YAML format)
      ! and stores it in fabm_model. After this the number of biogeochemical variables is fixed
      ! (access variable metadata in fabm_model%interior_state_variables, fabm_model%interior_diagnostic_variables)
      ! The model interacts with FABM and describes properties of all bgc variables and parameters
      self%fabm_model => fabm_create_model(fabm_cfg%fabm_config_file)
      if (.not. associated(self%fabm_model)) then
         call error("FABM model creation failed")
      end if

      ! Provide extents of the spatial domain (number of layers nz for a 1D column)
      ! Used to allocate memory for FABM-managed spatially explicit fields
      ! Because the entire extent is given some layers might be calculated that are not physical, be aware of that in the output
      call self%fabm_model%set_domain(grid%nz_grid)

      ! Allocate bottom index with one value and link FABM to it
      ! Set bottom location (pelagic-benthic interface) to 1 (default)
      allocate(self%bottom_index(1))
      self%bottom_index(1) = 1
      call self%fabm_model%set_bottom_index(self%bottom_index(1))

      ! Set kmax_bot to nz_grid if bottom_everywhere is true
      self%kmax_bot = 1
      if (fabm_cfg%bottom_everywhere) self%kmax_bot = grid%nz_grid

      ! At this point (after the call to fabm_create_model), memory should be
      ! allocated to hold the values of all size(fabm_model%*_state_variables) state variables
      ! All state variable values are combined in an array *_state with shape grid%nz_grid, state%n_fabm_*_state
      ! Interior state variables
      state%n_fabm_interior_state = size(self%fabm_model%interior_state_variables)
      if (state%n_fabm_interior_state > 0) then
         allocate(state%fabm_interior_state(grid%nz_grid, state%n_fabm_interior_state))
      end if
      ! Bottom state variables
      state%n_fabm_bottom_state = size(self%fabm_model%bottom_state_variables)
      if (state%n_fabm_bottom_state > 0) then
         allocate(state%fabm_bottom_state(self%kmax_bot, state%n_fabm_bottom_state))
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
         end do
      end if
      ! Bottom state variables: link for bottom-most layer to fulfill FABM requirements, link for other layers later
      if (state%n_fabm_bottom_state > 0) then
         do ivar = 1, state%n_fabm_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1,ivar))
            state%fabm_state_names(state%n_fabm_interior_state + ivar) = self%fabm_model%bottom_state_variables(ivar)%name
         end do
      end if
      ! Surface state variables
      if (state%n_fabm_surface_state > 0) then
         do ivar = 1, state%n_fabm_surface_state
            call self%fabm_model%link_surface_state_data(ivar, state%fabm_surface_state(ivar))
            state%fabm_state_names(state%n_fabm_interior_state + state%n_fabm_bottom_state + ivar) = self%fabm_model%surface_state_variables(ivar)%name
         end do
      end if
      
      ! Allocate interior tracer source terms [var_unit s-1]
      ! FABM increments rather than sets argument sms_int: initialize with 0 at every time step
      if (state%n_fabm_interior_state > 0) then
         allocate(self%sms_int(grid%nz_grid, state%n_fabm_interior_state))
      end if

      ! Allocate  fluxes over pelagic-benthic interface [var_unit m s-1] and bottom tracer source terms [var_unit s-1]
      if (state%n_fabm_interior_state > 0 .or. state%n_fabm_bottom_state > 0) then
         allocate(self%flux_bt(self%kmax_bot, state%n_fabm_interior_state))
         allocate(self%sms_bt(self%kmax_bot, state%n_fabm_bottom_state))
      end if

      ! Allocate retrieve fluxes over air-water surface [var_unit m s-1] and surface tracer source terms [var_unit s-1]
      if (state%n_fabm_interior_state > 0 .or. state%n_fabm_surface_state > 0) then
         allocate(self%flux_sf(state%n_fabm_interior_state))
         allocate(self%sms_sf(state%n_fabm_surface_state))
      end if

      ! Allocate local vertical velocities of pelagic state variables (sinking, floating, active movement) [m s-1]
      ! Movement through water, independent of water flow
      if (state%n_fabm_interior_state > 0) then
         allocate(self%velocity(grid%nz_grid, state%n_fabm_interior_state))
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
      
      if (fabm_cfg%bioshade_feedback) then
         do ivar = 1, size(self%fabm_model%interior_diagnostic_variables)
            if (self%fabm_model%interior_diagnostic_variables(ivar)%name == 'attenuation_coefficient_of_photosynthetic_radiative_flux') then
               self%fabm_model%interior_diagnostic_variables(ivar)%save = .true.
               self%att_index = ivar
            end if
         end do
      end if

      if (fabm_cfg%output_diagnostic_variables) then
         call set_fabm_diagnostic_vars(self, state, fabm_cfg)

         if (state%n_fabm_diagnostic_interior > 0) then
            allocate(state%fabm_diagnostic_interior(grid%nz_grid, state%n_fabm_diagnostic_interior))
         end if

         if (state%n_fabm_diagnostic_horizontal > 0) then
            allocate(state%fabm_diagnostic_horizontal(state%n_fabm_diagnostic_horizontal))
         end if
      end if
      
      ! Complete initialization and check whether FABM has all dependencies fulfilled
      ! (i.e., whether all required calls to fabm_model%link_*_data have been made and all required data have been provided)
      ! Stop with fatal error if not
      ! Selection of diagnostics that FABM will compute and store becomes frozen
      call self%fabm_model%start()

      ! Set FABM-provided initial values for state variables (tracers), typically space-independent.
      ! This sets the values of arrays sent to fabm_model%link_*_state_data,
      ! in this case those contained in *_state
      ! -> If model not initialized with custom or previously stored state (needs to be added as argument)
      ! Interior state variables
      call self%fabm_model%initialize_interior_state(1, grid%nz_grid)
      ! Bottom state variables
      call self%fabm_model%initialize_bottom_state()
      ! If bottom_everywhere is set, at every depth:
      ! FABM is pointed to location that holds state data for the current depth
      ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
      ! The bottom state at the current depth is initialized
      if (fabm_cfg%bottom_everywhere) then
         do k = 2, self%kmax_bot
            do ivar = 1, state%n_fabm_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k, ivar))
            end do
            self%bottom_index(1) = k
            call self%fabm_model%initialize_bottom_state()
         end do
         ! Reset Botom to 1
         do ivar = 1, state%n_fabm_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
         end do
         self%bottom_index(1) = 1
      end if
      ! Surface state variables
      call self%fabm_model%initialize_surface_state()

      ! Call the update function once as a first call to initialize fluxes, sources and vertical velocities
      call self%update(state, fabm_cfg, grid, .true.)
   end subroutine init

   ! The update function is called in the main loop of simstrat (in simstrat.f90) at every time step
   ! Particle atmospheric, pelagic and benthic fluxes and diffusion are computed to update bgc state variable values
   subroutine update(self, state, fabm_cfg, grid, first_call)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(StaggeredGrid) :: grid  
      logical, intent(in), optional :: first_call

      ! Local variables
      integer :: ivar, k, index
      logical :: first_call_local

      ! Default not first call
      if (present(first_call)) then
         first_call_local = first_call
      else
         first_call_local = .false.
      end if

      ! If it is not the first call: 1. calculate and 2. validate the new model state, stop if variables are not valid
      if (.not. first_call_local) then
         ! 1a. Time-integrate the advection-diffusion-reaction equations
         ! of all tracers, combining the Simstrat transport terms with the FABM biogeochemical source
         ! terms and fluxes (sms, flux) and vertical velocities (velocity). This results in an updated interior_state.
         do ivar = 1, state%n_fabm_interior_state  
            call diffusion_FABM_interior_state(self, state, fabm_cfg, grid, ivar)
            ! Check for negative values
            if (any(state%fabm_interior_state(1: grid%nz_occupied, ivar) < 0.0_RK)) then
               do k = 1, size(state%fabm_interior_state(1:grid%nz_occupied, ivar))
                  if (state%fabm_interior_state(k, ivar) < 0.0_RK) then
                     print *, 'FABM Interior Variable value is ', state%fabm_interior_state(k, ivar)
                     print *, 'at grid point ', k
                     print *, 'at time (days, seconds) = ', state%simulation_time
                  end if
               end do
               call error('FABM Variable '//self%fabm_model%interior_state_variables(ivar)%name//' below zero.')
            end if
         end do

         ! 1b. Direct time integration of source terms to update bottom_state
         do ivar = 1, state%n_fabm_bottom_state
            state%fabm_bottom_state(:, ivar) = state%fabm_bottom_state(:, ivar) + state%dt * self%sms_bt(:, ivar)
            ! Check for negative values
            if (any(state%fabm_bottom_state(:, ivar) + state%dt * self%sms_bt(:, ivar) < 0.0_RK)) then
               do k = 1, self%kmax_bot
                  if (state%fabm_bottom_state(k, ivar) + state%dt * self%sms_bt(k, ivar) < 0.0_RK) then
                     print *, 'FABM Bottom Variable value is ', state%fabm_bottom_state(k, ivar)
                     print *, 'at grid point ', k
                     print *, 'at time (days, seconds) = ', state%simulation_time
                  end if
                  call error('FABM Bottom Variable '//self%fabm_model%bottom_state_variables(ivar)%name//' below zero.')
               end do
            end if
         end do

         ! 1c. Direct time integration of source terms to update surface_state
         do ivar = 1, state%n_fabm_surface_state
            state%fabm_surface_state(ivar) = state%fabm_surface_state(ivar) + state%dt * self%sms_sf(ivar)
            ! Check for negative values
            if (state%fabm_surface_state(ivar) < 0.0_RK) then
               print *, 'FABM Surface Variable value is ', state%fabm_surface_state(ivar)
               print *, 'at time (days, seconds) = ', state%simulation_time
               call error('FABM Surface Variable '//self%fabm_model%surface_state_variables(ivar)%name//' below zero.')
            end if
         end do

         ! 2a.Interior state variables
         call self%fabm_model%check_interior_state(1, grid%nz_grid, fabm_cfg%repair_fabm, self%valid_int)
         
         ! 2b. Bottom state variables
         call self%fabm_model%check_bottom_state(fabm_cfg%repair_fabm, self%valid_bt)
         if (fabm_cfg%bottom_everywhere) then
            ! If bottom_everywhere is set, at every depth:
            ! FABM is pointed to location that holds state data for the current depth
            ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
            ! The bottom state at the current depth is validated
            do k = 2, self%kmax_bot
               do ivar = 1, state%n_fabm_bottom_state
                  call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k, ivar))
               end do
               self%bottom_index(1) = k
               call self%fabm_model%check_bottom_state(fabm_cfg%repair_fabm, self%valid_bt)
            end do
            ! Reset Bottom to 1
            do ivar = 1, state%n_fabm_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
            end do
            self%bottom_index(1) = 1
         end if

         ! 2c. Surface state variables
         call self%fabm_model%check_surface_state(fabm_cfg%repair_fabm, self%valid_sf)
         ! Error if out of bounds and not repaired
         if (.not. (self%valid_int .and. self%valid_bt .and. self%valid_sf)) then
            if (fabm_cfg%repair_fabm) then
               call warn("FABM Variable repaired")
            else
               call error("FABM Variable out of bounds")
            end if
         end if
      end if

      ! 1. Prepare all fields (e.g. light attenuation) FABM needs to compute fluxes and source terms
      ! Operates on entire active spatial domain
      call self%fabm_model%prepare_inputs()

      ! 2. Retrieve sources, fluxes and vertical velocities
      ! >0: flux into water
      ! <0: flux out of water
      
      ! 2a. Initialize with 0 and then retrieve interior tracer source terms
      ! FABM increments rather than sets argument sms_int: initialize with 0 at every time step
      if (state%n_fabm_interior_state > 0) then
         self%sms_int = 0.0
         call self%fabm_model%get_interior_sources(1, grid%nz_grid, self%sms_int)
         ! Set NaNs to 0
         if (any(ieee_is_nan(self%sms_int))) then
            call warn("FABM Interior Source contains NaNs, set to 0")
            where (ieee_is_nan(self%sms_int))
               self%sms_int = 0.0
            end where
         end if
      end if

      ! 2b. Initialize with 0 and then retrieve fluxes over pelagic-benthic interface and bottom tracer source terms
      ! FABM increments rather than sets argument flux_bt, sms_bt: initialize with 0 at every time step
      if (state%n_fabm_interior_state > 0 .or. state%n_fabm_bottom_state > 0) then
         if (state%n_fabm_interior_state > 0) self%flux_bt = 0.0
         if (state%n_fabm_bottom_state > 0) self%sms_bt = 0.0
         call self%fabm_model%get_bottom_sources(self%flux_bt(1, :), self%sms_bt(1, :))
         if (fabm_cfg%bottom_everywhere) then
            ! If bottom_everywhere is set, at every depth:
            ! FABM is pointed to location that holds state data for the current depth
            ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
            ! The fluxes and sources at the current depth are calculated
            do k = 2, self%kmax_bot
               do ivar = 1, state%n_fabm_bottom_state
                  call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k, ivar))
               end do
               self%bottom_index(1) = k
               call self%fabm_model%get_bottom_sources(self%flux_bt(k, :), self%sms_bt(k, :))
            end do
            ! Reset Botom to 1
            do ivar = 1, state%n_fabm_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
            end do
            self%bottom_index(1) = 1
         end if
         ! Set NaNs to 0
         if (state%n_fabm_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_bt))) then
               call warn("FABM Bottom Flux contains NaN, set to 0")
               where (ieee_is_nan(self%flux_bt))
                  self%flux_bt = 0.0
               end where
            end if
         end if
         if (state%n_fabm_bottom_state > 0) then
            if (any(ieee_is_nan(self%sms_bt))) then
               call warn("FABM Bottom Source contains NaN, set to 0")
               where (ieee_is_nan(self%sms_bt))
                  self%sms_bt = 0.0
               end where
            end if
         end if
      end if

      ! 2c. Initialize with 0 and then retrieve fluxes over air-water surface and surface tracer source terms
      ! FABM increments rather than sets argument flux_sf, sms_sf: initialize with 0 at every time step
      if (state%n_fabm_interior_state > 0 .or. state%n_fabm_surface_state > 0) then
         if (state%n_fabm_interior_state > 0) self%flux_sf = 0.0
         if (state%n_fabm_surface_state > 0) self%sms_sf = 0.0
         call self%fabm_model%get_surface_sources(self%flux_sf, self%sms_sf)
         ! Set NaNs to 0
         if (state%n_fabm_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_sf))) then
               call warn("FABM Surface Flux contains NaN, set to 0")
               where (ieee_is_nan(self%flux_sf))
                  self%flux_sf = 0.0
               end where
            end if
         end if
         if (state%n_fabm_surface_state > 0) then
            if (any(ieee_is_nan(self%sms_sf))) then
               call warn("FABM Surface Source contains NaN, set to 0")
               where (ieee_is_nan(self%sms_sf))
                  self%sms_sf = 0.0
               end where
            end if
         end if
      end if

      ! 3. Initialize with 0 and retrieve local vertical velocities of pelagic state variables
      ! >0: upward movement (floating, active movement)
      ! <0: downward movement (sinking, sedimentation, active movement)
      if (state%n_fabm_interior_state > 0) then
         if (first_call_local) self%velocity = 0.0
         call self%fabm_model%get_vertical_movement(1, grid%nz_grid, self%velocity)
         ! Set NaNs to 0
         if (any(ieee_is_nan(self%velocity))) then
            call warn("FABM Interior Velocity contains NaN, set to 0")
            where (ieee_is_nan(self%velocity))
               self%velocity = 0.0
            end where
         end if
      end if

      ! 4. Compute any remaining diagnostics not computed by preceding routines
      ! Operates on entire active spatial domain
      call self%fabm_model%finalize_outputs()

      ! 5. Save diagnostic variables
      if (fabm_cfg%output_diagnostic_variables) then
         do ivar = 1, state%n_fabm_diagnostic_interior
            index = self%diagnostic_index(ivar)
            state%fabm_diagnostic_interior(:, ivar) = self%fabm_model%get_interior_diagnostic_data(index)
         end do
         do ivar = 1, state%n_fabm_diagnostic_horizontal
            index = self%diagnostic_index(state%n_fabm_diagnostic_interior + ivar)
            state%fabm_diagnostic_horizontal(ivar) = self%fabm_model%get_horizontal_diagnostic_data(index)
         end do
      end if      
      
      ! 6. Retrieve attenuation coefficient after state variables have been calculated once
      if (fabm_cfg%bioshade_feedback .and. (.not. first_call_local)) then
         call absorption_update_fabm(self, state, grid)
      end if

      ! This is the ideal moment for the output: 
      ! All variables (enviromental, state, diagnostics, source, fluxes, vertical velocities) are in sync
   end subroutine update

   ! Diffusion algorithm for interior state variables: Simstrat transport terms and FABM biogeochemical terms integrated simultaneously
   ! Assuming small enough dt such that only fluxes between neighbouring layers are relevant
   subroutine diffusion_fabm_interior_state(self, state, fabm_cfg, grid, ivar)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(StaggeredGrid) :: grid
      integer, intent(in) :: ivar

      ! Local variables
      real(RK), dimension(grid%nz_occupied) :: velocity_up, velocity_down, sources, lower_diag_diff, main_diag_diff, upper_diag_diff, lower_diag, main_diag, upper_diag, rhs, AreaFactor_ext_1, AreaFactor_ext_2
      integer :: k

      ! Create linear system of equations with transport and source terms
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}

      ! Build diagonals of A with Simstrat transport terms (code from discretization module)
      ! With diffusivity for temperature (nuh)
      ! Upward flux from (1:grid%nz_occupied - 1) to (2:grid%nz_occupied): upper_diag(z) describes the upward flux into cell z
      upper_diag_diff(2:grid%nz_occupied) = state%dt*state%nuh(2:grid%nz_occupied)*grid%AreaFactor_1(2:grid%nz_occupied)
      upper_diag_diff(1) = 0.0_RK ! No upward flux into bottommost cell
      ! Downward flux from (2:grid%nz_occupied) to (1:grid%nz_occupied - 1): lower_diag(z) describes the downward flux into cell z
      lower_diag_diff(1:grid%nz_occupied-1) = state%dt*state%nuh(2:grid%nz_occupied)*grid%AreaFactor_2(1:grid%nz_occupied - 1)
      lower_diag_diff(grid%nz_occupied) = 0.0_RK ! No downward flux into uppermost cell
      ! 1 - downward flux out - upward flux out
      main_diag_diff(1:grid%nz_occupied) = 1.0_RK - upper_diag_diff(1:grid%nz_occupied) - lower_diag_diff(1:grid%nz_occupied)

      ! Add FABM fluxes to A as residual verical advection terms
      ! Benthic-pelagic and air-water fluxes are added as source terms below
      ! Area factors for external fluxes, minus sign to be closer to the AreaFactors from grid
      ! AreaFactor_ext_1 is for upward fluxes:
      ! Multiply by the layer through which the flux passes (z = i) and divide by the volume of the receiving cell (z = i)
      AreaFactor_ext_1(1:grid%nz_occupied) = -grid%Az(1:grid%nz_occupied) / (grid%h(1:grid%nz_occupied) * grid%Az_vol(1:grid%nz_occupied))
      ! AreaFactor_ext_2 is for downward fluxes: 
      ! Multiply by the layer through which the flux passes (z = i+1) and divide by the volume of the receiving cell (z = i)
      AreaFactor_ext_2(1:grid%nz_occupied) = -grid%Az(2:grid%nz_occupied+1) / (grid%h(1:grid%nz_occupied) * grid%Az_vol(1:grid%nz_occupied))
      ! Upward and downward movement rates for each layer
      ! if the variable moves upward velocity_up is positive and velocity_down zero
      ! vice versa if the variable moves downward
      do k = 1, grid%nz_occupied
         if (self%velocity(k, ivar) >= 0) then
            velocity_up(k) = self%velocity(k, ivar)
            velocity_down(k) = 0
         else
            velocity_up(k) = 0
            velocity_down(k) = -self%velocity(k, ivar)
         end if
      end do
      ! Add upward flux from (1:grid%nz_occupied-1) to (2:grid%nz_occupied): upper_diag(z) describes the upward flux into cell z
      ! upper_diag should be <= 0
      ! Benthic-pelagic flux is treated as source term below
      upper_diag(2:grid%nz_occupied) = upper_diag_diff(2:grid%nz_occupied) + state%dt*velocity_up(1:grid%nz_occupied-1)*AreaFactor_ext_1(2:grid%nz_occupied)
      upper_diag(1) = upper_diag_diff(1)
      ! Subtract upward flux from cell z to cell z+1
      ! main_diag should be >= 0
      ! Water-air flux is treated as sink term below
      main_diag(1:grid%nz_occupied-1) = main_diag_diff(1:grid%nz_occupied-1) - state%dt*velocity_up(1:grid%nz_occupied-1)*AreaFactor_ext_2(1:grid%nz_occupied-1)
      main_diag(grid%nz_occupied) = main_diag_diff(grid%nz_occupied)
      ! Add downward flux from (2:grid%nz_occupied) to (1:grid%nz_occupied - 1): lower_diag(z) describes the downward flux into cell z
      ! lower_diag should be <= 0
      ! Air-water flux is treated as source term below
      lower_diag(1:grid%nz_occupied-1) = lower_diag_diff(1:grid%nz_occupied-1) + state%dt*velocity_down(2:grid%nz_occupied)*AreaFactor_ext_2(1:grid%nz_occupied-1)
      lower_diag(grid%nz_occupied) = lower_diag_diff(grid%nz_occupied)
      ! Subtract downward flux from cell z to cell z-1
      ! main_diag should be >= 0
      ! Pelagic-benthic flux is treated as sink term below
      main_diag(2:grid%nz_occupied) = main_diag(2:grid%nz_occupied) - state%dt*velocity_down(2:grid%nz_occupied)*AreaFactor_ext_1(2:grid%nz_occupied)
      main_diag(1) = main_diag(1)

      ! Get source S^{n}
      ! Source at each layer
      sources = self%sms_int(1:grid%nz_occupied, ivar)
      ! Add pelagic-benthic and air-water flux [var_unit m s-1] as source [var_unit s-1]
      ! Convert to source by division by height of current bottom / uppermost layer [m]
      if (fabm_cfg%bottom_everywhere) then
         sources(:) = sources(:) + (self%flux_bt(1:grid%nz_occupied, ivar) / grid%h(1:grid%nz_occupied)) ! pelagic-benthic flux at every layer
      else
         sources(1) = sources(1) + (self%flux_bt(1, ivar) / grid%h(1)) ! pelagic-benthic flux only at bottommost layer
      end if
      sources(grid%nz_occupied) = sources(grid%nz_occupied) + (self%flux_sf(ivar) / grid%h(grid%nz_occupied))

      ! Calculate RHS (phi^{n}+dt*S^{n})
      rhs(:) = state%fabm_interior_state(1:grid%nz_occupied, ivar) + state%dt*sources(:)

      ! Solve LES to get phi^{n+1}
      call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%fabm_interior_state(1:grid%nz_occupied, ivar), grid%nz_occupied)
   end subroutine diffusion_fabm_interior_state

   ! Read names of diagnostic Vars in SetDiagnosticVars file
   subroutine set_fabm_diagnostic_vars(self, state, fabm_cfg)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg

      integer :: i, j, n, n_int, n_hor, unit, status
      character(len=100) :: line
      integer, parameter :: max_lines = 100
      character(len=80), dimension(max_lines) :: temp_names
      integer, dimension(max_lines) :: temp_index
      logical, dimension(max_lines) :: is_int

      open(newunit=unit, action='read', status='old', file=fabm_cfg%set_diag_vars, iostat = status)

      if (status .ne. 0) then
         call warn("No Output of FABM Diagnostic Variables")
      else
         write(6,*) 'Reading ', fabm_cfg%set_diag_vars
      end if

      n = 0
      n_int = 0
      n_hor = 0
      do
         read(unit, '(A)', iostat=status) line
         if (status .ne. 0) exit
         n = n + 1
         if (n > max_lines) then
            call error('Too many lines in '//fabm_cfg%set_diag_vars//', increase max_lines in set_fabm_diagnostic_vars')
         end if
         ! Trim line and assign to temp_names, cut to length 48
         temp_names(n) = adjustl(trim(line))
         do j = 1, size(self%fabm_model%interior_diagnostic_variables)
            if (self%fabm_model%interior_diagnostic_variables(j)%name == temp_names(n)) then
               self%fabm_model%interior_diagnostic_variables(j)%save = .true.
               temp_index(n) = j
               is_int(n) = .true.
               n_int = n_int + 1
            end if
         end do
         do j = 1, size(self%fabm_model%horizontal_diagnostic_variables)
            if (self%fabm_model%horizontal_diagnostic_variables(j)%name == temp_names(n)) then
               self%fabm_model%horizontal_diagnostic_variables(j)%save = .true.
               temp_index(n) = j
               is_int(n) = .false.
               n_hor = n_hor + 1
            end if
         end do
      end do

      state%n_fabm_diagnostic = n
      state%n_fabm_diagnostic_interior = n_int
      state%n_fabm_diagnostic_horizontal = n_hor

      ! Allocate array of proper size
      allocate(state%fabm_diagnostic_names(n))
      allocate(self%diagnostic_index(n))

      ! Copy names to array
      j = 1
      do i = 1, n
         if (is_int(i)) then
            state%fabm_diagnostic_names(j) = temp_names(i)
            self%diagnostic_index(j) = temp_index(i)
            j = j + 1
         end if
      end do
      do i = 1, n
         if (.not. is_int(i)) then
            state%fabm_diagnostic_names(j) = temp_names(i)
            self%diagnostic_index(j) = temp_index(i)
         end if
      end do
   end subroutine set_fabm_diagnostic_vars

   ! Light absorption feedback by FABM variables
   subroutine absorption_update_fabm(self, state, grid)
      ! Arguments
      class(SimstratFABM) :: self
      class(ModelState) :: state
      class(StaggeredGrid) :: grid

      ! Local variables
      real(RK), dimension(grid%nz_grid) :: attenuation_coefficient_of_photosynthetic_radiative_flux

      ! Retrieve attenuation_coefficient_of_photosynthetic_radiative_flux and set as absorb_vol
      attenuation_coefficient_of_photosynthetic_radiative_flux = self%fabm_model%get_interior_diagnostic_data(self%att_index)
      state%absorb_vol(1 : grid%nz_occupied) = attenuation_coefficient_of_photosynthetic_radiative_flux(1 : grid%nz_occupied)

      ! Interpolate to faces to be compatible with Simstrat temperature module
      call grid%interpolate_to_face(grid%z_volume, state%absorb_vol, grid%nz_occupied, state%absorb)
   end subroutine absorption_update_fabm

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
      if (associated(self%bottom_index)) deallocate(self%bottom_index)
      if (allocated(self%diagnostic_index)) deallocate(self%diagnostic_index)
   end subroutine deallocate_fabm

end module simstrat_fabm