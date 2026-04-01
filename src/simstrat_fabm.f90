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
   ! Main FABM module
   use fabm
   ! Definition of FABM types, add from lib/fabm/src/fabm_types.F90, added:
   ! type_interior_standard_variable: Interior standard variable symbol
   ! type_dependency_id: Interior environmental dependency identifier
   use fabm_types, only: type_interior_standard_variable, type_dependency_id
   ! Simstrat modules
   use strat_simdata
   use strat_kinds
   use strat_grid
   use strat_solver
   use utilities
   ! Intrinsic modules
   use, intrinsic :: ieee_arithmetic

   implicit none
   private

   ! Target variable for index of current bottom location
   ! Variable to enumerate over all layers if bottom_everywhere is true
   ! All associated code copied and adapted from https://gitlab.com/wateritech-public/waterecosystemstool/gotm
   integer, target :: k_bot
   ! Bottom size: 1 if bottom_everywhere is false, nz_grid if bottom_everywhere is true
   integer :: kmax_bot
   ! Index of attenuation_coefficient_of_photosynthetic_radiative_flux in FABM diagnostic variables
   integer :: att_index

   ! Necessary bools for evaluation of bottom everywhere
   logical :: interior_or_bottom_exist = .false.
   logical :: diagnostic_horizontal_exist =.false.

   ! Type for use of FABM
   ! Contains arrays accessed or set by FABM
   type, public :: SimstratFABM
      ! FABM model: holds metadata on all biogeochemical (bgc) models active within FABM
      class(type_fabm_model), pointer :: fabm_model => null()

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
      integer, pointer :: bottom_index => null()
      
      ! Define additional FABM variables, used by bgc models
      ! Simstrat contribution to attenuation_coefficient_of_shortwave_flux and attenuation_coefficient_of_photosynthetic_radiative_flux (background extinction)
      type(type_dependency_id) :: id_attenuation_coefficient_host
      ! Projection factor for benthic flux into horizontal layer volume
      type(type_interior_standard_variable) :: bot_pel_conv = type_interior_standard_variable(name='bot_pel_conv', units='-')
      type(type_fabm_interior_variable_id) :: id_bot_pel_conv
   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
      procedure, pass(self), public :: absorption_update_fabm
   end type SimstratFABM

contains

   ! Initialize: called once within the initialization of Simstrat
   ! Sets up memory, reads FABM configuration from fabm_cfg_file, links the external Simstrat variables and sets the initial conditions of FABM variables
   subroutine init(self, state, fabm_cfg, sim_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(SimConfig), intent(in) :: sim_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      integer :: ivar, index

      ! Make sure everything is deallocated
      call deallocate_fabm(self)

      ! Create the bgc models according to the configurations in FABMConfigFile
      ! The models interact with FABM and describe properties of all bgc variables and parameters
      self%fabm_model => fabm_create_model(fabm_cfg%config_file, initialize = .false.)
      if (.not. associated(self%fabm_model)) then
         call error('FABM model creation failed')
      end if

      ! Register bgc variables from Simstrat (<ID>, <NAME>, <UNITS>, <LONG_NAME>)
      ! Register simstrat contribution to attenuation coefficients (background extinction), if bioshade feedback is on
      if (fabm_cfg%bioshade_feedback) then
         call self%fabm_model%root%register_dependency(self%id_attenuation_coefficient_host, 'attenuation_coefficient_host', 'm-1', 'host contribution to attenuation coefficients')
      end if
      
      ! Add Simstrat bgc variables to FABM aggregate variables
      ! Add simstrat contribution to attenuation coefficients (background extinction), if bioshade feedback is on
      if (fabm_cfg%bioshade_feedback) then
         call self%fabm_model%root%add_to_aggregate_variable(fabm_standard_variables%attenuation_coefficient_of_shortwave_flux, self%id_attenuation_coefficient_host)
         call self%fabm_model%root%add_to_aggregate_variable(fabm_standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_attenuation_coefficient_host)
      end if

      ! Initialize the model: Reads run-time FABM model configuration from fabm_cfg%fabm_cfg_file (YAML format)
      ! and stores it in fabm_model. After this the number of biogeochemical variables is fixed
      ! (access variable metadata in fabm_model%interior_state_variables, fabm_model%interior_diagnostic_variables)
      call self%fabm_model%initialize()

      ! Provide extents of the spatial domain (number of layers nz for a 1D column)
      ! Used to allocate memory for FABM-managed spatially explicit fields
      ! Provide seconds_per_time_unit when time filtering functionality is used
      ! The product of seconds_per_time_unit and the t argument provided to prepare_inputs must have units of seconds
      call self%fabm_model%set_domain(grid%nz_grid, 1.0_RK)

      ! Set bottom location (pelagic-benthic interface) to 1 (default)
      ! Link FABM to bottom location
      k_bot = 1
      self%bottom_index => k_bot
      call self%fabm_model%set_bottom_index(self%bottom_index)

      ! Set kmax_bot to nz_grid if bottom_everywhere is true, 1 else
      kmax_bot = 1
      if (fabm_cfg%bottom_everywhere) kmax_bot = grid%nz_grid

      ! At this point (after the call to fabm_create_model), memory should be
      ! allocated to hold the values of all size(fabm_model%*_state_variables) state variables
      ! All state variable values are combined in an array *_state with shape grid%nz_grid, state%n_fabm_*_state
      ! Interior state variables
      state%n_fabm_interior_state = size(self%fabm_model%interior_state_variables)
      if (state%n_fabm_interior_state > 0) then
         interior_or_bottom_exist = .true.
         allocate(state%fabm_interior_state(grid%nz_grid, state%n_fabm_interior_state))
         allocate(state%min_int(state%n_fabm_interior_state))
         allocate(state%max_int(state%n_fabm_interior_state))
      end if
      ! Bottom state variables
      state%n_fabm_bottom_state = size(self%fabm_model%bottom_state_variables)
      if (state%n_fabm_bottom_state > 0) then
         interior_or_bottom_exist = .true.
         allocate(state%fabm_bottom_state(kmax_bot, state%n_fabm_bottom_state))
         allocate(state%min_bt(state%n_fabm_bottom_state))
         allocate(state%max_bt(state%n_fabm_bottom_state))
      end if
      ! Surface state variables
      state%n_fabm_surface_state = size(self%fabm_model%surface_state_variables)
      if (state%n_fabm_surface_state > 0) then
         allocate(state%fabm_surface_state(state%n_fabm_surface_state))
         allocate(state%min_sf(state%n_fabm_surface_state))
         allocate(state%max_sf(state%n_fabm_surface_state))
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
            state%min_int(ivar) = self%fabm_model%interior_state_variables(ivar)%minimum
            state%max_int(ivar) = self%fabm_model%interior_state_variables(ivar)%maximum
         end do
      end if
      ! Bottom state variables: link for bottom-most layer to fulfill FABM requirements, link for other layers later
      if (state%n_fabm_bottom_state > 0) then
         do ivar = 1, state%n_fabm_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1,ivar))
            state%fabm_state_names(state%n_fabm_interior_state + ivar) = self%fabm_model%bottom_state_variables(ivar)%name
            state%min_bt(ivar) = self%fabm_model%bottom_state_variables(ivar)%minimum
            state%max_bt(ivar) = self%fabm_model%bottom_state_variables(ivar)%maximum
         end do
      end if
      ! Surface state variables
      if (state%n_fabm_surface_state > 0) then
         do ivar = 1, state%n_fabm_surface_state
            call self%fabm_model%link_surface_state_data(ivar, state%fabm_surface_state(ivar))
            state%fabm_state_names(state%n_fabm_interior_state + state%n_fabm_bottom_state + ivar) = self%fabm_model%surface_state_variables(ivar)%name
            state%min_sf(ivar) = self%fabm_model%surface_state_variables(ivar)%minimum
            state%max_sf(ivar) = self%fabm_model%surface_state_variables(ivar)%maximum
         end do
      end if
      
      ! Allocate interior tracer source terms [var_unit s-1]
      if (state%n_fabm_interior_state > 0) then
         allocate(self%sms_int(grid%nz_grid, state%n_fabm_interior_state))
      end if

      ! Allocate fluxes over pelagic-benthic interface [var_unit m s-1] and bottom tracer source terms [var_unit s-1]
      if (state%n_fabm_interior_state > 0 .or. state%n_fabm_bottom_state > 0) then
         allocate(self%flux_bt(kmax_bot, state%n_fabm_interior_state))
         allocate(self%sms_bt(kmax_bot, state%n_fabm_bottom_state))
      end if

      ! Allocate fluxes over air-water surface [var_unit m s-1] and surface tracer source terms [var_unit s-1]
      if (state%n_fabm_interior_state > 0 .or. state%n_fabm_surface_state > 0) then
         allocate(self%flux_sf(state%n_fabm_interior_state))
         allocate(self%sms_sf(state%n_fabm_surface_state))
      end if

      ! Allocate local vertical velocities of pelagic state variables (sinking, floating, active movement) [m s-1]
      ! Movement through water, independent of water flow
      if (state%n_fabm_interior_state > 0) then
         allocate(self%velocity(grid%nz_grid, state%n_fabm_interior_state))
      end if

      ! Point FABM to fields that contain values for environmental data, all variables are assumed to be allocated
      ! Do this for all variables on FABM's standard variable list that the model can provide
      ! For this list, visit https://github.com/fabm-model/fabm/wiki/List-of-standard-variables
      ! Listed below are all non biogeochemical variables from that list
      
      ! Interior variables (arrays)
      ! Thickness of each layer [m]
      call self%fabm_model%link_interior_data(fabm_standard_variables%cell_thickness, grid%h(1:grid%nz_grid))
      ! Density of each layer [kg m-3]
      call self%fabm_model%link_interior_data(fabm_standard_variables%density, state%rho)
      ! Depth relative to surface [m]
      call self%fabm_model%link_interior_data(fabm_standard_variables%depth, grid%layer_depth)
      ! Shortwave radiative (SWR) flux [W m-2]
      call self%fabm_model%link_interior_data(fabm_standard_variables%downwelling_shortwave_flux, state%swr_vol)
      ! Photosynthetically active radiative (PAR) flux, a fraction of SWR flux [W m-2]
      call self%fabm_model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux, state%par_vol)
      ! Salinity at each layer [1e-3]
      call self%fabm_model%link_interior_data(fabm_standard_variables%practical_salinity, state%S)
      ! Pressure at each layer [dbar]: equal to layer depth, assuming one meter in depth is one dbar in pressure
      call self%fabm_model%link_interior_data(fabm_standard_variables%pressure, grid%layer_depth)
      ! Temperature at each layer [°C]
      call self%fabm_model%link_interior_data(fabm_standard_variables%temperature, state%T)
      ! Attenuation coefficient of SWR flux [m-1] and of PAR flux [m-1]
      ! If bioshade feedback is on, link the Simstrat contribution (background extinction)
      ! If bioshade feedback is off, Simstrat calculation of attenuation coefficient overwrites bgc model calculation
      ! Note: Passing the same attenuation coefficient for SWR and PAR is not exactly correct 
      ! since the attenuation coefficient for PAR can be higher than for SWR
      ! see for example: https://www.sciencedirect.com/science/article/pii/S1001074217317734
      if (fabm_cfg%bioshade_feedback) then
         call self%fabm_model%link_interior_data('attenuation_coefficient_host', state%background_extinction_vol)
      else
         call self%fabm_model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_shortwave_flux, state%absorb_vol, source=data_source_user)
         call self%fabm_model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, state%absorb_vol, source=data_source_user)
      end if
      ! Secchi depth [m], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%secchi_depth)
      ! Net rate of SWR energy absorption at each layer [W m-2], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%net_rate_of_absorption_of_shortwave_energy_in_layer)
      ! Projection factor for benthic flux into horizontal layer volume [m-1]
      self%id_bot_pel_conv = self%fabm_model%get_interior_variable_id(self%bot_pel_conv)
      call self%fabm_model%link_interior_data(self%id_bot_pel_conv, grid%dAz_norm)
      ! Vertical tracer diffusity [m2 s-1]: defined in GOTM, not defined in Simstrat
      ! Declaration would be in Simstrat_FABM type
      !type(type_interior_standard_variable) :: vertical_tracer_diffusivity = type_interior_standard_variable(name='vertical_tracer_diffusivity', units='m2 s-1')
      !call self%fabm_model%link_interior_data(self%vertical_tracer_diffusivity)

      ! Horizontal variables (scalars in a 1D model)
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
      ! Surface albedo [-], aggregate variable overwritten by Simstrat -> remove source=data_source_user to calculate by bgc models
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_albedo, state%wat_albedo, source=data_source_user)
      ! PAR flux at surface [W m-2]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux, state%par0)
      ! SWR flux at surface [W m-2]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, state%rad0)
      ! Surface drag coefficient [-], aggregate variable overwritten by Simstrat -> remove source=data_source_user to calculate by bgc models
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_drag_coefficient_in_air, state%C10, source=data_source_user)
      ! Surface specific humidity [-]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_specific_humidity, state%qa)
      ! Surface temperature [°C]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_temperature, state%T_atm)
      ! Total wind speed [m s-1]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%wind_speed, state%uv10)
      ! latitude [degree_north]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%latitude, state%Lat)
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

      ! Global variables (scalars)
      ! Number of days since start of the year [days]
      call self%fabm_model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, state%current_day_of_year)
      
      ! If bioshade feedback is on, set background_extinction_vol from FABMConfig and
      ! get index of attenuation coefficient FABM diagnostic variable
      if (fabm_cfg%bioshade_feedback) then        
         state%background_extinction_vol(:) = fabm_cfg%background_extinction
         att_index = 0
         do ivar = 1, size(self%fabm_model%interior_diagnostic_variables)
            if (self%fabm_model%interior_diagnostic_variables(ivar)%name == 'attenuation_coefficient_of_photosynthetic_radiative_flux') then
               self%fabm_model%interior_diagnostic_variables(ivar)%save = .true.
               att_index = ivar
            end if
         end do
         if (att_index == 0) then
            call error('Attenuation Coefficient not found in FABM Diagnostic Variables. Bioshade feedback not possible.')
         end if
      end if

      ! Write all diagnostic and aggregate variables to list_diagnostic
      ! Allocate arrays for diagnostic variables in SetDiagnosticVars and retrieve their index in the list of all interior and horizontal diagnostic variables
      if (fabm_cfg%output_diag_vars) then
         call list_diagnostic(self, fabm_cfg)
         call set_fabm_diagnostic_vars(self, state, fabm_cfg)
         if (state%n_fabm_diagnostic_interior > 0) then
            allocate(state%fabm_diagnostic_interior(grid%nz_grid, state%n_fabm_diagnostic_interior))
         end if
         if (state%n_fabm_diagnostic_horizontal > 0) then
            diagnostic_horizontal_exist = .true.
            allocate(state%fabm_diagnostic_horizontal(kmax_bot, state%n_fabm_diagnostic_horizontal))
         end if
      else
         call warn('Set OutputDiagnosticVars in FABMConfig to true to output FABM diagnostic variables.')
      end if

      ! Allocate arrays for repaired variables in list_repaired.dat and initialize them with the boundary value
      if (fabm_cfg%output_repaired_vars) then
         call set_fabm_repaired_vars(state, fabm_cfg)
         if (state%n_fabm_repaired_interior_min + state%n_fabm_repaired_interior_max > 0) then
            allocate(state%fabm_repaired_interior(grid%nz_grid, state%n_fabm_repaired_interior_min + state%n_fabm_repaired_interior_max))
            state%fabm_repaired_interior = ieee_value(state%fabm_repaired_interior, ieee_quiet_nan)
            do ivar = 1, state%n_fabm_interior_state 
               index = findloc(state%fabm_repaired_names, trim(self%fabm_model%interior_state_variables(ivar)%name)//'_minimum', dim = 1)
               if (index > 0) then
                  state%fabm_repaired_interior(:, index) = state%min_int(ivar)
                  continue
               end if
               index = findloc(state%fabm_repaired_names, trim(self%fabm_model%interior_state_variables(ivar)%name)//'_maximum', dim = 1)
               if (index > 0) then
                  state%fabm_repaired_interior(:, index) = state%max_int(ivar)
               end if
            end do
            if (any(ieee_is_nan(state%fabm_repaired_interior))) call error('Subroutine set_fabm_repaired_vars failed')
         end if
         if (state%n_fabm_repaired_bottom_min + state%n_fabm_repaired_bottom_max > 0) then
            allocate(state%fabm_repaired_bottom(kmax_bot, state%n_fabm_repaired_bottom_min + state%n_fabm_repaired_bottom_max))
            state%fabm_repaired_bottom = ieee_value(state%fabm_repaired_bottom, ieee_quiet_nan)
            do ivar = 1, state%n_fabm_bottom_state 
               index = findloc(state%fabm_repaired_names, trim(self%fabm_model%bottom_state_variables(ivar)%name)//'_minimum', dim = 1)
               if (index > 0) then
                  state%fabm_repaired_bottom(:, index - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min) = state%min_bt(ivar)
                  continue
               end if 
               index = findloc(state%fabm_repaired_names, trim(self%fabm_model%bottom_state_variables(ivar)%name)//'_maximum', dim = 1)
               if (index > 0) then
                  state%fabm_repaired_bottom(:, index - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min) = state%max_bt(ivar)
               end if
            end do
            if (any(ieee_is_nan(state%fabm_repaired_interior))) call error('Subroutine set_fabm_repaired_vars failed')
         end if
         if (state%n_fabm_repaired_surface_min + state%n_fabm_repaired_surface_max > 0) then
            allocate(state%fabm_repaired_surface(state%n_fabm_repaired_surface_min + state%n_fabm_repaired_surface_max))
            state%fabm_repaired_surface = ieee_value(state%fabm_repaired_surface, ieee_quiet_nan)
            do ivar = 1, state%n_fabm_surface_state 
               index = findloc(state%fabm_repaired_names, trim(self%fabm_model%surface_state_variables(ivar)%name)//'_minimum', dim = 1)
               if (index > 0) then
                  state%fabm_repaired_surface(index - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min - state%n_fabm_repaired_bottom_max - state%n_fabm_repaired_bottom_min) = state%min_sf(ivar)
                  continue
               end if 
               index = findloc(state%fabm_repaired_names, trim(self%fabm_model%surface_state_variables(ivar)%name)//'_maximum', dim = 1)
               if (index > 0) then
                  state%fabm_repaired_surface(index - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min - state%n_fabm_repaired_bottom_max - state%n_fabm_repaired_bottom_min) = state%max_sf(ivar)
               end if
            end do
            if (any(ieee_is_nan(state%fabm_repaired_interior))) call error('Subroutine set_fabm_repaired_vars failed')
         end if
      else
         call warn('FABM repaired variables not registered. Set OutputRepairedVars in FABMConfig to true to do so.')
      end if
      
      ! Complete initialization and check whether FABM has all dependencies fulfilled
      ! (i.e., whether all required calls to fabm_model%link_*_data have been made and all required data have been provided)
      ! Stop with fatal error if not
      ! Selection of diagnostics that FABM will compute and store becomes frozen
      call self%fabm_model%start()

      ! Set FABM-provided initial values for state variables (tracers), typically space-independent.
      ! This sets the values of arrays sent to fabm_model%link_*_state_data,
      ! in this case those contained in *_state
      ! If model is not initialized with custom or previously stored state
      if (.not. sim_cfg%continue_from_snapshot) then
         ! Interior state variables
         call self%fabm_model%initialize_interior_state(1, grid%nz_grid)
         ! Bottom state variables
         call self%fabm_model%initialize_bottom_state()
         ! If bottom_everywhere is set, at every depth:
         ! FABM is pointed to location that holds state data for the current depth
         ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
         ! The bottom state at the current depth is initialized
         if (fabm_cfg%bottom_everywhere) then
            do k_bot = 2, kmax_bot
               do ivar = 1, state%n_fabm_bottom_state
                  call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
               end do
               call self%fabm_model%initialize_bottom_state()
            end do
            ! Reset Botom to 1
            do ivar = 1, state%n_fabm_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
            end do
            k_bot = 1
         end if
         ! Surface state variables
         call self%fabm_model%initialize_surface_state()
      end if
   end subroutine init

   ! The update function is called in the main loop of simstrat (in simstrat.f90) at every time step
   ! Particle atmospheric, pelagic and benthic fluxes and diffusion are computed to update bgc state variable values
   subroutine update(self, state, fabm_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      integer :: ivar, ivar_diag, index, k

      ! If it is not the first call: 1. calculate and 2. validate the new model state, stop if variables are not valid
      if (.not. state%first_timestep) then
         ! 1a. Time-integrate the advection-diffusion-reaction equations
         ! of all tracers, combining the Simstrat transport terms with the FABM biogeochemical source
         ! terms and fluxes (sms, flux) and vertical velocities (velocity). This results in an updated interior_state.
         do ivar = 1, state%n_fabm_interior_state
            ! Special case for WET pelagic mirror variables (no diffusion)
            if (state%fabm_state_names(ivar)(len_trim(state%fabm_state_names(ivar))-2:) == '_PV') then
               if (fabm_cfg%bottom_everywhere) then
                  state%fabm_interior_state(:, ivar) = state%fabm_interior_state(:, ivar) + state%dt * self%flux_bt(:, ivar) * grid%dAz_norm(:)
               else
                  state%fabm_interior_state(1, ivar) = state%fabm_interior_state(1, ivar) + state%dt * self%flux_bt(1, ivar) * grid%dAz_norm(1)
               end if   
            else
               call diffusion_FABM_interior_state(self, state, fabm_cfg, grid, ivar)
            end if
         end do

         ! 1b. Direct time integration of source terms to update bottom_state
         do ivar = 1, state%n_fabm_bottom_state
            state%fabm_bottom_state(:, ivar) = state%fabm_bottom_state(:, ivar) + state%dt * self%sms_bt(:, ivar)
         end do

         ! 1c. Direct time integration of source terms to update surface_state
         do ivar = 1, state%n_fabm_surface_state
            state%fabm_surface_state(ivar) = state%fabm_surface_state(ivar) + state%dt * self%sms_sf(ivar)
         end do

         ! 2a. Simstrat check interior FABM bounds
         index = 1
         do ivar = 1, state%n_fabm_interior_state 
            if (any(state%fabm_interior_state(:, ivar) < state%min_int(ivar))) then
               if (fabm_cfg%output_repaired_vars) then
                  ! Register repaired variable and store out-of-bound value
                  index = findloc(state%fabm_repaired_names, trim(self%fabm_model%interior_state_variables(ivar)%name)//'_minimum', dim = 1)
                  if (index > 0) then
                     do k = 1, grid%nz_grid
                        if (state%fabm_interior_state(k, ivar) < state%min_int(ivar)) then
                           state%fabm_repaired_interior(k, index) = state%fabm_interior_state(k, ivar)
                        else
                           state%fabm_repaired_interior(k, index) = state%min_int(ivar)
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%interior_state_variables(ivar)%name, 'minimum', state%min_int(ivar))
                  end if
               end if
            end if  
            if (any(state%fabm_interior_state(:, ivar) > state%max_int(ivar))) then
               if (fabm_cfg%output_repaired_vars) then
                  ! Register repaired variable and store out-of-bound value
                  index = findloc(state%fabm_repaired_names, trim(self%fabm_model%interior_state_variables(ivar)%name)//'_maximum', dim = 1)
                  if (index > 0) then
                     do k = 1, grid%nz_grid
                        if (state%fabm_interior_state(k, ivar) > state%max_int(ivar)) then
                           state%fabm_repaired_interior(k, index) = state%fabm_interior_state(k, ivar)
                        else
                           state%fabm_repaired_interior(k, index) = state%min_int(ivar)
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%interior_state_variables(ivar)%name, 'maximum', state%max_int(ivar))
                  end if
               end if
            end if
         end do

         ! 2a. FABM check (and repair) interior state variables
         call self%fabm_model%check_interior_state(1, grid%nz_grid, fabm_cfg%repair_states, self%valid_int)

         ! 2a. Simstrat check for negative FABM interior state variables
         do ivar = 1, state%n_fabm_interior_state  
            if (any(state%fabm_interior_state(:, ivar) < 0.0_RK)) then
               do k = 1, grid%nz_grid
                  if (state%fabm_interior_state(k, ivar) < 0.0_RK) then
                     write (6, *) 'FABM Interior Variable value is ', state%fabm_interior_state(k, ivar)
                     write (6, *) 'at grid point ', k
                     write (6, *) 'at time (days, seconds) = ', state%simulation_time
                  end if
               end do
               call error('FABM Variable '//trim(self%fabm_model%interior_state_variables(ivar)%name)//' below zero.')
            end if
         end do

         ! 2b. Simstrat check bottom FABM bounds
         do ivar = 1, state%n_fabm_bottom_state  
            if (any(state%fabm_bottom_state(:, ivar) < state%min_bt(ivar))) then
               if (fabm_cfg%output_repaired_vars) then
                  ! Register repaired variable and store out-of-bound value
                  index = findloc(state%fabm_repaired_names, trim(self%fabm_model%bottom_state_variables(ivar)%name)//'_minimum', dim = 1)
                  if (index > 0) then
                     index = index - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min
                     do k = 1, kmax_bot
                        if (state%fabm_bottom_state(k, ivar) < state%min_bt(ivar)) then
                           state%fabm_repaired_bottom(k, index) = state%fabm_bottom_state(k, ivar)
                        else
                           state%fabm_repaired_bottom(k, index) = state%min_bt(ivar)
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%bottom_state_variables(ivar)%name, 'minimum', state%min_bt(ivar))
                  end if
               end if
            end if  
            if (any(state%fabm_bottom_state(:, ivar) > state%max_bt(ivar))) then
               if (fabm_cfg%output_repaired_vars) then
                  ! Register repaired variable and store out-of-bound value
                  index = findloc(state%fabm_repaired_names, trim(self%fabm_model%bottom_state_variables(ivar)%name)//'_maximum', dim = 1)
                  if (index > 0) then
                     index = index - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min
                     do k = 1, kmax_bot
                        if (state%fabm_bottom_state(k, ivar) > state%max_bt(ivar)) then
                           state%fabm_repaired_bottom(k, index) = state%fabm_bottom_state(k, ivar)
                        else
                           state%fabm_repaired_bottom(k, index) = state%max_bt(ivar)
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%bottom_state_variables(ivar)%name, 'maximum', state%max_bt(ivar))
                  end if
               end if
            end if
         end do
         
         ! 2b. FABM check (and repair) bottom state variables
         call self%fabm_model%check_bottom_state(fabm_cfg%repair_states, self%valid_bt)
         if (fabm_cfg%bottom_everywhere) then
            ! If bottom_everywhere is set, at every depth:
            ! FABM is pointed to location that holds state data for the current depth
            ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
            ! The bottom state at the current depth is validated
            do k_bot = 2, kmax_bot
               do ivar = 1, state%n_fabm_bottom_state
                  call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
               end do
               call self%fabm_model%check_bottom_state(fabm_cfg%repair_states, self%valid_bt)
            end do
            ! Reset Bottom to 1
            do ivar = 1, state%n_fabm_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
            end do
            k_bot = 1
         end if

         ! 2b. Simstrat check for negative FABM bottom state variables
         do ivar = 1, state%n_fabm_bottom_state
            if (any(state%fabm_bottom_state(:, ivar) < 0.0_RK)) then
               do k = 1, kmax_bot
                  if (state%fabm_bottom_state(k, ivar) < 0.0_RK) then
                     write (6, *) 'FABM Bottom Variable value is ', state%fabm_bottom_state(k, ivar)
                     write (6, *) 'at grid point ', k
                     write (6, *) 'at time (days, seconds) = ', state%simulation_time
                  end if
                  call error('FABM Bottom Variable '//trim(self%fabm_model%bottom_state_variables(ivar)%name)//' below zero.')
               end do
            end if
         end do

         ! 2c. Simstrat check surface FABM bounds
         do ivar = 1, state%n_fabm_surface_state  
            if (state%fabm_surface_state(ivar) < state%min_sf(ivar)) then
               if (fabm_cfg%output_repaired_vars) then
                  ! Register repaired variable and store out-of-bound value
                  index = findloc(state%fabm_repaired_names, trim(self%fabm_model%surface_state_variables(ivar)%name)//'_minimum', dim = 1)
                  if (index > 0) then
                     index = index - state%n_fabm_repaired_bottom_max - state%n_fabm_repaired_bottom_min - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min
                     state%fabm_repaired_surface(index) = state%min_sf(ivar)
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%surface_state_variables(ivar)%name, 'minimum', state%min_sf(ivar))
                  end if
               end if
            end if  
            if (state%fabm_surface_state(ivar) > state%max_sf(ivar)) then
               if (fabm_cfg%output_repaired_vars) then
                  ! Register repaired variable and store out-of-bound value
                  index = findloc(state%fabm_repaired_names, trim(self%fabm_model%surface_state_variables(ivar)%name)//'_maximum', dim = 1)
                  if (index > 0) then
                     index = index - state%n_fabm_repaired_bottom_max - state%n_fabm_repaired_bottom_min - state%n_fabm_repaired_interior_max - state%n_fabm_repaired_interior_min
                     state%fabm_repaired_surface(index) = state%max_sf(ivar)
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%surface_state_variables(ivar)%name, 'maximum', state%max_sf(ivar))
                  end if
               end if
            end if
         end do

         ! 2c. FABM check (and repair) surface state variables
         call self%fabm_model%check_surface_state(fabm_cfg%repair_states, self%valid_sf)

         ! 2c. Simstrat check for negative FABM surface state variables
         do ivar = 1, state%n_fabm_surface_state
            if (state%fabm_surface_state(ivar) < 0.0_RK) then
               write (6, *) 'FABM Surface Variable value is ', state%fabm_surface_state(ivar)
               write (6, *) 'at time (days, seconds) = ', state%simulation_time
               call error('FABM Surface Variable '//trim(self%fabm_model%surface_state_variables(ivar)%name)//' below zero.')
            end if
         end do

         ! 2. Error if FABM out of bounds and not repaired
         if (.not. (self%valid_int .and. self%valid_bt .and. self%valid_sf)) then
            if (.not. fabm_cfg%repair_states) then
               call error('FABM Variables out of bounds')
            end if
         end if
      end if

      ! Initialize (at first call) or update fluxes, sources and vertical movement
      if (state%first_timestep) then
         ivar = 2
         write (6, *) 'Initializing FABM calculations...'
      end if

      ! 1. Prepare all fields (e.g. light attenuation) FABM needs to compute fluxes and source terms
      ! Operates on entire active spatial domain
      ! To enable FABM's built-in time filters provide argument t that describes the model time
      call self%fabm_model%prepare_inputs(real(state%simulation_time(2), kind=RK))

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
            call error('FABM Interior Source contains NaNs, set to 0')
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
         ! Set NaNs to 0
         if (state%n_fabm_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_bt(1, :)))) then
               call warn('FABM Bottom Flux contains NaN, set to 0')
               where (ieee_is_nan(self%flux_bt(1, :)))
                  self%flux_bt(1, :) = 0.0
               end where
            end if
         end if
         if (state%n_fabm_bottom_state > 0) then
            if (any(ieee_is_nan(self%sms_bt(1, :)))) then
               call warn('FABM Bottom Source contains NaN, set to 0')
               where (ieee_is_nan(self%sms_bt(1, :)))
                  self%sms_bt(1, :) = 0.0
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
               call warn('FABM Surface Flux contains NaN, set to 0')
               where (ieee_is_nan(self%flux_sf))
                  self%flux_sf = 0.0
               end where
            end if
         end if
         if (state%n_fabm_surface_state > 0) then
            if (any(ieee_is_nan(self%sms_sf))) then
               call warn('FABM Surface Source contains NaN, set to 0')
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
         if (.not. state%first_timestep) self%velocity = 0.0
         call self%fabm_model%get_vertical_movement(1, grid%nz_grid, self%velocity)
         ! Set NaNs to 0
         if (any(ieee_is_nan(self%velocity))) then
            call warn('FABM Interior Velocity contains NaN, set to 0')
            where (ieee_is_nan(self%velocity))
               self%velocity = 0.0
            end where
         end if
      end if

      ! 4. Compute any remaining diagnostics not computed by preceding routines
      ! Operates on entire active spatial domain
      call self%fabm_model%finalize_outputs()

      ! 5. Retrieve values of diagnostic variables
      if (fabm_cfg%output_diag_vars) then
         if (state%n_fabm_diagnostic_interior > 0) then
            do ivar_diag = 1, state%n_fabm_diagnostic_interior
               index = state%diagnostic_index(ivar_diag)
               state%fabm_diagnostic_interior(:, ivar_diag) = self%fabm_model%get_interior_diagnostic_data(index)
            end do
         end if
         if (state%n_fabm_diagnostic_horizontal > 0) then
            do ivar_diag = 1, state%n_fabm_diagnostic_horizontal
               index = state%diagnostic_index(state%n_fabm_diagnostic_interior + ivar_diag)
               state%fabm_diagnostic_horizontal(1, ivar_diag) = self%fabm_model%get_horizontal_diagnostic_data(index)
            end do
         end if
      end if

      ! If bottom_everywhere is set, repeat relevant steps from 1. to 5. at every depth
      if (fabm_cfg%bottom_everywhere .and. (interior_or_bottom_exist .or. diagnostic_horizontal_exist)) then
         ! FABM is pointed to location that holds state data for the current depth
         ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
         ! Environmental data (dAz_norm) is updated
         ! The bottom fluxes and sources at the current depth are calculated
         ! Remaining diagnostics are calculated at the current depth
         ! The horizontal diagnostic variables at the current depth are retrieved
         ! Note: Since FABM does not distinguish between bottom and surface diagnostic variables, 
         !       surface diagnostic variables are also retrieved for every depth, with constant value
         do k_bot = 2, kmax_bot
            do ivar = 1, state%n_fabm_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
            end do
            call self%fabm_model%prepare_inputs(real(state%simulation_time(2), kind=RK))
            call self%fabm_model%get_bottom_sources(self%flux_bt(k_bot, :), self%sms_bt(k_bot, :))
            call self%fabm_model%finalize_outputs()
            if (fabm_cfg%output_diag_vars) then
               do ivar_diag = 1, state%n_fabm_diagnostic_horizontal
                  index = state%diagnostic_index(state%n_fabm_diagnostic_interior + ivar_diag)
                  state%fabm_diagnostic_horizontal(k_bot, ivar_diag) = self%fabm_model%get_horizontal_diagnostic_data(index)
               end do
            end if
         end do
         ! Set NaNs in bottom fluxes and sources to 0
         if (state%n_fabm_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_bt))) then
               call warn('FABM Bottom Flux contains NaN, set to 0')
               where (ieee_is_nan(self%flux_bt))
                  self%flux_bt = 0.0
               end where
            end if
         end if
         if (state%n_fabm_bottom_state > 0) then
            if (any(ieee_is_nan(self%sms_bt))) then
               call warn('FABM Bottom Source contains NaN, set to 0')
               where (ieee_is_nan(self%sms_bt))
                  self%sms_bt = 0.0
               end where
            end if
         end if
         ! Reset Bottom to 1
         do ivar = 1, state%n_fabm_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
         end do
         k_bot = 1
         call self%fabm_model%prepare_inputs(real(state%simulation_time(2), kind=RK))
         call self%fabm_model%finalize_outputs()
      end if

      ! This is the ideal moment for the output: 
      ! All variables (enviromental, state, diagnostics, source, fluxes, vertical velocities) are in sync
   end subroutine update

   ! Diffusion algorithm for interior state variables: Simstrat transport terms and FABM biogeochemical terms integrated simultaneously
   ! Assuming small enough dt such that only fluxes between neighbouring layers are relevant
   ! Also assuming that sides of the lake slope inward
   subroutine diffusion_fabm_interior_state(self, state, fabm_cfg, grid, ivar)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(StaggeredGrid), intent(in) :: grid
      integer, intent(in) :: ivar

      ! Local variables
      real(RK), dimension(grid%nz_grid) :: velocity_up, velocity_down, sources, lower_diag_diff, main_diag_diff, upper_diag_diff, lower_diag, main_diag, upper_diag, rhs, AreaFactor_ext_1, AreaFactor_ext_2
      integer :: k

      ! Create linear system of equations with transport and source terms
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}

      ! Build diagonals of A with Simstrat transport terms (code from discretization module)
      ! With diffusivity for temperature (nuh)
      ! Upward flux from (1:grid%nz_grid - 1) to (2:grid%nz_grid): upper_diag_diff(z) describes the upward flux into cell z
      upper_diag_diff(2:grid%nz_grid) = state%dt*state%nuh(2:grid%nz_grid)*grid%AreaFactor_1(2:grid%nz_grid)
      upper_diag_diff(1) = 0.0_RK ! No upward flux into bottommost cell
      ! Downward flux from (2:grid%nz_grid) to (1:grid%nz_grid - 1): lower_diag_diff(z) describes the downward flux into cell z
      lower_diag_diff(1:grid%nz_grid-1) = state%dt*state%nuh(2:grid%nz_grid)*grid%AreaFactor_2(1:grid%nz_grid - 1)
      lower_diag_diff(grid%nz_grid) = 0.0_RK ! No downward flux into uppermost cell
      ! 1 - downward flux out - upward flux out
      main_diag_diff(1:grid%nz_grid) = 1.0_RK - upper_diag_diff(1:grid%nz_grid) - lower_diag_diff(1:grid%nz_grid)

      ! Add FABM fluxes to A as residual verical advection terms
      ! Benthic-pelagic and air-water fluxes are added as source terms below
      ! Area factors for external fluxes, minus sign to be closer to the AreaFactors from grid
      ! AreaFactor_ext_1 is for upward fluxes:
      ! Multiply by the layer through which the flux passes (z = i) and divide by the volume of the receiving cell (z = i)
      AreaFactor_ext_1(1:grid%nz_grid) = -grid%Az(1:grid%nz_grid) / (grid%h(1:grid%nz_grid) * grid%Az_vol(1:grid%nz_grid))
      ! AreaFactor_ext_2 is for downward fluxes: 
      ! Multiply by the layer through which the flux passes (z = i+1) and divide by the volume of the receiving cell (z = i)
      AreaFactor_ext_2(1:grid%nz_grid) = -grid%Az(2:grid%nz_grid+1) / (grid%h(1:grid%nz_grid) * grid%Az_vol(1:grid%nz_grid))
      ! Upward and downward movement rates for each layer
      ! if the variable moves upward velocity_up is positive and velocity_down zero
      ! vice versa if the variable moves downward
      do k = 1, grid%nz_grid
         if (self%velocity(k, ivar) >= 0) then
            velocity_up(k) = self%velocity(k, ivar)
            velocity_down(k) = 0
         else
            velocity_up(k) = 0
            velocity_down(k) = -self%velocity(k, ivar)
         end if
      end do
      ! Add upward flux from (1:grid%nz_grid-1) to (2:grid%nz_grid): upper_diag(z) describes the upward flux into cell z
      ! upper_diag should be <= 0
      ! Benthic-pelagic flux is treated as source term below
      upper_diag(2:grid%nz_grid) = upper_diag_diff(2:grid%nz_grid) + state%dt*velocity_up(1:grid%nz_grid-1)*AreaFactor_ext_1(2:grid%nz_grid)
      upper_diag(1) = upper_diag_diff(1)
      ! Subtract upward flux from cell z to cell z+1
      ! main_diag should be >= 0
      ! Water-air flux is treated as sink term below
      main_diag(1:grid%nz_grid-1) = main_diag_diff(1:grid%nz_grid-1) - state%dt*velocity_up(1:grid%nz_grid-1)*AreaFactor_ext_2(1:grid%nz_grid-1)
      main_diag(grid%nz_grid) = main_diag_diff(grid%nz_grid)
      ! Add downward flux from (2:grid%nz_grid) to (1:grid%nz_grid - 1): lower_diag(z) describes the downward flux into cell z
      ! lower_diag should be <= 0
      ! Air-water flux is treated as source term below
      lower_diag(1:grid%nz_grid-1) = lower_diag_diff(1:grid%nz_grid-1) + state%dt*velocity_down(2:grid%nz_grid)*AreaFactor_ext_2(1:grid%nz_grid-1)
      lower_diag(grid%nz_grid) = lower_diag_diff(grid%nz_grid)
      ! Subtract downward flux from cell z to cell z-1
      ! main_diag should be >= 0
      ! Pelagic-benthic flux is treated as sink term below
      main_diag(2:grid%nz_grid) = main_diag(2:grid%nz_grid) - state%dt*velocity_down(2:grid%nz_grid)*AreaFactor_ext_1(2:grid%nz_grid)
      main_diag(1) = main_diag(1)

      ! Get source S^{n}
      ! Source at each layer
      sources = self%sms_int(:, ivar)
      ! Add pelagic-benthic and air-water flux [var_unit m s-1] as source [var_unit s-1]
      ! Convert bottom flux to source by multiplication by sediment area over layer volume (dAz_norm, [m])
      if (fabm_cfg%bottom_everywhere) then
         sources(:) = sources(:) + (self%flux_bt(:, ivar) * grid%dAz_norm(:)) ! pelagic-benthic flux at every layer
      else
         sources(1) = sources(1) + (self%flux_bt(1, ivar) * grid%dAz_norm(1)) ! pelagic-benthic flux at bottommost layer
      end if
      ! Convert surface flux to source by division by surface area over volume of uppermost layer [m]
      sources(grid%nz_grid) = sources(grid%nz_grid) + (self%flux_sf(ivar) * grid%Az(grid%nz_grid) / (grid%h(grid%nz_grid) * grid%Az_vol(grid%nz_grid)))

      ! Calculate RHS (phi^{n}+dt*S^{n})
      rhs(:) = state%fabm_interior_state(:, ivar) + state%dt*sources(:)

      ! Solve LES to get phi^{n+1}
      call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%fabm_interior_state(:, ivar), grid%nz_grid)
   end subroutine diffusion_fabm_interior_state

   ! Light absorption feedback by FABM variables
   subroutine absorption_update_fabm(self, state, fabm_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(in) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      real(RK), dimension(grid%nz_grid) :: attenuation_coefficient_of_photosynthetic_radiative_flux

      ! Retrieve attenuation_coefficient_of_photosynthetic_radiative_flux and set as absorb_vol
      if (state%first_timestep) then
         ! At the first timestep FABM variables have not yet been calculated -> set as background extinction from FABMConfig
         attenuation_coefficient_of_photosynthetic_radiative_flux(:) = fabm_cfg%background_extinction
      else
         attenuation_coefficient_of_photosynthetic_radiative_flux = self%fabm_model%get_interior_diagnostic_data(att_index)
      end if
      state%absorb_vol(:) = attenuation_coefficient_of_photosynthetic_radiative_flux(:)

      ! Interpolate to faces to be compatible with Simstrat temperature module
      call grid%interpolate_to_face(grid%z_volume, state%absorb_vol, grid%nz_grid, state%absorb)
   end subroutine absorption_update_fabm

   ! Output list of diagnostic variables
   subroutine list_diagnostic(self, fabm_cfg)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(FABMConfig), intent(in) :: fabm_cfg

      ! Local variables
      integer :: i, unit, status
      character(len=256) :: file_path

      ! Write the names of all interior diagnostic variables
      if (size(self%fabm_model%interior_diagnostic_variables) > 0) then
         ! Construct the file path
         file_path = trim(fabm_cfg%config_path)//'list_diagnostic_interior.dat'
         ! Create new file or overwrite already existing one
         open(newunit=unit, file=file_path, action='write', iostat = status)
         if (status .ne. 0) then
            call error('Failed to open or create file: ' // trim(file_path) // '. Please add FABM configurations folder.')
         end if
         ! Write the header
         write(unit, '(A)', advance = 'no') 'Short Name, '
         write(unit, '(A)', advance = 'no') 'Long Name, '
         write(unit, '(A)', advance = 'no') 'Units, '
         write(unit, '(A)') 'Output'
         ! Write the names and units of interior diagnostic variables
         do i = 1, size(self%fabm_model%interior_diagnostic_variables)
         write(unit, '(A)', advance='no') '"'//trim(self%fabm_model%interior_diagnostic_variables(i)%name)//'", '
         write(unit, '(A)', advance='no') '"'//trim(self%fabm_model%interior_diagnostic_variables(i)%long_name)//'", '
         write(unit, '(A)', advance='no') '"'//trim(self%fabm_model%interior_diagnostic_variables(i)%units)//'", '
         if (self%fabm_model%interior_diagnostic_variables(i)%output == 1) then
            write(unit, '(A)') '"Yes"'
         else
            write(unit, '(A)') '"No"'
         end if
         end do
         close(unit)
      end if

      ! Write the names of all horizontal diagnostic variables
      if (size(self%fabm_model%horizontal_diagnostic_variables) > 0) then
         ! Construct the file path for Horizontal Diagnostic Variable List
         file_path = trim(fabm_cfg%config_path)//'list_diagnostic_horizontal.dat'
         ! Creat new file or overwrite already existing one
         open(newunit=unit, file=file_path, action='write', iostat = status)
         ! Write the header
         write(unit, '(A)', advance = 'no') 'Short Name, '
         write(unit, '(A)', advance = 'no') 'Long Name, '
         write(unit, '(A)', advance = 'no') 'Units, '
         write(unit, '(A)') 'Output'
         ! Write the names and units of horizontal diagnostic variables
         do i = 1, size(self%fabm_model%horizontal_diagnostic_variables)
         write(unit, '(A)', advance='no') '"'//trim(self%fabm_model%horizontal_diagnostic_variables(i)%name)//'", '
         write(unit, '(A)', advance='no') '"'//trim(self%fabm_model%horizontal_diagnostic_variables(i)%long_name)//'", '
         write(unit, '(A)', advance='no') '"'//trim(self%fabm_model%horizontal_diagnostic_variables(i)%units)//'", '
         if (self%fabm_model%horizontal_diagnostic_variables(i)%output == 1) then
            write(unit, '(A)') '"Yes"'
         else
            write(unit, '(A)') '"No"'
         end if
         end do
         close(unit)
      end if
   end subroutine list_diagnostic

   ! Read names of diagnostic Vars in SetDiagnosticVars file
   subroutine set_fabm_diagnostic_vars(self, state, fabm_cfg)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg

      ! Local variables
      integer :: i, j, n, n_int, n_hor, unit, status, unit_list, status_list
      character(len=256) :: line, line_trim, short_name, short_name_trim, long_name, units, output, file_path
      integer, parameter :: max_lines = 10000 ! Maximum amount of diagnostic variables in SetDiagnosticVars file
      character(len=256), dimension(max_lines) :: temp_names ! Same len as fabm_diagnostic_names
      integer, dimension(max_lines) :: temp_index
      logical, dimension(max_lines) :: is_int

      ! Read from SetDiagVars
      open(newunit=unit, action='read', status='old', file=fabm_cfg%diag_vars, iostat = status)
      if (status .ne. 0) then
         call warn('No Output of FABM Diagnostic Variables')
         state%n_fabm_diagnostic = 0
         return
      else
         write(6,*) 'Reading ', trim(fabm_cfg%diag_vars)
      end if
      ! Skip header
      read(unit, '(A)', iostat=status)
      ! Initialize counts to 0
      n = 0
      n_int = 0
      n_hor = 0
      ! Read every diagnostic variable
      do
         read(unit, '(A)', iostat=status) line
         if (status .ne. 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:4) == '-Add') cycle
         if (line(1:4) == '----') cycle
         if (any(temp_names == trim(line))) cycle
         line_trim = trim(line)
         if ((line_trim(:11) == 'select_all_') .or. (line_trim(:14) == 'select_output_')) then
            ! Read from list_diagnostic_interior.dat
            file_path = trim(fabm_cfg%config_path)//'list_diagnostic_interior.dat'
            open(newunit=unit_list, action='read', status='old', file=file_path, iostat = status_list)
            read(unit_list, *, iostat=status_list)
            if (status_list .ne. 0) then
               call error('FABM Interior Diagnostic Variables have not been listed.')
            end if
            do
               read(unit_list, *, iostat=status_list) short_name, long_name, units, output
               if (status_list .ne. 0) exit
               if (len_trim(short_name) == 0) cycle
               short_name_trim = trim(short_name)
               if ((short_name_trim(:(len_trim(line)-11)) == line_trim(12:)) .or. (short_name_trim(:(len_trim(line)-14)) == line_trim(15:))) then
                  if (any(temp_names == trim(short_name))) cycle
                  if ((line_trim(:14) == 'select_output_') .and. output == 'No') cycle
                  n = n + 1
                  if (n > max_lines) then
                     call error('Too many lines in '//trim(fabm_cfg%diag_vars)//', increase max_lines in set_fabm_diagnostic_vars.')
                  end if
                  temp_names(n) = trim(short_name)
                  temp_index(n) = 0
               end if
            end do
            ! Close list_diagnostic_interior.dat
            close(unit_list)
            ! Read from list_diagnostic_horizontal.dat
            file_path = trim(fabm_cfg%config_path)//'list_diagnostic_horizontal.dat'
            open(newunit=unit_list, action='read', status='old', file=file_path, iostat = status_list)
            read(unit_list, *, iostat=status_list)
            if (status_list .ne. 0) then
               call error('FABM Horizontal Diagnostic Variables have not been listed.')
            end if
            do
               read(unit_list, *, iostat=status_list) short_name, long_name, units, output
               if (status_list .ne. 0) exit
               if (len_trim(short_name) == 0) cycle
               short_name_trim = trim(short_name)
               if ((short_name_trim(:(len_trim(line)-11)) == line_trim(12:)) .or. (short_name_trim(:(len_trim(line)-14)) == line_trim(15:))) then
                  if ((line_trim(:14) == 'select_output_') .and. output == 'No') cycle
                  if (any(temp_names == trim(short_name))) cycle
                  n = n + 1
                  if (n > max_lines) then
                     call error('Too many lines in '//trim(fabm_cfg%diag_vars)//', increase max_lines in set_fabm_diagnostic_vars.')
                  end if
                  temp_names(n) = trim(short_name)
                  temp_index(n) = 0
               end if
            end do
            ! Close list_diagnostic_horizontal.dat
            close(unit_list)
         else
            n = n + 1
            if (n > max_lines) then
               call error('Too many lines in '//trim(fabm_cfg%diag_vars)//', increase max_lines in set_fabm_diagnostic_vars.')
            end if
            temp_names(n) = trim(line)
            temp_index(n) = 0
         end if
      end do
      ! Close the file
      close(unit)

      ! Find diagnostic variable in FABM model
      do i = 1, n
         ! Set index at location in fabm_model%interior_diagnostic_variables
         do j = 1, size(self%fabm_model%interior_diagnostic_variables)
            if (self%fabm_model%interior_diagnostic_variables(j)%name == temp_names(i)) then
               self%fabm_model%interior_diagnostic_variables(j)%save = .true.
               temp_index(i) = j
               is_int(i) = .true.
               n_int = n_int + 1
            end if
         end do
         ! Set index at location in fabm_model%horizontal_diagnostic_variables
         do j = 1, size(self%fabm_model%horizontal_diagnostic_variables)
            if (self%fabm_model%horizontal_diagnostic_variables(j)%name == temp_names(i)) then
               self%fabm_model%horizontal_diagnostic_variables(j)%save = .true.
               temp_index(i) = j
               is_int(i) = .false.
               n_hor = n_hor + 1
            end if
         end do
         ! If not found in fabm_model
         if (temp_index(i) == 0) then
            call error('FABM diagnostic variable '//trim(temp_names(i))//' not found.')
         end if
      end do

      ! Allocate array of proper size
      state%n_fabm_diagnostic = n
      state%n_fabm_diagnostic_interior = n_int
      state%n_fabm_diagnostic_horizontal = n_hor
      allocate(state%fabm_diagnostic_names(n))
      allocate(state%diagnostic_index(n))

      ! Copy names and indices to array
      j = 1
      do i = 1, n
         if (is_int(i)) then
            state%fabm_diagnostic_names(j) = temp_names(i)
            state%diagnostic_index(j) = temp_index(i)
            j = j + 1
         end if
      end do
      do i = 1, n
         if (.not. is_int(i)) then
            state%fabm_diagnostic_names(j) = temp_names(i)
            state%diagnostic_index(j) = temp_index(i)
            j = j + 1
         end if
      end do
   end subroutine set_fabm_diagnostic_vars

   ! Add to list of repaired variables
   subroutine list_repaired(fabm_cfg, variable, boundary, boundary_value)
      ! Arguments
      class(FABMConfig), intent(in) :: fabm_cfg
      character(len=*), intent(in) :: variable
      character(len=*), intent(in) :: boundary
      real(RK), intent(in) :: boundary_value

      ! Local variables
      integer :: unit, status
      logical :: exists
      character(len=256) :: file_path, boundary_value_str
      character(len=256) :: line, new_line

      ! Construct the file path
      file_path = trim(fabm_cfg%config_path)//'list_repaired.dat'

      ! Convert boundary_value to string
      write(boundary_value_str, '(ES15.6)') boundary_value
      new_line = trim(variable)//', '//trim(boundary)//', '//trim(adjustl(boundary_value_str))

      ! Check if repaired variable case already present in RepairedVars
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      if (status .ne. 0) then
         call error('Failed to open or create file: ' // trim(file_path) // '. Please add FABM configurations folder.')
      end if
      ! Read the file line by line and check if the exact line exists
      rewind(unit)
      do
         read(unit, '(A)', iostat=status) line
         if (status .ne. 0) exit
         ! Check if the line already exists and exit subroutine if it does
         if (trim(line) == trim(new_line)) then
            close(unit)
            return
         end if
      end do
      close(unit)
      ! Reopen for appending
      open(newunit=unit, action='write', position='append', status='old', file=file_path, iostat=status)
      if (status .ne. 0) then
         call error('Failed to open or create file: ' // trim(file_path) // '. Please add FABM configurations folder.')
      else
         call warn('FABM variable '//trim(variable)//' added to '// trim(file_path)//'. Restart simulation to output '//trim(boundary)//' values.')
      end if

      ! Write variable name, its boundary and the boundary vale
      write(unit, '(A)') trim(new_line)

      ! Close the file after writing
      close(unit)
   end subroutine list_repaired

   ! Register repaired variables in RepairedVars file
   subroutine set_fabm_repaired_vars(state, fabm_cfg)
      ! Arguments
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg

      ! Local variables
      integer :: i, j, n, n_int_min, n_int_max, n_bt_min, n_bt_max, n_sf_min, n_sf_max, unit, status
      character(len=256) :: file_path
      character(len=256) :: name
      character(len=16) :: bound
      integer, parameter :: max_lines = 100 ! Maximum amount of repaired variables in RepairedVars file
      character(len=256), dimension(max_lines) :: temp_names
      character(len=16), dimension(max_lines) :: temp_bounds
      character(len=16), dimension(max_lines) :: temp_types

      ! Initialize counts to zero
      n = 0
      n_int_min = 0
      n_int_max = 0
      n_bt_min = 0
      n_bt_max = 0
      n_sf_min = 0
      n_sf_max = 0

      ! Construct the file path
      file_path = trim(fabm_cfg%config_path)//'list_repaired.dat'

      ! Read from RepairedVars
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      if (status .ne. 0) then
         call warn('No FABM Repaired Variables provided. Empty repaired variables file '// trim(file_path) //' created.')
         open(newunit=unit, action='write', status='new', file=file_path, iostat = status)
         if (status .ne. 0) then
            call error('Failed to open or create file: ' // trim(file_path) // '. Please add FABM configurations folder.')
         end if
         ! Write the header
         write(unit, '(A)', advance = 'no') 'Variable, '
         write(unit, '(A)', advance = 'no') 'Boundary reached, '
         write(unit, '(A)') 'Boundary value'
         close(unit)
         ! Write to state and allocate array of proper size
         state%n_fabm_repaired = n
         state%n_fabm_repaired_interior_min = n_int_min
         state%n_fabm_repaired_interior_max = n_int_max
         state%n_fabm_repaired_bottom_min = n_bt_min
         state%n_fabm_repaired_bottom_max = n_bt_max
         state%n_fabm_repaired_surface_min = n_sf_min
         state%n_fabm_repaired_surface_max = n_sf_max
         allocate(state%fabm_repaired_names(n))
         return
      else
         write(6,*) 'Reading ', trim(file_path)
      end if

      ! Skip header
      read(unit, *, iostat=status)

      ! Read every repaired variable and find in fabm_model
      do
         read(unit, *, iostat=status) name, bound
         if (status .ne. 0) exit
         if (len_trim(name) == 0) cycle  ! Skip empty lines, if any
         n = n + 1
         if (n > max_lines) then
            call error('Too many lines in '//trim(file_path)//', increase max_lines in set_fabm_repaired_vars')
         end if
         temp_names(n) = trim(name)
         temp_bounds(n) = trim(bound)
         temp_types(n) = 'undefined'
         ! Register as interior variable
         do j = 1, state%n_fabm_interior_state
            if (state%fabm_state_names(j) == temp_names(n)) then
               temp_types(n) = 'interior'
               if (temp_bounds(n) == 'minimum') then
                  n_int_min = n_int_min + 1
               else if  (temp_bounds(n) == 'maximum') then
                  n_int_max = n_int_max + 1
               else
                  call error('FABM repaired variable '//trim(temp_names(n))//' has invalid bound')
               end if
            end if
         end do
         ! Register as bottom variable
         do j = 1, state%n_fabm_bottom_state
            if (state%fabm_state_names(state%n_fabm_interior_state + j) == temp_names(n)) then
               if (temp_types(n) /= 'undefined') then
                  call error('Variable '//trim(temp_names(n))//' defined multiple times.')
               else
                  temp_types(n) = 'bottom'
               end if
               if (temp_bounds(n) == 'minimum') then
                  n_bt_min = n_bt_min + 1
               else if  (temp_bounds(n) == 'maximum') then
                  n_bt_max = n_bt_max + 1
               else
                  call error('FABM repaired variable '//trim(temp_names(n))//' has invalid bound')
               end if
            end if
         end do
         ! Register as surface variable
         do j = 1, state%n_fabm_surface_state
            if (state%fabm_state_names(state%n_fabm_interior_state + state%n_fabm_bottom_state + j) == temp_names(n)) then
               if (temp_types(n) /= 'undefined') then
                  call error('Variable '//trim(temp_names(n))//' defined multiple times.')
               else
                  temp_types(n) = 'surface'
               end if
               if (temp_bounds(n) == 'minimum') then
                  n_sf_min = n_sf_min + 1
               else if  (temp_bounds(n) == 'maximum') then
                  n_sf_max = n_sf_max + 1
               else
                  call error('FABM repaired variable '//trim(temp_names(n))//' has invalid bound')
               end if
            end if
         end do
         ! If not found in fabm_model
         if (temp_types(n) == 'undefined') then
            call error('FABM repaired variable '//trim(temp_names(n))//' not found. Control '//trim(file_path)//'.')
         end if
      end do
      ! Close the file
      close(unit)

      ! Add type and boundary to name
      do i = 1, n
         temp_names(i) = trim(temp_names(i))//'_'//trim(temp_bounds(i))
      end do

      ! Write to state and allocate array of proper size
      state%n_fabm_repaired = n
      state%n_fabm_repaired_interior_min = n_int_min
      state%n_fabm_repaired_interior_max = n_int_max
      state%n_fabm_repaired_bottom_min = n_bt_min
      state%n_fabm_repaired_bottom_max = n_bt_max
      state%n_fabm_repaired_surface_min = n_sf_min
      state%n_fabm_repaired_surface_max = n_sf_max
      allocate(state%fabm_repaired_names(n))

      ! Copy names to array for interior below minimum repaired variables
      j = 1
      do i = 1, n
         if (temp_types(i) == 'interior' .AND. temp_bounds(i) == 'minimum') then
            state%fabm_repaired_names(j) = temp_names(i)
            j = j + 1
         end if
      end do
      ! Copy names to array for interior above maximum repaired variables
      j = n_int_min + 1
      do i = 1, n
         if (temp_types(i) == 'interior' .AND. temp_bounds(i) == 'maximum') then
            state%fabm_repaired_names(j) = temp_names(i)
            j = j + 1
         end if
      end do
      ! Copy names to array for bottom below minimum repaired variables 
      j = n_int_min + n_int_max + 1
      do i = 1, n
         if (temp_types(i) == 'bottom' .AND. temp_bounds(i) == 'minimum') then
            state%fabm_repaired_names(j) = temp_names(i)
            j = j + 1
         end if
      end do
      ! Copy names to array for bottom above maximum repaired variables
      j = n_int_min + n_int_max + n_bt_min + 1
      do i = 1, n
         if (temp_types(i) == 'bottom' .AND. temp_bounds(i) == 'maximum') then
            state%fabm_repaired_names(j) = temp_names(i)
            j = j + 1
         end if
      end do
      ! Copy names to array for surface below minimum repaired variables
      j = n_int_min + n_int_max + n_bt_min + n_bt_max + 1
      do i = 1, n
         if (temp_types(i) == 'surface' .AND. temp_bounds(i) == 'minimum') then
            state%fabm_repaired_names(j) = temp_names(i)
            j = j + 1
         end if
      end do
      ! Copy names to array for surface above maximum repaired variables
      j = n_int_min + n_int_max + n_bt_min + n_bt_max + n_sf_min + 1
      do i = 1, n
         if (temp_types(i) == 'surface' .AND. temp_bounds(i) == 'maximum') then
            state%fabm_repaired_names(j) = temp_names(i)
            j = j + 1
         end if
      end do
   end subroutine set_fabm_repaired_vars

   ! Deallocate memory
   subroutine deallocate_fabm(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self

      ! Deallocate internal arrays
      if (allocated(self%sms_int)) deallocate(self%sms_int)
      if (allocated(self%flux_bt)) deallocate(self%flux_bt)
      if (allocated(self%sms_bt)) deallocate(self%sms_bt)
      if (allocated(self%flux_sf)) deallocate(self%flux_sf)
      if (allocated(self%sms_sf)) deallocate(self%sms_sf)
      if (allocated(self%velocity)) deallocate(self%velocity)
   end subroutine deallocate_fabm

   ! Finalize the coupling
   subroutine finalize_fabm(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      
      ! Finalize and deallocate the model
      if (associated(self%fabm_model)) then
         call self%fabm_model%finalize()
         deallocate(self%fabm_model)
         nullify(self%fabm_model)
         call fabm_finalize_library()
      end if

      ! Deallocate memory
      call deallocate_fabm(self)
      if (associated(self%bottom_index)) then
         deallocate(self%bottom_index)
         nullify(self%bottom_index)
      end if
   end subroutine finalize_fabm

end module simstrat_fabm