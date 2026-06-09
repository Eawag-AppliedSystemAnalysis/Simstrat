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
   ! Add types for additional FABM variables, used by specific bgc models
   ! Added from lib/fabm/src/fabm_types.F90:
   ! type_interior_standard_variable: Interior standard variable data type
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
   ! iterates over all layers if bottom_everywhere is true
   integer, target :: k_bot
   ! Bottom size: 1 if bottom_everywhere is false, nz_grid/nz_occupied if bottom_everywhere is true
   integer :: kmax_bot
   ! Validity of bottom variables at current bottom location
   logical :: valid_bot

   ! Index of attenuation_coefficient_of_photosynthetic_radiative_flux in FABM diagnostic variables
   integer :: att_index

   ! Timestep counter for FABM prepare_inputs function
   integer :: timestep_counter

   ! Necessary conditions for evaluation of bottom everywhere
   logical :: interior_or_bottom_exist = .false.
   logical :: diagnostic_horizontal_exist =.false.

   ! Definition of a manipulation action
   type, private :: Manipulation
      ! Variable name
      character(len=:), allocatable :: var_name
      ! Variable index
      integer :: var_index
      ! Start time (days)
      real(RK) :: start_time
      ! End time (days)
      ! Only relevant for action_type 2
      real(RK) :: end_time
      ! Action type (1: multiply, 2: add)
      integer :: action_type
      ! Action value
      ! Multiplied by timestep for action_type 2
      ! Divided by height for action_type 2 for depth varying variables
      real(RK) :: action_val
      ! Threshold: action_type 2 is not executed if all values are above
      real(RK) :: threshold
      ! Start depth for depth varying variables
      integer :: start_depth
      ! End depth for depth varying variables
      integer :: end_depth
   end type

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

      ! Arrays for interior state time-integration
      real(RK), dimension(:), allocatable :: upper_diag_disc, lower_diag_disc, main_diag_disc
      real(RK), dimension(:), allocatable :: AreaFactor_ext_1, AreaFactor_ext_2

      ! Variables for validity of state
      logical :: valid_int, valid_sf, valid_bt

      ! Index of current bottom location (pelagic-benthic interface)
      integer, pointer :: bottom_index => null()

      ! Arrays for diagnostic variables
      real(RK), dimension(:,:), pointer :: diagnostic_interior => null(), diagnostic_horizontal => null()
      ! Indices of diagnostic variables inside FABM
      integer, dimension(:), allocatable :: diagnostic_interior_index, diagnostic_horizontal_index

      ! Arrays for repaired variables
      real(RK), dimension(:,:), pointer :: repaired_interior => null(), repaired_bottom => null()
      real(RK), dimension(:), pointer :: repaired_surface => null()
      ! Names of repaired variables
      character(len=256), allocatable :: repaired_interior_names(:), repaired_bottom_names(:), repaired_surface_names(:)

      ! Manipulations of current state
      class(Manipulation), dimension(:), allocatable :: manipulations
      
      ! Define additional FABM variables, used by bgc models
      ! Simstrat contribution to attenuation_coefficient_of_shortwave_flux and attenuation_coefficient_of_photosynthetic_radiative_flux
      type(type_dependency_id) :: id_attenuation_coefficient_host
      ! Projection factor for benthic flux into horizontal layer volume
      type(type_interior_standard_variable) :: bot_pel_conv = type_interior_standard_variable(name='bot_pel_conv', units='-')
      type(type_fabm_interior_variable_id) :: id_bot_pel_conv
   contains
      ! Initialization and update
      procedure, pass :: init => fabm_init
      procedure, pass :: update => fabm_update

      ! Saving and loading
      procedure, pass :: save => fabm_save
      procedure, pass :: load => fabm_load

      ! Finalization and deallocation
      procedure, pass :: finalize => fabm_finalize
      procedure, pass :: deallocate => fabm_deallocate
   end type SimstratFABM

contains

   ! Initialize: called once within the initialization of Simstrat
   ! Sets up memory, reads FABM configuration from FABMConfigFile, links the external Simstrat variables and sets the initial conditions of FABM variables
   subroutine fabm_init(self, state, model_cfg, fabm_cfg, output_cfg, sim_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelConfig), intent(inout) :: model_cfg
      class(FABMConfig), intent(inout) :: fabm_cfg
      class(OutputConfig), intent(inout) :: output_cfg
      class(SimConfig), intent(in) :: sim_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      integer :: ivar, exitstat
      logical :: config_path_exists, initial_path_exists
      character(len=256) :: mkdirCmd

      ! Check if FABM configuration path exists
      inquire(file=fabm_cfg%config_path,exist=config_path_exists)
      ! Create FABM configuration folder if it does not exist
      if (.not. config_path_exists) then
         call warn('FABM Configuration path does not exist, create folder according to config file...')
         mkdirCmd = 'mkdir '//trim(fabm_cfg%config_path)
         call execute_command_line(mkdirCmd, exitstat = exitstat)
         ! mkdir does not seem to accept a path to a folder in execute_command_line, thus a default result folder "FABM_configurations" will be generated in this case.
         if (exitstat==1) then
            call warn('FABM Configuration path specified in config file could not be generated. Default result folder "FABM_configurations" was generated instead.')
            call execute_command_line('mkdir FABM_configurations')
            fabm_cfg%config_path = 'FABM_configurations'
         end if
      end if
      ! Normalize configuration path
      call normalize_path(fabm_cfg%config_path)

      ! Create the bgc models according to the configurations in FABMConfigFile
      ! The models interact with FABM and describe properties of all bgc variables and parameters
      self%fabm_model => fabm_create_model(fabm_cfg%config_file, initialize = .false.)
      if (.not. associated(self%fabm_model)) then
         call error('FABM model creation failed. Check FABMConfigFile.')
      end if

      ! Register bgc variables from Simstrat (<ID>, <NAME>, <UNITS>, <LONG_NAME>)

      ! Register Simstrat contribution to attenuation coefficients, if there is FABM contribution to light extinction
      if (fabm_cfg%input_extinction < 1.0_RK) then
         call self%fabm_model%root%register_dependency(self%id_attenuation_coefficient_host, 'attenuation_coefficient_host', 'm-1', 'host contribution to attenuation coefficients')
      end if
      
      ! Add Simstrat bgc variables to FABM aggregate variables

      ! Add Simstrat contribution to attenuation coefficients, if there is FABM contribution to light extinction
      if (fabm_cfg%input_extinction < 1.0_RK) then
         call self%fabm_model%root%add_to_aggregate_variable(fabm_standard_variables%attenuation_coefficient_of_shortwave_flux, self%id_attenuation_coefficient_host)
         call self%fabm_model%root%add_to_aggregate_variable(fabm_standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_attenuation_coefficient_host)
      end if

      ! Initialize the model: Reads run-time FABM model configuration from FABMConfigFile
      ! and stores it in fabm_model. After this the number of biogeochemical variables is fixed
      call self%fabm_model%initialize()

      ! Provide extents of the spatial domain (number of layers nz for a 1D column)
      ! Used to allocate memory for FABM-managed spatially explicit fields
      ! Provide seconds_per_time_unit when time filtering functionality is used
      ! The product of seconds_per_time_unit and the t argument provided to prepare_inputs must have units of seconds
      call self%fabm_model%set_domain(grid%nz_grid, seconds_per_time_unit=real(sim_cfg%timestep, kind=RK))

      ! Set bottom location (pelagic-benthic interface) to 1 (default)
      ! Link FABM to bottom location
      k_bot = 1
      self%bottom_index => k_bot
      call self%fabm_model%set_bottom_index(self%bottom_index)

      ! Set kmax_bot to nz_grid if bottom_everywhere is true, 1 else
      kmax_bot = 1
      if (fabm_cfg%bottom_everywhere) kmax_bot = grid%nz_grid

      ! Allocate memory to hold the values of all size(fabm_model%*_state_variables) state variables
      ! All state variable values are combined in an array *_state with shape grid%nz_grid, fabm_cfg%n_*_state

      ! Interior state variables
      fabm_cfg%n_interior_state = size(self%fabm_model%interior_state_variables)
      if (fabm_cfg%n_interior_state > 0) then
         interior_or_bottom_exist = .true.
         allocate(state%fabm_interior_state(grid%nz_grid, fabm_cfg%n_interior_state))
      end if
      ! Bottom state variables
      fabm_cfg%n_bottom_state = size(self%fabm_model%bottom_state_variables)
      if (fabm_cfg%n_bottom_state > 0) then
         interior_or_bottom_exist = .true.
         allocate(state%fabm_bottom_state(kmax_bot, fabm_cfg%n_bottom_state))
      end if
      ! Surface state variables
      fabm_cfg%n_surface_state = size(self%fabm_model%surface_state_variables)
      if (fabm_cfg%n_surface_state > 0) then
         allocate(state%fabm_surface_state(fabm_cfg%n_surface_state))
      end if

      ! Calculate total amount of states
      fabm_cfg%n_state = fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + fabm_cfg%n_surface_state

      ! Allocate memory for additional information (name, units, long_name, minimum, maximum, missing_value) on state variables
      if (fabm_cfg%n_state > 0) allocate(output_cfg%output_vars_fabm_state(fabm_cfg%n_state))

      ! Point FABM to fields that hold state variable data and set additional information

      ! Interior state variables
      if (fabm_cfg%n_interior_state > 0) then
         do ivar = 1, fabm_cfg%n_interior_state
            call self%fabm_model%link_interior_state_data(ivar, state%fabm_interior_state(:,ivar))
            output_cfg%output_vars_fabm_state(ivar)%name = self%fabm_model%interior_state_variables(ivar)%name
            output_cfg%output_vars_fabm_state(ivar)%long_name = self%fabm_model%interior_state_variables(ivar)%long_name
            output_cfg%output_vars_fabm_state(ivar)%units = self%fabm_model%interior_state_variables(ivar)%units
            output_cfg%output_vars_fabm_state(ivar)%minimum = self%fabm_model%interior_state_variables(ivar)%minimum
            output_cfg%output_vars_fabm_state(ivar)%maximum = self%fabm_model%interior_state_variables(ivar)%maximum
            ! Special case for WET pelagic mirror variables (treat as benthic variable)
            if (self%fabm_model%interior_state_variables(ivar)%name(len_trim(self%fabm_model%interior_state_variables(ivar)%name)-2:) == '_PV') then
               output_cfg%output_vars_fabm_state(ivar)%benthic = .true.
               if (fabm_cfg%bottom_everywhere) then 
                  output_cfg%output_vars_fabm_state(ivar)%values => state%fabm_interior_state(:, ivar)
                  output_cfg%output_vars_fabm_state(ivar)%volume_grid = .true.
               else
                  output_cfg%output_vars_fabm_state(ivar)%global_value => state%fabm_interior_state(1, ivar)
               end if
            else
               output_cfg%output_vars_fabm_state(ivar)%values => state%fabm_interior_state(:, ivar)
               output_cfg%output_vars_fabm_state(ivar)%volume_grid = .true.
            end if
         end do
      end if
      ! Bottom state variables: link for bottom-most layer to fulfill FABM requirements, link for other layers later
      if (fabm_cfg%n_bottom_state > 0) then
         do ivar = 1, fabm_cfg%n_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1,ivar))
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%name = self%fabm_model%bottom_state_variables(ivar)%name
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%long_name = self%fabm_model%bottom_state_variables(ivar)%long_name
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%units = self%fabm_model%bottom_state_variables(ivar)%units
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%minimum = self%fabm_model%bottom_state_variables(ivar)%minimum
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%maximum = self%fabm_model%bottom_state_variables(ivar)%maximum
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%benthic = .true.
            if (fabm_cfg%bottom_everywhere) then 
               output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%values => state%fabm_bottom_state(:, ivar)
               output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%volume_grid = .true.
            else
               output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + ivar)%global_value => state%fabm_bottom_state(1, ivar)
            end if
         end do
      end if
      ! Surface state variables
      if (fabm_cfg%n_surface_state > 0) then
         do ivar = 1, fabm_cfg%n_surface_state
            call self%fabm_model%link_surface_state_data(ivar, state%fabm_surface_state(ivar))
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%name = self%fabm_model%surface_state_variables(ivar)%name
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%long_name = self%fabm_model%surface_state_variables(ivar)%long_name
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%units = self%fabm_model%surface_state_variables(ivar)%units
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%minimum = self%fabm_model%surface_state_variables(ivar)%minimum
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%maximum = self%fabm_model%surface_state_variables(ivar)%maximum
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%global_value => state%fabm_surface_state(ivar)
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%volume_grid = .false.
            output_cfg%output_vars_fabm_state(fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state + ivar)%face_grid = .false.
         end do
      end if
      
      ! Allocate interior tracer source terms [var_unit s-1]
      if (fabm_cfg%n_interior_state > 0) then
         allocate(self%sms_int(grid%nz_grid, fabm_cfg%n_interior_state))
      end if

      ! Allocate fluxes over pelagic-benthic interface [var_unit m s-1] and bottom tracer source terms [var_unit s-1]
      if (fabm_cfg%n_interior_state > 0 .or. fabm_cfg%n_bottom_state > 0) then
         allocate(self%flux_bt(kmax_bot, fabm_cfg%n_interior_state))
         allocate(self%sms_bt(kmax_bot, fabm_cfg%n_bottom_state))
      end if

      ! Allocate fluxes over air-water surface [var_unit m s-1] and surface tracer source terms [var_unit s-1]
      if (fabm_cfg%n_interior_state > 0 .or. fabm_cfg%n_surface_state > 0) then
         allocate(self%flux_sf(fabm_cfg%n_interior_state))
         allocate(self%sms_sf(fabm_cfg%n_surface_state))
      end if

      ! Allocate local vertical velocities of pelagic state variables (sinking, floating, active movement) [m s-1]
      ! Movement through water, independent of water flow
      if (fabm_cfg%n_interior_state > 0) then
         allocate(self%velocity(grid%nz_grid, fabm_cfg%n_interior_state))
      end if

      ! Allocate and initialize diagonals and area factors for
      ! time-integration of advection-diffusion-reaction equations
      if (fabm_cfg%n_interior_state > 0) then
         allocate(self%upper_diag_disc(grid%nz_grid))
         allocate(self%lower_diag_disc(grid%nz_grid))
         allocate(self%main_diag_disc(grid%nz_grid))
         allocate(self%AreaFactor_ext_1(grid%nz_grid))
         allocate(self%AreaFactor_ext_2(grid%nz_grid))
         self%upper_diag_disc = 0.0_RK
         self%lower_diag_disc = 0.0_RK
         self%main_diag_disc = 0.0_RK
         self%AreaFactor_ext_1 = 0.0_RK
         self%AreaFactor_ext_2 = 0.0_RK
      end if

      ! Allocate and initialize absorption passed to FABM (absorb_to_fabm)
      ! and FABM contribution to absorption (absorb_from_fabm)
      allocate(state%absorb_to_fabm(grid%nz_grid))
      allocate(state%absorb_from_fabm(grid%nz_grid))
      state%absorb_to_fabm = 0.0_RK
      state%absorb_from_fabm = 0.0_RK
      
      ! If there is FABM contribution to light extinction, get index of attenuation coefficient FABM diagnostic variable
      if (fabm_cfg%input_extinction < 1.0_RK) then
         att_index = 0
         do ivar = 1, size(self%fabm_model%interior_diagnostic_variables)
            if (self%fabm_model%interior_diagnostic_variables(ivar)%name == 'attenuation_coefficient_of_photosynthetic_radiative_flux') then
               self%fabm_model%interior_diagnostic_variables(ivar)%save = .true.
               att_index = ivar
            end if
         end do
         if (att_index == 0) then
            call error('Attenuation Coefficient not found in FABM Diagnostic Variables. Bioshade feedback not possible: Set InputExtinctionFactor to 1.')
         end if
      end if

      ! Point FABM to fields that contain values for environmental data, all variables are assumed to be allocated
      call link_environmental_data(self, state, fabm_cfg, grid)

      ! Write all diagnostic and aggregate variables to list_diagnostic
      ! Allocate arrays for diagnostic variables in output_diagnostics
      call list_diagnostic(self, fabm_cfg)
      call set_diagnostic_vars(self, fabm_cfg, output_cfg, grid)

      ! Allocate arrays for repaired variables in list_repaired and initialize them with the boundary value
      if (fabm_cfg%output_repaired_vars) then
         call set_repaired_vars(self, fabm_cfg, output_cfg, grid)
      else
         fabm_cfg%n_repaired = 0
         fabm_cfg%n_repaired_interior = 0
         fabm_cfg%n_repaired_bottom = 0
         fabm_cfg%n_repaired_surface = 0
      end if

      ! Allocate and set manipulations structure
      call set_manipulations(self, fabm_cfg, output_cfg, sim_cfg, grid)
      
      ! Complete initialization and check whether FABM has all dependencies fulfilled
      ! (i.e., whether all required calls to fabm_model%link_*_data have been made and all required data have been provided)
      ! Stop with fatal error if not
      ! Selection of diagnostics that FABM will compute and store becomes frozen
      call self%fabm_model%start()

      ! If model is not initialized with previously stored state
      ! Set FABM-provided initial values for state variables
      ! This sets the values of arrays sent to fabm_model%link_*_state_data
      ! Custom initial states from FABM initial conditions path overwrite FABM-provided initial values
      if (.not. sim_cfg%continue_from_snapshot) then
         ! Check if FABM initial conditions path exists
         inquire(file=fabm_cfg%initial_path,exist=initial_path_exists)
         ! Create FABM intial conditions folder if it does not exist
         if (.not. initial_path_exists) then
            call warn('FABM initial conditions path does not exist, create folder according to config file...')
            mkdirCmd = 'mkdir '//trim(fabm_cfg%initial_path)
            call execute_command_line(mkdirCmd, exitstat = exitstat)
            if (exitstat==1) then
               call warn('FABM initial conditions specified in config file could not be generated. Default initial conditions folder "FABM_initial" was generated instead.')
               call execute_command_line('mkdir FABM_initial')
               fabm_cfg%initial_path = 'FABM_initial'
            end if
         end if
         ! Normalize initial path
         call normalize_path(fabm_cfg%initial_path)
         ! Set interior state variable initial values
         call self%fabm_model%initialize_interior_state(1, grid%nz_grid)
         ! Overwrite by values from FABM intial conditions folder if present
         do ivar = 1, fabm_cfg%n_interior_state
            call read_initial_data(state, model_cfg, fabm_cfg, output_cfg, grid, ivar)
         end do
         ! Set bottom state variable initial values
         call self%fabm_model%initialize_bottom_state()
         ! If bottom_everywhere is set, at every depth:
         ! FABM reads k_bot through bottom_index and is pointed to location that holds state data for the current depth
         ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
         ! The bottom state at the current depth is initialized
         if (fabm_cfg%bottom_everywhere) then
            do k_bot = 2, kmax_bot
               do ivar = 1, fabm_cfg%n_bottom_state
                  call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
               end do
               call self%fabm_model%initialize_bottom_state()
            end do
            ! Reset Botom to 1
            do ivar = 1, fabm_cfg%n_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
            end do
            k_bot = 1
         end if
         ! Overwrite by values from FABM intial conditions folder if present
         do ivar = 1, fabm_cfg%n_bottom_state
            call read_initial_data(state, model_cfg, fabm_cfg, output_cfg, grid, ivar, fabm_cfg%n_interior_state)
         end do
         ! Set surface state variable initial values
         call self%fabm_model%initialize_surface_state()
         ! Overwrite by values from FABM intial conditions folder if present
         do ivar = 1, fabm_cfg%n_surface_state
            call read_initial_data(state, model_cfg, fabm_cfg, output_cfg, grid, ivar, fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state)
         end do
      end if

      ! Check (and repair) initial state values
      ! Interior state variables
      call self%fabm_model%check_interior_state(1, grid%nz_grid, fabm_cfg%repair_states, self%valid_int)
      ! Bottom state variables
      call self%fabm_model%check_bottom_state(fabm_cfg%repair_states, self%valid_bt)
      ! If bottom_everywhere is set, at every depth:
      ! FABM reads k_bot through bottom_index and is pointed to location that holds state data for the current depth
      ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
      ! The bottom state at the current depth is checked
      if (fabm_cfg%bottom_everywhere) then
         do k_bot = 2, kmax_bot
            do ivar = 1, fabm_cfg%n_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
            end do
            call self%fabm_model%check_bottom_state(fabm_cfg%repair_states, valid_bot)
            self%valid_bt = self%valid_bt .and. valid_bot
         end do
         ! Reset Botom to 1
         do ivar = 1, fabm_cfg%n_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
         end do
         k_bot = 1
      end if
      ! Surface state variables
      call self%fabm_model%check_surface_state(fabm_cfg%repair_states, self%valid_sf)
      ! Error if FABM intial values out of bounds and not repaired
      if (.not. (self%valid_int .and. self%valid_bt .and. self%valid_sf)) then
         if (.not. fabm_cfg%repair_states) then
            call error('FABM initial values out of bounds')
         else
            call warn('FABM initial values out of bounds repaired.')
         end if
      end if

      ! Call the update function once as a first call to initialize fluxes, sources, vertical velocities and other diagnostic variables
      call fabm_update(self, state, fabm_cfg, output_cfg, grid, .true.)
   end subroutine fabm_init

   ! The update function is called in the main loop of simstrat (in simstrat.f90) at every time step
   ! Particle atmospheric, pelagic and benthic fluxes and diffusion are computed to update bgc state variable values
   subroutine fabm_update(self, state, fabm_cfg, output_cfg, grid, first_call)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(in) :: fabm_cfg
      class(OutputConfig), intent(in) :: output_cfg
      class(StaggeredGrid), intent(in) :: grid
      logical, intent(in), optional :: first_call

      ! Local variables
      real(RK) :: add
      integer :: ivar, ivar_diag, index, k
      logical :: first_call_local

      ! Default not first call
      if (present(first_call)) then
         first_call_local = first_call
      else
         first_call_local = .false.
      end if

      ! If it is not the first call: 1. calculate and 2. validate the new model state, stop if variables are not valid
      ! Diagnostic variables from the previous timestep are used for this calculation
      if (.not. first_call_local) then
         ! Update timestep_counter
         timestep_counter = timestep_counter + 1
         ! Update kmax_bot to only calculate on occupied grid cells
         if (fabm_cfg%bottom_everywhere) kmax_bot = grid%nz_occupied

         ! 1a. Time-integrate the advection-diffusion-reaction equations of all tracers, 
         ! combining the Simstrat transport terms with the FABM biogeochemical source terms and fluxes (sms, flux) and vertical velocities (velocity).
         ! This results in an updated interior_state.
         ! First get helper arrays for all variables
         call integrate_helper(self, state, grid)
         ! Separate time-integration for all variables
         do ivar = 1, fabm_cfg%n_interior_state
            ! No diffusion for benthic variables
            if (output_cfg%output_vars_fabm_state(ivar)%benthic) then
               state%fabm_interior_state(1:kmax_bot, ivar) = state%fabm_interior_state(1:kmax_bot, ivar) + state%dt * self%flux_bt(1:kmax_bot, ivar) * grid%dAz_norm(1:kmax_bot)  
            else
               call integrate_interior_state(self, state, grid, ivar)
            end if
         end do

         ! 1b. Direct time integration of source terms to update bottom_state
         do ivar = 1, fabm_cfg%n_bottom_state
            state%fabm_bottom_state(1:kmax_bot, ivar) = state%fabm_bottom_state(1:kmax_bot, ivar) + state%dt * self%sms_bt(1:kmax_bot, ivar)
         end do

         ! 1c. Direct time integration of source terms to update surface_state
         do ivar = 1, fabm_cfg%n_surface_state
            state%fabm_surface_state(ivar) = state%fabm_surface_state(ivar) + state%dt * self%sms_sf(ivar)
         end do

         ! 1d. Perform manipulations
         do ivar = 1, fabm_cfg%n_manipulations
            associate(man => self%manipulations(ivar), &
                      int_state => state%fabm_interior_state, &
                      bt_state => state%fabm_bottom_state, &
                      sf_state => state%fabm_surface_state)
               select case (man%action_type)
               ! Multiply by action_val from start_depth to end_depth (reversed since indexing is from bottom to top)
               ! if current datum is between start time and start time plus one timestep
               case(1)
                  if ((se_floats(man%start_time, state%datum)) &
                     .and. (man%start_time + state%dt / SECONDS_PER_DAY > state%datum)) then
                     index = man%var_index
                     if (index <= fabm_cfg%n_interior_state) then
                        ! No diffusion for benthic variables
                        if (output_cfg%output_vars_fabm_state(index)%benthic) then
                           if (fabm_cfg%bottom_everywhere) then
                              int_state(man%end_depth:man%start_depth,index) = int_state(man%end_depth:man%start_depth,index) * man%action_val
                           else
                              int_state(1,index) = int_state(1,index) * man%action_val
                           end if
                        else
                           int_state(man%end_depth:man%start_depth,index) = int_state(man%end_depth:man%start_depth,index) * man%action_val
                        end if
                     else if (index <= fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state) then
                        index = index - fabm_cfg%n_interior_state
                        if (fabm_cfg%bottom_everywhere) then
                           bt_state(man%end_depth:man%start_depth,index) = bt_state(man%end_depth:man%start_depth,index) * man%action_val
                        else
                           bt_state(1,index) = bt_state(1,index) * man%action_val
                        end if
                     else
                        index = index - fabm_cfg%n_interior_state - fabm_cfg%n_bottom_state
                        sf_state(index) = sf_state(index) * man%action_val
                     end if
                  end if
               ! Add action_val times timestep from start_depth to end_depth if current datum is between start time and end time
               ! Only if there is any value lower than threshold between start_depth and end_depth
               ! For variables on volume grid action_val is additionally divided by the height between start_depth and end_depth
               case(2)
                  if ((se_floats(man%start_time, state%datum)) &
                     .and. (ge_floats(man%end_time, state%datum))) then
                     index = man%var_index
                     if (index <= fabm_cfg%n_interior_state) then
                        if (output_cfg%output_vars_fabm_state(index)%benthic) then
                           if (fabm_cfg%bottom_everywhere) then
                              if (any(int_state(man%end_depth:man%start_depth,index) < man%threshold)) then
                                 add = man%action_val * state%dt / sum(grid%h(man%end_depth:man%start_depth))
                                 int_state(man%end_depth:man%start_depth,index) = int_state(man%end_depth:man%start_depth,index) + add
                              end if
                           else if (int_state(1,index) < man%threshold) then
                              int_state(1,index) = int_state(1,index) + state%dt * man%action_val
                           end if
                        else
                           if (any(int_state(man%end_depth:man%start_depth,index) < man%threshold)) then
                              add = man%action_val * state%dt / sum(grid%h(man%end_depth:man%start_depth))
                              int_state(man%end_depth:man%start_depth,index) = int_state(man%end_depth:man%start_depth,index) + add
                           end if
                        end if
                     else if (index <= fabm_cfg%n_interior_state + fabm_cfg%n_bottom_state) then
                        index = index - fabm_cfg%n_interior_state
                        if (fabm_cfg%bottom_everywhere) then
                           if (any(bt_state(man%end_depth:man%start_depth,index) < man%threshold)) then
                              add = man%action_val * state%dt / sum(grid%h(man%end_depth:man%start_depth))
                              bt_state(man%end_depth:man%start_depth,index) = bt_state(man%end_depth:man%start_depth,index) + add
                           end if
                        else if (bt_state(1,index) < man%threshold) then
                           bt_state(1,index) = bt_state(1,index) + state%dt * man%action_val
                        end if
                     else
                        index = index - fabm_cfg%n_interior_state - fabm_cfg%n_bottom_state
                        if (sf_state(index) < man%threshold) then
                           sf_state(index) = sf_state(index) + state%dt * man%action_val
                        end if
                     end if
                  end if
               end select
            end associate
         end do

         ! 2a. Simstrat register FABM interior state out of bounds
         if (fabm_cfg%output_repaired_vars) then
            do ivar = 1, fabm_cfg%n_interior_state
               if (any(state%fabm_interior_state(1:kmax_bot, ivar) < self%fabm_model%interior_state_variables(ivar)%minimum)) then
                  ! Register repaired variable or store out-of-bound value if already registered
                  if (fabm_cfg%n_repaired_interior > 0) then
                     index = findloc(self%repaired_interior_names, trim(self%fabm_model%interior_state_variables(ivar)%name)//'_minimum', dim = 1)
                  else
                     index = 0
                  end if
                  if (index > 0) then
                     do k = 1, kmax_bot
                        if (state%fabm_interior_state(k, ivar) < self%fabm_model%interior_state_variables(ivar)%minimum) then
                           self%repaired_interior(k, index) = state%fabm_interior_state(k, ivar)
                        else
                           self%repaired_interior(k, index) = self%fabm_model%interior_state_variables(ivar)%minimum
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%interior_state_variables(ivar)%name, 'minimum', self%fabm_model%interior_state_variables(ivar)%minimum)
                  end if
               end if
               if (any(state%fabm_interior_state(1:kmax_bot, ivar) > self%fabm_model%interior_state_variables(ivar)%maximum)) then
                  ! Register repaired variable or store out-of-bound value if already registered
                  if (fabm_cfg%n_repaired_interior > 0) then
                     index = findloc(self%repaired_interior_names, trim(self%fabm_model%interior_state_variables(ivar)%name)//'_maximum', dim = 1)
                  else
                     index = 0
                  end if
                  if (index > 0) then
                     do k = 1, kmax_bot
                        if (state%fabm_interior_state(k, ivar) > self%fabm_model%interior_state_variables(ivar)%maximum) then
                           self%repaired_interior(k, index) = state%fabm_interior_state(k, ivar)
                        else
                           self%repaired_interior(k, index) = self%fabm_model%interior_state_variables(ivar)%maximum
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%interior_state_variables(ivar)%name, 'maximum', self%fabm_model%interior_state_variables(ivar)%maximum)
                  end if
               end if
            end do
         end if

         ! 2a. FABM check (and repair) interior state variables
         call self%fabm_model%check_interior_state(1, kmax_bot, fabm_cfg%repair_states, self%valid_int)

         ! 2b. Simstrat register FABM bottom state out of bounds
         if (fabm_cfg%output_repaired_vars) then
            do ivar = 1, fabm_cfg%n_bottom_state  
               if (any(state%fabm_bottom_state(1:kmax_bot, ivar) < self%fabm_model%bottom_state_variables(ivar)%minimum)) then
                  ! Register repaired variable or store out-of-bound value if already registered
                  if (fabm_cfg%n_repaired_bottom > 0) then
                     index = findloc(self%repaired_bottom_names, trim(self%fabm_model%bottom_state_variables(ivar)%name)//'_minimum', dim = 1)
                  else
                     index = 0
                  end if
                  if (index > 0) then
                     do k = 1, kmax_bot
                        if (state%fabm_bottom_state(k, ivar) < self%fabm_model%bottom_state_variables(ivar)%minimum) then
                           self%repaired_bottom(k, index) = state%fabm_bottom_state(k, ivar)
                        else
                           self%repaired_bottom(k, index) = self%fabm_model%bottom_state_variables(ivar)%minimum
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%bottom_state_variables(ivar)%name, 'minimum', self%fabm_model%bottom_state_variables(ivar)%minimum)
                  end if
               end if
               if (any(state%fabm_bottom_state(1:kmax_bot, ivar) > self%fabm_model%bottom_state_variables(ivar)%maximum)) then
                  ! Register repaired variable or store out-of-bound value if already registered
                  if (fabm_cfg%n_repaired_bottom > 0) then
                     index = findloc(self%repaired_bottom_names, trim(self%fabm_model%bottom_state_variables(ivar)%name)//'_maximum', dim = 1)
                  else
                     index = 0
                  end if
                  if (index > 0) then
                     do k = 1, kmax_bot
                        if (state%fabm_bottom_state(k, ivar) > self%fabm_model%bottom_state_variables(ivar)%maximum) then
                           self%repaired_bottom(k, index) = state%fabm_bottom_state(k, ivar)
                        else
                           self%repaired_bottom(k, index) = self%fabm_model%bottom_state_variables(ivar)%maximum
                        end if
                     end do
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%bottom_state_variables(ivar)%name, 'maximum', self%fabm_model%bottom_state_variables(ivar)%maximum)
                  end if
               end if
            end do
         end if
         
         ! 2b. FABM check (and repair) bottom state variables
         call self%fabm_model%check_bottom_state(fabm_cfg%repair_states, self%valid_bt)
         if (fabm_cfg%bottom_everywhere) then
            ! If bottom_everywhere is set, at every depth:
            ! FABM reads k_bot through bottom_index and is pointed to location that holds state data for the current depth
            ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth
            ! The bottom state at the current depth is validated
            do k_bot = 2, kmax_bot
               do ivar = 1, fabm_cfg%n_bottom_state
                  call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
               end do
               call self%fabm_model%check_bottom_state(fabm_cfg%repair_states, valid_bot)
               self%valid_bt = self%valid_bt .and. valid_bot
            end do
            ! Reset Bottom to 1
            do ivar = 1, fabm_cfg%n_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
            end do
            k_bot = 1
         end if

         ! 2c. Simstrat register FABM surface state out of bounds
         if (fabm_cfg%output_repaired_vars) then
            do ivar = 1, fabm_cfg%n_surface_state
               if (state%fabm_surface_state(ivar) < self%fabm_model%surface_state_variables(ivar)%minimum) then
                  ! Register repaired variable or store out-of-bound value if already registered
                  if (fabm_cfg%n_repaired_surface > 0) then
                     index = findloc(self%repaired_surface_names, trim(self%fabm_model%surface_state_variables(ivar)%name)//'_minimum', dim = 1)
                  else
                     index = 0
                  end if
                  if (index > 0) then
                     self%repaired_surface(index) = self%fabm_model%surface_state_variables(ivar)%minimum
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%surface_state_variables(ivar)%name, 'minimum', self%fabm_model%surface_state_variables(ivar)%minimum)
                  end if
               end if
               if (state%fabm_surface_state(ivar) > self%fabm_model%surface_state_variables(ivar)%maximum) then
                  ! Register repaired variable or store out-of-bound value if already registered
                  if (fabm_cfg%n_repaired_surface > 0) then
                     index = findloc(self%repaired_surface_names, trim(self%fabm_model%surface_state_variables(ivar)%name)//'_maximum', dim = 1)
                  else
                     index = 0
                  end if
                  if (index > 0) then
                     self%repaired_surface(index) = self%fabm_model%surface_state_variables(ivar)%maximum
                  else
                     call list_repaired(fabm_cfg, self%fabm_model%surface_state_variables(ivar)%name, 'maximum', self%fabm_model%surface_state_variables(ivar)%maximum)
                  end if
               end if
            end do
         end if

         ! 2c. FABM check (and repair) surface state variables
         call self%fabm_model%check_surface_state(fabm_cfg%repair_states, self%valid_sf)

         ! 2d. Error if FABM out of bounds and not repaired
         if (.not. (self%valid_int .and. self%valid_bt .and. self%valid_sf)) then
            if (.not. fabm_cfg%repair_states) then
               call error('FABM variables out of bounds')
            end if
         end if
      else
         ! Initialize timestep_counter at first call
         timestep_counter = 0
         write (6, *) 'Initializing FABM diagnostics'
      end if
      
      ! Initialize (at first call) or update fluxes, sources, vertical movement and other diagnostic variables
      
      ! 1a. Prepare all fields (e.g. light attenuation) FABM needs to compute fluxes and source terms
      ! Operates on entire active spatial domain
      ! To enable FABM's built-in time filters provide argument t that describes the model time since start (in timesteps)
      call self%fabm_model%prepare_inputs(real(timestep_counter, kind=RK))

      ! 1b. Write FABM absorption to state, if there is FABM contribution to light extinction
      if (fabm_cfg%input_extinction < 1.0_RK) then
         ! Get total absorption
         state%absorb_from_fabm = self%fabm_model%get_interior_diagnostic_data(att_index)
         ! Check for negative light extinction inside FABM due to drastic changes in the biogeochemical calculation
         ! If there is negative light extinction, set absorption passed to FABM to background there 
         ! to avoid negative values and re-calculate all fields
         if (any(state%absorb_from_fabm < 0.0_RK)) then
            do k = 1, size(state%absorb_from_fabm)
               if (state%absorb_from_fabm(k) < 0.0_RK) then
                  state%absorb_to_fabm(k) = fabm_cfg%background_extinction
                  call warn('Negative light extinction value at height '//trim(str_int(k))//' inside FABM -- momentarily changed absorption passed to FABM to background there. Avoid by decreasing Timestep or changing InputExtinctionFactor.')
               end if
            end do
            call self%fabm_model%prepare_inputs(real(timestep_counter, kind=RK))
         end if
         ! Subtract absorption passed to FABM to get only biogeochemical contribution to extinction
         state%absorb_from_fabm = state%absorb_from_fabm - state%absorb_to_fabm
      end if

      ! 2. Retrieve sources, fluxes and vertical velocities
      ! >0: flux into water
      ! <0: flux out of water

      ! 2a. Initialize with 0 and then retrieve interior tracer source terms
      ! FABM increments rather than sets argument sms_int: initialize with 0 at every time step
      if (fabm_cfg%n_interior_state > 0) then
         self%sms_int(1:kmax_bot, :) = 0.0
         call self%fabm_model%get_interior_sources(1, kmax_bot, self%sms_int(1:kmax_bot, :))
         ! Set interior source NaNs to 0
         if (any(ieee_is_nan(self%sms_int(1:kmax_bot, :)))) then
            call warn('FABM Interior Source contains NaNs, set to 0.')
            where (ieee_is_nan(self%sms_int(1:kmax_bot, :)))
               self%sms_int(1:kmax_bot, :) = 0.0
            end where
         end if
      end if

      ! 2b. Initialize with 0 and then retrieve fluxes over pelagic-benthic interface and bottom tracer source terms
      ! FABM increments rather than sets argument flux_bt, sms_bt: initialize with 0 at every time step
      if (fabm_cfg%n_interior_state > 0 .or. fabm_cfg%n_bottom_state > 0) then
         if (fabm_cfg%n_interior_state > 0) self%flux_bt(1:kmax_bot, :) = 0.0
         if (fabm_cfg%n_bottom_state > 0) self%sms_bt(1:kmax_bot, :) = 0.0
         call self%fabm_model%get_bottom_sources(self%flux_bt(1, :), self%sms_bt(1, :))
         ! Set bottom flux NaNs to 0
         if (fabm_cfg%n_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_bt(1, :)))) then
               call warn('FABM Bottom Flux contains NaN, set to 0.')
               where (ieee_is_nan(self%flux_bt(1, :)))
                  self%flux_bt(1, :) = 0.0
               end where
            end if
         end if
         ! Set bottom source NaNs to 0
         if (fabm_cfg%n_bottom_state > 0) then
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
      if (fabm_cfg%n_interior_state > 0 .or. fabm_cfg%n_surface_state > 0) then
         if (fabm_cfg%n_interior_state > 0) self%flux_sf = 0.0
         if (fabm_cfg%n_surface_state > 0) self%sms_sf = 0.0
         call self%fabm_model%get_surface_sources(self%flux_sf, self%sms_sf)
         ! Set surface flux NaNs to 0
         if (fabm_cfg%n_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_sf))) then
               call warn('FABM Surface Flux contains NaN, set to 0')
               where (ieee_is_nan(self%flux_sf))
                  self%flux_sf = 0.0
               end where
            end if
         end if
         ! Set surface source NaNs to 0
         if (fabm_cfg%n_surface_state > 0) then
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
      if (fabm_cfg%n_interior_state > 0) then
         if (.not. state%first_timestep) self%velocity(1:kmax_bot, :) = 0.0
         call self%fabm_model%get_vertical_movement(1, kmax_bot, self%velocity(1:kmax_bot, :))
         ! Set interior velocity NaNs to 0
         if (any(ieee_is_nan(self%velocity(1:kmax_bot, :)))) then
            call warn('FABM Interior Velocity contains NaN, set to 0')
            where (ieee_is_nan(self%velocity(1:kmax_bot, :)))
               self%velocity(1:kmax_bot, :) = 0.0
            end where
         end if
      end if

      ! 4. Compute any remaining diagnostics not computed by preceding routines
      ! Operates on entire active spatial domain
      call self%fabm_model%finalize_outputs()

      ! 5. Retrieve values of diagnostic variables
      ! Interior diagnostics
      if (fabm_cfg%n_diagnostic_interior > 0) then
         do ivar_diag = 1, fabm_cfg%n_diagnostic_interior
            index = self%diagnostic_interior_index(ivar_diag)
            self%diagnostic_interior(:, ivar_diag) = self%fabm_model%get_interior_diagnostic_data(index)
         end do
      end if
      ! Horizontal diagnostics
      if (fabm_cfg%n_diagnostic_horizontal > 0) then
         do ivar_diag = 1, fabm_cfg%n_diagnostic_horizontal
            index = self%diagnostic_horizontal_index(ivar_diag)
            self%diagnostic_horizontal(1, ivar_diag) = self%fabm_model%get_horizontal_diagnostic_data(index)
         end do
      end if

      ! If bottom_everywhere is set and interior or bottom variables or horizontal diagnostics exits,
      ! repeat relevant steps from 1. to 5. at every depth
      if (fabm_cfg%bottom_everywhere .and. (interior_or_bottom_exist .or. diagnostic_horizontal_exist)) then
         ! The bottom (the location of the pelagic-benthic interface) is moved to the current depth (k_bot)
         ! FABM reads k_bot through bottom_index
         ! FABM is pointed to location that holds state data for the current depth (link_bottom_state_data)
         ! FABM fields are prepared (prepare_inputs) on entire domain
         ! The bottom fluxes and sources at the current depth are calculated (get_bottom_sources)
         ! Remaining diagnostics are calculated at the current depth (finalize_outputs) on entire domain
         ! The horizontal diagnostic variables at the current depth are retrieved (get_horizontal_diagnostic_data)
         ! Note: Since FABM does not distinguish between bottom and surface diagnostic variables, 
         ! surface diagnostic variables are also retrieved for every depth, with constant value
         do k_bot = 2, kmax_bot
            do ivar = 1, fabm_cfg%n_bottom_state
               call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(k_bot, ivar))
            end do
            if (interior_or_bottom_exist) then
               call self%fabm_model%prepare_inputs(real(timestep_counter, kind=RK))
               call self%fabm_model%get_bottom_sources(self%flux_bt(k_bot, :), self%sms_bt(k_bot, :))
            end if
            if (diagnostic_horizontal_exist) then
               call self%fabm_model%finalize_outputs()
               do ivar_diag = 1, fabm_cfg%n_diagnostic_horizontal
                  index = self%diagnostic_horizontal_index(ivar_diag)
                  self%diagnostic_horizontal(k_bot, ivar_diag) = self%fabm_model%get_horizontal_diagnostic_data(index)
               end do
            end if
         end do
         ! Set NaNs in bottom fluxes to 0
         if (fabm_cfg%n_interior_state > 0) then
            if (any(ieee_is_nan(self%flux_bt(1:kmax_bot, :)))) then
               call warn('FABM Bottom Flux contains NaN, set to 0')
               where (ieee_is_nan(self%flux_bt(1:kmax_bot, :)))
                  self%flux_bt(1:kmax_bot, :) = 0.0
               end where
            end if
         end if
         ! Set NaNs in bottom sources to 0
         if (fabm_cfg%n_bottom_state > 0) then
            if (any(ieee_is_nan(self%sms_bt(1:kmax_bot, :)))) then
               call warn('FABM Bottom Source contains NaN, set to 0')
               where (ieee_is_nan(self%sms_bt(1:kmax_bot, :)))
                  self%sms_bt(1:kmax_bot, :) = 0.0
               end where
            end if
         end if
         ! Reset Bottom to 1
         do ivar = 1, fabm_cfg%n_bottom_state
            call self%fabm_model%link_bottom_state_data(ivar, state%fabm_bottom_state(1, ivar))
         end do
         k_bot = 1
         ! Re-calculate all fields with bottom at 1
         call self%fabm_model%prepare_inputs(real(timestep_counter, kind=RK))
         call self%fabm_model%finalize_outputs()
      end if

      ! This is the ideal moment for the output: 
      ! All variables (enviromental, state, diagnostics, source, fluxes, vertical velocities) are in sync
   end subroutine fabm_update

   ! Point FABM to all variables on FABM's standard variable list that the model can provide
   ! For this list, visit https://github.com/fabm-model/fabm/wiki/List-of-standard-variables
   ! Update this subroutine in case of "FABM is lacking required data"
   subroutine link_environmental_data(self, state, fabm_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(FABMConfig), intent(inout) :: fabm_cfg
      class(StaggeredGrid), intent(in) :: grid
      
      ! Interior variables (arrays)
      ! Attenuation coefficient of SWR flux [m-1] and of PAR flux [m-1]
      ! If there is FABM contribution to light extinction, link the Simstrat contribution (absorb_to_fabm)
      ! If there is not, Simstrat calculation of attenuation coefficient overwrites bgc model calculation
      ! Note: Passing the same attenuation coefficient for SWR and PAR is not exactly correct 
      ! since the attenuation coefficient for PAR can be higher than for SWR
      ! see for example: https://www.sciencedirect.com/science/article/pii/S1001074217317734
      if (fabm_cfg%input_extinction < 1.0_RK) then
         call self%fabm_model%link_interior_data('attenuation_coefficient_host', state%absorb_to_fabm)
      else
         call self%fabm_model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_shortwave_flux, state%absorb_to_fabm, source=data_source_user)
         call self%fabm_model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, state%absorb_to_fabm, source=data_source_user)
      end if
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
      ! Alkalinity expressed as mole equivalent [mmol m-3], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%alkalinity_expressed_as_mole_equivalent)
      ! Fractional saturation of oxygen [1], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%fractional_saturation_of_oxygen)
      ! Mass concentration of suspended matter [g m-3], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%mass_concentration_of_suspended_matter)
      ! Mole concentration of * [mmol m-3], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%mole_concentration_of_*)
      ! Net rate of SWR energy absorption at each layer [W m-2], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%net_rate_of_absorption_of_shortwave_energy_in_layer)
      ! PH reported on total scale [1], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%ph_reported_on_total_scale)
      ! Secchi depth [m], not defined in Simstrat
      !call self%fabm_model%link_interior_data(fabm_standard_variables%secchi_depth)
      ! Projection factor for benthic flux into horizontal layer volume [m-1]
      self%id_bot_pel_conv = self%fabm_model%get_interior_variable_id(self%bot_pel_conv)
      call self%fabm_model%link_interior_data(self%id_bot_pel_conv, grid%dAz_norm)

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
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%wind_speed, state%uv10_gas)
      ! latitude [degree_north]
      call self%fabm_model%link_horizontal_data(fabm_standard_variables%latitude, state%Lat)
      ! Depth below geoid [m], only necessary when coupling with geodetic data, not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%bottom_depth_below_geoid)
      ! Bottom roughness [m], constant in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%bottom_roughness_length)
      ! Mole fraction of carbon dioxide in air [1e-6], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air)
      ! PAR flux in air [W m-2], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux_in_air)
      ! SWR flux in air [W m-2], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux_in_air)
      ! Longitude [degree_east], not defined in Simstrat
      !call self%fabm_model%link_horizontal_data(fabm_standard_variables%longitude)

      ! Global variables (scalars)
      ! Number of days since start of the year [days]
      call self%fabm_model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, state%current_day_of_year)
      ! Total * [mmol m-3], not defined in Simstrat
      !call self%fabm_model%link_scalar(fabm_standard_variables%total_*, state%current_day_of_year)
      ! Total iron [umol m-3], not defined in Simstrat
      !call self%fabm_model%link_scalar(fabm_standard_variables%total_iron, state%current_day_of_year)
   end subroutine link_environmental_data

   ! Read initial data and set in state
   subroutine read_initial_data(state, model_cfg, fabm_cfg, output_cfg, grid, ivar, ivar_offset)
      ! Arguments
      class(ModelState), intent(inout) :: state
      class(ModelConfig), intent(in) :: model_cfg
      class(FABMConfig), intent(in) :: fabm_cfg
      class(OutputConfig), intent(inout) :: output_cfg
      class(StaggeredGrid), intent(in) :: grid
      integer, intent(in) :: ivar
      integer, intent(in), optional :: ivar_offset

      ! Local variables
      integer :: n, unit, status, ivar_state
      character(len=256) :: file_path
      real(RK) :: depth, initial_val
      real(RK), allocatable :: temp_depths(:), temp_initials(:)
      real(RK), dimension(:), allocatable :: depths, initials
      logical :: is_int

      ! Allocation of temporary information arrays
      allocate(temp_depths(model_cfg%max_length_input_data))
      allocate(temp_initials(model_cfg%max_length_input_data))

      ! Set index including indexing offset
      ! Only interior state variables are passed without an offset
      if (present(ivar_offset)) then
         ivar_state = ivar + ivar_offset
         is_int = .false.
      else
         ivar_state = ivar
         is_int = .true.
      end if

      ! Check if initial data file exists, return if not
      file_path = trim(fabm_cfg%initial_path)//'/'//trim(output_cfg%output_vars_fabm_state(ivar_state)%name)//'_initial.dat'     
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      if (status .ne. 0) then
         close(unit)
         return
      else
         write(6,*) 'Reading ', trim(file_path)
      end if
      ! Skip header
      read(unit, *, iostat=status)
      ! Initialize counts to 0
      n = 0
      ! Read depth and initial value at every line
      ! Pass depths as positive values
      do
         read(unit, *, iostat=status) depth, initial_val
         if (status .ne. 0) exit
         n = n + 1
         if (n > model_cfg%max_length_input_data) then
            call error('Too many lines in '//trim(file_path)//', increase MaxLengthInputData in ModelConfig.')
         end if
         if (depth > 0.0_RK) then
            call error('One or several input depths of initial conditions file '//trim(file_path)//' are positive.')
         end if
         temp_depths(n) = abs(depth)
         temp_initials(n) = initial_val
      end do
      ! Return if no initial values provided
      ! Allocate depth and initial value array otherwise
      if (n==0) then
         call warn('File without initial values '//trim(file_path)//' ignored.')
         close(unit)
         return
      else
         allocate(depths(n))
         allocate(initials(n))
         depths(:) = temp_depths(:n)
         initials(:) = temp_initials(:n)
      end if
      ! Reverse order
      call reverse_in_place(depths)
      call reverse_in_place(initials)
      ! Set depths relative to absolute depth of lowest layer
      depths = grid%z_zero - depths
      ! Interior state variable
      ! Check if variable is on volume grid and interpolate to volume grid if yes
      ! Variable is not on volume grid for benthic variables if bottom everywhere is set to false
      if (is_int) then
         if (output_cfg%output_vars_fabm_state(ivar_state)%volume_grid) then
            if (n==1) then
               state%fabm_interior_state(1:grid%nz_grid, ivar) = initials(1)
            else
               call grid%interpolate_to_vol(depths, initials, n, state%fabm_interior_state(1:grid%nz_grid, ivar))
            end if
         else
            if (n==1) then
               state%fabm_interior_state(1, ivar) = initials(1)
            else
               call warn('File '//trim(file_path)//' has depth distributed initial values for global variable: only deepest initial value considered.')
               state%fabm_interior_state(1, ivar) = initials(1)
            end if
         end if
      ! Bottom state variable
      ! Check if variable is on volume grid and interpolate to volume grid if yes
      ! Variable is not on volume grid if bottom everywhere is set to false
      else if (output_cfg%output_vars_fabm_state(ivar_state)%benthic) then
         if (output_cfg%output_vars_fabm_state(ivar_state)%volume_grid) then
            if (n==1) then
               state%fabm_bottom_state(1:grid%nz_grid, ivar) = initials(1)
            else
               call grid%interpolate_to_vol(depths, initials, n, state%fabm_bottom_state(1:grid%nz_grid, ivar))
            end if
         else
            if (n==1) then
               state%fabm_bottom_state(1, ivar) = initials(1)
            else
               call warn('File '//trim(file_path)//' has depth distributed initial values for global variable: only deepest initial value considered.')
               state%fabm_bottom_state(1, ivar) = initials(1)
            end if
         end if
      ! Surface state variable
      ! Always a global variable
      else 
         if (n==1) then
            state%fabm_surface_state(ivar) = initials(1)
         else
            call warn('File '//trim(file_path)//' has depth distributed initial values for global variable: only highest initial value considered.')
            state%fabm_surface_state(ivar) = initials(n)
         end if
      end if
      ! Close the file
      call ok('Initial data file '//trim(file_path)//' successfully read.')
      close(unit)
   end subroutine read_initial_data

   ! Subroutine to help build transport matrix A in integrate_interior_state
   subroutine integrate_helper(self, state, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(StaggeredGrid), intent(in) :: grid

      ! Build diagonals of A with Simstrat transport terms (code from discretization module)
      ! With diffusivity for temperature (nuh)
      ! Upward flux from (1:grid%nz_occupied - 1) to (2:grid%nz_occupied): upper_diag_disc(z) describes the upward flux into cell z
      self%upper_diag_disc(2:grid%nz_occupied) = state%dt*state%nuh(2:grid%nz_occupied)*grid%AreaFactor_1(2:grid%nz_occupied)
      self%upper_diag_disc(1) = 0.0_RK ! No upward flux into bottommost cell
      ! Downward flux from (2:grid%nz_occupied) to (1:grid%nz_occupied - 1): lower_diag_disc(z) describes the downward flux into cell z
      self%lower_diag_disc(1:grid%nz_occupied-1) = state%dt*state%nuh(2:grid%nz_occupied)*grid%AreaFactor_2(1:grid%nz_occupied - 1)
      self%lower_diag_disc(grid%nz_occupied) = 0.0_RK ! No downward flux into uppermost cell
      ! 1 - downward flux out - upward flux out
      self%main_diag_disc(1:grid%nz_occupied) = 1.0_RK - self%upper_diag_disc(1:grid%nz_occupied) - self%lower_diag_disc(1:grid%nz_occupied)

      ! Area factors for external fluxes, minus sign to be closer to the AreaFactors from grid
      ! AreaFactor_ext_1 is for upward fluxes:
      ! Multiply by the layer through which the flux passes (z = i) and divide by the volume of the receiving cell (z = i)
      self%AreaFactor_ext_1(1:grid%nz_occupied) = -grid%Az(1:grid%nz_occupied) / (grid%h(1:grid%nz_occupied) * grid%Az_vol(1:grid%nz_occupied))
      ! AreaFactor_ext_2 is for downward fluxes: 
      ! Multiply by the layer through which the flux passes (z = i+1) and divide by the volume of the receiving cell (z = i)
      self%AreaFactor_ext_2(1:grid%nz_occupied) = -grid%Az(2:grid%nz_occupied+1) / (grid%h(1:grid%nz_occupied) * grid%Az_vol(1:grid%nz_occupied))
   end subroutine integrate_helper

   ! Integration algorithm for interior state variables: Simstrat transport terms and FABM biogeochemical terms integrated simultaneously
   ! Assuming small enough dt such that only fluxes between neighbouring layers are relevant
   ! Also assuming that sides of the lake slope inward
   subroutine integrate_interior_state(self, state, grid, ivar)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(StaggeredGrid), intent(in) :: grid
      integer, intent(in) :: ivar

      ! Local variables
      real(RK), dimension(grid%nz_occupied) :: velocity_up, velocity_down, sources, lower_diag, main_diag, upper_diag, rhs
      integer :: k

      ! Create linear system of equations with transport and source terms
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      ! Add FABM fluxes to A as residual verical advection terms
      ! Benthic-pelagic and air-water fluxes are added as source terms to S below

      ! Upward and downward movement rates for each layer
      ! if the variable moves upward velocity_up is positive and velocity_down zero
      ! vice versa if the variable moves downward
      do k = 1, grid%nz_occupied
         if (ge_floats(self%velocity(k, ivar), 0.0_RK)) then
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
      upper_diag(2:grid%nz_occupied) = self%upper_diag_disc(2:grid%nz_occupied) + state%dt*velocity_up(1:grid%nz_occupied-1)*self%AreaFactor_ext_1(2:grid%nz_occupied)
      upper_diag(1) = self%upper_diag_disc(1)
      ! Subtract upward flux from cell z to cell z+1
      ! main_diag should be >= 0
      ! Water-air flux is treated as sink term below
      main_diag(1:grid%nz_occupied-1) = self%main_diag_disc(1:grid%nz_occupied-1) - state%dt*velocity_up(1:grid%nz_occupied-1)*self%AreaFactor_ext_2(1:grid%nz_occupied-1)
      main_diag(grid%nz_occupied) = self%main_diag_disc(grid%nz_occupied)
      ! Add downward flux from (2:grid%nz_occupied) to (1:grid%nz_occupied - 1): lower_diag(z) describes the downward flux into cell z
      ! lower_diag should be <= 0
      ! Air-water flux is treated as source term below
      lower_diag(1:grid%nz_occupied-1) = self%lower_diag_disc(1:grid%nz_occupied-1) + state%dt*velocity_down(2:grid%nz_occupied)*self%AreaFactor_ext_2(1:grid%nz_occupied-1)
      lower_diag(grid%nz_occupied) = self%lower_diag_disc(grid%nz_occupied)
      ! Subtract downward flux from cell z to cell z-1
      ! main_diag should be >= 0
      ! Pelagic-benthic flux is treated as sink term below
      main_diag(2:grid%nz_occupied) = main_diag(2:grid%nz_occupied) - state%dt*velocity_down(2:grid%nz_occupied)*self%AreaFactor_ext_1(2:grid%nz_occupied)
      main_diag(1) = main_diag(1)

      ! Get source S^{n}
      ! Source at each layer
      sources = self%sms_int(1:grid%nz_occupied, ivar)
      ! Add pelagic-benthic and air-water flux [var_unit m s-1] as source [var_unit s-1]
      ! Convert bottom flux to source by multiplication by sediment area over layer volume (dAz_norm, [m])
      sources(1:kmax_bot) = sources(1:kmax_bot) + (self%flux_bt(1:kmax_bot, ivar) * grid%dAz_norm(1:kmax_bot))
      ! Convert surface flux to source by division by surface area over volume of uppermost layer [m]
      sources(grid%nz_occupied) = sources(grid%nz_occupied) + (self%flux_sf(ivar) * grid%Az(grid%nz_occupied) / (grid%h(grid%nz_occupied) * grid%Az_vol(grid%nz_occupied)))

      ! Calculate RHS (phi^{n}+dt*S^{n})
      rhs(1:grid%nz_occupied) = state%fabm_interior_state(1:grid%nz_occupied, ivar) + state%dt*sources(1:grid%nz_occupied)

      ! Solve LES to get phi^{n+1}
      call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%fabm_interior_state(1:grid%nz_occupied, ivar), grid%nz_occupied)
   end subroutine integrate_interior_state

   ! Output list of diagnostic variables in output_diagnostics
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
         file_path = trim(fabm_cfg%config_path)//'/list_diagnostic_interior.dat'
         ! Create new file or overwrite already existing one
         open(newunit=unit, file=file_path, action='write', iostat = status)
         if (status .ne. 0) then
            call error('Failed to open or create file: ' // trim(file_path) // '. Please check FABM configurations folder.')
         end if
         ! Write the header
         write(unit, '(A)', advance = 'no') 'Short Name, '
         write(unit, '(A)', advance = 'no') 'Long Name, '
         write(unit, '(A)', advance = 'no') 'Units, '
         write(unit, '(A)') 'Output'
         ! Write the names, units and output toggle of interior diagnostic variables
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
         ! Close the file
         close(unit)
      end if

      ! Write the names of all horizontal diagnostic variables
      if (size(self%fabm_model%horizontal_diagnostic_variables) > 0) then
         ! Construct the file path for Horizontal Diagnostic Variable List
         file_path = trim(fabm_cfg%config_path)//'/list_diagnostic_horizontal.dat'
         ! Creat new file or overwrite already existing one
         open(newunit=unit, file=file_path, action='write', iostat = status)
         ! Write the header
         write(unit, '(A)', advance = 'no') 'Short Name, '
         write(unit, '(A)', advance = 'no') 'Long Name, '
         write(unit, '(A)', advance = 'no') 'Units, '
         write(unit, '(A)') 'Output'
         ! Write the names, units and output toggle of horizontal diagnostic variables
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
         ! Close the file
         close(unit)
      end if
   end subroutine list_diagnostic

   ! Read names of diagnostic Vars in output_diagnostics file
   subroutine set_diagnostic_vars(self, fabm_cfg, output_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(FABMConfig), intent(inout) :: fabm_cfg
      class(OutputConfig), intent(inout) :: output_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      integer :: i, j, n, n_int, n_hor, n_off, unit, status, unit_list, status_list
      character(len=256) :: line, line_trim, short_name, short_name_trim, long_name, units, output, file_path
      integer, parameter :: max_lines = 2000 ! Maximum amount of diagnostic variables in output_diagnostics file
      character(len=256), allocatable :: temp_names(:)
      integer, dimension(max_lines) :: fabm_index
      logical, dimension(max_lines) :: is_int, is_hor

      ! Allocation of temporary information arrays
      allocate(temp_names(max_lines))

      ! Construct the file path
      file_path = trim(fabm_cfg%config_path)//'/output_diagnostics.dat'

      ! Read from output_diagnostics
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      ! Read the file or if it does not exist
      ! create an empty file and set amount of diagnostics to zero
      if (status .ne. 0) then
         call warn('No FABM diagnostic variables provided. Empty diagnostic variables file '// trim(file_path) //' created.')
         open(newunit=unit, action='write', status='new', file=file_path, iostat = status)
         if (status .ne. 0) then
            call error('Failed to open or create file: ' // trim(file_path) // '. Please check FABM configurations folder.')
         end if
         close(unit)
         fabm_cfg%n_diagnostic_interior = 0
         fabm_cfg%n_diagnostic_horizontal = 0
         fabm_cfg%n_diagnostic = 0
         return
      else
         write(6,*) 'Reading ', trim(file_path)
      end if
      ! Initialize counts to 0
      n = 0
      n_int = 0
      n_hor = 0
      n_off = 0
      ! Read every diagnostic variable
      do
         read(unit, '(A)', iostat=status) line
         if (status .ne. 0) exit
         if (len_trim(line) == 0) cycle
         if (any(temp_names == trim(line))) cycle
         line_trim = trim(line)
         ! Read from list of diagnostics in case of select diagnostics
         ! Otherwise just add the line
         if ((line_trim(:11) == 'select_all_') .or. (line_trim(:14) == 'select_output_')) then
            ! Read from list_diagnostic_interior.dat
            ! And add every line from there that has not been previously added
            file_path = trim(fabm_cfg%config_path)//'/list_diagnostic_interior.dat'
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
                     call error('Too many lines in '//trim(file_path)//', increase max_lines in set_diagnostic_vars.')
                  end if
                  temp_names(n) = trim(short_name)
                  fabm_index(n) = 0
               end if
            end do
            ! Close list_diagnostic_interior.dat
            close(unit_list)
            ! Read from list_diagnostic_horizontal.dat
            ! And add every line from there that has not been previously added
            file_path = trim(fabm_cfg%config_path)//'/list_diagnostic_horizontal.dat'
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
                     call error('Too many lines in '//trim(file_path)//', increase max_lines in set_diagnostic_vars.')
                  end if
                  temp_names(n) = trim(short_name)
                  fabm_index(n) = 0
               end if
            end do
            ! Close list_diagnostic_horizontal.dat
            close(unit_list)
         else
            n = n + 1
            if (n > max_lines) then
               call error('Too many lines in '//trim(file_path)//', increase max_lines in set_diagnostic_vars.')
            end if
            temp_names(n) = trim(line)
            fabm_index(n) = 0
         end if
      end do
      ! Close the file
      call ok('FABM diagnostic variables file '// trim(file_path) //' successfully read')
      close(unit)

      ! Find diagnostic variable in FABM model
      do i = 1, n
         ! Set index at location in fabm_model%interior_diagnostic_variables
         do j = 1, size(self%fabm_model%interior_diagnostic_variables)
            if (self%fabm_model%interior_diagnostic_variables(j)%name == temp_names(i)) then
               self%fabm_model%interior_diagnostic_variables(j)%save = .true.
               fabm_index(i) = j
               is_int(i) = .true.
               is_hor(i) = .false.
               n_int = n_int + 1
            end if
         end do
         ! Set index at location in fabm_model%horizontal_diagnostic_variables
         do j = 1, size(self%fabm_model%horizontal_diagnostic_variables)
            if (self%fabm_model%horizontal_diagnostic_variables(j)%name == temp_names(i)) then
               self%fabm_model%horizontal_diagnostic_variables(j)%save = .true.
               fabm_index(i) = j
               is_int(i) = .false.
               is_hor(i) = .true.
               n_hor = n_hor + 1
            end if
         end do
         ! Call warning if not found in fabm_model
         if (fabm_index(i) == 0) then
            call warn('FABM diagnostic variable '//trim(temp_names(i))//' not found.')
            is_int(i) = .false.
            is_hor(i) = .false.
            n_off = n_off + 1
         end if
      end do

      ! Set amount of diagnostic variables and allocate arrays of proper size
      fabm_cfg%n_diagnostic = n - n_off
      fabm_cfg%n_diagnostic_interior = n_int
      fabm_cfg%n_diagnostic_horizontal = n_hor
      if (n_int > 0) then
         allocate(self%diagnostic_interior(grid%nz_grid, n_int))
         allocate(self%diagnostic_interior_index(n_int))
      end if
      if (n_hor > 0) then
         diagnostic_horizontal_exist = .true.
         allocate(self%diagnostic_horizontal(kmax_bot, n_hor))
         allocate(self%diagnostic_horizontal_index(n_hor))
      end if

      ! Allocate array for and set output information (name, units, long_name, minimum, maximum, missing_value) on diagnostic variables
      ! Also store the index inside FABM
      if ((n - n_off) > 0) allocate(output_cfg%output_vars_fabm_diagnostic(n - n_off))
      ! j is the index in output_vars_fabm_diagnostic
      j = 1
      ! Interior diagnostic variables are set first
      do i = 1, n
         if (is_int(i)) then
            output_cfg%output_vars_fabm_diagnostic(j)%name = self%fabm_model%interior_diagnostic_variables(fabm_index(i))%name
            output_cfg%output_vars_fabm_diagnostic(j)%long_name = self%fabm_model%interior_diagnostic_variables(fabm_index(i))%long_name
            output_cfg%output_vars_fabm_diagnostic(j)%units = self%fabm_model%interior_diagnostic_variables(fabm_index(i))%units
            output_cfg%output_vars_fabm_diagnostic(j)%minimum = self%fabm_model%interior_diagnostic_variables(fabm_index(i))%minimum
            output_cfg%output_vars_fabm_diagnostic(j)%maximum = self%fabm_model%interior_diagnostic_variables(fabm_index(i))%maximum
            output_cfg%output_vars_fabm_diagnostic(j)%values => self%diagnostic_interior(:, j)
            output_cfg%output_vars_fabm_diagnostic(j)%volume_grid = .true.
            self%diagnostic_interior_index(j) = fabm_index(i)
            j = j + 1
         end if
      end do
      ! Horizontal diagnostic variables are set second
      ! Note: FABM does not differentiate between bottom and surface diagnostic variables: here all are treated as bottom diagnostic variables
      do i = 1, n
         if (is_hor(i)) then
            output_cfg%output_vars_fabm_diagnostic(j)%name = self%fabm_model%horizontal_diagnostic_variables(fabm_index(i))%name
            output_cfg%output_vars_fabm_diagnostic(j)%long_name = self%fabm_model%horizontal_diagnostic_variables(fabm_index(i))%long_name
            output_cfg%output_vars_fabm_diagnostic(j)%units = self%fabm_model%horizontal_diagnostic_variables(fabm_index(i))%units
            output_cfg%output_vars_fabm_diagnostic(j)%minimum = self%fabm_model%horizontal_diagnostic_variables(fabm_index(i))%minimum
            output_cfg%output_vars_fabm_diagnostic(j)%maximum = self%fabm_model%horizontal_diagnostic_variables(fabm_index(i))%maximum
            output_cfg%output_vars_fabm_diagnostic(j)%benthic = .true.
            if (fabm_cfg%bottom_everywhere) then
               output_cfg%output_vars_fabm_diagnostic(j)%values => self%diagnostic_horizontal(:, j - n_int)
               output_cfg%output_vars_fabm_diagnostic(j)%volume_grid = .true.
            else
               output_cfg%output_vars_fabm_diagnostic(j)%global_value => self%diagnostic_horizontal(1, j - n_int)
            end if
            self%diagnostic_horizontal_index(j - n_int) = fabm_index(i)
            j = j + 1
         end if
      end do
   end subroutine set_diagnostic_vars

   ! Add entry to list of repaired variables in list_repaired
   subroutine list_repaired(fabm_cfg, variable, boundary, boundary_value)
      ! Arguments
      class(FABMConfig), intent(in) :: fabm_cfg
      character(len=*), intent(in) :: variable
      character(len=*), intent(in) :: boundary
      real(RK), intent(in) :: boundary_value

      ! Local variables
      integer :: unit, status
      character(len=256) :: file_path, line, new_line

      ! Construct the file path
      file_path = trim(fabm_cfg%config_path)//'/list_repaired.dat'

      ! Construct the new repaired entry
      new_line = '"'//trim(variable)//'", "'//trim(boundary)//'", "'//trim(real_to_str(boundary_value, '(ES15.6)'))//'"'

      ! Check if repaired variable case already present in list_repaired
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      if (status .ne. 0) then
         call error('Failed to open or create file: ' // trim(file_path) // '. Please check FABM configurations folder.')
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
         call error('Failed to open or create file: ' // trim(file_path) // '. Please check FABM configurations folder.')
      else
         call warn('FABM variable '//trim(variable)//' added to '// trim(file_path)//'. Restart simulation to output '//trim(boundary)//' values.')
      end if

      ! Write variable name, its boundary and the boundary vale
      write(unit, '(A)') trim(new_line)

      ! Close the file after writing
      close(unit)
   end subroutine list_repaired

   ! Register repaired variables in list_repaired
   subroutine set_repaired_vars(self, fabm_cfg, output_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(FABMConfig), intent(inout) :: fabm_cfg
      class(OutputConfig), intent(inout) :: output_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      integer :: i, j, n, n_int, n_bt, n_sf, unit, status
      character(len=256) :: file_path
      character(len=256) :: name
      character(len=16) :: bound
      character(len=256), allocatable  :: temp_names(:)
      character(len=16), allocatable :: temp_bounds(:)
      character(len=16), allocatable :: temp_types(:)
      integer, parameter :: max_lines = 1000 ! Maximum amount of repaired variables in list_repaired
      integer, dimension(max_lines) :: fabm_index

      ! Allocation of temporary information arrays
      allocate(temp_names(max_lines))
      allocate(temp_bounds(max_lines))
      allocate(temp_types(max_lines))

      ! Initialize counts to zero
      n = 0
      n_int = 0
      n_bt = 0
      n_sf = 0

      ! Construct the file path
      file_path = trim(fabm_cfg%config_path)//'/list_repaired.dat'

      ! Read from list_repaired
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      ! Create empty file if no file is provided
      if (status .ne. 0) then
         call warn('No FABM Repaired Variables provided. Empty repaired variables file '// trim(file_path) //' created.')
         open(newunit=unit, action='write', status='new', file=file_path, iostat = status)
         if (status .ne. 0) then
            call error('Failed to open or create file: ' // trim(file_path) // '. Please check FABM configurations folder.')
         end if
         ! Write the header
         write(unit, '(A)', advance = 'no') 'Variable, '
         write(unit, '(A)', advance = 'no') 'Boundary reached, '
         write(unit, '(A)') 'Boundary value'
         close(unit)
         ! Set amount of repaired variables to 0
         fabm_cfg%n_repaired = 0
         fabm_cfg%n_repaired_interior = 0
         fabm_cfg%n_repaired_bottom = 0
         fabm_cfg%n_repaired_surface = 0
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
            call error('Too many lines in '//trim(file_path)//', increase max_lines in set_repaired_vars.')
         end if
         temp_names(n) = trim(name)
         if (.not. ((trim(bound) == 'minimum') .or. (trim(bound) == 'maximum'))) then
            call error('FABM repaired variable '//trim(temp_names(n))//' has invalid bound, set to "minimum" or "maximum".')
         end if
         temp_bounds(n) = trim(bound)
         temp_types(n) = 'undefined'
         fabm_index(n) = 0
         ! Search in interior variables
         do j = 1, fabm_cfg%n_interior_state
            if (self%fabm_model%interior_state_variables(j)%name == temp_names(n)) then
               temp_types(n) = 'interior'
               fabm_index(n) = j
               n_int = n_int + 1
            end if
         end do
         ! Search in bottom variables
         do j = 1, fabm_cfg%n_bottom_state
            if (self%fabm_model%bottom_state_variables(j)%name == temp_names(n)) then
               temp_types(n) = 'bottom'
               fabm_index(n) = j
               n_bt = n_bt + 1
            end if
         end do
         ! Search in surface variables
         do j = 1, fabm_cfg%n_surface_state
            if (self%fabm_model%surface_state_variables(j)%name == temp_names(n)) then
               temp_types(n) = 'surface'
               fabm_index(n) = j
               n_sf = n_sf + 1
            end if
         end do
         ! Call error if not found in fabm_model
         if (fabm_index(n) == 0) then
            call error('FABM repaired variable '//trim(temp_names(n))//' not found. Control '//trim(file_path)//'.')
         end if
      end do
      ! Close the file
      call ok('FABM repaired variables file '// trim(file_path) //' successfully read')
      close(unit)

      ! Set amount of repaired variables and allocate arrays of proper
      fabm_cfg%n_repaired = n
      fabm_cfg%n_repaired_interior = n_int
      fabm_cfg%n_repaired_bottom = n_bt
      fabm_cfg%n_repaired_surface = n_sf
      if (n_int > 0) then
         allocate(self%repaired_interior(grid%nz_grid, n_int))
         allocate(self%repaired_interior_names(n_int))
      end if
      if (n_bt > 0) then
         allocate(self%repaired_bottom(kmax_bot, n_bt))
         allocate(self%repaired_bottom_names(n_bt))
      end if
      if (n_sf > 0) then
         allocate(self%repaired_surface(n_sf))
         allocate(self%repaired_surface_names(n_sf))
      end if

      ! Allocate array for and set output information (name, units, long_name, minimum, maximum, missing_value) on repaired variables
      ! Also initialize state arrays to minimum / maximum value
      ! Also add names of repaired variables to arrays in self
      if (n > 0) allocate(output_cfg%output_vars_fabm_repaired(n))
      ! j is the index in output_vars_fabm_repaired
      j = 1
      ! Interior repaired variables are set first
      do i = 1, n
         if (temp_types(i) == 'interior') then
            output_cfg%output_vars_fabm_repaired(j)%name = trim(self%fabm_model%interior_state_variables(fabm_index(i))%name)//'_'//trim(temp_bounds(i))
            output_cfg%output_vars_fabm_repaired(j)%long_name = trim(self%fabm_model%interior_state_variables(fabm_index(i))%long_name)//' '//trim(temp_bounds(i))//' value'
            output_cfg%output_vars_fabm_repaired(j)%units = self%fabm_model%interior_state_variables(fabm_index(i))%units
            if (temp_bounds(i) == 'minimum') then
               self%repaired_interior(:, j) = self%fabm_model%interior_state_variables(fabm_index(i))%minimum
               output_cfg%output_vars_fabm_repaired(j)%maximum = self%fabm_model%interior_state_variables(fabm_index(i))%minimum
            else
               self%repaired_interior(:, j) = self%fabm_model%interior_state_variables(fabm_index(i))%maximum
               output_cfg%output_vars_fabm_repaired(j)%minimum = self%fabm_model%interior_state_variables(fabm_index(i))%maximum
            end if
            ! Special case for WET pelagic mirror variables (treat as benthic variable)
            if (self%fabm_model%interior_state_variables(fabm_index(i))%name(len_trim(self%fabm_model%interior_state_variables(fabm_index(i))%name)-2:) == '_PV') then
               output_cfg%output_vars_fabm_repaired(j)%benthic = .true.
               if (fabm_cfg%bottom_everywhere) then 
                  output_cfg%output_vars_fabm_repaired(j)%values => self%repaired_interior(:, j)
                  output_cfg%output_vars_fabm_repaired(j)%volume_grid = .true.
               else
                  output_cfg%output_vars_fabm_repaired(j)%global_value => self%repaired_interior(1, j)
               end if
            else
               output_cfg%output_vars_fabm_repaired(j)%values => self%repaired_interior(:, j)
               output_cfg%output_vars_fabm_repaired(j)%volume_grid = .true.
            end if
            self%repaired_interior_names(j) = trim(self%fabm_model%interior_state_variables(fabm_index(i))%name)//'_'//trim(temp_bounds(i))
            j = j + 1
         end if
      end do
      ! Bottom repaired variables are set second
      do i = 1, n
         if (temp_types(i) == 'bottom') then
            output_cfg%output_vars_fabm_repaired(j)%name = trim(self%fabm_model%bottom_state_variables(fabm_index(i))%name)//'_'//trim(temp_bounds(i))
            output_cfg%output_vars_fabm_repaired(j)%long_name = trim(self%fabm_model%bottom_state_variables(fabm_index(i))%long_name)//' '//trim(temp_bounds(i))//' value'
            output_cfg%output_vars_fabm_repaired(j)%units = self%fabm_model%bottom_state_variables(fabm_index(i))%units
            if (temp_bounds(i) == 'minimum') then
               self%repaired_bottom(:, j - n_int) = self%fabm_model%bottom_state_variables(fabm_index(i))%minimum
               output_cfg%output_vars_fabm_repaired(j)%maximum = self%fabm_model%bottom_state_variables(fabm_index(i))%minimum
            else
               self%repaired_bottom(:, j - n_int) = self%fabm_model%bottom_state_variables(fabm_index(i))%maximum
               output_cfg%output_vars_fabm_repaired(j)%minimum = self%fabm_model%bottom_state_variables(fabm_index(i))%maximum
            end if
            output_cfg%output_vars_fabm_repaired(j)%benthic = .true.
            if (fabm_cfg%bottom_everywhere) then 
               output_cfg%output_vars_fabm_repaired(j)%values => self%repaired_bottom(:, j - n_int)
               output_cfg%output_vars_fabm_repaired(j)%volume_grid = .true.
            else
               output_cfg%output_vars_fabm_repaired(j)%global_value => self%repaired_bottom(1, j - n_int)
            end if
            self%repaired_bottom_names(j - n_int) = trim(self%fabm_model%bottom_state_variables(fabm_index(i))%name)//'_'//trim(temp_bounds(i))
            j = j + 1
         end if
      end do
      ! Surface repaired variables are set last
      do i = 1, n
         if (temp_types(i) == 'surface') then
            output_cfg%output_vars_fabm_repaired(j)%name = trim(self%fabm_model%surface_state_variables(fabm_index(i))%name)//'_'//trim(temp_bounds(i))
            output_cfg%output_vars_fabm_repaired(j)%long_name = trim(self%fabm_model%surface_state_variables(fabm_index(i))%long_name)//' '//trim(temp_bounds(i))//' value'
            output_cfg%output_vars_fabm_repaired(j)%units = self%fabm_model%surface_state_variables(fabm_index(i))%units
            if (temp_bounds(i) == 'minimum') then
               self%repaired_surface(j - n_int - n_bt) = self%fabm_model%surface_state_variables(fabm_index(i))%minimum
               output_cfg%output_vars_fabm_repaired(j)%maximum = self%fabm_model%surface_state_variables(fabm_index(i))%minimum
            else
               self%repaired_surface(j - n_int - n_bt) = self%fabm_model%surface_state_variables(fabm_index(i))%maximum
               output_cfg%output_vars_fabm_repaired(j)%minimum = self%fabm_model%surface_state_variables(fabm_index(i))%maximum
            end if
            output_cfg%output_vars_fabm_repaired(j)%global_value => self%repaired_surface(j - n_int - n_bt)
            self%repaired_interior_names(j - n_int - n_bt) = trim(self%fabm_model%surface_state_variables(fabm_index(i))%name)//'_'//trim(temp_bounds(i))
            j = j + 1
         end if
      end do
   end subroutine set_repaired_vars

   ! Read manipulation namelist manipulations.nml
   subroutine set_manipulations(self, fabm_cfg, output_cfg, sim_cfg, grid)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      class(FABMConfig), intent(inout) :: fabm_cfg
      class(OutputConfig), intent(inout) :: output_cfg
      class(SimConfig), intent(in) :: sim_cfg
      class(StaggeredGrid), intent(in) :: grid

      ! Local variables
      character(len=256) :: file_path
      character(len=64) :: var_name
      integer :: i, ivar, n, unit, status, action_type
      real(RK) :: start_time, end_time, action_val, threshold, start_depth, end_depth, offset_above_level, depth_to_grid
      integer, parameter :: max_manipulations = 1000 ! Maximum amount of manipulations
      character(len=64), allocatable  :: temp_names(:)
      integer, dimension(max_manipulations) :: temp_indices, temp_types, temp_kstarts, temp_kstops
      real(RK), dimension(max_manipulations) :: temp_starts, temp_ends, temp_actions, temp_thresholds
      logical :: found
      namelist /manipulation/ var_name, start_time, end_time, action_type, action_val, threshold, start_depth, end_depth

      ! Allocation of temporary information arrays
      allocate(temp_names(max_manipulations))

      ! Calculate conversion to depth below lake level
      ! Difference between lowest depth and total amount of depths
      offset_above_level = grid%z_zero - grid%z_face(grid%nz_grid + 1)
      ! Calculate conversion from input depth to grid number
      ! Divide by difference of lowest depth to highest depth
      ! Multiply by difference of highest grid point to lowest (1)
      depth_to_grid = (grid%nz_grid - 1) / (grid%z_face(grid%nz_grid + 1) - grid%z_face(1))

      ! Construct the file path
      file_path = trim(fabm_cfg%config_path)//'/manipulations.nml'

      ! Open and read manipulations file
      open(newunit=unit, action='read', status='old', file=file_path, iostat = status)
      ! Create empty file and set amount of manipulations to zero if there is no file provided
      ! Otherwise read file
      if (status .ne. 0) then
         call warn('No FABM manipulations provided. Empty manipulations file '// trim(file_path) //' created.')
         open(newunit=unit, action='write', status='new', file=file_path, iostat = status)
         if (status .ne. 0) then
            call error('Failed to open or create file: ' // trim(file_path) // '. Please check FABM configurations folder.')
         end if
         ! Close the file and return
         close(unit)
         fabm_cfg%n_manipulations = 0
         return
      else
         write (6, *) 'Reading ', trim(file_path)
      end if
      ! Initialize count
      n = 0
      ! Read until end of file is reached
      do
         ! Initialize manipulation
         var_name = ''
         start_time = sim_cfg%start_datum
         end_time = sim_cfg%end_datum
         action_type = -1
         action_val = 0.0_RK
         threshold = huge(0.0_RK)
         start_depth = 0
         end_depth = -grid%z_zero
         ! In fortran the last namelist entry is already counted as EOF
         ! Exit to avoid error because trying to read an entry after the last
         if (status < 0) exit
         ! Read and check for invalid format or empty manipulation
         read (unit, nml=manipulation, iostat=status)
         if (status > 0) then
            call error('Invalid format in '//trim(file_path)//'.')
         end if
         if (action_type == -1) cycle
         ! Increase amount of manipulations
         n = n + 1
         ! Check if too many manipulations are provided
         if (n > max_manipulations) then
            call error('Too many manipulations in '//trim(file_path)//', increase max_manipulations in set_manipulations.')
         end if
         ! Search for var_name in state variables, skip manipulation if not found
         found = .false.
         do ivar = 1, fabm_cfg%n_state
            if (output_cfg%output_vars_fabm_state(ivar)%name == var_name) then
               temp_names(n) = trim(var_name)
               temp_indices(n) = ivar
               found = .true.
            end if
         end do
         if (.not. found) then
            call warn('Manipulation file '//trim(file_path)//' contains manipulation with variable '//trim(var_name)//' not found in FABM state: ignored.')
            n = n - 1   
            cycle
         end if
         ! Set start_time if it is before end_datum, skip manipulation otherwise
         if (se_floats(start_time, sim_cfg%end_datum)) then
            temp_starts(n) = start_time
         else
            call warn('Manipulation file '//trim(file_path)//' contains manipulation for variable '//trim(var_name)//' with start_time after end datum: ignored.')
            n = n - 1   
            cycle
         end if
         ! Set end_time if it is after start_datum and start_time, skip manipulation otherwise
         if ((ge_floats(end_time, sim_cfg%start_datum)) &
            .and. (ge_floats(end_time, start_time))) then
            temp_ends(n) = end_time
         else
            call warn('Manipulation file '//trim(file_path)//' contains manipulation for variable '//trim(var_name)//' with end_time before start_time or start datum: ignored.')
            n = n - 1   
            cycle
         end if    
         ! Set action_type if it is 1 or 2, skip manipulation otherwise
         if ((action_type == 1) .or. (action_type == 2)) then
            temp_types(n) = action_type
         else
            call warn('Manipulation file '//trim(file_path)//' contains manipulation for variable '//trim(var_name)//' with invalid action type (not 1 or 2): ignored.')
            n = n - 1   
            cycle
         end if
         if (compare_floats(action_val, 0.0_RK)) then
            call warn('Manipulation file '//trim(file_path)//' contains manipulation for variable '//trim(var_name)//' with action_val 0.0 or not provided.')
         end if
         ! Set action_val and threshold
         temp_actions(n) = action_val
         temp_thresholds(n) = threshold
         ! Depth input is rounded to nearest integer after conversion to grid number and substracted from amount of grid points
         ! Set kstart if start_depth is above or equal to lowest depth, skip manipulation otherwise
         ! Restrict kstart to amount of grid cells
         if (ge_floats(start_depth, -grid%z_zero, 1.0e-4_RK)) then
            temp_kstarts(n) = min(grid%nz_grid, grid%nz_grid + nint((start_depth + offset_above_level) * depth_to_grid))
         else
            call warn('Manipulation file '//trim(file_path)//' contains manipulation for variable '//trim(var_name)//' with start_depth below lowest depth: ignored.')
            n = n - 1   
            cycle
         end if
         ! Set kstop if end_depth is below or equal to 0, skip manipulation otherwise
         ! Minimum kstop is 1
         if ((se_floats(end_depth, -offset_above_level, 1.0e-4_RK)) &
            .and. (se_floats(end_depth, start_depth, 1.0e-4_RK))) then
            temp_kstops(n) = max(1, grid%nz_grid + nint((end_depth + offset_above_level) * depth_to_grid))
         else
            call warn('Manipulation file '//trim(file_path)//' contains manipulation for variable '//trim(var_name)//' with end_depth above 0 or start_depth: ignored.')
            n = n - 1   
            cycle
         end if
      end do
      ! Close the file and return if no valid manipulation has been provided
      call ok('FABM manipulations file '// trim(file_path) //' successfully read')
      close(unit)
      if (n == 0) then
         fabm_cfg%n_manipulations = 0
         return
      end if

      ! Allocate type of correct size
      allocate(self%manipulations(n))
      fabm_cfg%n_manipulations = n

      ! Assign values
      do i = 1, n
         self%manipulations(i)%var_name = trim(temp_names(i))
         self%manipulations(i)%var_index = temp_indices(i)
         self%manipulations(i)%start_time = temp_starts(i)
         self%manipulations(i)%end_time = temp_ends(i)
         self%manipulations(i)%action_type = temp_types(i)
         self%manipulations(i)%action_val = temp_actions(i)
         self%manipulations(i)%threshold = temp_thresholds(i)
         self%manipulations(i)%start_depth = temp_kstarts(i)
         self%manipulations(i)%end_depth = temp_kstops(i)
      end do
   end subroutine set_manipulations

   ! Save data stored in self
   subroutine fabm_save(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self

      ! Save arrays stored in self
      call save_matrix(80, self%sms_int)
      call save_matrix(80, self%velocity)
      call save_matrix(80, self%flux_bt)
      call save_matrix(80, self%sms_bt)
      call save_array(80, self%flux_sf)
      call save_array(80, self%sms_sf)
      call save_array(80, self%upper_diag_disc)
      call save_array(80, self%lower_diag_disc)
      call save_array(80, self%main_diag_disc)
      call save_array(80, self%AreaFactor_ext_1)
      call save_array(80, self%AreaFactor_ext_2)
      write(80) self%valid_int, self%valid_sf, self%valid_bt
      write(80) self%bottom_index
      call save_matrix_pointer(80, self%diagnostic_interior)
      call save_matrix_pointer(80, self%diagnostic_horizontal)
   end subroutine fabm_save

   ! Load data stored in self
   subroutine fabm_load(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self

      ! Load arrays stored in self
      call read_matrix(81, self%sms_int)
      call read_matrix(81, self%velocity)
      call read_matrix(81, self%flux_bt)
      call read_matrix(81, self%sms_bt)
      call read_array(81, self%flux_sf)
      call read_array(81, self%sms_sf)
      call read_array(81, self%upper_diag_disc)
      call read_array(81, self%lower_diag_disc)
      call read_array(81, self%main_diag_disc)
      call read_array(81, self%AreaFactor_ext_1)
      call read_array(81, self%AreaFactor_ext_2)
      read(81) self%valid_int, self%valid_sf, self%valid_bt
      read(81) self%bottom_index
      call read_matrix_pointer(81, self%diagnostic_interior)
      call read_matrix_pointer(81, self%diagnostic_horizontal)
   end subroutine fabm_load

   ! Finalize the coupling
   subroutine fabm_finalize(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self
      
      ! Nullify the model and finalize the library
      if (associated(self%fabm_model)) then
         nullify(self%fabm_model)
         call fabm_finalize_library()
      end if

      ! Deallocate memory allocated in self
      call fabm_deallocate(self)
   end subroutine fabm_finalize

   ! Deallocate all data allocated in self
   subroutine fabm_deallocate(self)
      ! Arguments
      class(SimstratFABM), intent(inout) :: self

      ! Deallocate arrays stored in self
      if (allocated(self%sms_int)) deallocate(self%sms_int)
      if (allocated(self%flux_bt)) deallocate(self%flux_bt)
      if (allocated(self%sms_bt)) deallocate(self%sms_bt)
      if (allocated(self%flux_sf)) deallocate(self%flux_sf)
      if (allocated(self%sms_sf)) deallocate(self%sms_sf)
      if (allocated(self%velocity)) deallocate(self%velocity)
      if (allocated(self%upper_diag_disc)) deallocate(self%upper_diag_disc)
      if (allocated(self%lower_diag_disc)) deallocate(self%lower_diag_disc)
      if (allocated(self%main_diag_disc)) deallocate(self%main_diag_disc)
      if (allocated(self%AreaFactor_ext_1)) deallocate(self%AreaFactor_ext_1)
      if (allocated(self%AreaFactor_ext_2)) deallocate(self%AreaFactor_ext_2)
      if (associated(self%bottom_index)) nullify(self%bottom_index)
      if (associated(self%diagnostic_interior)) then
         deallocate(self%diagnostic_interior)
         nullify(self%diagnostic_interior)
      end if
      if (associated(self%diagnostic_horizontal)) then
         deallocate(self%diagnostic_horizontal)
         nullify(self%diagnostic_horizontal)
      end if
      if (allocated(self%diagnostic_interior_index)) deallocate(self%diagnostic_interior_index)
      if (allocated(self%diagnostic_horizontal_index)) deallocate(self%diagnostic_horizontal_index)
      if (associated(self%repaired_interior)) then
         deallocate(self%repaired_interior)
         nullify(self%repaired_interior)
      end if
      if (associated(self%repaired_bottom)) then
         deallocate(self%repaired_bottom)
         nullify(self%repaired_bottom)
      end if
      if (associated(self%repaired_surface)) then
         deallocate(self%repaired_surface)
         nullify(self%repaired_surface)
      end if
      if (allocated(self%repaired_interior_names)) deallocate(self%repaired_interior_names)
      if (allocated(self%repaired_bottom_names)) deallocate(self%repaired_bottom_names)
      if (allocated(self%repaired_surface_names)) deallocate(self%repaired_surface_names)
      if (allocated(self%manipulations)) deallocate(self%manipulations)
   end subroutine fabm_deallocate

end module simstrat_fabm