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
!     |  Data structure definitions for simulation data
!<    +---------------------------------------------------------------+

module strat_simdata
   use strat_kinds
   use strat_grid
   use strat_consts
   use utilities
   implicit none
   private

   ! All Input files
   type, public :: InputConfig
      character(len=:), allocatable          :: MorphName
      character(len=:), allocatable          :: InitName
      character(len=:), allocatable          :: ForcingName
      character(len=:), allocatable          :: AbsorpName
      character(len=:), allocatable          :: GridName
      character(len=:), allocatable          :: QinpName
      character(len=:), allocatable          :: QoutName
      character(len=:), allocatable          :: TinpName
      character(len=:), allocatable          :: SinpName
      real(RK), dimension(:), allocatable    :: read_grid_array_from_json
      real(RK) :: read_grid_value_from_json
      integer :: grid_input_type
   end type

   ! Definition of a variable to log
   type, public :: LogVariable
      character(len=:), allocatable :: name
      real(RK), dimension(:), pointer :: values
      real(RK), pointer :: values_surf 
      logical :: volume_grid, face_grid
   end type

   ! Definition of a FABM variable to log
   type, public :: LogVariableFABM
      character(len=100), pointer, dimension(:) :: names
      real(RK), dimension(:,:), pointer :: values
   end type

   ! Definition of a FABM variable to log
   type, public :: LogVariableFABM_bound
      character(len=100), pointer, dimension(:) :: names
      real(RK), dimension(:), pointer :: values
   end type

   ! Logging configuration
   type, public :: OutputConfig
      character(len=:), allocatable :: PathOut
      character(len=:), allocatable :: zoutName
      character(len=:), allocatable :: toutName
      character(len=:), allocatable :: output_depth_reference
      real(RK), dimension(:), allocatable :: zout, zout_read
      real(RK), dimension(:), allocatable :: tout
      integer(8), dimension(:,:), allocatable :: simulation_times_for_output
      integer, dimension(:), allocatable :: n_timesteps_between_tout
      logical :: write_to_file, output_all
      integer :: number_output_vars
      character(len=20), dimension(:), allocatable :: output_var_names ! Names of output variables
      class(LogVariable), dimension(:), allocatable :: output_vars ! Output structure for Sismtrat variables
      class(LogVariableFABM), allocatable :: output_vars_fabm_interior_state ! Output structure for FABM interior state variables
      class(LogVariableFABM), allocatable :: output_vars_fabm_bottom_state ! Output structure for FABM bottom state variables
      class(LogVariableFABM_bound), allocatable :: output_vars_fabm_surface_state ! Output structure for FABM surface state variables
      class(LogVariableFABM), allocatable :: output_vars_fabm_diagnostic_interior ! Output structure for FABM interior diagnostic variables
      class(LogVariableFABM_bound), allocatable :: output_vars_fabm_diagnostic_horizontal ! Output structure for FABM horizontal diagnostic variables

      integer :: output_time_type, output_depth_type, thinning_interval
      real(RK) :: depth_interval, thinning_interval_read ! thinning_interval_read is a real to make sure that also values
      ! like 72.0 can be read (which are interpreted as a double)
   end type

   ! Simulation configuration
   type, public :: SimConfig
      integer :: timestep
      integer :: reference_year
      real(RK) :: start_datum
      real(RK) :: end_datum
      integer :: disp_simulation
      logical :: continue_from_snapshot = .false.
      logical :: save_text_restart = .false.
      logical :: use_text_restart = .false.
      logical :: show_bar = .true.
   end type

   ! Model configuration (read from file)
   type, public :: ModelConfig
      integer :: max_length_input_data
      logical :: couple_fabm
      integer :: turbulence_model
      logical :: split_a_seiche
      integer :: stability_func
      integer :: flux_condition
      integer :: forcing_mode
      logical :: user_defined_water_albedo
      logical :: use_filtered_wind
      integer :: seiche_normalization
      integer :: wind_drag_model
      integer :: inflow_mode
      logical :: bottom_friction
      integer :: ice_model
      integer :: snow_model
   end type

   ! FABM configuration (read from file)
   type, public :: FABMConfig
      ! Directory of YAML file with all biogeochemical configuration
      character(len=:), allocatable :: fabm_config_file
      ! Path to Inflow files
      character(len=:), allocatable :: path_fabm_inflow
      ! Directory of file with names of diagnostic variables to output
      character(len=:), allocatable :: set_diag_vars
      ! Whether to output diagnostic variables
      logical :: output_diagnostic_variables
      ! Whether to clip all state variables to valid range from bgc models when update is called
      logical :: repair_fabm
      ! Whether there is a pelagic-benthic interface at every depth
      logical :: bottom_everywhere
      ! Whether to calculate the attenuation coefficient with FABM
      logical :: bioshade_feedback
   end type

   ! Model params (read from file)
   type, public :: ModelParam
      real(RK) :: Lat
      real(RK) :: p_air
      real(RK) :: a_seiche
      real(RK) :: a_seiche_w
      real(RK) :: strat_sumr
      real(RK) :: q_NN
      real(RK) :: f_wind
      real(RK) :: C10_constant
      real(RK) :: CD
      real(RK) :: fgeo
      real(RK) :: p_sw_water
      real(RK) :: p_lw
      real(RK) :: p_windf
      real(RK) :: p_absorb
      real(RK) :: beta_sol
      real(RK) :: wat_albedo
      real(RK) :: p_sw_ice
      real(RK) :: freez_temp
      real(RK) :: snow_temp
      real(RK) :: seiche_ini
      real(RK) :: w_ice_ini
      real(RK) :: b_ice_ini
      real(RK) :: snow_ini
      !real(RK) :: k_min
   end type

   ! Model state (self is actually the simulation data!!!)
   type, public :: ModelState
      ! Iteration variables
      integer :: current_year ! Current year of simulation, used for zenith angle dependent water albedo
      integer :: current_month ! Current month of simulation, used for zenith angle dependent water albedo
      real(RK) :: current_day ! Current day of simulation, used for zenith angle dependent water albedo
      real(RK) :: datum, dt
      integer(8), dimension(2) :: simulation_time, simulation_time_old
      logical :: first_timestep = .true.

      ! Variables located on z_cent grid
      ! Note that for these variables the value at 0 z.b. U(0) is not used
      real(RK), dimension(:), allocatable :: U, V ! Water velocities [m/s]
      real(RK), dimension(:), pointer :: T, S ! Temperature [°C], Salinity [‰], FABM needs pointer attribute
      real(RK), dimension(:), allocatable :: dS ! Source/sink for salinity
      real(RK), dimension(:, :), allocatable :: Q_inp ! Horizontal inflow [m^3/s]
      real(RK), dimension(:), allocatable :: Q_inp_bound, Q_inp_bound_con ! Bottom- / Surface-bound horizontal inflow [m^2/s] (absolute and concentration dependent)
      real(RK), dimension(:), pointer :: rho ! Water density [kg/m^3], FABM needs pointer attribute
      integer :: n_pH
      
      ! All FABM biogeochemical state variable values in an array *_state
      ! In Simstrat_FABM allocated with shape (grid%nz_grid, size(n_fabm_*_state))
      real(RK), dimension(:,:), pointer :: fabm_interior_state, fabm_bottom_state
      real(RK), dimension(:), pointer :: fabm_surface_state
      integer :: n_fabm_state, n_fabm_interior_state, n_fabm_bottom_state, n_fabm_surface_state, n_fabm_diagnostic, n_fabm_diagnostic_interior, n_fabm_diagnostic_horizontal
      character(len=100), dimension(:), pointer :: fabm_state_names ! Names of FABM state variables used in the simulation
      real(RK), dimension(:,:), pointer :: fabm_diagnostic_interior ! State matrix of FABM diagnostic variables
      real(RK), dimension(:), pointer :: fabm_diagnostic_horizontal ! State matrix of FABM diagnostic variables
      character(len=100), dimension(:), pointer :: fabm_diagnostic_names ! Names of FABM diagnostic variables in output
   
      ! Variables located on z_upp grid
      real(RK), dimension(:), allocatable :: k, ko ! Turbulent kinetic energy (TKE) [J/kg]
      real(RK), dimension(:), allocatable :: avh
      real(RK), dimension(:), allocatable :: eps ! TKE dissipation rate [W/kg]
      real(RK), dimension(:), allocatable :: num, nuh ! Turbulent viscosity (momentum) and diffusivity(temperature)
      real(RK), dimension(:), allocatable :: P, B ! Shear stress production [W/kg], buoyancy production [W/kg]
      real(RK), dimension(:), allocatable :: NN ! Brunt-Väisälä frequency [s-2]
      real(RK), dimension(:), allocatable :: cmue1, cmue2 ! Model constants
      real(RK), dimension(:), allocatable :: P_Seiche ! Production of TKE [W/kg] and seiche energy [J]
      real(RK) :: E_Seiche
      real(RK) :: gamma ! Proportionality constant for loss of seiche energy

      real(RK), dimension(:), allocatable :: absorb ! Absorption coeff [m-1]
      real(RK), dimension(:), pointer :: absorb_vol ! Absorption coeff on vol grid [m-1], FABM needs pointer attribute
      real(RK) :: u10, v10, Wf ! Wind speeds, wind factor
      real(RK), pointer :: uv10 ! pointer attribute needed for FABM
      real(RK) :: drag, u_taus ! Drag
      real(RK), pointer :: u_taub ! Bottom stress, FABM needs pointer attribute
      real(RK) :: tx, ty ! Shear stress
      real(RK), pointer :: C10 ! Wind drag coefficient, FABM needs pointer attribute
      real(RK) :: SST, heat, heat_snow, heat_ice, heat_snowice! Sea surface temperature and heat flux

      real(RK), pointer :: T_atm, qa ! Air temp and specific humidity at surface
      real(RK), pointer :: Cloud ! Cloud area fraction, FABM needs ponter attribute
      real(RK), dimension(:), allocatable :: rad, rad_vol ! Solar radiation (in water)
      real(RK), dimension(:), pointer :: swr_vol ! Shortwave radiation [J/s/m2] (used for FABM: needs pointer attribute)
      real(RK), dimension(:), pointer :: par_vol ! Photosynthetically active radiation (fraction of swr, used for FABM: needs pointer attribute)
      real(RK), dimension(:), allocatable :: Q_vert ! Vertical exchange between boxes
      real(RK), dimension(9,12) :: albedo_data  ! Experimental monthly albedo data for determination of current water albedo
      real(RK) :: albedo_water   ! Current water albedo
      integer :: lat_number ! Latitude band (used for determination of albedo)

      ! Snow and Ice
      real(RK), allocatable :: snow_h ! Snow layer height [m]
      real(RK), allocatable :: total_ice_h ! Total ice layer height [m]
      real(RK), allocatable :: black_ice_h ! Black ice layer height [m]
      real(RK), allocatable :: white_ice_h ! Snowice layer height [m]
      real(RK) :: snow_dens ! Snow density [kg m-3]
      real(RK) :: ice_temp ! Ice temperature [°C]
      real(RK) :: precip ! Precipiation in water eqvivalent hight [m]
      real(RK), pointer :: ice_area_fraction ! Ice area fraction, FABM needs pointer attribute

      !For saving heatflux
      real(RK), allocatable :: ha ! Incoming long wave [W m-2]
      real(RK), allocatable :: hw ! Outgoing long wave [W m-2]
      real(RK), allocatable :: hk ! Sensible flux [W m-2]
      real(RK), allocatable :: hv ! Latent heat [W m-2]
      real(RK), pointer :: rad0 !  Solar radiation at surface  [W m-2], FABM needs pointer attribute
      real(RK), pointer :: par0 ! Photosynthetically active radiation [W m-2], FABM needs pointer attribute

      real(RK) :: cde, cm0
      real(RK) ::  fsed
      real(RK), dimension(:), allocatable     :: fgeo_add

      ! Pointers to parameters for FABM (needs pointer attribute)
      real(RK), pointer :: Lat ! latitude [degree_north]
      real(RK), pointer :: p_air ! Surface air pressure [Pa]
      real(RK), pointer :: wat_albedo ! Surface albedo [-]

   contains
      procedure, pass :: init => model_state_init
      procedure, pass :: save => save_model_state
      procedure, pass :: load => load_model_state
   end type

   ! Structure that encapsulates a full program state
   type, public :: SimulationData
      type(InputConfig), public   :: input_cfg
      type(OutputConfig), public  :: output_cfg
      type(SimConfig), public     :: sim_cfg
      type(ModelConfig), public   :: model_cfg
      type(FABMConfig), public    :: fabm_cfg
      type(ModelParam), public    :: model_param
      type(ModelState), public    :: model
      type(StaggeredGrid), public :: grid
   contains
      procedure, pass :: init => simulation_data_init
   end type

contains
   subroutine simulation_data_init(self, param, state_size)
      class(SimulationData), intent(inout) :: self
      class(ModelParam), intent(in) :: param
      integer, intent(in) :: state_size
      ! Init model data structures
      call self%model%init(param, state_size)
   end subroutine

   ! Allocates all arrays of the model state in the correct size
   subroutine model_state_init(self, param, state_size)
      class(ModelState), intent(inout) :: self
      class(ModelParam), intent(in) :: param
      integer, intent(in) :: state_size

      ! Values on volume grid
      ! Important: Size is smaller than vars on upper grid.
      !            https://en.wikipedia.org/wiki/Off-by-one_error#Fencepost_error ;-)
      allocate (self%U(state_size))
      allocate (self%V(state_size))
      allocate (self%T(state_size))
      allocate (self%S(state_size))
      allocate (self%dS(state_size))
      allocate (self%rho(state_size))
      allocate (self%avh(state_size))

      ! Values on z_upp grid
      allocate (self%k(state_size + 1))
      allocate (self%ko(state_size + 1))
      allocate (self%eps(state_size + 1))
      allocate (self%num(state_size + 1))
      allocate (self%nuh(state_size + 1))
      allocate (self%P(state_size + 1))
      allocate (self%B(state_size + 1))
      allocate (self%NN(state_size + 1))
      allocate (self%cmue1(state_size + 1))
      allocate (self%cmue2(state_size + 1))
      allocate (self%P_Seiche(state_size + 1))

      allocate (self%absorb(state_size + 1))
      allocate (self%absorb_vol(state_size))
      allocate (self%rad(state_size + 1))
      allocate (self%rad_vol(state_size))
      allocate (self%swr_vol(state_size))
      allocate (self%par_vol(state_size))
      allocate (self%Q_vert(state_size + 1))

      allocate (self%snow_h)
      allocate (self%total_ice_h)
      allocate (self%black_ice_h)
      allocate (self%white_ice_h)
      allocate (self%ice_area_fraction)

      allocate (self%ha)
      allocate (self%hw)
      allocate (self%hk)
      allocate (self%hv)
      allocate (self%rad0)
      allocate (self%par0)

      ! Init to zero
      self%U = 0.0_RK
      self%V = 0.0_RK
      self%T = 0.0_RK
      self%S = 0.0_RK
      self%dS = 0.0_RK
      self%rho = 0.0_RK

      self%k = 0.0_RK
      self%ko = 0.0_RK
      self%eps = 0.0_RK
      self%num = 0.0_RK
      self%nuh = 0.0_RK
      self%P = 0.0_RK
      self%B = 0.0_RK
      self%NN = 0.0_RK
      self%cmue1 = 0.0_RK
      self%cmue2 = 0.0_RK
      self%P_Seiche = 0.0_RK
      self%E_Seiche = 0.0_RK

      self%absorb = 0.0_RK
      self%absorb_vol = 0.0_RK
      self%rad = 0.0_RK
      self%rad_vol = 0.0_RK
      self%swr_vol = 0.0_RK
      self%par_vol = 0.0_RK
      self%Q_vert = 0.0_RK

      self%snow_h = 0.0_RK
      self%total_ice_h = 0.0_RK
      self%black_ice_h = 0.0_RK
      self%white_ice_h = 0.0_RK
      self%ice_temp = 0.0_RK
      self%snow_dens = rho_s_0
      self%precip = 0.0_RK 
      self%ice_area_fraction = 0.0_RK 
   
      self%ha = 0.0_RK
      self%hw = 0.0_RK
      self%hk = 0.0_RK 
      self%hv = 0.0_RK
      self%rad0 = 0.0_RK
      self%par0 = 0.0_RK
      self%n_pH = 0

      ! init pointers
      allocate(self%uv10)
      self%uv10 = 0.0_RK
      allocate(self%u_taub)
      self%u_taub = 0.0_RK      
      allocate(self%C10)
      self%C10 = 0.0_RK
      allocate(self%T_atm)
      self%T_atm = 0.0_RK
      allocate(self%qa)
      self%qa = 0.0_RK
      allocate(self%Cloud)
      self%Cloud = 0.0_RK

      ! Get initial values from parameters
      allocate(self%Lat)
      self%Lat = param%Lat
      allocate(self%p_air)
      self%p_air = param%p_air
      allocate(self%wat_albedo)
      self%wat_albedo = param%wat_albedo

      self%simulation_time_old = 0

   end subroutine

   ! save model state unformatted
   subroutine save_model_state(self, couple_fabm, inflow_mode)
      implicit none
      class(ModelState), intent(in) :: self
      logical, intent(in) :: couple_fabm
      integer, intent(in) :: inflow_mode

      !write(80) self%current_year, self%current_month, self%current_day, self%datum
      !write(80) self%simulation_time(1), self%simulation_time(2), self%simulation_time_old(1), self%simulation_time_old(2)
      call save_array(80, self%U)
      call save_array(80, self%V)
      call save_array_pointer(80, self%T)
      call save_array_pointer(80, self%S)
      call save_array(80, self%dS)
      call save_array_pointer(80, self%rho)
      call save_array(80, self%k)
      call save_array(80, self%ko)
      call save_array(80, self%avh)
      call save_array(80, self%eps)
      call save_array(80, self%num)
      call save_array(80, self%nuh)
      call save_array(80, self%P)
      call save_array(80, self%B)
      call save_array(80, self%NN)
      call save_array(80, self%cmue1)
      call save_array(80, self%cmue2)
      call save_array(80, self%P_Seiche)
      write(80) self%E_Seiche, self%gamma
      call save_array(80, self%absorb)
      call save_array_pointer(80, self%absorb_vol)
      write(80) self%u10, self%v10, self%uv10, self%Wf
      write(80) self%u_taub, self%drag, self%u_taus
      write(80) self%tx, self%ty
      write(80) self%C10
      write(80) self%SST, self%heat, self%heat_snow, self%heat_ice, self%heat_snowice
      write(80) self%T_atm
      write(80) self%Cloud, self%qa, self%Lat, self%p_air, self%wat_albedo
      call save_array(80, self%rad)
      call save_array(80, self%rad_vol)
      call save_array_pointer(80, self%swr_vol)
      call save_array_pointer(80, self%par_vol)
      write(80) self%albedo_data
      write(80) self%albedo_water
      write(80) self%lat_number
      write(80) self%snow_h
      write(80) self%total_ice_h
      write(80) self%black_ice_h
      write(80) self%white_ice_h
      write(80) self%snow_dens
      write(80) self%ice_temp
      write(80) self%precip
      write(80) self%ice_area_fraction
      write(80) self%ha
      write(80) self%hw
      write(80) self%hk
      write(80) self%hv
      write(80) self%rad0
      write(80) self%par0
      write(80) self%cde, self%cm0
      write(80) self%fsed
      call save_array(80, self%fgeo_add)
      if (couple_fabm) then
         call save_matrix_pointer(80, self%fabm_interior_state)
         call save_matrix_pointer(80, self%fabm_bottom_state)
         call save_array_pointer(80, self%fabm_surface_state)
         call save_matrix_pointer(80, self%fabm_diagnostic_interior)
         call save_array_pointer(80, self%fabm_diagnostic_horizontal)
      end if
      if (inflow_mode > 0) then
         call save_matrix(80, self%Q_inp)
         call save_array(80, self%Q_vert)
         if (couple_fabm) then
            call save_array(80, self%Q_inp_bound)
            call save_array(80, self%Q_inp_bound_con)
         end if
      end if
   end subroutine

   ! load model state unformatted
   subroutine load_model_state(self, couple_fabm, inflow_mode)
      implicit none
      class(ModelState), intent(inout) :: self
      logical, intent(in) :: couple_fabm
      integer, intent(in) :: inflow_mode

      read(81) self%current_year, self%current_month, self%current_day, self%datum
      read(81) self%simulation_time(1), self%simulation_time(2), self%simulation_time_old(1), self%simulation_time_old(2)
      call read_array(81, self%U)
      call read_array(81, self%V)
      call read_array_pointer(81, self%T)
      call read_array_pointer(81, self%S)
      call read_array(81, self%dS)
      call read_array_pointer(81, self%rho)
      call read_array(81, self%k)
      call read_array(81, self%ko)
      call read_array(81, self%avh)
      call read_array(81, self%eps)
      call read_array(81, self%num)
      call read_array(81, self%nuh)
      call read_array(81, self%P)
      call read_array(81, self%B)
      call read_array(81, self%NN)
      call read_array(81, self%cmue1)
      call read_array(81, self%cmue2)
      call read_array(81, self%P_Seiche)
      read(81) self%E_Seiche, self%gamma
      call read_array(81, self%absorb)
      call read_array_pointer(81, self%absorb_vol)
      read(81) self%u10, self%v10, self%uv10, self%Wf
      read(81) self%u_taub, self%drag, self%u_taus
      read(81) self%tx, self%ty
      read(81) self%C10
      read(81) self%SST, self%heat, self%heat_snow, self%heat_ice, self%heat_snowice
      read(81) self%T_atm
      read(81) self%Cloud, self%qa, self%Lat, self%p_air, self%wat_albedo
      call read_array(81, self%rad)
      call read_array(81, self%rad_vol)
      call read_array_pointer(81, self%swr_vol)
      call read_array_pointer(81, self%par_vol)
      read(81) self%albedo_data
      read(81) self%albedo_water
      read(81) self%lat_number
      read(81) self%snow_h
      read(81) self%total_ice_h
      read(81) self%black_ice_h
      read(81) self%white_ice_h
      read(81) self%snow_dens
      read(81) self%ice_temp
      read(81) self%precip
      read(81) self%ice_area_fraction
      read(81) self%ha
      read(81) self%hw
      read(81) self%hk
      read(81) self%hv
      read(81) self%rad0
      read(81) self%par0
      read(81) self%cde, self%cm0
      read(81) self%fsed
      call read_array(81, self%fgeo_add)
      if (couple_fabm) then
         call read_matrix_pointer(81, self%fabm_interior_state)
         call read_matrix_pointer(81, self%fabm_bottom_state)
         call read_array_pointer(81, self%fabm_surface_state)
         call read_matrix_pointer(80, self%fabm_diagnostic_interior)
         call read_array_pointer(80, self%fabm_diagnostic_horizontal)
      end if
      if (inflow_mode > 0) then
         call read_matrix(81, self%Q_inp)
         call read_array(81, self%Q_vert)
         if (couple_fabm) then
            call read_array(81, self%Q_inp_bound)
            call read_array(81, self%Q_inp_bound_con)
         end if
      end if
   end subroutine

end module strat_simdata
