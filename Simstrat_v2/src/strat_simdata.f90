!     +---------------------------------------------------------------+
!     |  Data structure definitions for simulation data
!     +---------------------------------------------------------------+

module strat_simdata
   use strat_kinds
   use strat_grid
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

   ! Logging configuration
   type, public :: OutputConfig
      character(len=:), allocatable :: PathOut
      character(len=:), allocatable :: zoutName
      character(len=:), allocatable :: toutName
      real(RK), dimension(:), allocatable :: zout
      real(RK), dimension(:), allocatable :: tout
      real(RK), dimension(:), allocatable :: n_timesteps_between_tout
      real(RK), dimension(:), allocatable :: adjusted_timestep
      logical :: write_to_file
      class(LogVariable), dimension(:), allocatable :: output_vars

      integer :: output_time_type, output_depth_type, thinning_interval
      real(RK) :: depth_interval
   end type

   ! Simulation configuration
   type, public :: SimConfig
      integer :: timestep
      real(RK) :: start_datum
      real(RK) :: end_datum
   end type

   ! Model configuration (read from file)
   type, public :: ModelConfig
      integer :: max_length_input_data
      logical :: couple_aed2
      integer :: turbulence_model
      integer :: stability_func
      integer :: flux_condition
      integer :: forcing_mode
      logical :: use_filtered_wind
      integer :: seiche_normalization
      integer :: wind_drag_model
      integer :: inflow_placement
      integer :: pressure_gradients
      logical :: salinity_transport
      integer :: disp_simulation
      integer :: ice_model
      integer :: snow_model
   end type

   ! Model params (read from file)
   type, public :: ModelParam
      real(RK) :: Lat
      real(RK) :: p_air
      real(RK) :: a_seiche
      real(RK) :: q_NN
      real(RK) :: f_wind
      real(RK) :: C10_constant
      real(RK) :: CD
      real(RK) :: fgeo
      real(RK) :: k_min
      real(RK) :: p_radin
      real(RK) :: p_windf
      real(RK) :: beta_sol
      real(RK) :: beta_snow_ice
      real(RK) :: albsw
      real(RK) :: ice_albedo ! Ice albedo
      real(RK) :: snow_albedo ! Snow albedo   
      real(RK) :: freez_temp ! Water freez point [kg m-3]     
   end type

   ! Model state (this is actually the simulation data!!!)
   type, public :: ModelState
      ! Iteration variables
      integer :: i, j, output_counter, model_step_counter
      real(RK) :: datum, dt

      ! Variables located on z_cent grid
      ! Note that for these variables the value at 0 z.b. U(0) is not used
      real(RK), dimension(:), allocatable :: U, V ! Water velocities [m/s]
      real(RK), dimension(:), allocatable :: T, S ! Temperature [°C], Salinity [‰]
      real(RK), dimension(:), allocatable :: dS ! Source/sink for salinity
      real(RK), dimension(:, :), allocatable :: Q_inp ! Horizontal inflow [m^3/s]
      real(RK), dimension(:), allocatable :: rho ! Water density [kg/m^3]
   
      ! Variables located on z_upp grid
      real(RK), dimension(:), allocatable :: k, ko ! Turbulent kinetic energy (TKE) [J/kg]
      real(RK), dimension(:), allocatable :: avh
      real(RK), dimension(:), allocatable :: eps ! TKE dissipation rate [W/kg]
      real(RK), dimension(:), allocatable :: num, nuh ! Turbulent viscosity (momentum) and diffusivity (temperature)
      real(RK), dimension(:), allocatable :: P, B ! Shear stress production [W/kg], buoyancy production [W/kg]
      real(RK), dimension(:), allocatable :: NN ! Brunt-Väisälä frequency [s-2]
      real(RK), dimension(:), allocatable :: cmue1, cmue2 ! Model constants
      real(RK), dimension(:), allocatable :: P_Seiche ! Production of TKE [W/kg] and seiche energy [J]
      real(RK) :: E_Seiche
      real(RK) :: gamma ! Proportionality constant for loss of seiche energy

      real(RK), dimension(:), allocatable :: absorb ! Absorption coeff [m-1]
      real(RK) :: u10, v10, uv10, Wf ! Wind speeds, wind factor
      real(RK) :: u_taub, drag, u_taus ! Drag
      real(RK) :: tx, ty ! Shear stress
      real(RK) :: C10 ! Wind drag coefficient
      real(RK) :: SST, heat , heat_snow_ice! Sea surface temperature and heat flux
      !real(RK) :: rad0 ! Solar radiation at surface
      real(RK) :: T_atm ! Air temp at surface   
      real(RK), dimension(:), allocatable :: rad ! Solar radiation (in water)
      real(RK), dimension(:), allocatable :: Q_vert ! Vertical exchange between boxes

      ! Snow and Ice
      real(RK), allocatable :: snow_h ! Snow layer height [m]
      real(RK), allocatable :: ice_h ! Ice layer height [m]   
      real(RK) :: snow_dens ! On ice snow density [kg m-3]   
      real(RK) :: ice_temp ! Ice density [kg m-3]
      real(RK) :: snow_temp ! Ice density [kg m-3]
      real(RK) :: precip ! precipiation in water eqvivalent hight [m] 
   
      !For saving heatflux 
      real(RK), allocatable :: ha ! Incoming long wave [W m-2]
      real(RK), allocatable :: hw ! Outgoing long wave [W m-2]
      real(RK), allocatable :: hk ! Sensible flux [W m-2]
      real(RK), allocatable :: hv ! Latent heat [W m-2]
      real(RK), allocatable :: rad0 !  Solar radiation at surface  [W m-2]
   
      real(RK) :: cde, cm0
      real(RK) ::  fsed
      real(RK), dimension(:), allocatable     :: fgeo_add
      logical :: has_advection
      integer :: nz_input
      logical :: has_salinity_grad, has_salinity

   contains
      procedure, pass :: init => model_state_init
   end type

   ! Structure that encapsulates a full program state
   type, public :: SimulationData
      type(InputConfig), public   :: input_cfg
      type(OutputConfig), public  :: output_cfg
      type(SimConfig), public     :: sim_cfg
      type(ModelConfig), public   :: model_cfg
      type(ModelParam), public    :: model_param
      type(ModelState), public    :: model
      type(StaggeredGrid), public :: grid
   contains
      procedure, pass :: init => simulation_data_init
   end type

contains
   subroutine simulation_data_init(this, state_size)
      class(SimulationData), intent(inout) :: this
      integer, intent(in) :: state_size
      ! init model data structures
      call this%model%init(state_size)

   end subroutine

   ! Allocates all arrays of the model state in the correct size
   subroutine model_state_init(this, state_size)
      class(ModelState), intent(inout) :: this
      integer, intent(in) :: state_size

      ! Values on volume grid
      ! Important: Size is smaller than vars on upper grid.
      !            https://en.wikipedia.org/wiki/Off-by-one_error#Fencepost_error ;-)
      allocate (this%U(state_size))
      allocate (this%V(state_size))
      allocate (this%T(state_size))
      allocate (this%S(state_size))
      allocate (this%dS(state_size))
      allocate (this%Q_inp(1:4, state_size + 1))
      allocate (this%rho(state_size))
      allocate (this%avh(state_size))

      ! Values on z_upp grid
      allocate (this%k(state_size + 1))
      allocate (this%ko(state_size + 1))
      allocate (this%eps(state_size + 1))
      allocate (this%num(state_size + 1))
      allocate (this%nuh(state_size + 1))
      allocate (this%P(state_size + 1))
      allocate (this%B(state_size + 1))
      allocate (this%NN(state_size + 1))
      allocate (this%cmue1(state_size + 1))
      allocate (this%cmue2(state_size + 1))
      allocate (this%P_Seiche(state_size + 1))

      allocate (this%absorb(state_size + 1))
      allocate (this%rad(state_size + 1))
      allocate (this%Q_vert(state_size + 1))

      allocate (this%snow_h) 
      allocate (this%ice_h)  
   
      allocate (this%ha) 
      allocate (this%hw) 
      allocate (this%hk) 
      allocate (this%hv) 
      allocate (this%rad0) 
   
      ! init to zero
      this%U = 0.0_RK
      this%V = 0.0_RK
      this%T = 0.0_RK
      this%S = 0.0_RK
      this%dS = 0.0_RK
      this%Q_inp = 0.0_RK
      this%rho = 0.0_RK

      this%k = 0.0_RK
      this%ko = 0.0_RK
      this%eps = 0.0_RK
      this%num = 0.0_RK
      this%nuh = 0.0_RK
      this%P = 0.0_RK
      this%B = 0.0_RK
      this%NN = 0.0_RK
      this%cmue1 = 0.0_RK
      this%cmue2 = 0.0_RK
      this%P_Seiche = 0.0_RK

      this%absorb = 0.0_RK
      this%rad = 0.0_RK
      this%Q_vert = 0.0_RK
    
      this%snow_h = 0.0_RK
      this%ice_h = 0.0_RK
      this%ice_temp = 0.0_RK 
      this%snow_temp = 0.0_RK
   
      this%ha = 0.0_RK
      this%hw = 0.0_RK
      this%hk = 0.0_RK 
      this%hv = 0.0_RK 
      this%rad0 = 0.0_RK   
   
   end subroutine

end module strat_simdata
