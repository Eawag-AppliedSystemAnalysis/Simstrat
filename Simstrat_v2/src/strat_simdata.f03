module strat_simdata
  use strat_kinds
  implicit none
  private

  ! Common Types
  type, public :: InputConfig
    character(len=:), allocatable          ::MorphName
    character(len=:), allocatable          :: InitName
    character(len=:), allocatable          :: ForcingName
    character(len=:), allocatable          :: AbsorpName
    character(len=:), allocatable          :: GridName
    character(len=:), allocatable          :: QinpName
    character(len=:), allocatable          :: QoutName
    character(len=:), allocatable          :: TinpName
    character(len=:), allocatable          :: SinpName
  end type

  type, public :: OutputConfig
    character(len=:), allocatable :: PathOut
    character(len=:), allocatable :: zoutName
    logical :: write_on_the_fly
    integer :: thinning_interval
  end type

  type, public :: SimConfig
    integer :: timestep
    integer :: start_datum
    integer :: end_datum
  end type

  type, public :: ModelConfig
    integer :: max_nr_grid_cells
    logical :: couple_aed2
    logical :: use_buffered_forcing
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
    integer :: disp_diagnostic
    integer :: data_averaging
  end type

  type, public :: ModelParam
    real(RK) :: Lat
    real(RK) :: p_air
    real(RK) :: a_seiche
    real(RK) :: q_NN
    real(RK) :: f_wind
    real(RK) :: C10
    real(RK) :: CD
    real(RK) :: fgeo
    real(RK) :: k_min
    real(RK) :: p_radin
    real(RK) :: p_windf
    real(RK) :: beta_sol
    real(RK) :: albsw
  end type

  type, public :: ModelState
    ! Iteration variables
    integer :: step,itera,i,j,std
    real(RK) :: datum

    ! Variables located on z_cent grid
    ! Note that for these variables the value at 0 z.b. U(0) is not used
    real(RK), dimension(:), allocatable :: U,V       ! Water velocities [m/s]
    real(RK), dimension(:), allocatable :: T,S       ! Temperature [°C], Salinity [‰]
    real(RK), dimension(:), allocatable :: dS        ! Source/sink for salinity
    real(RK), dimension(:), allocatable :: Q_inp     ! Horizontal inflow [m^3/s]
    real(RK), dimension(:), allocatable :: rho       ! Water density [kg/m^3]

    ! Variables located on z_upp grid
    real(RK), dimension(:), allocatable :: k,ko      ! Turbulent kinetic energy (TKE) [J/kg]
    real(RK), dimension(:), allocatable :: eps       ! TKE dissipation rate [W/kg]
    real(RK), dimension(:), allocatable :: num,nuh   ! Turbulent viscosity (momentum) and diffusivity (temperature)
    real(RK), dimension(:), allocatable :: P,B       ! Shear stress production [W/kg], buoyancy production [W/kg]
    real(RK), dimension(:), allocatable :: NN        ! Brunt-Väisälä frequency [s-2]
    real(RK), dimension(:), allocatable :: cmue1,cmue2 ! Model constants
    real(RK), dimension(:), allocatable :: P_Seiche    ! Production of TKE [W/kg] and seiche energy [J]
    real(RK) :: E_Seiche
    real(RK) :: gamma     ! Proportionality constant for loss of seiche energy

    real(RK), dimension(:), allocatable :: absorb    ! Absorption coeff [m-1]
    real(RK) :: u10, v10, uv10, Wf   ! Wind speeds, wind factor
    real(RK) :: u_taub,drag,u_taus   ! Drag
    real(RK) :: tx,ty     ! Shear stress
    real(RK) :: SST, heat ! Sea surface temperature and heat flux
    real(RK) :: rad0 ! Solar radiation at surface
    real(RK), dimension(:), allocatable :: rad  ! Solar radiation (in water)
    real(RK), dimension(:), allocatable :: Q_vert    ! Vertical exchange between boxes

  contains
    procedure, pass :: init => model_state_init
  end type

  type, public :: SimulationData
        type(InputConfig), public   :: input_cfg
        type(OutputConfig), public  :: output_cfg
        type(SimConfig), public     :: sim_cfg
        type(ModelConfig), public   :: model_cfg
        type(ModelParam),  public   :: model_param
        type(ModelState),  public   :: model

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

      ! Values on z_cent grid
      ! Important: Size is smaller than vars on upper grid.
      !            https://en.wikipedia.org/wiki/Off-by-one_error#Fencepost_error ;-)
      allocate(this%U(state_size - 1))
      allocate(this%V(state_size - 1))
      allocate(this%T(state_size - 1))
      allocate(this%S(state_size - 1))
      allocate(this%dS(state_size - 1))
      allocate(this%Q_inp(state_size - 1))
      allocate(this%rho(state_size - 1))

      ! Values on z_upp grid
      allocate(this%k(state_size))
      allocate(this%ko(state_size))
      allocate(this%eps(state_size))
      allocate(this%num(state_size))
      allocate(this%nuh(state_size))
      allocate(this%P(state_size))
      allocate(this%B(state_size))
      allocate(this%NN(state_size))
      allocate(this%cmue1(state_size))
      allocate(this%cmue2(state_size))
      allocate(this%P_Seiche(state_size))

      allocate(this%absorb(state_size))
      allocate(this%rad(state_size))
      allocate(this%Q_vert(state_size))

      write(*,*) size(this%T)

    end subroutine




end module strat_simdata
