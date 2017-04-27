module strat_simdata
  use strat_kinds
  implicit none
  private

  ! Common Types
  type, public :: InputConfig
  end type

  type, public :: OutputConfig
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
    real(RK), dimension(:), allocatable :: P_Seiche,E_Seiche    ! Production of TKE [W/kg] and seiche energy [J]
    real(RK), dimension(:), allocatable :: gamma     ! Proportionality constant for loss of seiche energy

    real(RK), dimension(:), allocatable :: absorb    ! Absorption coeff [m-1]
    real(RK), dimension(:), allocatable :: u10, v10, uv10, Wf   ! Wind speeds, wind factor
    real(RK), dimension(:), allocatable :: u_taub,drag,u_taus   ! Drag
    real(RK), dimension(:), allocatable :: tx,ty     ! Shear stress
    real(RK), dimension(:), allocatable :: SST, heat ! Sea surface temperature and heat flux
    real(RK), dimension(:), allocatable :: rad0,rad  ! Solar radiation (at surface and in water)
    real(RK), dimension(:), allocatable :: Q_vert    ! Vertical exchange between boxes
  end type

  type, public :: Simulation
        class(ModelState), allocatable, public :: model
  end type



end module strat_simdata
