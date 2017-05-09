!     +---------------------------------------------------------------+
!     |  Simstrat model for simulation of                             |
!     |  vertical transport in lakes and reservoirs                   |
!     +---------------------------------------------------------------+

program simstrat_main
  use strat_kinds
  use strat_inputfile, only : SimstratSimulationFactory
  use strat_outputfile
  use strat_simdata, only : SimulationData
  use strat_forcing
  use utilities
  use strat_stability, only : StabilityModule
  use strat_windshear
  use strat_statevar
  use strat_temp
  use strat_solver
  use strat_discretization
  use strat_keps
  use strat_turbulence
  use strat_transport
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(SimstratSimulationFactory) :: factory
  class(SimulationData), pointer :: simdata
  type(ThomasAlgSolver) :: solver
  type(EulerIDiscretizationMFQ) :: euler_i_disc
  type(EulerIDiscretizationKEPS) :: euler_i_disc_keps
  type(ForcingModule) :: mod_forcing
  type(StabilityModule) :: mod_stability
  type(SimpleLogger) :: logger
  type(TempModelVar) :: mod_temperature
  type(UVModelVar) :: mod_u, mod_v
  type(KModelVar) :: mod_k
  type(TransportModVar) :: mod_s
  type(TurbulenceModule) :: mod_turbulence

  character(len=100) :: arg
  character(len=:), allocatable :: ParName

  !print some information
  write(*,*) 'Simstrat version '//version
  write(*,*) 'This software has been developed at eawag - Swiss Federal Institute of Aquatic Science and Technology'
  write(*,*) ''

  !get first cli argument
  call get_command_argument(1,arg)
  ParName = trim(arg)
  if(ParName=='') ParName='simstrat.par'

  !initialize model from inputfiles
  call factory%initialize_model(ParName, simdata)

  ! initialize Discretization
  call euler_i_disc%init(simdata%grid)
  call euler_i_disc_keps%init(simdata%grid)

  !initialize forcing module
  call mod_forcing%init(simdata%model_cfg, &
                        simdata%model_param, &
                        simdata%input_cfg%ForcingName, &
                        simdata%grid)

  ! Setup logger
  call logger%initialize(simdata%output_cfg, simdata%grid)

  ! initialize simulation modules
  call mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)
  call mod_turbulence%init(simdata%grid, simdata%model_cfg, simdata%model_param)

  call mod_temperature%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%T)

  call mod_u%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%U)
  call mod_u%assign_shear_stress(simdata%model%tx)

  call mod_v%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%V)
  call mod_v%assign_shear_stress(simdata%model%ty)

  call mod_s%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%S)
  call mod_s%assign_external_source(simdata%model%dS)

  call mod_k%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%K)


  call run_simulation()

  call logger%close()

contains

  subroutine run_simulation()
    integer :: i
    call logger%log(0.0_RK) ! Write initial conditions

    ! todo: time control!
    simdata%model%dS(simdata%grid%ubnd_vol) = 0.1
    simdata%model%dt = 0.5

    do i=1,50
    ! Read forcing file
    call mod_forcing%update(simdata%model)

    ! Update physics
    call mod_stability%update(simdata%model)
    !advection%update()
    call mod_forcing%update_coriolis(simdata%model)

    ! Update and solve U and V - terms
    call mod_u%update(simdata%model, simdata%model_param)
    call mod_v%update(simdata%model, simdata%model_param)

    ! Update and solve t - terms
    write(*,*) "TTTTTTTTTTTTTTTTTT"
    call mod_temperature%update(simdata%model, simdata%model_param)

    ! Update and solve transportation terms (here: Salinity S only)
    !call mod_S%update(simdata%model, simdata%model_param)

    ! update turbulence states
    call mod_turbulence%update(simdata%model, simdata%model_param)

    ! Solve k & eps
  !  call mod_k%update(simdata%model, simdata%model_param)
    !k%update_and_solve()
    !eps%update_and_solve()

  !  call mod_turbulence%update_post_eps(simdata%model)

    call logger%log(simdata%model%datum)

    simdata%model%datum = simdata%model%datum + simdata%model%dt
  end do

  end subroutine

end program simstrat_main
