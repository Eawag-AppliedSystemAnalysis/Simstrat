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
  use strat_windshear, only: WindShearModule
  use strat_statevar
  use strat_temp
  use strat_solver
  use strat_discretization
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(SimstratSimulationFactory) :: factory
  class(SimulationData), pointer :: simdata
  type(ThomasAlgSolver) :: solver
  type(EulerIDiscretization) :: euler_i_disc
  type(ForcingModule) :: mod_forcing
  type(StabilityModule) :: mod_stability
  type(SimpleLogger) :: logger
  type(TempModelVar) :: mod_temperature

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

  !initialize forcing module
  call mod_forcing%init(simdata%model_cfg, &
                        simdata%model_param, &
                        simdata%input_cfg%ForcingName, &
                        simdata%grid)

  ! Setup logger
  call logger%initialize(simdata%output_cfg, simdata%grid)

  ! initialize simulation modules
  call mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)
  call mod_temperature%init(simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%T)

  call run_simulation()

  call logger%close()

contains

  subroutine run_simulation()

    call logger%log(0.0_RK) ! Write initial conditions

    ! todo: time control!
    simdata%model%dt = 0.5

    ! Read forcing file
    call mod_forcing%update(simdata%model)

    ! Update physics
    call mod_stability%update(simdata%model)
    !advection%update()
    call mod_forcing%update_corriolis(simdata%model)

    ! Update and solve U and V - terms
  !  call uv%update()
  !  call uv%solve()

    ! Update and solve t - terms
    call mod_temperature%update(simdata%model, simdata%model_param)

    !s%update_and_solve()
    !turbulence%update()
    !k%update_and_solve()
    !eps%update_and_solve()
    !dissipation%update()

    write(*,*) simdata%model%T

    call logger%log(simdata%model%datum)


  end subroutine

end program simstrat_main
