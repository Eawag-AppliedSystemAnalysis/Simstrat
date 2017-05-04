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
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(SimstratSimulationFactory) :: factory
  class(SimulationData), pointer :: simdata
  type(ForcingModule) :: mod_forcing
  type(StabilityModule) :: mod_stability
  type(SimpleLogger) :: logger

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

  !initialize forcing module
  call mod_forcing%init(simdata%model_cfg, simdata%input_cfg%ForcingName, simdata%grid)

  ! Setup logger
  call logger%initialize(simdata%output_cfg, simdata%grid)

  ! initialize simulation modules
  call mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)

  call run_simulation()

  call logger%close()

contains

  subroutine run_simulation()


    ! Read forcing file
    call mod_forcing%update(simdata%model, simdata%model_param)

    call mod_stability%update(simdata%model)
    !advection%update()

    ! Update and solve U and V - terms
  !  call uv%update()
  !  call uv%solve()

    ! Update and solve t - terms
    ! call t%update()
    ! call t%solve()
    ! and so on

    !s%update_and_solve()
    !turbulence%update()
    !k%update_and_solve()
    !eps%update_and_solve()
    !dissipation%update()

    call logger%log(simdata%model%datum)


  end subroutine

end program simstrat_main
