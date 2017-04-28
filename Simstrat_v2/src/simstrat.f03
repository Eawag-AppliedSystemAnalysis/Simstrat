!     +---------------------------------------------------------------+
!     |  Simstrat model for simulation of                             |
!     |  vertical transport in lakes and reservoirs                   |
!     +---------------------------------------------------------------+

program simstrat_main
  use strat_kinds
  use simstrat_inputfile_module, only : SimstratSimulationFactory
  use strat_simdata, only : SimulationData
  use strat_windshear, only: WindShearModule
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(SimstratSimulationFactory) :: factory
  class(SimulationData), pointer :: simdata

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
 !simdata%model_cfg%max_nr_grid_cells = 2016
 write(*,*) simdata%model_cfg

!  call run_simulation()

contains

  subroutine run_simulation()

    class(WindShearModule), allocatable :: uv
    type(SimulationData) :: sim_data

    write(*,*) "Hallo welt"

    !forcing%update()
    !stability%update()
    !advection%update()

    ! Update and solve U and V - terms
    call uv%update()
    call uv%solve()

    ! Update and solve t - terms
    ! call t%update()
    ! call t%solve()
    ! and so on

    !s%update_and_solve()
    !turbulence%update()
    !k%update_and_solve()
    !eps%update_and_solve()
    !dissipation%update()

    !model%log()

  end subroutine

end program simstrat_main
