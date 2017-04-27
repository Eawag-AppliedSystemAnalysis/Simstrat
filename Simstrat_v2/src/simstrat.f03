!     +---------------------------------------------------------------+
!     |  Simstrat model for simulation of                             |
!     |  vertical transport in lakes and reservoirs                   |
!     +---------------------------------------------------------------+

program simstrat_main
  use strat_kinds
  use strat_simdata, only : Simulation
  use strat_windshear, only: WindShearModule
  use, intrinsic :: ieee_arithmetic

  implicit none

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

  call run_simulation()

contains

  subroutine run_simulation()

    class(WindShearModule), allocatable :: uv
    class(Simulation), allocatable :: sim_data

    allocate(sim_data)
    allocate(sim_data%model)

    sim_data%model%step = 26

    write(*,*) sim_data%model%step
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
