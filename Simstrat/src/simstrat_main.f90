!     +---------------------------------------------------------------+
!     |  Simstrat model for simulation of                             |
!     |  vertical transport in lakes and reservoirs                   |
!     +---------------------------------------------------------------+

program simstrat_main
  use simstrat_kinds
  use simstrat_inputfile_module
  use simstrat_model_module
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(SimstratModelFactory) :: factory
  type(SimstratModel), pointer :: model
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
  model => factory%initialize_simstrat_model_from_inputfiles(ParName)

  !run the model
  call model%run

end program simstrat_main
