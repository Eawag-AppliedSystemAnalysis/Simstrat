module simstrat_kinds
  implicit none
  private

  integer, parameter, public :: RK = kind(0.d0) !Real kind
  character(len=3), parameter, public :: version = '2.0'

end module simstrat_kinds
