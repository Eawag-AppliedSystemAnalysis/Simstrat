module strat_windshear
  use strat_kinds
  implicit none

  private

  type, extends(StateVariable), public :: WindShearModule
     contains
       procedure, nopass, public :: update => windshear_update
       procedure, nopass, public :: solve => windshear_solve
   end type
contains

  subroutine windshear_solve()
      implicit none

      write(*,*) "Windshear SOLVE!"
  end subroutine

  subroutine windshear_update()
      implicit none

      write(*,*) "Windshear UPDATE!"
  end subroutine

end module strat_windshear
