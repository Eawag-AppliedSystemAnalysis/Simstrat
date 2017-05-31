module Coriolis
  use SimstratModel
  implicit none

  private
  public doCoriolis

contains
  !####################################################################
  subroutine doCoriolis(U,V)
  !####################################################################

      implicit none

      ! Global variables
	    real(dp), intent(inout) :: U(0:),V(0:)

      ! Local variables
	    real(dp) :: U_old(0:nz)

      ! Calculation
      U_old(1:nz) = U(1:nz)
      U(1:nz) =  U(1:nz)*cos(Cori*dt) + V(1:nz)*sin(Cori*dt)
      V(1:nz) =-U_old(1:nz)*sin(Cori*dt) + V(1:nz)*cos(Cori*dt)

      return
  end subroutine doCoriolis

end module Coriolis
