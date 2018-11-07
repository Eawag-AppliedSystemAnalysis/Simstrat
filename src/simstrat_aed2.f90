!     +---------------------------------------------------------------+
!     |  Simstrat - AED2 interface
!     +---------------------------------------------------------------+

module simstrat_aed2
   use strat_simdata

   implicit none
   private

   type, public :: SimstratAED2
      class(AED2Config), pointer :: cfg

   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
   end type SimstratAED2

contains

   subroutine init(self, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg

   end subroutine


   subroutine update(self, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg


   end subroutine



end module simstrat_aed2