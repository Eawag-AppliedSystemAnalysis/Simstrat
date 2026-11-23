!     +---------------------------------------------------------------+
!     |  Stub of the module simstrat_aed2 for Simstrat without AED2
!     +---------------------------------------------------------------+
!     Used by the FoBiS build mode release-gnu-noaed2 instead of simstrat_aed2.f90 (issue #55):
!     Simstrat then compiles without libaed2. If CoupleAED2 is true in the par-file, the simulation stops.

module simstrat_aed2
   use strat_simdata
   use strat_grid
   use utilities

   implicit none
   private

   type, public :: SimstratAED2
   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
      procedure, pass(self), public :: deallocate => aed2_deallocate
   end type SimstratAED2

contains

   subroutine init(self, state, grid, model_cfg, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(AED2Config), target :: aed2_cfg

      call error('Simstrat was compiled without AED2 (build mode release-gnu-noaed2). Set ModelConfig.CoupleAED2 to false.')
   end subroutine

   subroutine update(self, state)
      implicit none
      class(SimstratAED2) :: self
      class(ModelState) :: state
   end subroutine

   subroutine aed2_deallocate(self)
      implicit none
      class(SimstratAED2) :: self
   end subroutine

end module simstrat_aed2
