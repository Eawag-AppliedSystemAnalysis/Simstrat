!     +---------------------------------------------------------------+
!     |  Interface and implementation for discretization methods
!     |  - Currently, there are two impelmentations:
!     |  - EulerIDiscretizationMFQ: Euler implicit discretization for Mean Flow Quantities
!     |  - EulerIDiscretizationKEPS: Euler Implicit discretization for k and epsilon
!     +---------------------------------------------------------------+

module strat_discretization
   use strat_kinds
   use strat_simdata
   use strat_grid

   implicit none
   private

   ! Interface
   type, abstract, public :: Discretization
      class(StaggeredGrid), pointer :: grid
   contains
      procedure, pass(self), public :: init => generic_init
      procedure(generic_create_LES), deferred, pass(self), public ::create_LES
   end type

   ! Implementation for Mean flow Quantities
   type, extends(Discretization), public :: EulerIDiscretizationMFQ
   contains
      procedure, pass, public :: create_LES => euleri_create_LES_MFQ
   end type

   ! Implementation for K and eps
   type, extends(Discretization), public :: EulerIDiscretizationKEPS
   contains
      procedure, pass, public :: create_LES => euleri_create_LES_KEPS
   end type

contains
   subroutine generic_init(self, grid)
      class(Discretization), intent(inout) :: self
      type(StaggeredGrid), target :: grid
      self%grid => grid
   end subroutine

   ! Abstract method definition for creating a linear equation system (LES)
   subroutine generic_create_LES(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(Discretization), intent(inout) :: self
      real(RK), dimension(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
   end subroutine

   ! Implementation for Mean Flow Quantities
   subroutine euleri_create_LES_MFQ(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(EulerIDiscretizationMFQ), intent(inout) :: self
      real(RK), dimension(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
      integer :: n
      n = size(main_diag)

      ! Build diagonals
      upper_diag(1) = 0.0_RK
      upper_diag(2:n) = dt*nu(1:n - 1)*self%grid%AreaFactor_1(2:n)
      lower_diag(1:n - 1) = dt*nu(1:n - 1)*self%grid%AreaFactor_2(1:n - 1)
      lower_diag(n) = 0.0_RK
      main_diag(1:n) = 1.0_RK - upper_diag(1:n) - lower_diag(1:n) + boundaries(1:n)*dt

      ! Calculate RHS
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      rhs(1:n) = var(1:n) + dt*sources(1:n)
   end subroutine

! Implementation for K and EPS (note different usage of nu)
   subroutine euleri_create_LES_KEPS(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(EulerIDiscretizationKEPS), intent(inout) :: self
      real(RK), dimension(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
      integer :: n
      n = size(main_diag)

      ! Build diagonals
      upper_diag(1) = 0
      upper_diag(n) = 0
      lower_diag(1) = 0
      lower_diag(n) = 0
      upper_diag(2:n - 1) = dt*nu(1:n - 2)*self%grid%AreaFactor_k1(1:n - 2) !todo:check indices of nu!
      lower_diag(2:n - 1) = dt*nu(2:n - 1)*self%grid%AreaFactor_k2(1:n - 2)
      main_diag(1:n) = 1.0_RK - upper_diag(1:n) - lower_diag(1:n) + boundaries(1:n)*dt

      ! Calculate RHS
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      rhs(1:n) = var(1:n) + dt*sources(1:n)
      rhs(1) = 0
      rhs(n) = 0

   end subroutine

end module strat_discretization
