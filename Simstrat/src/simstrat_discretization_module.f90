module simstrat_discretization_module
  use simstrat_kinds
  implicit none

  !############################################################################
  !# Abstract Types
  !############################################################################

  !#####################################################
  type, abstract, public :: SimstratDiscretizationScheme
  !#####################################################
    ! Subroutines a discretization-scheme / grid for simstrat has to provide
    integer, public ::                              nz_mfq, nz_tq                    !Number of (active) grid cells
    contains
      procedure(generic_init), deferred, pass(self), public :: initialize                             !initialize the discretization scheme
      procedure(generic_precondition), deferred, pass(self), public :: precondition                   !pre-calculate constants that are needed for the calculation of the LES
      procedure(generic_interpolate), deferred, pass(self), public :: interpolateMFQ
      procedure(generic_export), deferred, pass(self), public :: exportMFQ
      procedure(generic_export), deferred, pass(self), public :: exportKEPS
      procedure(generic_create_les), deferred, pass(self), public :: createLinearEquationSystem_MFQ   !create diagonals of the linear equation system for mean flow quantities
      procedure(generic_create_les), deferred, pass(self), public :: createLinearEquationSystem_k     !create diagonals of LES for k
      procedure(generic_create_les), deferred, pass(self), public :: createLinearEquationSystem_eps   !create diagonals of LES for eps
      procedure(generic_volume), deferred, pass(self), public :: totalvolume
  end type

  ! generic subroutine definitions of a discretization scheme and stencil for simstrat 
  abstract interface
    subroutine generic_init(self,&                      !reference to SimstratDiscretizationScheme object
                            nz_max,&                    !maximum of memory that should be allocated (grid might grow/shrink due to water level changes)
                            z_faces,&                   !volume faces in height above sediment
                            z_area,&                    !positions where lake bathymetry is defined
                            area)                       !areas/lake bathymetry
      import SimstratDiscretizationScheme, RK
      implicit none

      class(SimstratDiscretizationScheme), intent(inout) :: self
      integer, intent(in) :: nz_max
      real(RK), dimension(:), intent(in) :: z_faces
      real(RK), dimension(:), intent(in) :: z_area, area

    end subroutine

    subroutine generic_precondition(self,&
                                  dt)                 !constant time step
      import SimstratDiscretizationScheme, RK
      implicit none

      class(SimstratDiscretizationScheme), intent(inout) :: self
      real(RK), intent(in) :: dt

    end subroutine

    pure function generic_interpolate(self, z_inp, var_inp) result(var)
      import SimstratDiscretizationScheme, RK
      implicit none

      class(SimstratDiscretizationScheme), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z_inp, var_inp
      real(RK), dimension(:), allocatable :: var

    end function

    pure function generic_export(self, z_exp, var_exp) result(var)
      import SimstratDiscretizationScheme, RK
      implicit none

      class(SimstratDiscretizationScheme), intent(in) :: self
      real(RK), dimension(:), intent(in) :: z_exp, var_exp
      real(RK), dimension(:), allocatable :: var

    end function

    subroutine generic_create_les(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
      import SimstratDiscretizationScheme, RK
      implicit none

      class(SimstratDiscretizationScheme), intent(in) :: self
      real(RK), dimension(:), intent(in) :: var, sources, boundary_cond
      real(RK), dimension(:), intent(in) :: nu, fluxes
      real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

    end subroutine

    pure function generic_volume(self) result(volume)
      import SimstratDiscretizationScheme, RK
      implicit none

      class(SimstratDiscretizationScheme), intent(in) :: self
      real(RK) :: volume
    end function
  end interface

end module simstrat_discretization_module
