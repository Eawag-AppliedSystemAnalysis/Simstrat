module simstrat_finite_volume_implementation
  use simstrat_kinds
  use utilities
  use simstrat_discretization_module
  implicit none

  private

  !##################################################################################
  !# Implementation of discretization scheme / grid
  !##################################################################################

  type, abstract, extends(SimstratDiscretizationScheme), public :: StaggeredFiniteVolumeDiscretization
    private

      !definition of the finite volumes
      real(RK), dimension(:), allocatable, public ::  z_centres, z_faces,&  !position of the volume centres and volume faces
                                                      h_centres, h_faces,&  !hight of the volumes, distance between the volumes
                                                      a_centres, a_faces,&  !area at the volume centres, area at the volume faces
                                                      a_extended,&
                                                      dAdz, dAdz_faces

      !precondition
      real(RK), dimension(:), allocatable, public ::  volumes, volumes_faces, form_base, form_1, form_2, form_base_k, form_k1, form_k2, form_eps
      real(RK), public :: dt
    contains
      procedure, pass(self), public :: initialize => initialize_staggered_finite_volume_discretization
      procedure, pass(self), public :: interpolateMFQ => interpolate_on_volume_centres
      procedure, pass(self), public :: exportMFQ => export_from_volume_centres
      procedure, pass(self), public :: exportKEPS => export_from_volume_faces
      procedure, pass(self), public :: interpolateOnVolumeCentres => interpolate_on_volume_centres
      procedure, pass(self), public :: interpolateOnVolumeFaces => interpolate_on_volume_faces
      procedure, pass(self), public :: precondition => precondition_staggered_finite_volume_discretization
      procedure, pass(self), public :: totalvolume => get_total_volume
  end type StaggeredFiniteVolumeDiscretization

  type, extends(StaggeredFiniteVolumeDiscretization), public :: StaggeredFVImplicitEulerScheme
    contains
      procedure, pass(self), public :: createLinearEquationSystem_MFQ => les_MFQ_ie
      procedure, pass(self), public :: createLinearEquationSystem_k => les_k_ie
      procedure, pass(self), public :: createLinearEquationSystem_eps => les_eps_ie
  end type StaggeredFVImplicitEulerScheme

  type, extends(StaggeredFiniteVolumeDiscretization), public :: StaggeredFVCrankNicolsonScheme
    contains
      procedure, pass(self), public :: createLinearEquationSystem_MFQ => les_MFQ_cn
      procedure, pass(self), public :: createLinearEquationSystem_k => les_k_cn
      procedure, pass(self), public :: createLinearEquationSystem_eps => les_eps_cn
  end type StaggeredFVCrankNicolsonScheme

contains


  !############################################################################
  !# Initialize
  !############################################################################
  subroutine initialize_staggered_finite_volume_discretization( self,&
                                                                nz_max,&
                                                                z_faces,&
                                                                z_area,&
                                                                area)
    implicit none

    !arguments
    class(StaggeredFiniteVolumeDiscretization), intent(inout) :: self
    integer, intent(in) :: nz_max
    real(RK), dimension(:), intent(in) :: z_faces
    real(RK), dimension(:), intent(in) :: z_area, area
    real(RK) :: test

    !local variables
    integer :: nz, nz_faces

    nz_faces = size(z_faces)
    nz = nz_faces-1

    !allocate memory
    allocate(self%z_centres(nz_max), self%z_faces(nz_max), &
             self%h_centres(nz_max), self%h_faces(nz_max), &
             self%a_centres(nz_max), self%a_faces(nz_max),self%a_extended(nz_max),&
             self%dAdz(nz_max),self%dAdz_faces(nz_max))

    allocate(self%volumes(nz_max), self%volumes_faces(nz_max),&
             self%form_base(nz_max), self%form_base_k(nz_max),&
             self%form_1(nz_max), self%form_2(nz_max),&
             self%form_k1(nz_max), self%form_k2(nz_max),&
             self%form_eps(nz_max))

    !initialize discretization
    self%nz_mfq = nz                                                                             ! number of volumes
    self%nz_tq = nz+1

    self%z_faces(1:nz_faces) = z_faces
    self%z_centres(1:nz)     = 0.5_RK*(self%z_faces(1:nz_faces-1)+self%z_faces(2:nz_faces))  ! position of the volume centres

    self%h_centres(1:nz)     = self%z_faces(2:nz_faces)-self%z_faces(1:nz_faces-1)           ! distance between volume centres

    self%h_faces(2:nz) = self%z_centres(2:nz)-self%z_centres(1:nz-1)                   ! height of the volumes (distance between volume faces)
    self%h_faces(1) = self%h_faces(2)
    self%h_faces(nz_faces) = self%h_faces(nz)

    self%a_faces(1:nz_faces) = self%interpolateOnVolumeFaces(z_area, area)                   ! Interpolated areas at the positions of the volume faces
    self%a_centres(1:nz)     = (self%a_faces(1:nz_faces-1)+self%a_faces(2:nz_faces))/2.0_RK  ! area at the volume centres
    self%a_extended(1) = self%a_faces(1)
    self%a_extended(nz+2) = self%a_faces(nz)
    self%a_extended(2:nz+1) = self%a_centres(1:nz)

    self%dAdz(1:nz)          = (self%a_faces(2:nz_faces)-self%a_faces(1:nz_faces-1))/self%h_centres(1:nz) !gradient at volume centres

    self%dAdz_faces(2:nz) = (self%a_centres(2:nz)-self%a_centres(1:nz-1))/self%h_faces(2:nz)!approxiamte dAdz_faces
    self%dAdz_faces(1) = self%dAdz_faces(2)
    self%dAdz_faces(nz+1) = self%dAdz_faces(nz)
  end subroutine initialize_staggered_finite_volume_discretization

  pure function get_total_volume(self) result(volume)
    implicit none

    class(StaggeredFiniteVolumeDiscretization), intent(in) :: self
    real(RK) :: volume
    volume = sum(self%h_centres(1:self%nz_mfq)*self%a_centres(1:self%nz_mfq))
  end function

  !############################################################################
  !# Grid interpolation
  !############################################################################
  pure function interpolate_on_volume_centres(self, z_inp, var_inp) result(var)
    implicit none

    class(StaggeredFiniteVolumeDiscretization), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z_inp, var_inp

    real(RK), dimension(:), allocatable :: var
    allocate(var(self%nz_mfq))

    call Interp(z_inp, var_inp,size(z_inp),self%z_centres,var,self%nz_mfq)
  end function interpolate_on_volume_centres

  pure function interpolate_on_volume_faces(self, z_inp, var_inp) result(var)
    implicit none

    class(StaggeredFiniteVolumeDiscretization), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z_inp, var_inp

    real(RK), dimension(:), allocatable :: var
    allocate(var(self%nz_tq))

    call Interp(z_inp, var_inp,size(z_inp),self%z_faces,var,self%nz_tq)
  end function interpolate_on_volume_faces

  pure function export_from_volume_centres(self, z_exp, var_exp) result(var)
    implicit none

    class(StaggeredFiniteVolumeDiscretization), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z_exp, var_exp

    real(RK), dimension(:), allocatable :: var
    allocate(var(size(z_exp)))

    call Interp(self%z_centres(1:self%nz_mfq), var_exp(1:self%nz_mfq), self%nz_mfq, z_exp,var, size(z_exp))
  end function

  pure function export_from_volume_faces(self, z_exp, var_exp) result(var)
    implicit none

    class(StaggeredFiniteVolumeDiscretization), intent(in) :: self
    real(RK), dimension(:), intent(in) :: z_exp, var_exp

    real(RK), dimension(:), allocatable :: var
    allocate(var(size(z_exp)))

    call Interp(self%z_faces(1:self%nz_tq), var_exp(1:self%nz_tq), self%nz_tq, z_exp,var, size(z_exp))
  end function

  !############################################################################
  !# Precondition
  !############################################################################
  subroutine precondition_staggered_finite_volume_discretization(self, dt)
    implicit none

    class(StaggeredFiniteVolumeDiscretization), intent(inout) :: self
    real(RK), intent(in) :: dt

    ! To increase readability
    associate(nz => self%nz_mfq, nz_faces => self%nz_tq, a_faces => self%a_faces, a_centres => self%a_centres, a_extended=>self%a_extended, h_centres => self%h_centres, h_faces=>self%h_faces, form_base => self%form_base, form_base_k=>self%form_base_k)

      self%volumes(1:nz) = a_centres(1:nz)*h_centres(1:nz)
      self%volumes_faces(1:nz_faces) = a_faces(1:nz_faces)*h_faces(1:nz_faces)
      self%dt = dt

      ! constants used to build the tridiagonals
      form_base(1:nz-1) = -4.0_RK*dt*a_faces(2:nz_faces-1)/(h_centres(2:nz)+h_centres(1:nz-1))

      self%form_1(1) = 0.0_RK !no-flux boundary
      self%form_1(2:nz) = form_base(1:nz-1)/(a_faces(3:nz_faces)+a_faces(2:nz_faces-1))/h_centres(2:nz)

      self%form_2(nz) = 0.0_RK !no-flux boundary
      self%form_2(1:nz-1) = form_base(1:nz-1)/(a_faces(2:nz_faces-1)+a_faces(1:nz_faces-2))/h_centres(1:nz-1)

      form_base_k(1:nz_faces-1) = -2.0_RK*dt*a_extended(2:nz_faces)/h_centres(1:nz)

      self%form_k1(1) = 0.0_RK
      self%form_k1(2:nz_faces) = form_base_k(1:nz_faces-1)/(a_extended(3:nz_faces+1)+a_extended(2:nz_faces))/h_faces(2:nz_faces)

      self%form_k2(nz_faces) = 0.0_RK
      self%form_k2(1:nz_faces-1) = form_base_k(1:nz_faces-1)/(a_extended(2:nz_faces)+a_extended(1:nz_faces-1))/h_faces(1:nz_faces-1)

      self%form_eps(1:nz_faces) = (a_extended(1:nz_faces)-a_extended(2:nz_faces+1))/h_faces(1:nz_faces)/a_extended(2:nz_faces+1)
    end associate

  end subroutine

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals
  !############################################################################
  ! Core routine: construct the basic form of the diagonals that can be evolved to implicit euler or crank-nicolson scheme
  subroutine staggered_finite_volume_discretization_les_core(dt, form_param1, form_param2, volumes, a_faces, nz_centres, nu, fluxes, sources, ld, md, ud, rhs)
    implicit none

    !arguments
    integer, intent(in) :: nz_centres
    real(RK), intent(in) :: dt
    real(RK), dimension(:), intent(in) :: form_param1, form_param2, volumes
    real(RK), dimension(:), intent(in) :: sources !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes, a_faces !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

    !local variables
    integer :: nz_faces
    real(RK), dimension(nz_centres) :: nu1, nu2, Q
    real(RK), dimension(nz_centres+1) :: f

    nz_faces = nz_centres + 1
    nu1 = nu(1:nz_faces-1)
    nu2 = nu(2:nz_faces)

    !calculate fluxes from specific fluxes
    f(1:nz_faces) = fluxes(1:nz_faces)*a_faces(1:nz_faces)
    Q(1:nz_centres) = (f(1:nz_faces-1)-f(2:nz_faces))/volumes(1:nz_centres)

    !# upper and lower diagonals are both the same for implicit euler
    !# and crank-nicolson
    !MP: Todo: Check with NU is correct!
    ld(1:nz_centres) = form_param1(1:nz_centres)*nu2(1:nz_centres)
    ud(1:nz_centres) = form_param2(1:nz_centres)*nu2(1:nz_centres)
    !# right hand side
    !MP: Correction: Q should not be multiplied with dt
    rhs(1:nz_centres) = dt*(sources(1:nz_centres))+Q(1:nz_centres)

    return
  end subroutine staggered_finite_volume_discretization_les_core

  !############################################################################
  !# Implicit Euler
  !############################################################################

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals for mean flow quantities
  !############################################################################
  subroutine les_MFQ_ie(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
    implicit none

    class(StaggeredFVImplicitEulerScheme), intent(in) :: self
    real(RK), dimension(:), intent(in) :: var, sources, boundary_cond !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

    call staggered_finite_volume_discretization_les_core(self%dt, self%form_1, self%form_2, self%volumes, self%a_faces, self%nz_mfq, nu, fluxes, sources, ld, md, ud, rhs)

    md(1:self%nz_mfq) = (1.0_RK - ld(1:self%nz_mfq) - ud(1:self%nz_mfq)) + boundary_cond(1:self%nz_mfq)*self%dt
    rhs(1:self%nz_mfq) = rhs(1:self%nz_mfq) + var(1:self%nz_mfq)

    return
  end subroutine les_MFQ_ie

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals for turbulent quantity k
  !############################################################################
  subroutine les_k_ie(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
    implicit none

    class(StaggeredFVImplicitEulerScheme), intent(in) :: self
    real(RK), dimension(:), intent(in) :: var, sources, boundary_cond !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

    call staggered_finite_volume_discretization_les_core(self%dt, self%form_k1, self%form_k2, self%volumes_faces, self%a_extended, self%nz_tq, nu, fluxes, sources, ld, md, ud, rhs)

    md(1:self%nz_tq) = (1.0_RK - ld(1:self%nz_tq) - ud(1:self%nz_tq)) + boundary_cond(1:self%nz_tq)*self%dt
    rhs(1:self%nz_tq) = rhs(1:self%nz_tq) + var(1:self%nz_tq)
  end subroutine les_k_ie

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals for eps
  !############################################################################
  subroutine les_eps_ie(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
    implicit none

    class(StaggeredFVImplicitEulerScheme), intent(in) :: self
    real(RK), dimension(:), intent(in) :: var, sources, boundary_cond !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

    associate(nz_tq=>self%nz_tq)

      call staggered_finite_volume_discretization_les_core(self%dt, self%form_k1, self%form_k2, self%volumes_faces, self%a_extended, nz_tq, nu, fluxes, sources, ld, md, ud, rhs)

      md(1:nz_tq) = (1.0_RK - ld(1:nz_tq) - ud(1:nz_tq)) + boundary_cond(1:nz_tq)*self%dt
      rhs(1:nz_tq) = rhs(1:nz_tq) + var(1:nz_tq)
    end associate
  end subroutine les_eps_ie

  !############################################################################
  !# Crank Nicolson
  !############################################################################

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals for mean flow quantities
  !############################################################################
  subroutine les_MFQ_cn(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
    implicit none

    class(StaggeredFVCrankNicolsonScheme), intent(in) :: self
    real(RK), dimension(:), intent(in) :: var, sources, boundary_cond !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

    !call staggered_finite_volume_discretization_les_core(self%dt, self%form_1, self%form_2, self%volumes, self%a_faces, self%nz_mfq, nu, fluxes, sources, ld, md, ud, rhs)

    !!Implicit Part
    !md(1:self%nz_mfq) = (1.0_RK - 0.5_RK*(ld(1:self%nz_mfq) + ud(1:self%nz_mfq)))
    !!Explicit Part
    !!...
    !!RHS
    !rhs(1:self%nz_mfq) = rhs(1:self%nz_mfq) + var(1:self%nz_mfq)

    return
  end subroutine les_MFQ_cn

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals for turbulent quantity k
  !############################################################################
  pure subroutine les_k_cn(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
    implicit none

    class(StaggeredFVCrankNicolsonScheme), intent(in) :: self
    real(RK), dimension(:), intent(in) :: var, sources, boundary_cond !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

  end subroutine les_k_cn

  !############################################################################
  !# Create Linear Equation System / Matrix Diagonals for eps
  !############################################################################
  pure subroutine les_eps_cn(self, var, nu, boundary_cond, fluxes, sources, ld, md, ud, rhs)
    implicit none

    class(StaggeredFVCrankNicolsonScheme), intent(in) :: self
    real(RK), dimension(:), intent(in) :: var, sources, boundary_cond !centres
    real(RK), dimension(:), intent(in) :: nu, fluxes !faces
    real(RK), dimension(:), intent(inout) :: ld, md, ud, rhs

  end subroutine les_eps_cn
end module simstrat_finite_volume_implementation
