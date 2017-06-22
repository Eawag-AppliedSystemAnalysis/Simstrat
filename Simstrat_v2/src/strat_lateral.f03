module strat_lateral
  use strat_kinds
  use strat_simdata
  use strat_consts
  use strat_grid
  use utilities
  implicit none
  private

type, abstract, public :: GenericLateralModule
  class(ModelConfig), pointer :: cfg
  class(StaggeredGrid), pointer :: grid
  class(ModelParam), pointer :: param

  ! Variables that where either marked with "save" before, or that have been
  ! global, but only used in the lateral environment:
  real(RK), dimension(:,:), allocatable   :: z_Inp, Q_start, Q_end, depth_surfaceFlow
  real(RK) :: tb_start(1:4), tb_end(1:4)        ! Input depths, start time, end time
  integer :: nz_old, eof(1:4)
  integer :: nval(1:4), nval_deep(1:4), nval_surface(1:4)      ! Number of values

contains
  procedure, pass :: init => lateral_generic_init
  procedure(lateral_generic_update), deferred, pass :: update
  procedure, pass :: surface_flow => lateral_generic_surface_flow
end type

type, extends(GenericLateralModule), public :: LateralRhoModule
contains
  procedure, pass, public :: update => lateral_rho_update
end type

type, extends(GenericLateralModule), public:: LateralModule
contains
  procedure, pass, public :: update => lateral_update
end type

contains
  subroutine lateral_generic_update(self, state)
    implicit none
    class(GenericLateralModule) :: self
    class(ModelState) :: state
  end subroutine

  subroutine lateral_generic_init(self, model_config, model_param, grid)
    implicit none
    class(GenericLateralModule) :: self
    class(StaggeredGrid), target :: grid
    class(ModelConfig), target :: model_config
    class(ModelParam), target :: model_param

    self%cfg => model_config
    self%param => model_param
    self%grid => grid
  end subroutine

  subroutine lateral_generic_surface_flow(self, Q_surface, Q_total, depth_surface,i)
        implicit none
        class(GenericLateralModule) :: self
        ! Global declarations
        real(RK), intent(in) :: Q_surface, depth_surface
        real(RK), intent(inout) :: Q_total(1:)
        integer, intent(in) :: i

        ! Local variables
        real(RK) :: Q_surf_integral, add_surf
        integer :: j

        Q_surf_integral = Q_surface

        do j=1,self%grid%ubnd_fce

          if (((Q_surf_integral>0) .and. (.not. i==2)) .or. ((Q_surf_integral<0) .and. (i==2))) then
            Q_total(self%grid%ubnd_fce+1-j) = Q_total(self%grid%ubnd_fce+1-j) + Q_surf_integral
          else
            exit
          end if

          add_surf = Q_surface*self%grid%h(self%grid%ubnd_fce+1-j)/depth_surface
          Q_surf_integral = Q_surf_integral - add_surf
        end do
        return
  end subroutine

  subroutine lateral_rho_update(self, state)
    implicit none
    class(LateralRhoModule) :: self
    class(ModelState) :: state
  end subroutine

  subroutine lateral_update(self, state)
    implicit none
    class(LateralModule) :: self
    class(ModelState) :: state

    associate(datum => state%datum, &
              idx => state%std, &
              Q_inp => state%Q_inp, &     ! Q_inp is the input at each depth for each time step, Q_vert is the integrated net water input
              Q_vert => state%Q_vert)


  end associate
  end subroutine

end module
