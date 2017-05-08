module strat_turbulence
  use strat_kinds
  use strat_consts
  use strat_grid
  use strat_simdata
  implicit none
  private

  ! Common Types
  type, public :: TurbulenceModule
    class(StaggeredGrid), pointer :: grid
    class(ModelConfig), pointer :: model_cfg
    class(ModelParam), pointer :: model_param

  contains
    procedure, pass :: init => turbulence_module_init
    procedure, pass :: update => turbulence_module_update

    procedure, pass :: do_production => turbulence_module_do_production

  end type
contains

  subroutine turbulence_module_init(self, grid, model_cfg, model_param)
    implicit none
    class(TurbulenceModule) :: self
    class(StaggeredGrid), target :: grid
    class(ModelConfig), target :: model_cfg
    class(ModelParam), target :: model_param

    self%grid => grid
    self%model_cfg => model_cfg
    self%model_param => model_param
  end subroutine

  subroutine turbulence_module_update(self, state)
    implicit none
    class(TurbulenceModule) :: self
    class(ModelState) :: state
    real(RK), dimension(self%grid%l_fce) :: beta

    call self%do_production(state)

  end subroutine

  subroutine turbulence_module_do_production(self, state)
  !####################################################################
      implicit none
      class(TurbulenceModule) :: self
      class(ModelState) :: state

      associate(grid=>self%grid, &
                ubnd_vol => self%grid%ubnd_vol, &
                ubnd_fce =>self%grid%ubnd_fce)

      ! Equation 5 (left) of Goudsmit, 2002
      ! P is defined on the inner faces
      state%P = 0
      state%P(2:ubnd_fce-1) = (state%U(2:ubnd_vol)-state%U(1:ubnd_vol-1))**2+(state%V(2:ubnd_vol)-state%V(1:ubnd_vol-1))**2
      state%P(2:ubnd_fce-1) = state%P(2:ubnd_fce-1)*state%num(2:ubnd_fce-1)*grid%meanint(2:ubnd_fce-1)**2

      ! Equation 5 (right) of Goudsmit, 2002
      state%B = 0
      state%B(2:ubnd_fce-1) = -state%nuh(2:ubnd_fce-1)*state%NN(2:ubnd_fce-1)

      return
    end associate
  end subroutine

end module
