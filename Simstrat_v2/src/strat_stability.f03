module strat_stability
  use strat_kinds
  use strat_consts
  use strat_grid
  use strat_simdata
  implicit none
  private

  ! Common Types
  type, public :: StabilityModule
    class(StaggeredGrid), pointer :: grid
    class(ModelConfig), pointer :: model_cfg
    class(ModelParam), pointer :: model_param

  contains
    procedure, pass :: init => stability_module_init
    procedure, pass :: update => stability_module_update
    procedure, pass :: update_cmue_cn => stability_module_update_cmue_cn
    procedure, pass :: update_cmue_qe => stability_module_update_cmue_qe
    !procedure, pass :: update_NN_no_sal => stability_module_update_NN_no_sal
    !procedure, pass :: update_NN_no_sal_grad => stability_module_update_NN_no_sal_grad
    procedure, pass :: update_NN => stability_module_update_NN

  end type
contains

  subroutine stability_module_init(self, grid, model_cfg, model_param)
    implicit none
    class(StabilityModule) :: self
    class(StaggeredGrid), target :: grid
    class(ModelConfig), target :: model_cfg
    class(ModelParam), target :: model_param

    self%grid => grid
    self%model_cfg => model_cfg
    self%model_param => model_param
  end subroutine

  subroutine stability_module_update(self, state)
    implicit none
    class(StabilityModule) :: self
    class(ModelState) :: state
    real(RK), dimension(self%grid%l_fce) :: beta

    !Do buoyancy update (update NN)
    if(state%has_salinity) then
        if(state%has_salinity_grad) then
          call self%update_NN(state%T, state%S, state%rho, state%NN)
        else
          write(*,*) "!!!! BUOYANCY NOT IMPLEMENTED FOR THIS CASE"
      !    call self%update_NN_no_sal_grad(state%T, state%S, state%rho, state%NN)
        end if
    else ! salinity zero everywhere
      write(*,*) "!!!! BUOYANCY NOT IMPLEMENTED FOR THIS CASE"
      !call self%update_NN_no_sal(state%T, state%S, state%rho, state%NN)
    end if

    !update cmue depending on selected stabilty function
    if(self%model_cfg%stability_func == 1) then
      call self%update_cmue_cn(state%cmue1, state%cmue2)
    else if(self%model_cfg%stability_func == 2) then

      beta = state%NN*(state%k/state%eps)**2
      call self%update_cmue_qe(beta, state%cmue1, state%cmue2, state%cde)
    end if

  end subroutine

  subroutine stability_module_update_NN(self,T,S,rho,NN)
    implicit none
    class(StabilityModule) :: self
    ! Global variables
    real(RK), dimension(:), intent(in) :: T, S
    real(RK), dimension(:), intent(inout) :: NN, rho

    ! Local variables
    real(RK) :: a(self%grid%l_fce),buoy(self%grid%l_fce)
    real(RK) :: rho0t(self%grid%l_fce),rho0st(self%grid%l_fce)
    integer :: i

    associate(grd => self%grid)

    do i=1,grd%ubnd_fce-1
        rho0t(i)= 0.9998395_RK+T(i)*(6.7914e-5_RK+T(i)*(-9.0894e-6_RK+T(i)*&
            (1.0171e-7_RK+T(i)*(-1.2846e-9_RK+T(i)*(1.1592e-11_RK+T(i)*(-5.0125e-14_RK))))))
        rho0st(i)= (8.181e-4_RK+T(i)*(-3.85e-6_RK+T(i)*(4.96e-8_RK)))*S(i)
        rho(i)= rho_0*(rho0t(i)+rho0st(i))
        !if (fc/=0) then
        !    rho(i) = rho0t(i)*(1-fc)+fc*rho(i)
        !end
        buoy(i)= -g*(rho(i)-rho_0)/rho_0
    end do

    NN(2:grd%ubnd_fce-1) = grd%meanint(1:grd%ubnd_vol-1)*(buoy(2:grd%ubnd_fce-1)-buoy(1:grd%ubnd_fce-2))
    NN(1)= NN(2)
    NN(grd%ubnd_fce)= NN(grd%ubnd_fce-1)
    
  end associate
  end subroutine



  subroutine stability_module_update_cmue_cn(self, cmue1, cmue2)
    implicit none

    ! Global variables
    class(StabilityModule) :: self
    real(RK), dimension(:), intent(inout) :: cmue1,cmue2

    !Standard version of k-eps model
    cmue1 = cmue
    cmue2 = cmue/Prndtl

    ! set boundaries
    cmue1(1) = cmue1(2)
    cmue2(1) = cmue2(2)
    cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce-1)
    cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce-1)
  end subroutine

  subroutine stability_module_update_cmue_qe(self, beta,  cmue1, cmue2, cde)
    implicit none
    class(StabilityModule) :: self
    real(RK), dimension(:), intent(inout) ::beta,cmue1,cmue2
    real(RK) :: cde
    real(RK) :: gh,sm,sh
    integer :: i

    do i=2,self%grid%ubnd_fce-1
      gh = -cde**2*0.5_RK*beta(i)
      if(gh> 0.02) gh = gh-(gh-0.02_RK)**2/(gh+0.0233_RK-2*0.02_RK)
      if(gh<-0.28) gh = -0.28_RK

      sm = 1.0_RK-3*c1-6*a1/b1-3*a2*gh*((b2-3*a2)*(1.0_RK-6*a1/b1)-3*c1*(b2+6*a1))
      sm = a1*sm/((1.0_RK-3*a2*gh*(6*a1+b2))*(1.0_RK-9*a1*a2*gh))
      sh = a2*(1.0_RK-6*a1/b1)/(1.0_RK-3*a2*gh*(6*a1+b2))

      cmue1(i) = sqrt(2.0_RK)*cde*sm
      cmue2(i) = sqrt(2.0_RK)*cde*sh
    end do

    ! set boundaries
    cmue1(1) = cmue1(2)
    cmue2(1) = cmue2(2)
    cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce-1)
    cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce-1)
  end subroutine

end module
