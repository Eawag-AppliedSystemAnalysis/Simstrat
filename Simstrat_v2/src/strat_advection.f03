module strat_advection
  use strat_kinds
  use strat_simdata
  use strat_consts
  use strat_grid
  use utilities
  implicit none
  private

type, public :: AdvectionModule
  class(ModelConfig), pointer :: cfg
  class(StaggeredGrid), pointer :: grid
  class(ModelParam), pointer :: param


contains
  procedure, pass :: init => advection_init
  procedure, pass :: update => advection_update
end type

contains

  subroutine advection_init(self, model_config, model_param, grid)
    implicit none
    class(AdvectionModule) :: self
    class(StaggeredGrid), target :: grid
    class(ModelConfig), target :: model_config
    class(ModelParam), target :: model_param

    self%cfg => model_config
    self%param => model_param
    self%grid => grid

  end subroutine


  subroutine advection_update(self, state)
    implicit none
    class(AdvectionModule) :: self
    class(ModelState) :: state

    real(RK) :: top_z, top_h, curr_depth
    real(RK) :: top
    real(RK) :: dh, dh_i(1:2), h_div_2, h_mult_2      ! depth differences
    real(RK) :: dt_i(1:2)          ! first and second time step
    real(RK) :: AreaFactor_adv(1:self%grid%nz_grid_max)
    integer :: i, t_i

    associate(grid => self%grid, nz_occupied => self%grid%nz_occupied, dt => state%dt)
    !todo: grid%lake_level_old = z_upp(nz)

    !Depth difference compared to previous timestep
    top_z = grid%z_face(grid%nz_occupied)
    top_h = grid%h(grid%nz_occupied)
    curr_depth = grid%depth

    dh = state%Q_vert(grid%nz_occupied)/grid%Az(grid%nz_occupied)*state%dt
    h_div_2 = 0.5_RK*grid%h(grid%nz_occupied-1)    ! Take second highest box since the top box might not be at the full height
    h_mult_2 = 2_RK*grid%h(grid%nz_occupied-1)

    ! Calculate timestep splitting
    !Split timestep depending on situation
    if (dh==0.) then                          ! If volume does not change, take one normal time step
        dt_i(1) = dt
    else if ((dh+top_z)>=curr_depth) then        ! If surface level reached, take a step until surface
        dt_i(1) = (curr_depth - top_z)/dh*dt
    else if (((dh+top_h)>h_div_2) .and.&   ! If top box>0.5*lower box and <2*lower box, take one time step
            ((dh+top_h)<h_mult_2)) then
        dt_i(1) = dt
    else if ((dh+top_h)<=h_div_2) then     ! If top box<=0.5*lower box, first step until top box=0.5*lower box
        dt_i(1) = abs((top_h-h_div_2)/dh)*dt
    else                                     ! If top box>=2*lower box, first step until top box = 2*lower box
        dt_i(1) = abs((2*grid%h(grid%nz_occupied-1)-top_h)/dh)*dt
    end if
    dt_i(2) = dt-dt_i(1)                       ! Rest of timestep



    ! FB 2016: Revision
   do t_i=1,2 !First and (if needed) second timestep
       AreaFactor_adv(1:nz_occupied) = dt_i(t_i)/(self%grid%Az(1:nz_occupied)*self%grid%h(1:nz_occupied))     ! Area factor for dt(t_i)
       dh_i(t_i) = dh*dt_i(t_i)/dt                         ! Depth difference for dt(t_i)


   end do      !end do t_i=1,2
  end associate
  end subroutine

  subroutine advection_merge_box(self, state, dh)
    implicit none
    class(AdvectionModule) :: self
    class(ModelState) :: state
    real(RK) :: dh
    real(RK) :: w_a,w_b,k_norm
    associate(ubnd_fce => self%grid%ubnd_fce, ubnd_vol => self%grid%ubnd_vol)

   ! New values of the state variables are weighted averages
    !determine weighting an normalization connstant
    w_a = 0.5_RK*self%grid%Az(ubnd_fce)
    w_b = self%grid%Az(ubnd_fce- 1)

    ! shrink grid by one (this also updates ubnd_fce/vol)
    call self%grid%shrink(dh)

    ! update quantities in new top box (based on former top box and current value)
    state%U(ubnd_vol) = (w_a * state%U(ubnd_vol+1) + w_b * state%U(ubnd_vol)) / (w_a+w_b)
    state%V(ubnd_vol) = (w_a * state%V(ubnd_vol+1) + w_b * state%V(ubnd_vol)) / (w_a+w_b)
    state%T(ubnd_vol) = (w_a * state%T(ubnd_vol+1) + w_b * state%T(ubnd_vol)) / (w_a+w_b)
    state%S(ubnd_vol) = (w_a * state%S(ubnd_vol+1) + w_b * state%S(ubnd_vol)) / (w_a+w_b)

    state%k(ubnd_fce) = (w_a * state%k(ubnd_fce+1) + w_b * state%k(ubnd_fce)) / (w_a+w_b)
    state%eps(ubnd_fce) = (w_a * state%eps(ubnd_fce+1) + w_b * state%eps(ubnd_fce)) / (w_a+w_b)
    state%Q_vert(ubnd_fce) = (w_a * state%Q_vert(ubnd_fce+1) + w_b * state%Q_vert(ubnd_fce)) / (w_a+w_b)

    ! update area factors
    call self%grid%update_area_factors()

  end  associate
  end subroutine

  subroutine advection_add_box(self, state, dh)
    implicit none
    class(AdvectionModule) :: self
    class(ModelState) :: state
    real(RK) :: dh
    associate(ubnd_fce => self%grid%ubnd_fce, ubnd_vol => self%grid%ubnd_vol)

    ! extend grid by one (also updates ubnd_vol etc)
    call self%grid%grow(dh)


    ! Update quantities in new grid element
    state%U(ubnd_vol)   = state%U(ubnd_vol - 1)
    state%V(ubnd_vol)   = state%V(ubnd_vol - 1)
    state%T(ubnd_vol)   = state%T(ubnd_vol - 1)
    state%S(ubnd_vol)   = state%S(ubnd_vol - 1)
    state%Q_vert(ubnd_fce) = state%Q_vert(ubnd_fce -1 )       ! Vertical discharge of new box

    state%k(ubnd_fce)   = state%k(ubnd_fce - 1)
    state%eps(ubnd_fce) = state%eps(ubnd_fce - 1)

    call self%grid%update_area_factors()

  end  associate
  end subroutine

end module
