module strat_keps
use strat_kinds
use strat_consts
use strat_simdata
use strat_statevar
use strat_grid
use strat_solver
implicit none
private

type, extends(ModelVariable), public :: KModelVar
contains
  procedure, pass(self), public :: calc_terms => k_var_calc_terms
  procedure, pass(self), public :: post_solve => k_var_post_solve
end type

type, extends(ModelVariable), public :: EpsModelVar
contains
  procedure, pass(self), public :: calc_terms => eps_var_calc_terms
  procedure, pass(self), public :: post_solve => eps_var_post_solve
end type


contains
  subroutine k_var_calc_terms(self, state, param, sources, boundaries)
    class(KModelVar), intent(inout) :: self
    class(ModelState), intent(inout) :: state
    class(ModelParam), intent(inout) :: param
    real(RK),  dimension(:) ::  sources, boundaries
    real(RK) :: rhs_0, rhs_ubnd
    real(RK) :: pminus(2:self%grid%ubnd_fce-1), pplus(2:self%grid%ubnd_fce-1), Prod, Buoy, Diss

    integer :: i
    associate(grid => self%grid, &
              ubnd_fce => self%grid%ubnd_fce, &
              ubnd_vol => self%grid%ubnd_vol)


    !!!!!!!! Precalculations !!!!!!!!
    sources = 0
    boundaries = 0
    state%ko(1:ubnd_fce) = state%k(1:ubnd_fce) ! ko = TKE at old time step

    ! indices of avh not clear! do they need to be changed to 3:ubnd_fcd-2 ??
    state%avh(2:ubnd_fce-1) = 0.5_RK/sig_k*(state%num(1:ubnd_fce-1)+state%num(2:ubnd_fce)) ! average num for TKE

    if (self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
        state%avh(1) = 0.0_RK
        state%avh(ubnd_fce) = 0.0_RK
    else
        state%avh(1)=2*state%u_taub**4/(state%eps(1)+state%eps(2))        ! = 0 for no shear stress
        state%avh(ubnd_fce)=2*state%u_taus**4/(state%eps(ubnd_fce)+state%eps(ubnd_fce-1))   ! = 0 for no shear stress
    end if

    do i=2,ubnd_fce-1
        Prod = state%P(i)+state%P_Seiche(i)                   ! Add seiche energy
        Buoy=state%B(i)
        Diss=state%eps(i)
        if (Prod+Buoy>0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
        else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
        end if
    end do

    !!!!!!!! Define sources !!!!!!!!
    sources(2:ubnd_fce-1) = pplus(2:ubnd_fce-1)

    !!!!!! Define boundary conditions !!!!
    boundaries(2:ubnd_fce-1) = pminus(2:ubnd_fce-1)/state%k(2:ubnd_fce-1)

    if(self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
      ! K(0) and K(ubnd_fce) are assigned in the post processing function
    else ! no fluxes, unity A-matrix + condition on RHS
      rhs_0    = state%u_taub**2/sqrt(state%cm0*state%cde)
      rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

      ! Trick to have rhs(1) = rhs_0 and rhs(ubnd_fce) = rhs_ubnd
      ! in discretization, rhs is calculated as rhs(1) = var(1) + sources(1)*dt
      ! Given the equation below, the following results:
      ! rhs(1) = var(1) + (-var(1)/dt + rhs_0/dt)*dt = var(1) -var(1) + rhs_0 = rhs_0
      sources(1) = -(self%var(1)/state%dt)+(rhs_0/state%dt)
      sources(ubnd_fce) =  -(self%var(ubnd_fce)/state%dt)+(rhs_ubnd/state%dt)
    end if
  end associate
  end subroutine

  subroutine k_var_post_solve(self, state, param)
    class(KModelVar), intent(inout) :: self
    class(ModelState), intent(inout) :: state
    class(ModelParam), intent(inout) :: param
    integer :: i
    associate(grid => self%grid, &
              ubnd_fce => self%grid%ubnd_fce, &
              ubnd_vol => self%grid%ubnd_vol)

    if(self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
    ! Define TKE at boundary (no flux)
      self%var(1) = self%var(2)
      self%var(self%grid%ubnd_fce) = self%var(self%grid%ubnd_fce - 1)
    end if

    ! check lower limit of k
    do i=1,ubnd_fce-1
        if(self%var(i)<k_min) self%var(i)=k_min             ! Lower limit of TKE
    end do

    end associate
  end subroutine





  subroutine eps_var_calc_terms(self, state, param, sources, boundaries)
    class(EpsModelVar), intent(inout) :: self
    class(ModelState), intent(inout) :: state
    class(ModelParam), intent(inout) :: param
    real(RK),  dimension(:) ::  sources, boundaries
    real(RK) :: cee3, rhs_0, rhs_ubnd
    real(RK) :: flux(1:self%grid%ubnd_fce)
    real(RK) :: pminus(2:self%grid%ubnd_fce-1), pplus(2:self%grid%ubnd_fce-1), Prod, Buoy, Diss

    integer :: i
    associate(grid => self%grid, &
              ubnd_fce => self%grid%ubnd_fce, &
              ubnd_vol => self%grid%ubnd_vol)


    !!!!!!!! Precalculations !!!!!!!!
    sources = 0
    boundaries = 0

    ! indices of avh not clear! do they need to be changed to 3:ubnd_fcd-2 ??
    ! indices of avh not clear! do they need to be changed to 3:ubnd_fcd-2 ??
    state%avh(1:ubnd_fce-1) = 0.5_RK/sig_e*(state%num(1:ubnd_fce-1)+state%num(2:ubnd_fce)) ! average num for TKE

    if (self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
        flux(1) = state%avh(1)*(state%cde*((state%ko(1))**1.5_RK))/(kappa*(K_s+0.5_RK*grid%h(1)))**2
        flux(ubnd_fce) = state%avh(ubnd_fce-1)*(state%cde*((state%ko(ubnd_fce-1))**1.5_RK))/(kappa*(z0+0.5_RK*grid%h(ubnd_vol)))**2
        do i=2,ubnd_fce -1
            flux(i) = state%num(i)/sig_e*(state%cde*((state%ko(i))**1.5_RK)/(kappa*(z0+0.25_RK*(grid%h(i-1)+grid%h(i))))**2)
        end do
        state%avh(1) = 0.0_RK
        state%avh(ubnd_fce-1) = 0.0_RK
    else
        state%avh(1)=2*state%u_taub**4/sig_e/(state%eps(1)+state%eps(2))        ! = 0 for no shear stress
        state%avh(ubnd_fce-1)=2*state%u_taus**4/sig_e/(state%eps(ubnd_fce)+state%eps(ubnd_fce-1))   ! = 0 for no shear stress
    end if

    do i=2,ubnd_fce-1
        if (state%B(i)>0) then
            cee3=1.
        else
            cee3=ce3
        end if
        Prod=ce1*state%eps(i)/state%ko(i)*(state%P(i)+state%P_Seiche(i))        ! New code plus seiche
          !!DEBUG!!
          Prod=ce1*state%eps(i)/state%ko(i)*(state%P(i))        ! New code plus seiche
          !! END DEBUG!!
        Buoy=cee3*state%eps(i)/state%ko(i)*state%B(i)
        Diss=ce2*state%eps(i)*state%eps(i)/state%ko(i)

        if (Prod+Buoy>0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
        else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
        end if
    end do

    !!!!!!!! Define sources !!!!!!!!

    sources(2:ubnd_fce-1) = pplus(2:ubnd_fce-1)

    !!!!!! Define boundary conditions !!!!
    boundaries(2:ubnd_fce-1) = pminus(2:ubnd_fce-1)/state%eps(2:ubnd_fce-1)
    if(self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
      sources(2:ubnd_fce-1) = sources(2:ubnd_fce-1)+ flux(2:ubnd_fce-1)*grid%AreaFactor_eps(2:ubnd_fce-1) ! AreaFactor_eps = 1/A * dA/dz (at epsilon posions)
      if (grid%Az(1)/=0) then                          ! Flux from bottom only!
          sources(2)= sources(2)+flux(1)*(grid%Az(1)+grid%Az(2))/(grid%Az(2)*(grid%h(1)+grid%h(2)))
      end if
        sources(ubnd_fce-1)= sources(ubnd_fce-1)+flux(ubnd_fce)*(grid%Az(ubnd_fce)+grid%Az(ubnd_fce-1))/(grid%Az(ubnd_fce-1)*(grid%h(ubnd_fce-1)+grid%h(ubnd_fce-2)))
      ! eps(0) and eps(ubnd_fce) are assigned in the post processing function
    else ! no fluxes, unity A-matrix + condition on RHS
      rhs_0    = state%u_taub**2/sqrt(state%cm0*state%cde)
      rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

      ! Trick to have rhs(1) = rhs_0 and rhs(ubnd_fce) = rhs_ubnd
      ! in discretization, rhs is calculated as rhs(1) = var(1) + sources(1)*dt
      ! Given the equation below, the following results:
      ! rhs(1) = var(1) + (-var(1)/dt + rhs_0/dt)*dt = var(1) -var(1) + rhs_0 = rhs_0
      sources(1) = -(self%var(1)/state%dt)+(rhs_0/state%dt)
      sources(ubnd_fce) =  -(self%var(ubnd_fce)/state%dt)+(rhs_ubnd/state%dt)
    end if

  end associate
  end subroutine

  subroutine eps_var_post_solve(self, state, param)
    class(EpsModelVar), intent(inout) :: self
    class(ModelState), intent(inout) :: state
    class(ModelParam), intent(inout) :: param
    integer :: i
    real(RK) ::epslim
    associate(grid => self%grid, &
              ubnd_fce => self%grid%ubnd_fce, &
              ubnd_vol => self%grid%ubnd_vol)


    if(self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
    ! Define eps at boundary (no flux)
      self%var(1)= self%var(2) +(state%cde*((state%ko(2   ))**1.5_RK)/(kappa*(K_s+grid%h(1 )))**2)*grid%h(1 )
      self%var(ubnd_fce)=   self%var(ubnd_fce-1)+(state%cde*((state%ko(ubnd_fce-1))**1.5_RK)/(kappa*(z0 +grid%h(ubnd_vol)))**2)*grid%h(ubnd_vol)
    end if

    ! check lower limit of eps / update num
    do i=1,ubnd_fce
        if (state%NN(i)>0) then
            epslim = 0.212_RK*state%k(i)*sqrt(state%NN(i))
        else
            epslim = eps_min
        end if
        if(state%eps(i)<epslim) then
          state%eps(i)=epslim
        end if
        if (state%eps(i)<0) then
            write(6,*) 'Dissipation negative'
        end if


        state%num(i)= state%cmue1(i)*state%k(i)*state%k(i)/state%eps(i)+1.5e-6_RK
        state%nuh(i)= state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i)+1.5e-7_RK

    end do

    state%num(1)= kappa*state%u_taub*K_s+avh_min
    state%num(ubnd_fce)= kappa*state%u_taus*z0 +avh_min
    end associate
  end subroutine





end module strat_keps
