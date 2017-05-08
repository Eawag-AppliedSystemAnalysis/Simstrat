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

contains
  subroutine k_var_calc_terms(self, state, param, sources, boundaries)
    class(KModelVar), intent(inout) :: self
    class(ModelState), intent(inout) :: state
    class(ModelParam), intent(inout) :: param
    real(RK),  dimension(:) ::  sources, boundaries
    real(RK) :: rhs_0, rhs_ubnd
    real(RK) :: pminus(self%grid%ubnd_fce), pplus(self%grid%ubnd_fce), Prod, Buoy, Diss

    integer :: i
    associate(grid => self%grid, &
              ubnd_fce => self%grid%ubnd_fce, &
              ubnd_vol => self%grid%ubnd_vol)


    !!!!!!!! Precalculations !!!!!!!!
    sources = 0
    boundaries = 0
    state%ko(1:ubnd_fce) = state%k(1:ubnd_fce) ! ko = TKE at old time step

    ! indices of avh not clear! do they need to be changed to 3:ubnd_fcd-2 ??
    state%avh(3:ubnd_fce-2) = 0.5_RK/sig_k*(state%num(1:ubnd_fce-2)+state%num(2:ubnd_fce-1)) ! average num for TKE

    if (self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
        state%avh(2) = 0.0_RK
        state%avh(ubnd_fce-1) = 0.0_RK
    else
        state%avh(2)=2*state%u_taub**4/(state%eps(0)+state%eps(1))        ! = 0 for no shear stress
        state%avh(ubnd_fce-1)=2*state%u_taus**4/(state%eps(ubnd_fce)+state%eps(ubnd_fce-1))   ! = 0 for no shear stress
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
    boundaries(2:ubnd_fce-1) = pminus(2:ubnd_fce-1)/state%eps(2:ubnd_fce-1)

    if(self%cfg%flux_condition==1 .and. self%cfg%turbulence_model == 1) then
      ! K(0) and K(ubnd_fce) are assigned in the post processing function
    else ! no fluxes, unity A-matrix + condition on RHS
      rhs_0    = state%u_taub**2/sqrt(state%cm0*state%cde)
      rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

      !Setting avh and boundaries to zero, results in a unit A - matrix
      ! as lower_diag/upper_diag are calculated using avh as multiplier
      ! and the main_diag is simply 1-lower_diag-upper_diag +boundaries*dt
      boundaries = 0
      state%avh = 0

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
    do i=1,ubnd_fce
        if(self%var(i)<k_min) self%var(i)=k_min             ! Lower limit of TKE
    end do

    end associate
  end subroutine


end module strat_keps
