module simstrat_model_module
  use simstrat_kinds
  use simstrat_model_constants
  use simstrat_solver
  use simstrat_discretization_module, only: SimstratDiscretizationScheme
  use simstrat_finite_volume_implementation
  use utilities
  implicit none

  private

  !SimstratSimulationChain describes a sequence of procedures that has to be performed in every time step
  !###########################################
  type, public :: SimstratSimulationChain
  !###########################################
    procedure(generic_step_procedure), nopass, pointer :: step_procedure => null() !This is the actual procedure
    class(SimstratSimulationChain), pointer :: next => null() !next element in the chain
    contains
      procedure, pass(self) :: perform_simulation_step => link_step_procedure
  end type

  !A SimstratStateVariable is a special case of a SimulationChain step
  type, extends(SimstratSimulationChain), public :: SimstratStateVariable
    character(len=:), allocatable :: state_name
    real(RK), dimension(:), pointer :: state_value
    procedure(generic_source_term), nopass, pointer :: calculate_source_terms => null()
    contains
      procedure, pass(self) :: perform_simulation_step => chain_solve_state_variable
  end type

  !SimstratStateVariable_k is a special case of a SimstratStateVariable (because we have to use a specific interface to the discretization scheme)
  type, extends(SimstratStateVariable), public :: SimstratStateVariable_k
    contains
      procedure, pass(self) :: perform_simulation_step => chain_solve_state_variable_k
  end type

  !SimstratStateVariable_k is a special case of a SimstratStateVariable (because we have to use a specific interface to the discretization scheme)
  type, extends(SimstratStateVariable), public :: SimstratStateVariable_eps
    contains
      procedure, pass(self) :: perform_simulation_step => chain_solve_state_variable_eps
  end type

  type, public :: StateVariablePointer
    class(SimstratStateVariable), pointer :: ptr
  end type

  type, public :: SimstratOutputVariablePointer
    character(len=:), allocatable :: name
    real(RK), dimension(:), pointer :: value
  end type

  !Interface to forcing inputfiles
  !###########################################
  type, abstract, public :: SimstratForcing
  !###########################################
    !Forcing inputfile
    integer :: n_forcing
    integer :: forcing_type
    logical :: use_wind_filt
    real(RK), dimension(:), allocatable :: datum_forcing, u, v, Tair, Fsol, vap, cloud_coverage, Wf, SST, Hnet
    real(RK) :: u_datum, v_datum, Tair_datum, Fsol_datum, vap_datum, cloud_coverage_datum, Wf_datum, SST_datum, Hnet_datum
    !type(SimstratStateVariable), pointer :: T
    !Absorption inputfile
    integer :: n_absorption, n_depth_absorption
    real(RK), dimension(:), allocatable :: datum_absorption, depth_absorption, time_absorption
    real(RK), dimension(:,:), allocatable :: absorption
    real(RK), dimension(:), allocatable :: absorption_datum

    contains
      procedure(generic_forcing_init), deferred, pass(self) :: initialize
      procedure(generic_forcing_get), deferred, nopass :: interpolate_forcing_at_date
      procedure(generic_forcing_process), deferred, nopass :: process_forcing_at_date
  end type

  !###########################################
  type, abstract, public :: SimstratOutputLogger
  !###########################################
    contains
      procedure(generic_log_init), deferred, pass(self), public :: initialize
      procedure(generic_log), deferred, pass(self), public :: log
      procedure(generic_log_close), deferred, pass(self), public :: close
  end type

  !###########################################
  type, public :: SimstratModel
  !###########################################
    integer :: nz_max

    !Names of inputfiles
    character(len=:), allocatable          :: ParName,MorphName,InitName,ForcingName,AbsorpName
    character(len=:), allocatable          :: GridName,zoutName,PathOut
    character(len=:), allocatable          :: QinpName,QoutName,TinpName,SinpName
    character(len=:), allocatable          :: stencil_type

    !Model Configuration
    integer :: Mod, Stab, ModFlux, ForcingType, ModSNorm, ModC10, ModInflow, Pgrad, ModSal, disp_sim, disp_dgn, igoal
    logical :: UseWindFilt
    logical :: couple_aed2

    !Model parameter
    real(RK) :: Lat, p_air, a_seiche, q_NN, f_wind, C10, CD, fgeo, k_min, p_radin, p_windf, beta_sol, albsw
    real(RK) :: cori
    real(RK) :: cde, cm0
    integer :: salctr, delsal

    !Morphology transformation
    real(RK) :: z_zero
    real(RK) :: A_surf, volume, depth

    !Simulation parameter
    real(RK) :: dt, t_start, t_end

    class(SimstratDiscretizationScheme), pointer :: discretization
    !class(AbstractSimstratModelIntegrator), pointer :: integrator

    class(SimstratForcing), pointer :: forcing
    logical :: use_buffered_forcing

    !State variables
    type(StateVariablePointer), dimension(:), allocatable :: state_vars
    real(RK), dimension(:), pointer :: T, S, u, v, k, eps
    !type(SimstratStateVariable), dimension(:), allocatable :: state_vars

    !additional variables needed during simulation
    real(RK) :: fsed
    real(RK), dimension(:), allocatable :: fgeo_add, fsed_add
    !related to state variables
    real(RK), dimension(:), allocatable :: L, ko !TKE dissipation rate [W/kg]
    real(RK), dimension(:), pointer :: num, nuh, num_eps !Turbulent viscosity (momentum) and diffusivity (temperature)
    real(RK), dimension(:), pointer :: P, B, NN !Shear stress production [W/kg], buoyancy production [W/kg], Brunt-Väisälä frequency [s-2]

    real(RK), dimension(:), allocatable :: dS !Source/sink for salinity
    real(RK), dimension(:), allocatable :: cmue1,cmue2 !Model constants
    real(RK), dimension(:), allocatable :: rho !density
    real(RK), dimension(:), allocatable :: P_Seiche !Production of TKE [W/kg]
    real(RK) :: E_Seiche  !seiche energy [J]
    real(RK) :: gamma !Constant of proportionality for loss of seiche energy

    !Absorption
    real(RK), dimension(:), allocatable :: ga1 !Absorption coeff [m-1]

    !Forcing
    real(RK) :: u10,v10,uv10,Wf !Absorption coeff [m-1], wind drag, wind speeds
    real(RK) :: u_taub,drag,u_taus
    real(RK) :: tx,ty,SST,heat !Shear stress, sea surface temperature and heat flux
    real(RK) :: rad0 !Solar radiation at surface
    !real(RK), dimension(:), allocatable :: rad !Solar radiation in water

    !In- and outflow
    real(RK), dimension(:), allocatable :: Qvert, Q_inp !Vertical and horizontal flows

    !Simulation chain
    class(SimstratSimulationChain), pointer :: simulation_chain

    class(SimstratOutputLogger), pointer :: output_logger
    type(SimstratOutputVariablePointer), dimension(:), allocatable :: output_vars
    logical :: write_on_the_fly
    integer :: thinning_interval

  contains
    procedure, pass(model) :: initialize => initializeSimstratModel
    procedure, pass(model) :: run => runSimstratModel
    procedure, pass(model) :: destroy => destroySimstratModel
  end type


  abstract interface
    subroutine generic_step_procedure(datum, model)
      import SimstratModel, RK
      implicit none

      class(SimstratModel), intent(inout) :: model
      real(RK), intent(in) :: datum
    end subroutine

    subroutine generic_source_term(datum, model, boudary_cond, fluxes, sourcessource_terms)
      import SimstratModel, RK
      implicit none

      class(SimstratModel), intent(in) :: model
      real(RK), dimension(:), intent(inout) :: boudary_cond, fluxes, sourcessource_terms
      real(RK), intent(in) :: datum
    end subroutine

    subroutine generic_forcing_init(self, ForcingName, forcing_type, use_wind_filt, disp_dgn, AbsorpName, z_zero)
      import SimstratForcing, SimstratStateVariable, RK
      implicit none

      class(SimstratForcing), intent(inout) :: self
      character(len=:), intent(in), allocatable :: ForcingName, AbsorpName
      integer, intent(in) :: forcing_type
      logical, intent(in) :: use_wind_filt
      integer, intent(in) :: disp_dgn
      real(RK), intent(in) :: z_zero
    end subroutine

    subroutine generic_forcing_get(datum, model)
      import RK, SimstratModel
      implicit none

      real(RK), intent(in) :: datum
      class(SimstratModel), intent(inout) :: model
    end subroutine

    subroutine generic_forcing_process(datum, model)
      import RK, SimstratModel
      implicit none

      real(RK), intent(in) :: datum
      class(SimstratModel), intent(inout) :: model
    end subroutine

    subroutine generic_log_init(self, model)
      import SimstratOutputLogger, SimstratModel
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(SimstratModel), intent(in) :: model
    end subroutine

    subroutine generic_log(self, datum, model)
      import SimstratOutputLogger, RK, SimstratModel
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
      class(SimstratModel), intent(in) :: model
      real(RK), intent(in) :: datum
    end subroutine

    subroutine generic_log_close(self)
      import SimstratOutputLogger, SimstratModel
      implicit none

      class(SimstratOutputLogger), intent(inout) :: self
    end subroutine
  end interface

contains
  subroutine initializeSimstratModel(model, nz_max)
    implicit none
    class(SimstratModel), intent(inout) :: model
    integer, optional, intent(in) :: nz_max

    if(present(nz_max)) model%nz_max = nz_max

    associate(nz => model%nz_max)
      allocate(model%L(nz), model%ko(nz),&
               model%num(nz), model%nuh(nz), model%num_eps(nz),&
               model%P(nz), model%B(nz), model%NN(nz),&
               model%ga1(nz),&
               model%dS(nz),&
               model%cmue1(nz), model%cmue2(nz),&
               model%rho(nz),&
               model%P_Seiche(nz),&
   !            model%rad(nz),&
               model%Qvert(nz), model%Q_inp(nz),&
               model%fgeo_add(nz), model%fsed_add(nz))
   end associate

   model%L = 0.0_RK
   model%ko = 0.0_RK
   model%num = 0.0_RK
   model%nuh = 0.0_RK
   model%num_eps = 0.0_RK
   model%NN = 0.0_RK
   model%P = 0.0_RK
   model%B = 0.0_RK
   model%ga1 = 0.0_RK
   model%dS = 0.0_RK
   model%cmue1 = 0.0_RK
   model%cmue2 = 0.0_RK
   model%rho = 0.0_RK
   model%P_Seiche = 0.0_RK
   !model%rad = 0.0_RK
   model%Qvert = 0.0_RK
   model%Q_inp = 0.0_RK
   model%fgeo_add = 0.0_RK
   model%fsed_add = 0.0_RK
  end subroutine

  subroutine runSimstratModel(model)
    implicit none
    class(SimstratModel), target, intent(inout) :: model

    real(RK) :: datum
    real(RK) :: cpustart, cpufinish
    class(SimstratSimulationChain), pointer :: chain_step

    class(SimstratModel), pointer :: modelptr

    modelptr => model

    write(*,*) 'Run simulation...'
    call cpu_time(cpustart)

    !current time of the simulation
    datum = model%t_start
    !precondition the discretization scheme
    call model%discretization%precondition(model%dt)

    !initialize output logger
    call model%output_logger%initialize(model)

    !loop through time steps
    do while (datum<model%t_end)
      !loop through simulation chain
      chain_step => model%simulation_chain
      do while(associated(chain_step))
        call chain_step%perform_simulation_step(datum, model)
        chain_step => chain_step%next
      end do

      datum = datum + model%dt/86400
      call model%output_logger%log(datum, model)
    end do
    call cpu_time(cpufinish)
    write(*,fmt="(1X,A,F8.2,A)") 'Simulation completed after ',cpufinish-cpustart,' seconds.'

    !close output logger
    call model%output_logger%close()
  end subroutine

  subroutine destroySimstratModel(model)
    implicit none
    class(SimstratModel), intent(inout) :: model
  end subroutine

  subroutine link_step_procedure(self, datum, model)
    implicit none
    class(SimstratSimulationChain), intent(inout) :: self
    class(SimstratModel), intent(inout) :: model
    real(RK), intent(in) :: datum

    call self%step_procedure(datum, model)

  end subroutine

  subroutine chain_solve_state_variable(self, datum, model)
    implicit none
    class(SimstratStateVariable), intent(inout) :: self
    class(SimstratModel), pointer, intent(inout) :: model
    real(RK), intent(in) :: datum

    real(RK), dimension(model%discretization%nz_mfq) :: md, rhs
    real(RK), dimension(model%discretization%nz_mfq) :: ld, ud
    real(RK), dimension(model%discretization%nz_mfq) :: boundary_cond, sources
    real(RK), dimension(model%discretization%nz_tq) :: fluxes
    md = 0.0_RK
    ld = 0.0_RK
    ud = 0.0_RK
    rhs = 0.0_RK
    boundary_cond = 0.0_RK
    fluxes = 0.0_RK
    sources = 0.0_RK

    !get source terms
    call self%calculate_source_terms(datum, model, boundary_cond, fluxes, sources)

    !create diagonals
    call model%discretization%createLinearEquationSystem_MFQ(self%state_value, model%nuh, boundary_cond, fluxes, sources, ld, md, ud, rhs)

    !solve linear equation system
    call solve_tridiag_thomas(ld, md, ud, rhs, self%state_value, model%discretization%nz_mfq)

  end subroutine

  subroutine chain_solve_state_variable_k(self, datum, model)
    implicit none
    class(SimstratStateVariable_k), intent(inout) :: self
    class(SimstratModel), intent(inout) :: model
    real(RK), intent(in) :: datum

    real(RK), dimension(model%discretization%nz_tq) :: md, rhs
    real(RK), dimension(model%discretization%nz_tq) :: ld, ud
    real(RK), dimension(model%discretization%nz_tq) :: boundary_cond, sources
    real(RK), dimension(model%discretization%nz_tq+1) :: fluxes
    real(RK), dimension(model%discretization%nz_tq+1) :: num_k
    md = 0.0_RK
    ld = 0.0_RK
    ud = 0.0_RK
    rhs = 0.0_RK
    num_k = 0.0_RK
    boundary_cond = 0.0_RK
    fluxes = 0.0_RK
    sources = 0.0_RK

    !get source terms
    !call self%calculate_source_terms(datum, model, boundary_cond, fluxes, sources)

    associate(nz_tq=>model%discretization%nz_tq, eps=>model%eps, u_taub=>model%u_taub, u_taus=>model%u_taus, Mod=>model%Mod, ModFlux=>model%ModFlux, num=>model%num, cm0=>model%cm0, cde=>model%cde, ko=>model%ko)
      !This is very specific to k and would not work for other variables that are initialized as SimstratStateVariable_k
      ko(1:nz_tq) = self%state_value(1:nz_tq) ! ko = TKE at old time step

      !Interpolate num for TKE
      num_k(2:nz_tq) = 0.5_RK/sig_k*(num(1:nz_tq-1)+num(2:nz_tq))
      if ((ModFlux==1).and.(Mod==1)) then
          num_k(1) = 0.0_RK
          num_k(nz_tq+1) = 0.0_RK
      else
          num_k(1)=2*u_taub**4/(eps(1)+eps(2))        ! = 0 for no shear stress
          num_k(nz_tq+1)=2*u_taus**4/(eps(nz_tq)+eps(nz_tq-1))   ! = 0 for no shear stress
      end if

      !create diagonals
      call model%discretization%createLinearEquationSystem_k(self%state_value, num_k, boundary_cond, fluxes, sources, ld, md, ud, rhs)

      !In the future: find a better solution for this. This should somehow go into the k_terms
      if ((ModFlux==1).and.(Mod==1)) then
        !solve linear equation system
        call solve_tridiag_thomas(ld, md, ud, rhs, self%state_value, nz_tq)
        !call solve_tridiag_thomas(ld(2:nz_tq-1), md(2:nz_tq-1), ud(2:nz_tq-1), rhs(2:nz_tq-1), self%state_value(2:nz_tq-1), nz_tq-2)
        self%state_value(1)= self%state_value(2)                                ! Define TKE at boundary (no flux)
        self%state_value(model%discretization%nz_tq)= self%state_value(model%discretization%nz_tq-1)
      else
        ud(1)= 0.0_RK
        md(1)= 1.0_RK
        rhs(1)= u_taub**2/sqrt(cm0*cde)

        md(nz_tq)= 1.0_RK
        ld(nz_tq)= 0.0_RK
        rhs(nz_tq)= u_taus**2/sqrt(cm0*cde)

        call solve_tridiag_thomas(ld, md, ud, rhs, self%state_value, nz_tq)
      end if
      where(self%state_value < model%k_min) self%state_value = model%k_min
    end associate
  end subroutine

  subroutine chain_solve_state_variable_eps(self, datum, model)
    implicit none
    class(SimstratStateVariable_eps), intent(inout) :: self
    class(SimstratModel), intent(inout) :: model
    real(RK), intent(in) :: datum

    real(RK), dimension(model%discretization%nz_tq) :: md, rhs
    real(RK), dimension(model%discretization%nz_tq) :: ld, ud, NN_faces
    real(RK), dimension(model%discretization%nz_tq) :: boundary_cond, sources
    real(RK), dimension(model%discretization%nz_tq+1) :: fluxes
    real(RK), dimension(model%discretization%nz_tq-2) :: eps_sub
    real(RK) :: epslim
    integer :: i
    md = 0.0_RK
    ld = 0.0_RK
    ud = 0.0_RK
    rhs = 0.0_RK
    boundary_cond = 0.0_RK
    fluxes = 0.0_RK
    sources = 0.0_RK

    select type(discretization=>model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)

        associate(nz_tq=>discretization%nz_tq, eps=>self%state_value, h_centres=>discretization%h_centres, u_taub=>model%u_taub, u_taus=>model%u_taus, Mod=>model%Mod, ModFlux=>model%ModFlux, num=>model%num, cm0=>model%cm0, cde=>model%cde, ko=>model%ko, num_eps=>model%num_eps, NN=>model%NN, k=>model%k, cmue1=>model%cmue1, cmue2=>model%cmue2, nuh=>model%nuh, L=>model%L)
          num_eps(2:nz_tq) = 0.5_RK/sig_e*(num(1:nz_tq-1)+num(2:nz_tq))
          if ((ModFlux==1).and.(Mod==1)) then
              num_eps(2) = 0.0_RK
              num_eps(nz_tq+1) = 0.0_RK
          else
              num_eps(2) = 2*u_taub**4/sig_e/(eps(1)+eps(2))   ! = 0 for no shear stress
              num_eps(nz_tq+1) = 2*u_taus**4/sig_e/(eps(nz_tq)+eps(nz_tq-1))   ! = 0 for no shear stress
          end if

          call self%calculate_source_terms(datum, model, boundary_cond, fluxes, sources)

          !create diagonals
          call model%discretization%createLinearEquationSystem_eps(eps, num_eps, boundary_cond, fluxes, sources, ld, md, ud, rhs)

          if ((ModFlux==1).and.(Mod==1)) then
            call solve_tridiag_thomas(ld, md, ud, rhs, eps, nz_tq)
            call solve_tridiag_thomas(ld(2:nz_tq-1), md(2:nz_tq-1), ud(2:nz_tq-1), rhs(2:nz_tq-1), eps(2:nz_tq-1), nz_tq-2)
            !In the future: find a better solution for this. This should somehow go into the eps_terms
            !its too specific for simstrat_model_modle
            eps(1 )= eps(2   )+(cde*((ko(2   ))**1.5)/(kappa*(K_s+h_centres(1 )))**2)*h_centres(1 )
            eps(nz_tq)= eps(nz_tq-1)+(cde*((ko(nz_tq-1))**1.5)/(kappa*(z0 +h_centres(nz_tq-1)))**2)*h_centres(nz_tq-1)
          else
            ud(1)= 0.0_RK
            md(1)= 1.0_RK
            rhs(1)= cde*k(1)**1.5_RK/kappa/K_s

            md(nz_tq)= 1.0_RK
            ld(nz_tq)= 0.0_RK
            rhs(nz_tq)= cde*k(nz_tq)**1.5_RK/kappa/z0

            call solve_tridiag_thomas(ld, md, ud, rhs, eps, nz_tq)
          end if

          NN_faces(2:nz_tq) = NN(1:nz_tq-1)
          NN_faces(1) = NN_faces(2)

          do i=1,nz_tq
              if (NN_faces(i)>0) then
                  epslim= 0.212*k(i)*sqrt(NN_faces(i))
              else
                  epslim= eps_min
              end if
              if(eps(i)<epslim) then
                eps(i)=epslim
              end if
              if (eps(i)<0) then
                  write(*,*) 'Dissipation negative'
              end if

              num(i)= cmue1(i)*k(i)**2/eps(i)+1.5e-6
              nuh(i)= cmue2(i)*k(i)**2/eps(i)+1.5e-7
              L(i)= cde*sqrt(ko(i)*ko(i)*ko(i))/eps(i)
          end do

          num(2 )= kappa*u_taub*K_s+avh_min
          num(nz_tq)= kappa*u_taus*z0 +avh_min
        end associate
      class default
        call error('chain not implemeted for this discretization')
    end select

  end subroutine
end module simstrat_model_module
