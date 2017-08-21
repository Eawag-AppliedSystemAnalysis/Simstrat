!     +---------------------------------------------------------------+
!     |  Simstrat model for simulation of                             |
!     |  vertical transport in lakes and reservoirs                   |
!     +---------------------------------------------------------------+

program simstrat_main
   use strat_kinds
   use strat_inputfile, only: SimstratSimulationFactory
   use strat_outputfile
   use strat_simdata, only: SimulationData
   use strat_forcing
   use utilities
   use strat_stability, only: StabilityModule
   use strat_windshear
   use strat_statevar
   use strat_temp
   use strat_solver
   use strat_discretization
   use strat_keps
   use strat_turbulence
   use strat_transport
   use strat_absorption
   use strat_advection
   use strat_lateral
   use, intrinsic :: ieee_arithmetic

   implicit none

   ! Instantiate all modules
   ! note that some are pointers/targets for polymorphism reasons
   type(SimstratSimulationFactory) :: factory
   class(SimulationData), pointer :: simdata
   type(ThomasAlgSolver) :: solver
   type(EulerIDiscretizationMFQ) :: euler_i_disc
   type(EulerIDiscretizationKEPS) :: euler_i_disc_keps
   type(ForcingModule) :: mod_forcing
   type(StabilityModule) :: mod_stability
   type(SimpleLogger) :: logger
   type(TempModelVar) :: mod_temperature
   type(UVModelVar) :: mod_u, mod_v
   type(KModelVar) :: mod_k
   type(EpsModelVar) :: mod_eps
   type(TransportModVar) :: mod_s
   type(TurbulenceModule) :: mod_turbulence
   type(AbsorptionModule) :: mod_absorption
   type(AdvectionModule) :: mod_advection
   type(LateralModule), target :: mod_lateral_normal
   type(LateralRhoModule), target :: mod_lateral_rho
   class(GenericLateralModule), pointer :: mod_lateral

   character(len=100) :: arg
   character(len=:), allocatable :: ParName

   !print some information
   write (*, *) 'Simstrat version '//version
   write (*, *) 'This software has been developed at eawag - Swiss Federal Institute of Aquatic Science and Technology'
   write (*, *) ''

   !get first cli argument
   call get_command_argument(1, arg)
   ParName = trim(arg)
   if (ParName == '') ParName = 'simstrat.par'

   !initialize model from inputfiles
   call factory%initialize_model(ParName, simdata)

   ! initialize Discretization
   call euler_i_disc%init(simdata%grid)
   call euler_i_disc_keps%init(simdata%grid)

   !initialize forcing module
   call mod_forcing%init(simdata%model_cfg, &
                         simdata%model_param, &
                         simdata%input_cfg%ForcingName, &
                         simdata%grid)

   ! initialize absorption module
   call mod_absorption%init(simdata%model_cfg, &
                            simdata%model_param, &
                            simdata%input_cfg%AbsorpName, &
                            simdata%grid)

   ! initialize advection module
   call mod_advection%init(simdata%model_cfg, &
                           simdata%model_param, &
                           simdata%grid)

   ! initliaze lateral module based on configuration
   if (simdata%model_cfg%inflow_placement == 1) then
      ! Gravity based inflow
      mod_lateral => mod_lateral_rho
   else
      mod_lateral => mod_lateral_normal
   end if
   call mod_lateral%init(simdata%model_cfg, &
                         simdata%model_param, &
                         simdata%grid)

   ! Setup logger
   call logger%initialize(simdata%output_cfg, simdata%grid)

   ! initialize simulation modules
   call mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)
   call mod_turbulence%init(simdata%grid, simdata%model_cfg, simdata%model_param)

   ! Set temperature state var to have nu_h as nu and T as model variable
   call mod_temperature%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%T)

   ! Set U and V var to have num as nu and U reps V as model variable
   ! also, assign shear stress in model for this variable
   call mod_u%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%U)
   call mod_u%assign_shear_stress(simdata%model%tx)

   call mod_v%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%V)
   call mod_v%assign_shear_stress(simdata%model%ty)

   ! Set mod_s (transport module) to have nuh as nu and to manipulate S based on dS
   call mod_s%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%S)
   call mod_s%assign_external_source(simdata%model%dS)

   ! Set up K and eps state vars with keps discretization and avh as nu
   call mod_k%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%K)
   call mod_eps%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%eps)

   call run_simulation()

   ! close logger files after simulation
   call logger%close()

contains

   subroutine run_simulation()
      integer :: i
      call logger%log(0.0_RK) ! Write initial conditions

      ! Time step (currently fixed - could be changed in future)
      simdata%model%dt = 0.1

      do i = 1, 10000
         simdata%model%std = i
         if (simdata%model%datum >= simdata%sim_cfg%end_datum) then
            exit
         end if

         ! Read forcing file
         call mod_forcing%update(simdata%model)

         ! Update absorption
         call mod_absorption%update(simdata%model)

         ! Update physics
         call mod_stability%update(simdata%model)
         call mod_lateral%update(simdata%model)
         call mod_advection%update(simdata%model)
         call mod_forcing%update_coriolis(simdata%model)

         ! Update and solve U and V - terms
         call mod_u%update(simdata%model, simdata%model_param)
         call mod_v%update(simdata%model, simdata%model_param)

         ! Update and solve t - terms
         call mod_temperature%update(simdata%model, simdata%model_param)

         ! Update and solve transportation terms (here: Salinity S only)
         call mod_S%update(simdata%model, simdata%model_param)

         ! update turbulence states
         call mod_turbulence%update(simdata%model, simdata%model_param)

         ! Solve k & eps
         call mod_k%update(simdata%model, simdata%model_param)
         call mod_eps%update(simdata%model, simdata%model_param)

         ! Call logger to write files
         call logger%log(simdata%model%datum)

         !increase datum
         simdata%model%datum = simdata%model%datum + simdata%model%dt
      end do
   end subroutine

end program simstrat_main
