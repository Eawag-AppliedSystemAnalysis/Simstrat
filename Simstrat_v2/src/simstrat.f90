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
   use strat_ice     
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
   type(IceModule) :: mod_ice   
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

   if (simdata%model%has_advection) then
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
   end if

   ! Setup logger
   call logger%initialize(simdata%output_cfg, simdata%grid)

   ! initialize simulation modules
   call mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)
   call mod_turbulence%init(simdata%grid, simdata%model_cfg, simdata%model_param)
   call mod_ice%init(simdata%model_cfg, simdata%model_param, simdata%grid)

   ! Set temperature state var to have nu_h as nu and T as model variable
   call mod_temperature%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%T, simdata%grid%ubnd_vol)

   ! Set U and V var to have num as nu and U reps V as model variable
   ! also, assign shear stress in model for this variable
   call mod_u%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%U, simdata%grid%ubnd_vol)
   call mod_u%assign_shear_stress(simdata%model%tx)

   call mod_v%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%V, simdata%grid%ubnd_vol)
   call mod_v%assign_shear_stress(simdata%model%ty)

   ! Set mod_s (transport module) to have nuh as nu and to manipulate S based on dS
   call mod_s%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%S, simdata%grid%ubnd_vol)
   call mod_s%assign_external_source(simdata%model%dS)

   ! Set up K and eps state vars with keps discretization and avh as nu
   call mod_k%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%K, simdata%grid%ubnd_fce)
   call mod_eps%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%eps, simdata%grid%ubnd_fce)

   call run_simulation()

   ! close logger files after simulation
   call logger%close()

contains

   subroutine run_simulation()
      integer :: i
      !call logger%log(0.0_RK) ! Write initial conditions

      ! Run simulation until end datum
      do while (simdata%model%datum<simdata%sim_cfg%end_datum)

         !increase datum and step
         simdata%model%datum = simdata%model%datum + simdata%model%dt/86400
         simdata%model%std = simdata%model%std + 1
         
         ! Read forcing file
         call mod_forcing%update(simdata%model)

         ! Update absorption
         call mod_absorption%update(simdata%model)

         ! Update physics
         call mod_stability%update(simdata%model)
         ! If there is inflow/outflow do advection part
         if (simdata%model%has_advection) then
            call mod_lateral%update(simdata%model)
            call mod_advection%update(simdata%model)
         end if
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
		 
         ! Update ice 
         if (simdata%model_cfg%ice_model == 1) then 
            call mod_ice%update(simdata%model, simdata%model_param)   
         end if      
		 
         ! Call logger to write files
         if (mod(simdata%model%std,simdata%output_cfg%thinning_interval)==0) then
            call logger%log(simdata%model%datum)
         end if

         ! Display simulation (datum, lake surface, temperature at bottom, temperature at surface)

         ! Standard display: display when logged
         if (simdata%model_cfg%disp_simulation==1) then
            if (mod(simdata%model%std,simdata%output_cfg%thinning_interval)==0) then
               write(6,'(F10.4,F10.5,F10.5,F10.5)') simdata%model%datum, simdata%grid%z_face(simdata%grid%ubnd_fce), &
               simdata%model%T(simdata%grid%nz_occupied), simdata%model%T(1)
            end if
         ! Extra display: display every iteration
         else if (simdata%model_cfg%disp_simulation==2) then
            write(6,'(F10.4,F10.4,F10.4,F10.4)') simdata%model%datum, simdata%grid%z_face(simdata%grid%ubnd_fce), &
            simdata%model%T(simdata%grid%ubnd_vol), simdata%model%T(1)
         end if
         
      end do
   end subroutine

end program simstrat_main
