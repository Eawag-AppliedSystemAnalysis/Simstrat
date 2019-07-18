!<  +---------------------------------------------------------------+
!     Simstrat model for simulation of
!     vertical transport in lakes and reservoirs
!<  +---------------------------------------------------------------+

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
   type(InterpolatingLogger) :: logger
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
   write (6, *) 'Simstrat version '//version
   write (6, *) 'This software has been developed at eawag - Swiss Federal Institute of Aquatic Science and Technology'
   write (6, *) ''

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

      ! initialize lateral module based on configuration
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
   call logger%initialize(simdata%sim_cfg, simdata%output_cfg, simdata%grid)

   ! Initialize simulation modules
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

   ! Close logger files after simulation
   call logger%close()

contains

   subroutine run_simulation()
      !! run the marching time loop

      call logger%log(simdata%model%datum) ! Write initial conditions

      ! Run simulation until end datum or until no more results are required by the output time file
      do while (simdata%model%datum<simdata%sim_cfg%end_datum .and. simdata%model%output_counter<=size(simdata%output_cfg%tout))

         ! ****************************************
         ! ***** Update counters and timestep *****
         ! ****************************************

         ! Increase iteration counter
         simdata%model%model_step_counter = simdata%model%model_step_counter + 1

         ! If output times are specified in file
         if(simdata%output_cfg%thinning_interval==0) then
            ! At the first iteration or always after the model output was logged
            if (simdata%model%model_step_counter==1) then
               ! Adjust timestep so that the next output time is on the "time grid"
               simdata%model%dt = simdata%output_cfg%adjusted_timestep(simdata%model%output_counter)
               ! Don't log anymore until the next output time is reached
               simdata%output_cfg%write_to_file = .false.
               ! Log initial state if first output time corresponds to simulation start
               if (simdata%model%dt < 1e-6) then
                  ! Log initial conditions and display
                  call logger%log(simdata%model%datum)
                  if (simdata%sim_cfg%disp_simulation > 0) then
                     write(6,'(F12.4,F20.4,F15.4,F15.4)') simdata%model%datum, simdata%grid%lake_level, &
                     simdata%model%T(simdata%grid%ubnd_vol), simdata%model%T(1)
                  end if
                  ! Update counters and timestep
                  simdata%model%model_step_counter = simdata%model%model_step_counter + 1
                  simdata%model%output_counter = simdata%model%output_counter + 1
                  simdata%model%dt = simdata%output_cfg%adjusted_timestep(simdata%model%output_counter)
               end if

            else if (simdata%output_cfg%write_to_file) then
               ! Adjust timestep so that the next output time is on the "time grid"
               simdata%model%dt = simdata%output_cfg%adjusted_timestep(simdata%model%output_counter)
               ! Don't log anymore until the next output time is reached
               simdata%output_cfg%write_to_file = .false.
            end if

            ! If the next output time is reached
            if(simdata%model%model_step_counter==sum(int(simdata%output_cfg%n_timesteps_between_tout(1:simdata%model%output_counter)))) then
               ! Increase output counter
               simdata%model%output_counter = simdata%model%output_counter + 1
               ! Log next model state
               simdata%output_cfg%write_to_file = .true.
            end if
         else ! If output frequency is given in file
            ! If the next output time is reached
            if(mod(simdata%model%model_step_counter,simdata%output_cfg%thinning_interval)==0) then
               ! Log next model state
               simdata%output_cfg%write_to_file = .true.
            else
               ! Don't log
               simdata%output_cfg%write_to_file = .false.
            end if
         end if

         ! Advance to the next timestep
         simdata%model%datum = simdata%model%datum + simdata%model%dt/86400

         ! ************************************
         ! ***** Compute next model state *****
         ! ************************************

         ! Read forcing file
         call mod_forcing%update(simdata%model)

         ! Update absorption
         call mod_absorption%update(simdata%model)

         ! Update physics
         call mod_stability%update(simdata%model)

         ! If there is inflow/outflow do advection part
         if (simdata%model%has_advection) then
            ! Treat inflow/outflow
            call mod_lateral%update(simdata%model)
            ! Set old lake level (before it is changed by advection module)
            simdata%grid%lake_level_old = simdata%grid%z_face(simdata%grid%ubnd_fce)
            ! Update lake advection using the inflow/outflow data
            call mod_advection%update(simdata%model)
            ! Update lake level
            simdata%grid%lake_level = simdata%grid%z_face(simdata%grid%ubnd_fce)
         end if

         ! Display to screen
         if (simdata%model%model_step_counter==1 .and. simdata%sim_cfg%disp_simulation/=0) then
            write(6,*)
            write(6,*) ' -------------------------- '
            write(6,*) '   SIMULATION IN PROGRESS   '
            write(6,*) ' -------------------------- '
            write(6,*)
            if(simdata%sim_cfg%disp_simulation/=0) write(6,'(A12, A20, A20, A20)') 'Time [d]','Surface level [m]','T_surf [degC]','T_bottom [degC]'
         end if

         ! Update Coriolis
         call mod_forcing%update_coriolis(simdata%model)

         ! Update and solve U and V - terms
         call mod_u%update(simdata%model, simdata%model_param)
         call mod_v%update(simdata%model, simdata%model_param)

         ! Update and solve T - terms
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
         if (simdata%output_cfg%write_to_file) then
            call logger%log(simdata%model%datum)
         end if

         ! ***********************************
         ! ***** Log to file and display *****
         ! ***********************************

         ! Standard display: display when logged: datum, lake surface, T(1), T(surf)
         if (simdata%sim_cfg%disp_simulation==1 .and. simdata%output_cfg%write_to_file) then
            write(6,'(F12.4,F16.4,F20.4,F20.4)') simdata%model%datum, simdata%grid%lake_level, &
            simdata%model%T(simdata%grid%nz_occupied), simdata%model%T(1)

         ! Extra display: display every iteration: datum, lake surface, T(1), T(surf)
         else if (simdata%sim_cfg%disp_simulation==2) then
            write(6,'(F12.4,F20.4,F15.4,F15.4)') simdata%model%datum, simdata%grid%lake_level, &
            simdata%model%T(simdata%grid%ubnd_vol), simdata%model%T(1)
         end if

         ! This logical is used to do some allocation in the forcing, absorption and lateral subroutines during the first timestep
         simdata%model%first_timestep = .false.
      end do

      if (simdata%sim_cfg%disp_simulation/=0) then
         write(6,*)
         write(6,*) ' -------------------------- '
         write(6,*) '    SIMULATION COMPLETED    '
         write(6,*) ' -------------------------- '
         write(6,*)
      end if
   end subroutine

end program simstrat_main
