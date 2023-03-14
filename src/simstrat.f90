! ---------------------------------------------------------------------------------
!     Simstrat a physical 1D model for lakes and reservoirs
!
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     Copyright (C) 2020, Eawag
!
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>. 
! ---------------------------------------------------------------------------------
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
   use simstrat_aed2
   use strat_lateral
   use forbear
   use csv_module
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
   type(SimstratAED2) :: mod_aed2
   type(LateralModule), target :: mod_lateral_normal
   type(LateralRhoModule), target :: mod_lateral_rho
   class(GenericLateralModule), pointer :: mod_lateral
   ! Instantiate progress bar object
   type(bar_object):: bar
   type(csv_file), dimension(:), allocatable :: csv_files

   character(len=100) :: arg
   character(len=:), allocatable :: ParName
   character(len=:), allocatable :: snapshot_file_path
   character(len=:), allocatable :: file_text_restart
   character(len=:), allocatable :: file_text_restart2
   logical :: continue_from_snapshot = .false.
   integer(8),dimension(2) :: simulation_end_time
   real(RK) :: new_start_datum

   ! Print some information
   write (6, *) 'Simstrat version '//version
   write (6, *) 'Coupled with the biogeochemical library AED2'
   write (6, *) 'This software has been developed at Eawag - Swiss Federal Institute of Aquatic Science and Technology'
   write (6, *) ''

   ! Get first cli argument
   call get_command_argument(1, arg)
   ParName = trim(arg)
   if (ParName == '') ParName = 'simstrat.par'

   ! Initialize model from input files
   call factory%initialize_model(ParName, simdata)

   ! Initialize Discretization
   call euler_i_disc%init(simdata%grid)
   call euler_i_disc_keps%init(simdata%grid)

   ! Initialize forcing module
   call mod_forcing%init(simdata%model_cfg, &
                         simdata%model_param, &
                         simdata%input_cfg%ForcingName, &
                         simdata%grid)

   ! Initialize albedo data used for water albedo calculation, if switch is off
   if (simdata%model_cfg%user_defined_water_albedo) then
      simdata%model%albedo_water = simdata%model_param%wat_albedo
   else
      call mod_forcing%init_albedo(simdata%model, simdata%sim_cfg)
   end if

   ! Initialize absorption module
   call mod_absorption%init(simdata%model_cfg, &
                            simdata%model_param, &
                            simdata%input_cfg%AbsorpName, &
                            simdata%grid)

   ! Initialize biochemical model "AED2" if used
   if (simdata%model_cfg%couple_aed2) then
      call mod_aed2%init(simdata%model, simdata%grid, simdata%model_cfg, simdata%aed2_cfg)
   end if

   ! If there is advection (due to inflow)
   if (simdata%model_cfg%inflow_mode > 0) then
      ! initialize advection module
      call mod_advection%init(simdata%model, simdata%model_cfg, simdata%model_param, simdata%grid)

      ! initialize lateral module based on configuration
      if (simdata%model_cfg%inflow_mode == 1) then

         ! Gravity based inflow
         mod_lateral => mod_lateral_normal
      else if (simdata%model_cfg%inflow_mode == 2) then
         ! User defined inflow depths
         mod_lateral => mod_lateral_rho
      end if
      call mod_lateral%init(simdata%model, simdata%model_cfg, simdata%input_cfg, simdata%aed2_cfg, simdata%model_param, simdata%grid)
   else
      call warn('Lake in-/outflow is turned off')
   end if

   ! Binary simulation snapshot file
   snapshot_file_path = simdata%output_cfg%PathOut//'/simulation-snapshot.dat'

   ! Text output files for simulation restart
   file_text_restart = simdata%output_cfg%PathOut//'/initial_conditions_for_restart.dat'
   file_text_restart2 = simdata%output_cfg%PathOut//'/seiche_ice_for_restart.dat'

   if (simdata%sim_cfg%continue_from_snapshot) then
      inquire (file=snapshot_file_path, exist=continue_from_snapshot)
      print *,"Snapshot is available and used",continue_from_snapshot
   end if

   ! Setup logger
   call logger%initialize(simdata%model, simdata%sim_cfg, simdata%model_cfg, simdata%aed2_cfg, simdata%output_cfg, simdata%grid, continue_from_snapshot)

   ! Calculate simulation_end_time, which is a tuple of integers (days, seconds)

   ! If output times are at regular intervals
   if (simdata%output_cfg%thinning_interval > 0) then
      ! Compute number of simulation days
      simulation_end_time(1) = int(floor(simdata%sim_cfg%end_datum - simdata%sim_cfg%start_datum))
      ! Compute number of simulation seconds (in addition to the days calculated above)
      simulation_end_time(2) = (simdata%sim_cfg%end_datum - simdata%sim_cfg%start_datum - real(simulation_end_time(1), RK)) * SECONDS_PER_DAY + 0.5
   
   ! If output times are user defined
   else
      simulation_end_time = simdata%output_cfg%simulation_times_for_output(:, &
            size(simdata%output_cfg%simulation_times_for_output,2))
      
      if (simdata%sim_cfg%end_datum < simdata%sim_cfg%start_datum + real(simulation_end_time(1)) + real(simulation_end_time(2))/SECONDS_PER_DAY) then
         call error('Some of the output times are larger than the simulation duration')
      end if
   end if

   ! Initialize simulation modules
   call mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)
   call mod_turbulence%init(simdata%model, simdata%grid, simdata%model_cfg, simdata%model_param)
   call mod_ice%init(simdata%model, simdata%model_cfg, simdata%model_param, simdata%grid)

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
      call ok("Start day: "//real_to_str(simdata%sim_cfg%start_datum, '(F7.1)'))
      new_start_datum = simdata%sim_cfg%start_datum
      if (continue_from_snapshot) then
         call load_snapshot(snapshot_file_path, simdata%model_cfg%couple_aed2)
         call ok("Simulation snapshot successfully read. Snapshot day: "//real_to_str(simdata%model%datum, '(F7.1)'))
         call logger%calculate_simulation_time_for_next_output(simdata%model%simulation_time)
         new_start_datum = simdata%model%datum
      else
         call logger%log(simdata)
      end if
      call ok("End day: "//real_to_str(simdata%sim_cfg%end_datum, '(F7.1)'))
      call logger%start()

      ! initialize a bar with the progress percentage counter
      if (simdata%sim_cfg%show_bar) then
         call bar%initialize(filled_char_string='#', &
            prefix_string=' Simulation progress |',  &
            suffix_string='| ', add_progress_percent=.true., &
            add_date_time=.true., &
            max_value=(simdata%sim_cfg%end_datum-new_start_datum))

            ! start the progress bar
         call bar%start
      end if

      ! Run the simulation loop
      ! Run simulation until end datum or until no more results are required by the output time file
      ! Simulation time is a tuple of 2 integers (days, seconds)
      do while (simdata%model%simulation_time(1) < simulation_end_time(1) .or. &
         simdata%model%simulation_time(1) == simulation_end_time(1) .and. simdata%model%simulation_time(2) < simulation_end_time(2))

         ! Advance to the next timestep
         simdata%model%simulation_time_old = simdata%model%simulation_time

         simdata%model%simulation_time(2) = simdata%model%simulation_time(2) + simdata%sim_cfg%timestep

         ! If second counter is larger than 86400, change day
         if (simdata%model%simulation_time(2) >= SECONDS_PER_DAY) then
            simdata%model%simulation_time(1) = simdata%model%simulation_time(1) + 1
            simdata%model%simulation_time(2) = simdata%model%simulation_time(2) - SECONDS_PER_DAY
         end if

         simdata%model%datum = datum(simdata%sim_cfg%start_datum, simdata%model%simulation_time)

         ! ************************************
         ! ***** Compute next model state *****
         ! ************************************

         ! Update water albedo
         if (.not. simdata%model_cfg%user_defined_water_albedo) then
            call mod_forcing%update_albedo(simdata%model)
         end if

         ! Update forcing
         call mod_forcing%update(simdata%model)

         ! Update absorption (except if AED2 is off or if AED2 is on but bioshade feedback is off)
         if (simdata%model_cfg%couple_aed2) then
            if (.not. simdata%aed2_cfg%bioshade_feedback) then
               call mod_absorption%update(simdata%model)
            end if
         else
            call mod_absorption%update(simdata%model)
         end if
!write(6,*) 'before stab',simdata%model%k(simdata%grid%ubnd_fce),simdata%model%eps(simdata%grid%ubnd_fce),simdata%model%num(simdata%grid%ubnd_fce),simdata%model%nuh(simdata%grid%ubnd_fce),&
!simdata%model%NN(simdata%grid%ubnd_fce), simdata%model%rho(simdata%grid%ubnd_vol)
         ! Update physics
         call mod_stability%update(simdata%model)
!write(6,*) 'after stab',simdata%model%k(simdata%grid%ubnd_fce),simdata%model%eps(simdata%grid%ubnd_fce),simdata%model%num(simdata%grid%ubnd_fce),simdata%model%nuh(simdata%grid%ubnd_fce),&
!simdata%model%NN(simdata%grid%ubnd_fce), simdata%model%rho(simdata%grid%ubnd_vol)

         ! If there is inflow/outflow do advection part
         if (simdata%model_cfg%inflow_mode > 0) then
            ! Treat inflow/outflow
            call mod_lateral%update(simdata%model)
            ! Set old lake level (before it is changed by advection module)
            simdata%grid%lake_level_old = simdata%grid%z_face(simdata%grid%ubnd_fce)

            ! Update lake advection using the inflow/outflow data
            call mod_advection%update(simdata%model)
            ! Update lake level
            simdata%grid%lake_level = simdata%grid%z_face(simdata%grid%ubnd_fce)
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

         ! Update biogeochemistry
         if (simdata%model_cfg%couple_aed2) then
            call mod_aed2%update(simdata%model)
         end if

         ! Call logger to write files
         call logger%log(simdata)

         ! This logical is used to do some allocation in the forcing, absorption and lateral subroutines during the first timestep
         simdata%model%first_timestep = .false.

         !update the progress bar
         if (simdata%sim_cfg%show_bar) then
            call bar%update(current=(simdata%model%datum-new_start_datum))
         end if

      end do
      if (simdata%sim_cfg%continue_from_snapshot) call save_snapshot(snapshot_file_path, simdata%model_cfg%couple_aed2)
      if (simdata%sim_cfg%save_text_restart) call save_restart(file_text_restart, file_text_restart2)
   end subroutine

   ! Function to save the necessary variables for a restart from a text file (similar to initial conditions)
   ! The depth grid includes both grid components (volume and surface) of the staggered grid
   subroutine save_restart(file_path_1, file_path_2)
      character(len=*), intent(in) :: file_path_1
      character(len=*), intent(in) :: file_path_2

      logical :: status_ok
      type(csv_file), dimension(2):: save_files
      real(RK) :: total_grid(size(simdata%grid%z_face) + size(simdata%grid%z_volume))
      real(RK),dimension(size(simdata%grid%z_face(1:simdata%grid%ubnd_fce)) + size(simdata%grid%z_volume(1:simdata%grid%ubnd_vol))) :: depth, U, V, T, S, k, eps, num, nuh
      integer :: i
      real(RK),dimension(9) :: row1
      real(RK),dimension(4) :: row2

      call save_files(1)%open(file_path_1, n_cols=9,status_ok=status_ok)
      call save_files(1)%add(["depth (m) ","u (m/s)   ","v (m/s)   ","T (°C)   ","S (g/kg)  ","k (J/kg)  ","eps (W/kg)","num (m2/s)","nuh (m2/s)"])
      call save_files(1)%next_row()
      depth(1) = simdata%grid%z_face(1)

      do i=2,2*size(simdata%grid%z_volume(1:simdata%grid%ubnd_vol)),2
         depth(i+1) = simdata%grid%z_face(i/2+1)
         depth(i) = simdata%grid%z_volume(i/2)
      end do

      call Interp(simdata%grid%z_volume(1:simdata%grid%ubnd_vol), simdata%model%U(1:simdata%grid%ubnd_vol), size(simdata%model%U), depth, U, size(U))
      call Interp(simdata%grid%z_volume(1:simdata%grid%ubnd_vol), simdata%model%V(1:simdata%grid%ubnd_vol), size(simdata%model%V), depth, V, size(V))
      call Interp(simdata%grid%z_volume(1:simdata%grid%ubnd_vol), simdata%model%T(1:simdata%grid%ubnd_vol), size(simdata%model%T), depth, T, size(T))
      call Interp(simdata%grid%z_volume(1:simdata%grid%ubnd_vol), simdata%model%S(1:simdata%grid%ubnd_vol), size(simdata%model%S), depth, S, size(S))
      call Interp(simdata%grid%z_face(1:simdata%grid%ubnd_fce), simdata%model%k(1:simdata%grid%ubnd_vol), size(simdata%model%k), depth, k, size(k))
      call Interp(simdata%grid%z_face(1:simdata%grid%ubnd_fce), simdata%model%eps(1:simdata%grid%ubnd_vol), size(simdata%model%eps), depth, eps, size(eps))
      call Interp(simdata%grid%z_face(1:simdata%grid%ubnd_fce), simdata%model%num(1:simdata%grid%ubnd_vol), size(simdata%model%num), depth, num, size(num))
      call Interp(simdata%grid%z_face(1:simdata%grid%ubnd_fce), simdata%model%nuh(1:simdata%grid%ubnd_vol), size(simdata%model%nuh), depth, nuh, size(nuh))

      call reverse_in_place(depth)
      depth = depth - depth(1)
      call reverse_in_place(U)
      call reverse_in_place(V)
      call reverse_in_place(T)
      call reverse_in_place(S)
      call reverse_in_place(k)
      call reverse_in_place(eps)
      call reverse_in_place(num)
      call reverse_in_place(nuh)

      do i=1,size(simdata%grid%z_face(1:simdata%grid%ubnd_fce)) + size(simdata%grid%z_volume(1:simdata%grid%ubnd_vol))
         row1 = [depth(i),U(i),V(i),T(i),S(i),k(i),eps(i),num(i),nuh(i)]
         call save_files(1)%add(row1, real_fmt="(ES24.15)")
         call save_files(1)%next_row()
      end do
      call save_files(1)%close(status_ok)

      call save_files(2)%open(file_path_2,n_cols=4,status_ok=status_ok)
      call save_files(2)%add(['E_seiche (J) ', 'BlackIceH (m)', 'WhiteIceH (m)', 'SnowH (m)    '])

      call save_Files(2)%next_row()
      row2 = [simdata%model%E_seiche,simdata%model%black_ice_h,simdata%model%white_ice_h,simdata%model%snow_h]
      call save_files(2)%add(row2,real_fmt="(ES24.12)")

      call save_files(2)%close(status_ok)
   end subroutine

   subroutine save_snapshot(file_path, couple_aed2)
      implicit none
      character(len=*), intent(in) :: file_path
      logical, intent(in) :: couple_aed2

      open(80, file=file_path, Form='unformatted', Action='Write')
      call simdata%model%save(couple_aed2)
      call simdata%grid%save()
      call mod_absorption%save()
      if (simdata%model_cfg%inflow_mode > 0) then
         call mod_lateral%save()
      end if
      call logger%save()
      close(80)
   end subroutine

   subroutine load_snapshot(file_path, couple_aed2)
      implicit none
      character(len=*), intent(in) :: file_path
      logical, intent(in) :: couple_aed2

      open(81, file=file_path, Form='unformatted', Action='Read')
      call simdata%model%load(couple_aed2)
      call simdata%grid%load()
      call mod_absorption%load()
      if (simdata%model_cfg%inflow_mode > 0) then
         call mod_lateral%load()
      end if
      call logger%load()
      close(81)
   end subroutine

end program simstrat_main
