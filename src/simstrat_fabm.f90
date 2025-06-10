! ---------------------------------------------------------------------------------
!     Simstrat a physical 1D model for lakes and reservoirs
!
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     Copyright (C) 2020, Eawag
!     FABM: Copyright (C) 2014, Bolding & Bruggeman ApS
!     GOTM: Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!           Original author(s): Jorn Bruggeman (rem.: for every subroutine)
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
!<    +---------------------------------------------------------------+
!     |  Simstrat - FABM interface
!<    +---------------------------------------------------------------+

!#include "cppdefs.h"
!#include "fabm_version.h"

module simstrat_fabm
   use strat_simdata
   use strat_grid
   use strat_solver
   use utilities
   use fabm
   use fabm_types
   use fabm_expressions
   use fabm_config
   use fabm_driver
   if _FABM_API_VERSION_ > 0
      use fabm_v0_compatibility
   end if
   use yaml_settings

   ! use aed2_common; aed2_core

   ! use input,only: register_input, type_scalar_input, type_profile_input
   ! use settings
   ! use field_manager

   implicit none
   private

   !Public member functions
   public configure_gotm_fabm, gotm_fabm_create_model, init_gotm_fabm, init_gotm_fabm_state, start_gotm_fabm
   public set_env_gotm_fabm,do_gotm_fabm
   public clean_gotm_fabm
   public fabm_calc
   public freshwater_impact
   public register_observation
   public calculate_conserved_quantities, total0

   ! Passed through from fabm_types, used by hosts to provide additional inputs:
   public standard_variables

   ! Variables below must be accessible for gotm_fabm_output
   public dt,h,save_inputs

   ! Variables below must be accessible for getm_fabm
   public cc_transport

   ! Optional additional forcing
   public fabm_airp          ! air pressure in Pa
   public fabm_julianday     ! Julian day
   public fabm_calendar_date ! Subroutine that computes y/m/d from Julian day

   ! Error handling
   ! type,extends(type_base_driver) :: type_gotm_driver
   ! contains
   !     procedure :: fatal_error => gotm_driver_fatal_error
   !     procedure :: log_message => gotm_driver_log_message
   ! end type

   ! Private data members
   real(rk) :: cnpar
   integer :: w_adv_method,w_adv_discr,ode_method,split_factor,configuration_method
   logical :: fabm_calc,repair_state, &
      bioshade_feedback,bioalbedo_feedback,biodrag_feedback, &
      freshwater_impact,salinity_relaxation_to_freshwater_flux, &
      save_inputs, no_surface

   ! Arrays for observations, relaxation times and FABM variable identifiers associated with the observations.
   type type_1d_state_info
      real(RK),pointer,dimension(:) :: obs => null()
      real(RK),pointer,dimension(:) :: relax_tau => null()
      real(RK),allocatable,dimension(:) :: diff_flux
   end type

   type type_0d_state_info
      real(RK),pointer :: obs => null()
      real(RK),pointer :: relax_tau => null()
   end type
   
   real(RK),allocatable,dimension(:),target :: horizontal_expression_data

   ! Observation indices (from obs_0d, obs_1d) for pelagic and benthic state variables.
   type(type_1d_state_info),allocatable :: cc_info(:)
   type(type_0d_state_info),allocatable :: cc_ben_info(:)

   interface register_observation
      module procedure register_bulk_observation
      module procedure register_horizontal_observation
      module procedure register_scalar_observation
   end interface

   type(type_bulk_variable_id),save :: temp_id,salt_id,rho_id,h_id,swr_id,par_id,pres_id,nuh_id
   type(type_horizontal_variable_id),save :: lon_id,lat_id,windspeed_id,par_sf_id,cloud_id,taub_id,swr_sf_id

   ! Some parameters
   logical,save :: save_diag = .true.
   logical,save :: compute_light = .false.

   integer,parameter :: gotmrk = kind(_ONE_)

   ! Calendar date interface
   interface 
   subroutine calendar_date_interface(julian,yyyy,mm,dd)
      integer :: julian
      integer :: yyyy,mm,dd
   end subroutine
   end interface   

   procedure(calendar_date_interface),pointer :: fabm_calendar_date

   ! Input (all from GOTM)
   public configure_simstrat_fabm_input, init_simstrat_fabm_input
   public type_input_variable, first_input_variable
   integer,parameter :: maxpathlen = 256
   integer,parameter :: max_variable_count_per_file = 256
   ! Information on an observed variable
   type type_input_variable
      !type (type_scalar_input) :: scalar_input
      !type (type_profile_input) :: profile_input
      type (type_bulk_variable_id) :: interior_id          ! FABM identifier of pelagic variable (not associated if variable is not pelagic)
      type (type_horizontal_variable_id) :: horizontal_id        ! FABM identifier of horizontal variable (not associated if variable is not horizontal)
      type (type_scalar_variable_id) :: scalar_id            ! FABM identifier of scalar variable (not associated if variable is not scalar)
      integer :: ncid = -1 ! NetCDF id in output file (only used if this variable is included in output)
      real(RK) :: relax_tau ! Relaxation times
      real(RK) :: relax_tau_bot ! Relaxation times for bottom layer (depth-dependent variables only)
      real(RK) :: relax_tau_surf ! Relaxation times for surface layer (depth-dependent variables only)
      real(RK) :: h_bot, h_surf ! Thickness of bottom and surface layers (for relaxation rates that vary per layer)
      real(RK),allocatable,dimension(:) :: relax_tau_1d ! Relaxation times for profiles (depth-dependent variables)
      type(type_input_variable),pointer :: next => null() ! Next variable in current input file
   end type
   type(type_input_variable),pointer,save :: first_input_variable => null()
   type(type_input_variable),pointer,save :: last_input_variable  => null()
   type,extends(type_dictionary_populator) :: type_fabm_input_populator
   contains
      procedure :: create => fabm_input_create
   end type
   type (type_fabm_input_populator) :: fabm_input_populator

   ! Type for use of FABM
   type,public :: SimstratFABM
      class(fabm_config),pointer :: fabm_cfg
      class(StaggeredGrid),pointer :: grid
      class(type_fabm_model),pointer :: model
      !class(type_fabm_model),pointer,save,public :: model => null()

      ! Arrays for state and diagnostic variables
      real(RK),pointer,dimension(:,:) :: cc ! water quality array: nlayers, nvars
      !real(RK),allocatable,dimension(:,:),target :: cc
      real(RK),allocatable,dimension(:,:) :: cc_diag
      real(RK),allocatable,dimension(:) :: cc_diag_hz
      !real(RK),allocatable,dimension(:) :: cc_diag_sheet
      !real(RK),pointer,dimension(:) :: tss
      !real(RK),pointer,dimension(:) :: sed_zones

      ! Variables to hold time spent on advection, diffusion, sink/source terms.
      integer(8) :: clock_adv,clock_diff,clock_source

      ! Arrays for fluxes of state variables
      real(RK),pointer,dimension(:) :: flux_atm
      real(RK),pointer,dimension(:) :: flux_ben
      real(RK),pointer,dimension(:,:) :: flux_pel
      real(RK),pointer,dimension(:,:) :: flux_zone

      ! Arrays for work, vertical movement, and cross-boundary fluxes
      real(RK),allocatable,dimension(:,:) :: ws
      real(RK),allocatable,dimension(:)   :: total,local, &
         sfl,bfl,Qsour,Lsour,DefaultRelaxTau,cc_old,curh,curnuh,iweights 
      logical,allocatable, dimension(:)   :: cc_transport
      integer,allocatable, dimension(:)   :: posconc

      ! Arrays for environmental variables not supplied externally.
      real(RK),pointer,dimension(:) :: par,pres,swr,k_par,z,nuh_ct, &
         uva,uvb,nir
      !real(RK),allocatable,dimension(:),target :: par,pres,swr,k_par,z,nuh_ct,uva,uvb,nir

      ! Column pointers; type (aed2_column_t),pointer,dimension(:) :: column, column_sed

      ! External variables
      real(RK) :: dt,dt_eff   ! External and internal time steps
      integer  :: w_adv_ctr   ! Scheme for vertical advection (0 if not used)
      real(RK),pointer,dimension(:) :: nuh,h,bioshade,w,rho, &
         SRelaxTau,sProf,salt
      real(RK),pointer :: precip,evap,bio_drag_scale,bio_albedo

      ! Some more variables
      real(RK),pointer :: I_0,A,g1,g2, &
         decimal_yearday
      !real(RK), target :: decimal_yearday
      integer,pointer :: yearday,secondsofday
      logical :: fabm_ready

      real(RK),pointer :: fabm_airp
      integer,pointer :: fabm_julianday

      logical :: check_conservation
      real(RK),allocatable :: local(:,:)
      real(RK),allocatable :: total0(:)
      real(RK),allocatable :: change_in_total(:)

      integer :: repair_interior_count
      integer :: repair_surface_count
      integer :: repair_bottom_count

      character(len=256),public :: yaml_file = 'fabm.yaml'

      character(len=48),pointer :: names(:)
      character(len=48),pointer :: bennames(:)
      character(len=48),pointer :: diagnames(:)
      character(len=48),pointer :: diagnames_sheet(:)

      integer,allocatable,dimension(:) :: externalid

      real(RK),allocatable,dimension(:) :: min_, max_

      integer :: n_AED2_state_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
      integer :: zone_var = 0

   contains
      procedure,pass(self),public :: init
      procedure,pass(self),public :: update
   end type SimstratFABM

contains
   !include 'simstrat_aed2_subroutines.f90'
   !include 'simstrat_aed2_physics.f90'
   !#include "../src/util/ode_solvers_template.F90"

   ! Initialize (sets up memory, reads FABM configuration from fabm.yaml,
   ! links the external Simstrat variables and sets the initial conditions of FABM variables)
   ! After this the number of biogeochemical variables is fixed.
   ! (access variable metadata in model%interior_state_variables, model%interior_diagnostic_variables)
   subroutine init(self, state, grid, model_cfg, fabm_cfg)
      use yaml_settings
      !use util, only: UPSTREAM, P2, MUSCL, Superbee, P2_PDM
      
      ! Arguments
      class(SimstratFABM) :: self
      class(ModelState) :: state
      class(StaggeredGrid),target :: grid
      class(ModelConfig),target :: model_cfg
      class(fabm_config),target :: fabm_cfg
      class(type_settings),intent(inout) :: cfg
      type(type_settings),pointer :: branch

      ! Local variables
      !integer :: i, output_level
      !logical :: file_exists, in_output
      !character(len=80) :: fname
      !character(len=64) :: models(64)
      !namelist /aed2_models/ models
      !type(aed2_variable_t),pointer :: tvar
      !integer i, status, av, v, sv

      ! Configuration
      associate (n_AED2_state_vars => self%n_AED2_state_vars, &
             n_vars => self%n_vars, &
             n_vars_ben => self%n_vars_ben, &
             n_vars_diag => self%n_vars_diag, &
             n_vars_diag_sheet => self%n_vars_diag_sheet)

         ! ! AED2 config file
         ! fname = aed2_cfg%aed2_config_file

         ! if ( aed2_init_core('.') /= 0 ) call error("Initialisation of aed2_core failed")
         ! call aed2_print_version

         ! ! Create model tree
         ! write (6,*) "     Processing aed2_models config from ", trim(fname)
         ! open(50,file=fname,action='read',status='old',iostat=status)
         ! if ( status /= 0 ) then
         !    call error("Cannot open file " // trim(fname))
         !    stop
         ! end if

         ! models = ''
         ! read(50, nml=aed2_models, iostat=status)
         ! if ( status /= 0 ) then
         !    call error("Cannot read namelist entry aed2_models")
         !    stop
         ! end if

         ! do i=1,size(models)
         !    if (models(i)=='') exit
         !    call aed2_define_model(models(i), 50)
         ! end do

         ! ! Finished reading AED2 config
         ! close(50)
         ! write (6,*) "      AED2 file parsing completed."

         ! ! Assign number of different variables
         ! n_AED2_state_vars = aed2_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

         ! ! Print variable information to screen
         ! print "(/,5X,'AED2 : n_AED2_state_vars = ',I3,' ; MaxLayers         = ',I4)",n_AED2_state_vars,self%grid%nz_grid
         ! print "(  5X,'AED2 : n_vars      = ',I3,' ; n_vars_ben        = ',I3)",n_vars,n_vars_ben
         ! print "(  5X,'AED2 : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I3,/)",n_vars_diag,n_vars_diag_sheet

         ! ! Check variable dependencies
         ! call check_data(self)

         ! ! Allocate space for the allocatables/pointers of this module
         ! call allocate_memory(self)

         ! ! Allocate memory for AED2 state and inflow matrix used by Simstrat
         ! state%AED2_state => self%cc
         ! state%AED2_diagnostic => self%cc_diag
         ! state%AED2_diagnostic_sheet => self%cc_diag_sheet
         ! state%n_AED2_state = n_vars  + n_vars_ben
         ! state%n_AED2_diagnostic = n_vars_diag
         ! state%n_AED2_diagnostic_sheet = n_vars_diag_sheet

         ! ! Define column pointer (which is the object that is handed over to AED2 at every timestep)
         ! ! It containes external (Simstrat) variables like T and S, but also the variables of this (SimstratAED2) module
         ! call define_column(self, state)
         ! !if (benthic_mode .GT. 1) call define_sed_column(column_sed, n_zones, flux, flux_atm, flux_ben)

         ! ! Assign name, min and max values of variables, print names to screen
         ! call assign_var_names(self)
         ! allocate(state%AED2_state_names(n_vars))
         ! state%AED2_state_names => self%names

         ! allocate(state%AED2_diagnostic_names(n_vars_diag))
         ! state%AED2_diagnostic_names => self%diagnames

         ! allocate(state%AED2_diagnostic_names_sheet(n_vars_diag_sheet))
         ! state%AED2_diagnostic_names_sheet => self%diagnames_sheet

         ! ! Now set initial values of AED2 variables
         ! v = 0 ; sv = 0;
         ! do av=1,self%n_AED2_state_vars
         !    if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"
         !    if ( .not. ( tvar%extern .or. tvar%diag) ) then  ! neither global nor diagnostic variable
         !       if ( tvar%sheet ) then
         !          sv = sv + 1
         !          call AED2_InitCondition(self, self%cc(:, n_vars + sv), tvar%name, tvar%initial)
         !       else
         !          v = v + 1
         !          call AED2_InitCondition(self, self%cc(:, v), tvar%name, tvar%initial)
         !       end if
         !    end if
         ! end do

         ! write(*,"(/,5X,'----------  AED2 config : end  ----------',/)")
      end associate

      ! Add grid and fabm_cfg to SimstratFABM object
      self%grid => grid
      self%fabm_cfg => fabm_cfg

      ! send message
      !LEVEL1 'init_gotm_fabm_yaml'
      ! Initialize all configuration variables with the get subroutine from fabm/yaml/yaml_settings.F90
      call cfg%get(fabm_calc, 'use', 'enable FABM', &
                  default=.false.)
      call cfg%get(yaml_file, 'yaml_file', 'FABM configuration file', &
                  default='fabm.yaml', display=display_advanced)
      call cfg%get(freshwater_impact, 'freshwater_impact', 'enable dilution/concentration by precipitation/evaporation', &
                  default=.true.) ! disable to check mass conservation
      branch => cfg%get_child('feedbacks', 'feedbacks to physics')
      call branch%get(bioshade_feedback, 'shade', 'interior light absorption', &
                  default=.false.)
      call branch%get(bioalbedo_feedback, 'albedo', 'surface albedo', &
                  default=.false.)
      call branch%get(biodrag_feedback, 'surface_drag', 'surface drag', &
                  default=.false.)
      call cfg%get(repair_state, 'repair_state', 'clip state to minimum/maximum boundaries', &
                  default=.false.)
      branch => cfg%get_child('numerics', display=display_advanced)
      call branch%get(ode_method, 'ode_method', 'time integration scheme applied to source terms', &
                  options=(/ option(1, 'Forward Euler', 'FE'), option(2, 'Runge-Kutta 2', 'RK2'), option(3, 'Runge-Kutta 4', 'RK4'), &
                  option(4, 'first-order Patanker', 'Patankar1'), option(5, 'second-order Patanker', 'Patankar2'), option(7, 'first-order modified Patanker', 'MP1'), &
                  option(8, 'second-order modified Patanker', 'MP2'), option(10, 'first-order extended modified Patanker', 'EMP1'), &
                  option(11, 'second-order extended modified Patankar', 'EMP2') /), default=1)
      call branch%get(split_factor, 'split_factor', 'number of substeps used for source integration', &
                  minimum=1,maximum=100,default=1)
      call branch%get(w_adv_discr, 'w_adv_discr', 'vertical advection scheme for settling/rising', options=&
               (/ option(UPSTREAM, 'first-order upstream', 'upstream'), option(P2, 'third-order upstream-biased polynomial', 'P2'), &
                  option(Superbee, 'third-order TVD with Superbee limiter', 'Superbee'), option(MUSCL, 'third-order TVD with MUSCL limiter', 'MUSCL'), &
                  option(P2_PDM, 'third-order TVD with ULTIMATE QUICKEST limiter', 'P2_PDM') /), default=P2_PDM)
      call branch%get(cnpar, 'cnpar', '"implicitness" of diffusion scheme', '1', &
                  minimum=_ZERO_,default=_ONE_)
      if 0 then
         call cfg%get(salinity_relaxation_to_freshwater_flux, 'salinity_relaxation_to_freshwater_flux', '', &
                     default=.false.)
      else
         salinity_relaxation_to_freshwater_flux = .false.
      end if
      branch => cfg%get_child('debug', display=display_advanced)
      call branch%get(no_surface, 'no_surface', 'disable surface processes', &
                  default=.false.) ! disables surface exchange; useful to check mass conservation
      call branch%get(save_inputs, 'save_inputs', 'include additional forcing fields in output', &
                  default=.false.)
      call cfg%get(configuration_method, 'configuration_method', 'configuration file', &
                  options=(/option(-1, 'auto-detect (prefer fabm.yaml)', 'auto'), option(0, 'fabm.nml', 'nml'), option(1, 'fabm.yaml', 'yaml')/), &
                  default=-1, display=display_advanced)
      ! Send message that reading is done
      !if (fabm_calc) then
      !   LEVEL2 'Reading configuration from:'
      !   LEVEL3 trim(yaml_file)
      !end if
      !LEVEL2 'done.'

      ! Provide FABM with an object for communication with host, abort if FABM is not enabled
      allocate(type_gotm_driver::driver)
      if (.not. fabm_calc) return
      fabm_ready = .false.

      ! Initialize the model, handled by fabm. or:
      model => fabm_create_model()
      ! create from YAML file with subroutine from fabm/fabm_v0_compatibility
      if (_FABM_API_VERSION_ > 0) then
         allocate(model)
         call fabm_create_model_from_yaml_file(model,trim(yaml_file))
      else
         ! Create model tree
         if (configuration_method==-1) then
            configuration_method = 1
            inquire(file=trim(yaml_file),exist=file_exists)
            if (.not.file_exists) then
               inquire(file='fabm.nml',exist=file_exists)
               if (file_exists) configuration_method = 0
            end if
         end if
         select case (configuration_method)
         case (0)
            model => fabm_create_model_from_file(namlst)
         case (1)
            allocate(model)
            call fabm_create_model_from_yaml_file(model,path=trim(yaml_file))
         end select   
      end if

      ! Initialize the Simstrat-FABM driver by reaing settings from fabm/yaml
      ! Input parameters
      integer,intent(in) :: nlev
      !character(len=*),intent(in) :: fname
      real(rk),optional,intent(in) :: dt
      class(type_field_manager),intent(inout),optional :: field_manager
      ! Local variables
      integer :: i
      logical :: in_output
      ! Send message
      !LEVEL1 'post_init_gotm_fabm'
      if (fabm_calc) then
         clock_adv    = 0
         clock_diff   = 0
         clock_source = 0
         repair_interior_count = 0
         repair_surface_count = 0
         repair_bottom_count = 0
      end if

      ! Provide extents of the spatial domain (number of layers nz for a 1D column) or:
      call model%set_domain(nz)
      ! initialize model tree (creates metadata and assigns variable identifiers)
      ! fabm v0 compatibiity subroutines
      call fabm_set_domain(model,nlev,dt)
      if (_FABM_API_VERSION_ == 0) then
            call model%set_bottom_index(1)
            call model%set_surface_index(nlev)
      end if

      ! Point FABM to environmental data
      ! Do this for all variables on FABM's standard variable list that the model can provide.
      ! In Simstrat: ModelVariable contains pointer var to variable data

      ! allocate in simstrat: real, dimension(nz), target :: temp, salt
      call model%link_interior_data(fabm_standard_variables%temperature, state%model%T)
      call model%link_interior_data(fabm_standard_variables%practical_salinity, state%model%S)

      ! allocate in simstrat: real, target :: wind, yearday
      call model%link_horizontal_data(fabm_standard_variables%wind_speed, wind)
      call model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, yearday)

      ! Allocate memory to hold the values of all size(model%*_state_variables) state variables.
      ! Combine all state variable values in an array *_state with shape (nz,)size(model%*_state_variables).
      real, dimension(:,:), target, allocatable :: state
      real, dimension(:), target, allocatable :: bottom_state, surface_state
      allocate(state(nz, size(model%interior_state_variables)))
      allocate(bottom_state (size(model%bottom_state_variables)))
      allocate(surface_state(size(model%surface_state_variables)))

      ! Point FABM to your state variable data
      do ivar=1,size(model%interior_state_variables)
         call model%link_interior_state_data(ivar, state(:,:,:,ivar))
      end do
      do ivar=1,size(model%bottom_state_variables)
         call model%link_bottom_state_data(ivar, bottom_state(:,:,ivar))
      end do
      do ivar=1,size(model%surface_state_variables)
         call model%link_surface_state_data(ivar, surface_state(:,:,ivar))
      end do

      ! Complete initialization and check whether FABM has all dependencies fulfilled
      ! (i.e., whether all required calls to model%link_*_data have been made)
      call model%start()

      ! Initialize the tracers
      ! This sets the values of arrays sent to model%link_*_state_data,
      ! in this case those contained in *_state.
      call model%initialize_interior_state(1, nz)
      call model%initialize_surface_state() 
      call model%initialize_bottom_state()
   end subroutine

   ! The update function is called in the main loop of simstrat (in simstrat.f90) at each timestep
   ! Particle mobility (sedimentation), light absorption feedback by AED2 variables, atmospheric,
   ! pelagic and benthic fluxes, advection and diffusion are computed.
   subroutine update(self, state)
      use,intrinsic :: ieee_arithmetic

      implicit none
      class(SimstratFABM) :: self
      class(ModelState) :: state

      ! Local variables
      !type(aed2_variable_t),pointer :: tvar
      !real(RK) :: min_C
      !integer :: v, i, lev, r
      !real(RK), dimension(self%grid%ubnd_vol) :: tmp

      ! Calculate local pressure
      self%pres(1:self%grid%ubnd_vol) = -self%grid%z_volume(1:self%grid%ubnd_vol)

      self%cc_diag = 0.
      self%cc_diag_sheet = 0.

      ! (3) Calculate source/sink terms due to settling rising of state
      ! variables in the water column (note that settling into benthos
      ! is done in aed2_do_benthos)
      if (self%aed2_cfg%particle_mobility) then
         v = 0
         do i = 1,self%n_AED2_state_vars
            if ( aed2_get_var(i, tvar) ) then
               if ( .not. (tvar%sheet .or. tvar%diag .or. tvar%extern) ) then
               v = v + 1
               ! only for state_vars that are not sheet
                  if ( .not. ieee_is_nan(tvar%mobility) ) then
                     self%ws(:, v) = tvar%mobility
                     min_C = tvar%minimum
                     call Mobility(self, state, min_C, self%ws(:, v), self%cc(:, v))
                  end if
               end if
            end if
         end do
      end if

      call check_states(self)

      ! Compute shading of AED2 variables. If bioshade feedback is off, then the "normal" absorption is computed in the main loop of Simstrat
      if (self%aed2_cfg%bioshade_feedback) then
         call absorption_updateAED2(self, state)
      end if

      self%par(:) = state%rad_vol(:)*rho_0*cp

      ! Calculate irradiance spectrum from par (factors from GLM)
      self%nir(:) = (self%par(:)/0.45) * 0.51
      self%uva(:) = (self%par(:)/0.45) * 0.035
      self%uvb(:) = (self%par(:)/0.45) * 0.005

      call calculate_fluxes(self, state)

      ! Update the water column layers using the biochemical reaction of AED2
      do v = 1, self%n_vars
         do lev = 1, self%grid%nz_occupied
            self%cc(lev, v) = self%cc(lev, v) + state%dt*self%flux_pel(lev, v)
         end do
      end do

      ! Now update benthic variables, depending on whether zones are simulated
      if ( self%aed2_cfg%benthic_mode .gt. 1 ) then
         call error("The use of sediment zones is currently not implemented in Simstrat-AED2")
           ! Loop through benthic state variables to update their mass
           !do v = n_vars+1, n_vars+n_vars_ben
              ! Loop through each sediment zone
              !do lev = 1, n_zones
                 ! Update the main cc_sed data array with the
                 !z_cc(lev, v) = z_cc(lev, v)+ dt_eff*flux_zone(lev, v)
              !end do
           !end do
      else
         do v = self% n_vars + 1, self%n_vars + self%n_vars_ben
            self%cc(1, v) = self%cc(1, v) + state%dt*self%flux_ben(v)
         end do
      end if

      ! If simulating sediment zones, distribute cc-sed benthic properties back
      ! into main cc array, mainly for plotting
      !if ( benthic_mode .GT. 1 ) call copy_from_zone(cc, cc_diag, cc_diag_sheet, wlev)

      ! Diffusive transport of AED2 variables (advective transport of AED2 variables is done in the usual Simstrat routines (lateral/lateral_rho))
      do v=1, self%n_vars
         call diffusion_AED2_state(self, state, v)
      end do

      ! Validate the model state, stop if variables are not valid.
      ! Only run after first time step has been run.
      logical, parameter :: repair = .false.
      logical :: valid_int, valid_sf, valid_bt
      call model%check_interior_state(1, nz, repair, valid_int)
      call model%check_surface_state(repair, valid_sf)
      call model%check_bottom_state(repair, valid_bt)
      if (.not. (valid_int .and. valid_sf .and. valid_bt) .and. .not. repair) stop

      ! Prepare all fields FABM needs to compute source terms (e.g., light)
      call model%prepare_inputs()

      ! Retrieve fluxes and tracer source terms at surface.
      real(rk) :: flux_sf(size(model%interior_state_variables))
      real(rk) :: sms_sf(size(model%surface_state_variables))
      flux_sf = 0
      sms_sf = 0
      call model%get_surface_sources(flux_sf, sms_sf)

      ! Retrieve fluxes and tracer source terms at bottom.
      ! Consider that there is a bottom at every layer.
      real(rk) :: flux_bt(size(model%interior_state_variables))
      real(rk) :: sms_bt(size(model%bottom_state_variables))
      flux_bt = 0
      sms_bt = 0
      call model%get_bottom_sources(flux_bt, sms_bt)
      
      ! Retrieve interior tracer source terms (tracer units s-1).
      real(rk) :: sms(nz,size(model%interior_state_variables))
      sms = 0
      call model%get_interior_sources(1, nz, sms)

      ! Retrieve vertical velocities (sinking, floating, active movement) in m s-1.
      real(rk) :: velocity(nz,size(model%interior_state_variables))
      call model%get_vertical_movement(1, nz, velocity)

      ! Compute any remaining diagnostics
      call model%finalize_outputs()

      ! Here you would time-integrate the advection-diffusion-reaction equations
      ! of all tracers, combining the transport terms with the biogeochemical source
      ! terms dy and vertical velocities w. This should result in an updated interior_state.
   end subroutine

   ! Calculate fluxes (with AED2 methods)
   subroutine calculate_fluxes(self, state)
      use,intrinsic :: ieee_arithmetic

      ! Arguments
      class(SimstratAED2),intent(inout) :: self
      class(ModelState),intent(in) :: state

      ! Local variables
      integer :: lev,zon,v_start,v_end,av,sv,sd
      real(RK) :: scale
      !real(RK), dimension(self%grid%nz_occupied, self%n_vars)    :: flux_pel_pre
      !real(RK), dimension(self%aed2_cfg%n_zones, self%n_vars) :: flux_pel_z
      logical :: splitZone
      type(aed2_variable_t),pointer :: tvar
      
      ! Begin
      self%flux_pel = zero_
      self%flux_atm = zero_
      self%flux_ben = zero_

      ! Start with calculating all flux terms for rhs in mass/m3/s
      ! Includes (1) benthic flux, (2) surface exchange and (3) water column kinetics
      ! as calculated by glm

      ! (1) BENTHIC FLUXES
      if ( self%aed2_cfg%benthic_mode .gt. 1 ) then
         call error("The use of sediment zones is currently not implemented in Simstrat-AED2")
            ! Multiple static sediment zones are simulated, and therfore overlying
            ! water conditions need to be aggregated from multiple cells/layers, and output flux
            ! needs disaggregating from each zone back to the overlying cells/layers

            !do zon=1,self%aed2_cfg%n_zones
               ! Reinitialise flux_ben to be repopulated for this zone
               !flux_ben = zero_
               !flux_pel_pre = zero_

               ! If multiple benthic zones, we must update the benthic variable pointer for the new zone
               !if ( self%zone_var .ge. 1 ) then
                  !column_sed(zone_var)%cell_sheet => z_sed_zones(zon)
                  !MH WE NEED A COLUMN TO CC VAR MAP FOR BENTHIC GUYS
                  !CAB Yes, a map (or 2 maps) would be better, but QnD since this all needs reworking
                  !sv = 0 ; sd = 0
                  !do av=1,self%n_AED2_state_vars
                  !   if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"
                  !   if ( .not. tvar%extern .and. tvar%sheet ) then
                  !      if ( tvar%diag ) then
                  !         sd = sd + 1
                  !         column(av)%cell_sheet => z_diag_hz(zon, sd)
                  !      else
                  !         sv = sv + 1
                  !         column(av)%cell_sheet => z_cc(zon, self%n_vars + sv)
                  !      end if
                  !   end if
                  !end do
                  !print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
               !end if
               !if ( self%aed2_cfg%benthic_mode .eq. 3 ) then
                  ! Zone is able to operated on by riparian and dry methods
                  !call aed2_calculate_riparian(column_sed, zon, z_pc_wet(zon))
                  !if (z_pc_wet(zon) .eq. 0. ) call aed2_calculate_dry(column_sed, zon)
               !end if
               ! Calculate temporal derivatives due to benthic processes.
               ! They are stored in flux_ben (benthic vars) and flux_pel (water vars)
               !flux_pel_pre = flux_pel

               !print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
               !call aed2_calculate_benthic(column_sed, zon)

               ! Record benthic fluxes in the zone array
               !flux_zon(zon, :) = flux_ben(:)

               ! Now we have to find out the water column flux that occured and
               ! disaggregate it to relevant layers
               !flux_pel_z(zon,:) = flux_pel(zon,:)-flux_pel_pre(zon,:)
            !end do

            ! Disaggregation of zone induced fluxes to overlying layers
            !v_start = 1 ; v_end = self%n_vars
            !zon = self%aed2_cfg%n_zones
            !do lev=self%grid%nz_occupied,1,-1
               !if ( zon .ne. 1 ) then
                  !splitZone = zz(lev-1) < zone_heights(zon-1)
               !else
                  !splitZone = .FALSE.
               !end if

               !if (splitZone) then
                  !scale = (zone_heights(zon-1) - zz(lev-1)) / (zz(lev) - zz(lev-1))
                  !flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end) * scale

                  !zon = zon - 1

                  !flux_pel(lev,v_start:v_end) = flux_pel(lev,v_start:v_end) + &
                  !                              flux_pel_z(zon,v_start:v_end) * (1.0 - scale)
               !else
                  !flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end)
               !end if
            !end do
            ! Limit flux out of bottom waters to concentration of that layer
            ! i.e. don't flux out more than is there & distribute
            ! bottom flux into pelagic over bottom box (i.e., divide by layer height).
            ! scaled to proportion of area that is "bottom"
            !do lev=1,self%grid%nz_occupied
               !if(lev>1)flux_pel(lev, :) = flux_pel(lev, :) * (self%grid%Az_vol(lev) - self%grid%Az_vol(lev - 1))/self%grid%Az_vol(lev)
               !flux_pel(lev, :) = max(-1.0 * self%cc(lev, :), flux_pel(lev, :)/self%grid%h(lev))
            !end do
      else
         ! Sediment zones are not simulated and therefore just operate on the bottom-most
         ! GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
         ! applied on flanks of the remaining layers, but note this is not suitable for
         ! model configurations where mass balance of benthic variables is required.

         ! Calculate temporal derivatives due to exchanges at the sediment/water interface
         !if ( self%zone_var .GE. 1 ) column(self%zone_var)%cell_sheet => z_sed_zones(1)
         call aed2_calculate_benthic(self%column, 1)

         ! Limit flux out of bottom layers to concentration of that layer
         ! i.e. don't flux out more than is there
         ! & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
         self%flux_pel(1, :) = max(-1.0 * self%cc(1, :), self%flux_pel(1, :)/self%grid%h(1))

         if ( self%aed2_cfg%benthic_mode .EQ. 1 ) then
            do lev=2,self%grid%nz_occupied
               ! Calculate temporal derivatives due to benthic fluxes.
               call aed2_calculate_benthic(self%column, lev)

               ! Limit flux out of bottom layers to concentration of that layer
               ! i.e. don't flux out more than is there
               ! & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
               ! scaled to proportion of area that is "bottom"
               self%flux_pel(lev, :) = max(-1.0 * self%cc(lev, :), self%flux_pel(lev, :)/self%grid%h(lev))
               self%flux_pel(lev, :) = self%flux_pel(lev, :) * (self%grid%Az(lev) - self%grid%Az(lev - 1))/self%grid%Az(lev)
            end do
         end if
      end if

      ! (2) SURFACE FLUXES
      ! Calculate temporal derivatives due to air-water exchange.
      if (.not. (state%total_ice_h > 0)) then ! no surface exchange under ice cover
         call aed2_calculate_surface(self%column, self%grid%nz_occupied)

         ! Distribute the fluxes into pelagic surface layer
         self%flux_pel(self%grid%nz_occupied, :) = self%flux_pel(self%grid%nz_occupied, :) + self%flux_atm(:)/self%grid%h(self%grid%nz_occupied)
      end if

      ! (3) WATER COLUMN KINETICS
      ! Add pelagic sink and soustatuse terms for all depth levels.
      do lev=1,self%grid%nz_occupied
         call aed2_calculate(self%column, lev)
      end do
   end subroutine calculate_fluxes

   ! Copy of diffusion algorithm used for Simstrat state variables
   subroutine diffusion_AED2_state(self, state, var_index)
      ! Arguments
      class(SimstratAED2) :: self
      class(ModelState) :: state
      integer :: var_index

      ! Local variables
      real(RK), dimension(self%grid%ubnd_vol) :: boundaries, sources, lower_diag, main_diag, upper_diag, rhs

      boundaries = 0.
      sources = 0.

      if (var_index == state%n_pH) state%AED2_state(:,state%n_pH) = 10.**(-state%AED2_state(:,state%n_pH))
      call euleri_create_LES_MFQ_AED2(self, state%AED2_state(:,var_index), state%nuh, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, state%dt)
      call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%AED2_state(:,var_index), self%grid%ubnd_vol)
      if (var_index == state%n_pH) state%AED2_state(:,state%n_pH) = -log10(state%AED2_state(:,state%n_pH))
   end subroutine

   ! Copy of disretization of Simstrat mean quantities
   subroutine euleri_create_LES_MFQ_AED2(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(SimstratAED2), intent(inout) :: self
      real(RK), dimension(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
      integer :: n

      n=self%grid%ubnd_vol

      ! Build diagonals
      upper_diag(1) = 0.0_RK
      upper_diag(2:n) = dt*nu(2:n)*self%grid%AreaFactor_1(2:n)
      lower_diag(1:n - 1) = dt*nu(2:n)*self%grid%AreaFactor_2(1:n-1)
      lower_diag(n) = 0.0_RK
      main_diag(1:n) = 1.0_RK - upper_diag(1:n) - lower_diag(1:n) + boundaries(1:n)*dt

      ! Calculate RHS
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      rhs(1:n) = var(1:n) + dt*sources(1:n)
   end subroutine

   ! Allocate memory for all variables (done by FABM)
   subroutine allocate_memory(self)
      implicit none

      ! Arguments
      class(SimstratAED2) :: self

      ! Local variables
      integer status
      allocate(self%column(self%n_AED2_state_vars),stat=status)
      allocate(self%column_sed(self%n_AED2_state_vars),stat=status)

      ! names = grab the names from info
      allocate(self%names(self%n_vars),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (names)'
      allocate(self%bennames(self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (bennames)'
      allocate(self%diagnames(self%n_vars_diag),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (diagnames)'
      allocate(self%diagnames_sheet(self%n_vars_diag_sheet),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (diagnames_sheet)'

      ! Now that we know how many vars we need, we can allocate space for them
      allocate(self%cc(self%grid%nz_grid, (self%n_vars + self%n_vars_ben)),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (CC)'
      self%cc = 0.         !# initialise to zeroFastatus

      ! Allocate memory for fluxes
      allocate(self%flux_atm(self%n_vars + self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (flux_atm)'

      allocate(self%flux_ben(self%n_vars + self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (flux_ben)'

      allocate(self%flux_pel(self%grid%nz_occupied, self%n_vars + self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (flux_pel)'

      !allocate(self%flux_zone(self%aed2_cfg%n_zones, self%n_vars + self%n_vars_ben),stat=status)
      !if (status /= 0) stop 'allocate_memory(): Error allocating (flux_zone)'

      ! Min, max values
      allocate(self%min_(self%n_vars + self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (min_)'

      allocate(self%max_(self%n_vars + self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (max_)'


      !# Allocate diagnostic variable array and set all values to zero.
      !# (needed because time-integrated/averaged variables will increment rather than set the array)
      allocate(self%cc_diag(self%grid%nz_grid, self%n_vars_diag),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (cc_diag)'
      self%cc_diag = zero_

      !# Allocate diagnostic variable array and set all values to zero.
      !# (needed because time-integrated/averaged variables will increment rather than set the array)
      allocate(self%cc_diag_sheet(self%n_vars_diag_sheet),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (cc_diag_sheet)'
      self%cc_diag_sheet = zero_

      !# Allocate array with vertical movement rates (m/s, positive for upwards),
      !# and set these to the values provided by the model.
      allocate(self%ws(self%grid%nz_grid, self%n_AED2_state_vars),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (ws)'
      self%ws = zero_

      !# Allocate array for photosynthetically active radiation (PAR).
      !# This will be calculated internally during each time step.
      allocate(self%par(self%grid%nz_grid),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (par)'
      self%par = zero_

      allocate(self%nir(self%grid%nz_grid),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (nir)'
      self%nir = zero_
      allocate(self%uva(self%grid%nz_grid),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (uva)'
      self%uva = zero_
      allocate(self%uvb(self%grid%nz_grid),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (uvb)'
      self%uvb = zero_

      allocate(self%sed_zones(self%grid%nz_grid))
      !# Allocate array for local pressure.
      !# This will be calculated [approximated] from layer depths internally
      !# during each time step.
      allocate(self%pres(self%grid%nz_grid),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (pres)'
      self%pres = zero_

      allocate(self%tss(self%grid%nz_grid),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (tss)'
      self%tss = zero_

      allocate(self%externalid(self%n_AED2_state_vars))
   end subroutine

   ! Assign names to all variables (done by FABM)
   subroutine assign_var_names(self)
      implicit none

      ! Arguments
      class(SimstratAED2) :: self

      ! Local variables
      type(aed2_variable_t),pointer :: tvar
      integer i, j

      print "(5X,'Configured variables to simulate:')"

      j = 0
      do i=1,self%n_AED2_state_vars
         if ( aed2_get_var(i, tvar) ) then
            if ( .not. (tvar%sheet .or. tvar%diag .or. tvar%extern) ) then
               j = j + 1
               self%names(j) = trim(tvar%name)
               self%min_(j) = tvar%minimum
               self%max_(j) = tvar%maximum
               print *,"     S(",j,") AED2 pelagic(3D) variable: ", trim(self%names(j))
            end if
         end if
      end do

      j = 0
      do i=1,self%n_AED2_state_vars
         if ( aed2_get_var(i, tvar) ) then
            if ( tvar%sheet .and. .not. (tvar%diag .or. tvar%extern) ) then
               j = j + 1
               self%bennames(j) = trim(tvar%name)
               self%min_(self%n_vars+j) = tvar%minimum
               self%max_(self%n_vars+j) = tvar%maximum
               print *,"     B(",j,") AED2 benthic(2D) variable: ", trim(self%bennames(j))
            end if
         end if
      end do

      j = 0
      do i=1,self%n_AED2_state_vars
         if ( aed2_get_var(i, tvar) ) then
            if ( tvar%diag ) then
               if ( .not.  tvar%sheet ) then
                  j = j + 1
                  self%diagnames(j) = trim(tvar%name)
                  print *,"     D(",j,") AED2 diagnostic 3Dvariable: ", trim(tvar%name)
               end if
            end if
         end if
      end do

      j = 0
      do i=1,self%n_AED2_state_vars
         if ( aed2_get_var(i, tvar) ) then
            if ( tvar%diag ) then
               if (tvar%sheet ) then
                  j = j + 1
                  self%diagnames_sheet(j) = trim(tvar%name)
                  print *,"     D(",j,") AED2 diagnostic 2Dvariable: ", trim(tvar%name)
               end if
            end if
         end if
      end do      

   end subroutine

   ! Pass variables to FABM
   subroutine define_column(self, state)
      !-------------------------------------------------------------------------------
      ! Set up the current column pointers
      !-------------------------------------------------------------------------------
      ! Arguments
      class(SimstratAED2) :: self
      class(ModelState) :: state

      ! Local variables
      integer :: av !, i
      integer :: v, d, sv, sd, ev
      type(aed2_variable_t), pointer :: tvar
      !-------------------------------------------------------------------------------
      ! Begin
      associate(column => self%column)

         v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
         do av=1,self%n_AED2_state_vars
            if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"

            if ( tvar%extern ) then !# global variable
               ev = ev + 1
               select case (tvar%name)
                  case ( 'temperature' ) ; column(av)%cell => state%T
                  case ( 'salinity' )    ; column(av)%cell => state%S
                  case ( 'density' )     ; column(av)%cell => state%rho
                  case ( 'layer_ht' )    ; column(av)%cell => self%grid%h(1:self%grid%nz_occupied)
                  case ( 'extc_coef' )   ; column(av)%cell => state%absorb_vol
                  case ( 'tss' )         ; column(av)%cell => self%tss
                  case ( 'par' )         ; column(av)%cell => self%par
                  case ( 'nir' )         ; column(av)%cell => self%nir
                  case ( 'uva' )         ; column(av)%cell => self%uva
                  case ( 'uvb' )         ; column(av)%cell => self%uvb
                  case ( 'pressure' )    ; column(av)%cell => self%pres
                  case ( 'depth' )       ; column(av)%cell => self%grid%layer_depth
                  case ( 'sed_zone' )    ; column(av)%cell_sheet => self%sed_zones(1)
                  case ( 'wind_speed' )  ; column(av)%cell_sheet => state%uv10
                  case ( 'rain')         ; column(av)%cell_sheet => state%rain
                  case ( 'par_sf' )      ; column(av)%cell_sheet => state%rad0
                  case ( 'taub' )        ; column(av)%cell_sheet => state%u_taub
                  case ( 'lake_depth' )  ; column(av)%cell_sheet => self%grid%lake_level
                  case ( 'layer_area' )  ; column(av)%cell => self%grid%Az_vol
                  case default ; call error("External variable "//TRIM(tvar%name)//" not found.")
               end select
            elseif ( tvar%diag ) then  !# Diagnostic variable
               if ( tvar%sheet ) then
                  sd = sd + 1
                  column(av)%cell_sheet => self%cc_diag_sheet(sd)
               else
                  d = d + 1
                  column(av)%cell => self%cc_diag(:,d)
               end if
            else    !# state variable
               if ( tvar%sheet ) then
                  sv = sv + 1
                  if ( tvar%bot ) then
                     column(av)%cell_sheet => self%cc(1, self%n_vars + sv)
      !            print *,'av',av,sv
                  elseif ( tvar%top ) then
                     column(av)%cell_sheet => self%cc(self%grid%nz_occupied, self%n_vars + sv)
                  end if

                  column(av)%flux_ben => self%flux_ben(self%n_vars + sv)
                  column(av)%flux_atm => self%flux_atm(self%n_vars + sv)
               else
                  v = v + 1
                  column(av)%cell => self%cc(:,v)
                  column(av)%flux_atm => self%flux_atm(v)
                  column(av)%flux_pel => self%flux_pel(:,v)
                  column(av)%flux_ben => self%flux_ben(v)
               end if
            end if
         end do
      end associate
   end subroutine define_column

   ! Check that all variable dependencies have been met (done by FABM)
   subroutine check_data(self)
      !-------------------------------------------------------------------------------
      ! Check that all variable dependencies have been met
      !-------------------------------------------------------------------------------
      ! Arguments
      class(SimstratAED2) :: self

      ! Local variables
      integer :: av
      integer :: v, d, sv, sd, ev, err_count
      type(aed2_variable_t),pointer :: tvar
      !-------------------------------------------------------------------------------
      ! Begin
      v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
      err_count = 0

      do av=1,self%n_AED2_state_vars
         if ( .not.  aed2_get_var(av, tvar) ) then
            call error("Error getting variable info")
            stop
         end if

         if ( tvar%extern ) then !# global variable
            ev = ev + 1
            select case (tvar%name)
               case ( 'temperature' ) ; tvar%found = .true.
               case ( 'salinity' )    ; tvar%found = .true.
               case ( 'density' )     ; tvar%found = .true.
               case ( 'layer_ht' )    ; tvar%found = .true.
               case ( 'extc_coef' )   ; tvar%found = .true.
               case ( 'tss' )         ; tvar%found = .true.
               case ( 'par' )         ; tvar%found = .true.
               case ( 'nir' )         ; tvar%found = .true.
               case ( 'uva' )         ; tvar%found = .true.
               case ( 'uvb' )         ; tvar%found = .true.
               case ( 'pressure' )    ; tvar%found = .true.
               case ( 'depth' )       ; tvar%found = .true.
               case ( 'sed_zone' )    ; tvar%found = .true.
               case ( 'wind_speed' )  ; tvar%found = .true.
               case ( 'rain' )        ; tvar%found = .true.
               case ( 'par_sf' )      ; tvar%found = .true.
               case ( 'taub' )        ; tvar%found = .true.
               case ( 'lake_depth' )  ; tvar%found = .true.
               case ( 'layer_area' )  ; tvar%found = .true.
               case default ; call error("ERROR: external variable "//trim(tvar%name)//" not found.")
            end select
         elseif ( tvar%diag ) then  !# Diagnostic variable
            if ( tvar%sheet ) then
               sd = sd + 1
            else
               d = d + 1
            end if
         else    !# state variable
            if ( tvar%sheet ) then
               sv = sv + 1
            else
               v = v + 1
            end if
         end if
         if ( .not. tvar%found ) then
            call error("Undefined variable " //trim(tvar%name))
            err_count = err_count + 1
         end if
      enddo

      if ( self%n_vars < v ) print *,"More vars than expected",v,self%n_vars
      if ( self%n_vars_ben < sv ) print *,"More sheet vars than expected"
      if ( self%n_vars_diag < d ) print *,"More diag vars than expected"
      if ( self%n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

      if ( err_count > 0 ) then
         call error("In AED2 configuration")
         stop
      end if
   end subroutine check_data

   ! Check whether states make sense (done by FABM)
   subroutine check_states(self)
      use,intrinsic :: ieee_arithmetic

      implicit none

      ! Arguments
      class(SimstratAED2) :: self

      ! Local variables
         type(aed2_variable_t), pointer :: tv
         integer :: i,v,lev
      !
      !-------------------------------------------------------------------------------
      ! Begin
         do lev=1, self%grid%nz_occupied
            call aed2_equilibrate(self%column, lev) ! Should probably moved to the update routine for clarity
            v = 0
            do i=1,self%n_AED2_state_vars
               if ( aed2_get_var(i, tv) ) then
                  if ( .not. (tv%diag .or. tv%extern) ) then
                     v = v + 1
                     if ( .not. ieee_is_nan(self%min_(v)) ) then
                     if ( self%cc(lev, v) < self%min_(v) ) self%cc(lev, v) = self%min_(v);
                     end if
                     if ( .not. ieee_is_nan(self%max_(v)) ) then
                     if ( self%cc(lev, v) > self%max_(v) ) self%cc(lev, v) = self%max_(v)
                     end if
                  end if
               end if
         end do
      end do
   end subroutine check_states

   ! Initialize input
   subroutine configure_simstrat_fabm_input()

      class (type_settings), pointer :: cfg

      cfg => settings_store%get_child('fabm/input', populator=fabm_input_populator)
   end subroutine configure_simstrat_fabm_input

   ! Add inputs
   subroutine append_input(input_variable)
      type (type_input_variable), target :: input_variable

      if (.not. associated(first_input_variable)) then
         first_input_variable => input_variable
      else
         last_input_variable%next => input_variable
      end if
      last_input_variable => input_variable
   end subroutine

   ! Allocate input variable, search in interior/horizontal/global variables
   subroutine fabm_input_create(self, pair)
      class (type_fabm_input_populator), intent(inout) :: self
      type (type_key_value_pair),        intent(inout) :: pair

      type (type_input_variable), pointer :: input_variable
      class (type_simstrat_settings), pointer :: branch
      character(len=attribute_length)     :: fabm_name
      integer :: i

      allocate(input_variable)

      ! First search in interior variables
      input_variable%interior_id = model%get_bulk_variable_id(pair%name)

      if (fabm_is_variable_used(input_variable%interior_id)) then
         fabm_name = fabm_get_variable_name(model, input_variable%interior_id)
         call type_input_create(pair, input_variable%profile_input, trim(input_variable%interior_id%variable%long_name), trim(input_variable%interior_id%variable%units), default=0._rk, pchild=branch)
         do i = 1, size(model%state_variables)
            if (fabm_name == model%state_variables(i)%name) then
               call branch%get(input_variable%relax_tau, 'relax_tau', 'relaxation time scale', 's', minimum=0._rk, default=1.e15_rk)
               call branch%get(input_variable%relax_tau_bot, 'relax_tau_bot', 'relaxation time scale for bottom layer', 's', minimum=0._rk, default=1.e15_rk)
               call branch%get(input_variable%relax_tau_surf, 'relax_tau_surf', 'relaxation time scale for surface layer', 's', minimum=0._rk, default=1.e15_rk)
               call branch%get(input_variable%h_bot, 'thickness_bot', 'thickness of bottom relaxation layer', 'm', minimum=0._rk, default=0._rk)
               call branch%get(input_variable%h_surf, 'thickness_surf', 'thickness of surface relaxation layer', 'm', minimum=0._rk, default=0._rk)
               exit
            end if
         end do
      else
         ! Variable was not found among interior variables. Try variables defined on horizontal slice of model domain (e.g., benthos)
         input_variable%horizontal_id = model%get_horizontal_variable_id(pair%name)
         if (fabm_is_variable_used(input_variable%horizontal_id)) then
            fabm_name = fabm_get_variable_name(model, input_variable%horizontal_id)
            call type_input_create(pair, input_variable%scalar_input, trim(input_variable%horizontal_id%variable%long_name), trim(input_variable%horizontal_id%variable%units), default=0._rk, pchild=branch)
            do i = 1, size(model%bottom_state_variables)
               if (fabm_name == model%bottom_state_variables(i)%name) then
                  call branch%get(input_variable%relax_tau, 'relax_tau', 'relaxation time scale', 's', minimum=0._rk, default=1.e15_rk)
                  exit
               end if
            end do
         else
            ! Variable was not found among interior or horizontal variables. Try global scalars.
            input_variable%scalar_id = model%get_scalar_variable_id(pair%name)
            if (.not. fabm_is_variable_used(input_variable%scalar_id)) then
               FATAL 'Variable '//pair%name//', referenced among FABM inputs was not found in model.'
               stop 'simstrat_fabm_input:init_simstrat_fabm_input'
            end if
            call type_input_create(pair, input_variable%scalar_input, trim(input_variable%scalar_id%variable%long_name), trim(input_variable%scalar_id%variable%units), default=0._rk, pchild=branch)
         end if
      end if
      call append_input(input_variable)
   end subroutine fabm_input_create

   ! Register input variables with the Simstrat-FABM driver
   subroutine init_simstrat_fabm_input(nlev,h)

      ! DESCRIPTION:
      ! Initialize files with observations on FABM variables.

      ! USES:
      use settings

      ! INPUT PARAMETERS:
      integer,          intent(in) :: nlev
      real(RK),         intent(in) :: h(1:nlev)

      ! LOCAL VARIABLES:
      type (type_input_variable), pointer :: curvariable
      integer                             :: k
      real(RK)                            :: db,ds,depth

      ! Calculate depth (used to determine whether in surface/bottom/bulk for relaxation times)
      depth = sum(h)

      curvariable => first_input_variable
      do while (associated(curvariable))
         if (fabm_is_variable_used(curvariable%interior_id)) then
            call register_input(curvariable%profile_input)

            allocate(curvariable%relax_tau_1d(0:nlev))
            curvariable%relax_tau_1d = curvariable%relax_tau

            ! Apply separate relaxation times for bottom and surface layer, if specified.
            db = _ZERO_
            ds = depth
            do k=1,nlev
               db = db+0.5*h(k)
               ds = ds-0.5*h(k)
               if (db<=curvariable%h_bot) curvariable%relax_tau_1d(k) = curvariable%relax_tau_bot
               if (ds<=curvariable%h_surf) curvariable%relax_tau_1d(k) = curvariable%relax_tau_surf
               db = db+0.5*h(k)
               ds = ds-0.5*h(k)
            end do

         ! Register observed variable with the simstrat-FABM driver.
            call register_observation(curvariable%interior_id, curvariable%profile_input%data, curvariable%relax_tau_1d)
         else
            call register_input(curvariable%scalar_input)
            if (fabm_is_variable_used(curvariable%horizontal_id)) then
               ! Horizontal variable
               call register_observation(curvariable%horizontal_id, curvariable%scalar_input%value, curvariable%relax_tau)
            else
               ! Scalar variable
               call register_observation(curvariable%scalar_id, curvariable%scalar_input%value)
            end if
         end if
         curvariable => curvariable%next
      end do
   end subroutine init_simstrat_fabm_input

   ! Process FABM initial conditions to Simstrat
   ! Should be done by the 4 above (all from GOTM)
   subroutine AED2_InitCondition(self, var, varname, default_val)
      !#################################### written/copied by A. Gaudard, 2015
         implicit none

         class(SimstratAED2) :: self
         real(RK), intent(inout) :: var(1:self%grid%nz_grid) ! Vector of initial conditions
         real(RK), intent(in) :: default_val ! Depth-independent value (default from aed2.nml)
         character(len=*), intent(in) :: varname ! Identifying the variable

         real(RK) :: z_read(self%grid%max_length_input_data), var_read(self%grid%max_length_input_data)
         real(RK) :: z_read_depth
         character(len=100) :: fname
         integer :: i,nval

         fname = trim(self%aed2_cfg%path_aed2_initial)//trim(varname)//'_ini.dat'
         open(14,action='read',status='unknown',err=1,file=fname)       ! Opens initial conditions file
         write(6,*) 'reading initial conditions of ', trim(varname)
         read(14,*)                                ! Skip header
         do i=1,self%grid%max_length_input_data                             ! Read initial values
               read(14,*,end=9) z_read(i),var_read(i)
         end do
      9   nval = i - 1                               ! Number of values
         if (nval<0) then
               write(6,*) 'Error reading ', trim(varname), ' initial conditions file (no data found).'
               stop
         end if
         close(14)
         do i=1,nval
               z_read(i) = abs(z_read(i))               ! Make depths positive
         end do
         z_read_depth = z_read(1)                     ! Initial depth (top-most)

         call reverse_in_place(z_read(1:nval))
         z_read(1:nval) = self%grid%z_zero - z_read(1:nval)
         call reverse_in_place(var_read(1:nval))

         if (nval==1) then
               write(6,*) '      Only one row! Water column will be initially homogeneous.'
               var_read(1:self%grid%nz_grid) = var_read(1)
         else
               call Interp(z_read(1:nval), var_read(1:nval), nval, self%grid%z_volume, var, self%grid%nz_grid)
         end if
         return

      1   write(6,*) '   File ''',trim(fname),''' not found. Initial conditions set to default value from file ''AED2.nml''.'
         var(1:self%grid%nz_grid) = default_val !File not found: value from fabm.nml (constant)
         return
   end subroutine AED2_InitCondition

   ! Light absorption feedback by AED2 variables
   subroutine absorption_updateAED2(self, state)

      ! Arguments
      class(SimstratAED2) :: self
      class(ModelState) :: state

      ! Local variables
      integer :: i
      real(RK) :: bio_extinction

      do i=self%grid%nz_occupied, 1, -1
         bio_extinction = 0.0_RK
         call aed2_light_extinction(self%column, i, bio_extinction)
         state%absorb_vol(i) = self%aed2_cfg%background_extinction + bio_extinction

      end do
      ! Interpolate to faces to be compatible with Simstrat temperature module
      call self%grid%interpolate_to_face(self%grid%z_volume, state%absorb_vol, self%grid%nz_occupied, state%absorb)

   end subroutine

   ! The mobility algorithm is taken from the GLM code (http://aed.see.uwa.edu.au/research/models/GLM/)
   !
   ! Assumptions:
   ! 1) movement direction has at most one change down the layers              *
   ! 2) sides of the lake slope inward (ie bottom is narrower than top)
   subroutine mobility(self, state, min_C, settling_v, conc)
      ! Arguments
      class(SimstratAED2) :: self
      class(ModelState) :: state
      real(RK), intent(in) :: min_C
      real(RK), dimension(:), intent(in) :: settling_v
      real(RK), dimension(:), intent(inout) :: conc

      ! Local variables
      real(RK) :: dtMax, tdt, tmp
      real(RK), dimension(self%grid%ubnd_vol) :: mins, vols, Y
      integer :: dirChng, signum, i, count


      ! determine mobility timestep i.e. maximum time step that particles
      ! will not pass through more than one layer in a time step 
      ! (probably not used in Simstrat as the timestep needs to be small)

      dtMax = state%dt
      dirChng = 0   ! this represents the layer at which direction switches from sinking to rising or visa versa
      signum = sign(1.0_RK,settling_v(1))  ! positive for rising, negative for sinking

      do i = 1,self%grid%ubnd_vol
         ! for convenience
         vols(i) = self%grid%Az_vol(i)*self%grid%h(i)
         mins(i) = min_C*vols(i)
         Y(i) = conc(i)*vols(i)

         !look for the change of direction
         if (signum .ne. sign(1.0_RK,settling_v(i))) then
            signum = -signum
            dirChng = i-1
         end if

         ! check if all movement can be from within this cell
         if (abs(settling_v(i)*state%dt) > self%grid%h(i)) then
            tdt = self%grid%h(i)/abs(settling_v(i))
            if (tdt < dtMax) dtMax = tdt
         end if

         ! check if movement can all be into the next cell.
         ! if movement is settling, next is below, otherwise next is above.
         ! check also for top or bottom in case of oopsies

         if (settling_v(i) > 0.) then
            if ((i < self%grid%ubnd_vol) .and. (abs(settling_v(i))*state%dt) > self%grid%h(i+1)) then
               tdt = self%grid%h(i+1)/abs(settling_v(i))
               if(tdt < dtMax) dtMax = tdt
            end if
         else if ((i > 1) .and. (abs(settling_v(i))*state%dt) > self%grid%h(i-1)) then
            tdt = self%grid%h(i-1)/abs(settling_v(i))
            if(tdt < dtMax) dtMax = tdt
         end if
      end do ! end find maximum time step dtMax
      if (dirChng == 0 .and. settling_v(1) > 0. ) dirChng = self%grid%ubnd_vol ! all rising
      if (dirChng == 0 .and. settling_v(self%grid%ubnd_vol) < 0.) dirChng = self%grid%ubnd_vol ! all sinking

      tdt = dtMax
      count = 0
      do
         ! do this in steps of dtMax, but at least once
         ! each time tdt is dtMax, except, possibly, the last which is whatever was left.
         count = count + 1    ! counter
         if (count*dtMax > state%dt) tdt = state%dt - (count - 1)*dtMax ! Last timestep

         ! 2 possibilities
         ! 1) lower levels rising, upper levels sinking
         ! 2) lower levels sinking, upper levels rising
         if (settling_v(1) > 0. ) then ! lower levels rising
            if (settling_v(self%grid%ubnd_vol) < 0.) then !top levels are sinking
               call Sinking(self, Y, conc, settling_v, vols, mins, tdt, self%grid%ubnd_vol, dirChng, tmp)
               Y(dirChng) = Y(dirChng) + tmp
               conc(dirChng) = Y(dirChng)/vols(dirChng)
            end if
            call Rising(self, Y, conc, settling_v, vols, mins, tdt, 1, dirChng)
         else ! lower levels sinking
            call Sinking(self, Y, conc, settling_v, vols, mins, tdt, dirChng, 1, tmp)
            if ( settling_v(self%grid%ubnd_vol) > 0.) then !top levels are rising
               call Rising(self, Y, conc, settling_v, vols, mins, tdt, dirChng, self%grid%ubnd_vol)
            end if
            if (state%dt > 0.) exit
         end if
      end do

   end subroutine

   ! Rising is the easier on the two since the slope means we dont need to look
   ! at relative areas (the cell above will always be >= to the current cell)
   ! all matter is moved to the next cell      
   ! for each cell except the top :  
   ! 1) calculate how much is going to move  
   ! 2) subtract amount that must now move 
   ! 3) add the amount moved from previous cell to current cell
   ! 4) fix concentration 
   ! for the top cell :
   ! 1) add the amount moved from previous cell to current cell
   ! 2) fix concentration
   subroutine Rising(self, Y, conc, settling_v, vols, mins, dt, start_i, end_i)
      ! Arguments
      class(SimstratAED2) :: self
      real(RK), dimension(:), intent(inout) :: conc, Y
      real(RK), dimension(:), intent(in) :: settling_v, vols, mins
      real(RK) :: dt
      integer :: start_i, end_i

      ! Local variables
      real(RK) :: mov, moved
      integer i

      mov = 0.
      moved = 0.

      do i = start_i,end_i
         ! speed times time (=h) time area * concen = mass to move
         mov = (settling_v(i) * dt) * self%grid%Az_vol(i)*conc(i)
         ! if removing that much would bring it below min conc
         if ((Y(i) + moved - mov) < mins(i) ) mov = Y(i) + moved - mins(i)

         Y(i) = Y(i) + moved - mov;
         conc(i) = Y(i) / vols(i) ! return it to a concentration
         moved = mov ! for the next step
      end do
      ! nothing rises out of the end cell, but we still add that which came from below
      Y(end_i) = Y(end_i) + moved
      conc(end_i) = Y(end_i) / vols(end_i)

   end subroutine

   ! for each cell :
   ! 1) calculate how much is going to move
   ! 2) subtract amount that must now move
   ! 3) add the amount moved from previous cell to current cell
   ! 4) fix concentration
   ! 5) compute the amount that will go to the next cell
   subroutine Sinking(self, Y, conc, settling_v, vols, mins, dt, start_i, end_i, moved)
      ! Arguments
      class(SimstratAED2) :: self
      real(RK), dimension(:), intent(inout) :: conc, Y
      real(RK), dimension(:), intent(in) :: settling_v, vols, mins
      real(RK) :: dt
      real(RK), intent(out) :: moved
      integer :: start_i, end_i

      ! Local variables
      real(RK) :: mov
      integer :: i

      mov = 0.
      moved = 0.

      do i = start_i,end_i, -1
         ! speed times time (=h) time area * concen = mass to move
         mov = (abs(settling_v(i)) * dt) * self%grid%Az_vol(i) * conc(i)
         ! if removing that much would bring it below min conc
         if ((Y(i) + moved - mov) < mins(i)) mov = Y(i) + moved - mins(i)
         Y(i) = Y(i) + moved - mov
         conc(i) = Y(i) / vols(i) ! return it to a concentration

         ! so now mov has how much has moved out of the cell, but not all
         ! of that will go into the next cell (FB: part of it sediments on the flancs)
         if ( i > 1 )  then
            moved = mov * (self%grid%Az_vol(i-1) / self%grid%Az_vol(i)) ! for the next step
         else
            moved = mov ! we are about to exit anyway.
         end if
      end do
   end subroutine
end module simstrat_fabm