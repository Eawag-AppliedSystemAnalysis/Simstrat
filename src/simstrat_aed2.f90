!     +---------------------------------------------------------------+
!     |  Simstrat - AED2 interface
!     +---------------------------------------------------------------+

module simstrat_aed2
   use strat_simdata
   use strat_grid
   use strat_solver
   use utilities
   use aed2_common
   use aed2_core

   implicit none
   private

   type, public :: SimstratAED2
      class(AED2Config), pointer :: aed2_cfg
      class(StaggeredGrid), pointer :: grid

      !# Arrays for state and diagnostic variables
      real(RK),pointer,dimension(:,:) :: cc !# water quality array: nlayers, nvars
      real(RK),pointer,dimension(:,:) :: cc_diag
      real(RK),pointer,dimension(:) :: cc_diag_hz
      real(RK),pointer,dimension(:) :: tss
      real(RK),pointer,dimension(:) :: sed_zones

      ! Arrays for fluxes of state variables
      real(RK),pointer,dimension(:) :: flux_atm
      real(RK),pointer,dimension(:) :: flux_ben
      real(RK),pointer,dimension(:,:) :: flux_pel
      real(RK),pointer,dimension(:,:) :: flux_zone

      !# Arrays for work, vertical movement, and cross-boundary fluxes
      real(RK),allocatable,dimension(:,:) :: ws
      real(RK),allocatable,dimension(:)   :: total
      real(RK),allocatable,dimension(:)   :: local

      !# Arrays for environmental variables not supplied externally.
      real(RK),pointer,dimension(:) :: par, pres
      real(RK),pointer,dimension(:) :: uva, uvb, nir

      !# External variables
      integer  :: w_adv_ctr    ! Scheme for vertical advection (0 if not used)

      character(len=48),allocatable :: names(:)
      character(len=48),allocatable :: bennames(:)

      integer,allocatable,dimension(:) :: externalid

      real(RK),allocatable,dimension(:) :: min_, max_

      integer :: n_aed2_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
      integer :: zone_var = 0

      ! Variables for in-/outflow
      real(RK), dimension(:, :), allocatable   :: z_Inp, Q_start, Qs_start, Q_end, Qs_end, Q_read_start, Q_read_end
      real(RK), dimension(:, :), allocatable   :: Inp_read_start, Inp_read_end, Qs_read_start, Qs_read_end, Q_inp
      real(RK), dimension(:), allocatable :: tb_start, tb_end ! Input depths, start time, end time
      integer, dimension(:), allocatable :: eof, nval, nval_deep, nval_surface

   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
      procedure, pass(self), public :: absorption_updateAED2
   end type SimstratAED2

contains

   subroutine init(self, state, grid, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(ModelState) :: state
      class(AED2Config), target :: aed2_cfg
      class(StaggeredGrid), target :: grid

      ! Local variables
      character(len=80) :: fname
      character(len=64) :: models(64)
      namelist /aed2_models/ models
      type(aed2_variable_t),pointer :: tvar
      integer i, status, av, v, sv

      self%grid => grid
      self%aed2_cfg => aed2_cfg

      associate (n_aed2_vars => self%n_aed2_vars, &
                 n_vars => self%n_vars, &
                 n_vars_ben => self%n_vars_ben, &
                 n_vars_diag => self%n_vars_diag, &
                 n_vars_diag_sheet => self%n_vars_diag_sheet)

         ! AED2 config file
         fname = 'aed2.nml'

         if ( aed2_init_core('.') /= 0 ) call error("Initialisation of aed2_core failed")
         call aed2_print_version

         ! Create model tree
         write (6,*) "     Processing aed2_models config from ", trim(fname)
         open(50,file=fname,action='read',status='old',iostat=status)
         if ( status /= 0 ) then
            call error("Cannot open file " // trim(fname))
            stop
         end if

         models = ''
         read(50, nml=aed2_models, iostat=status)
         if ( status /= 0 ) then
            call error("Cannot read namelist entry aed2_models")
            stop
         end if

         do i=1,size(models)
            if (models(i)=='') exit
            call aed2_define_model(models(i), 50)
         end do

         !# should be finished with this file
         close(50)
         write (6,*) "      AED2 file parsing completed."

         ! Assign number of different variables
         n_aed2_vars = aed2_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

         ! Print variable information to screen
         print "(/,5X,'AED2 : n_aed2_vars = ',I3,' ; MaxLayers         = ',I4)",n_aed2_vars,self%grid%nz_grid
         print "(  5X,'AED2 : n_vars      = ',I3,' ; n_vars_ben        = ',I3)",n_vars,n_vars_ben
         print "(  5X,'AED2 : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I3,/)",n_vars_diag,n_vars_diag_sheet

         ! Check variable dependencies
         call check_data(self)

         ! Allocate space for the allocatables/pointers of this module
         call allocate_memory(self)

         ! Allocate memory for AED2 state and inflow matrix used by Simstrat
         !allocate(state%AED2_state(self%grid%nz_grid, n_vars + n_vars_ben))
         state%AED2_state => self%cc
         allocate(state%AED2_inflow(self%grid%nz_grid, n_vars + n_vars_ben))

         ! Assign name, min and max values of variables, print names to screen
         call assign_var_names(self)

         !# Now set initial values
         v = 0 ; sv = 0;
         do av=1,self%n_aed2_vars
            if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"
            if ( .not. ( tvar%extern .or. tvar%diag) ) then  !# neither global nor diagnostic variable
               if ( tvar%sheet ) then
                  sv = sv + 1
                  call AED2_InitCondition(self, self%cc(:, n_vars + sv), tvar%name, tvar%initial)
               else
                  v = v + 1
                  call AED2_InitCondition(self, self%cc(:, v), tvar%name, tvar%initial)
               end if
            end if
         end do

         write(*,"(/,5X,'----------  AED2 config : end  ----------',/)")
      end associate
   end subroutine


   subroutine update(self, state)
      implicit none
      class(SimstratAED2) :: self
      class(ModelState) :: state

      ! Local variables
      type (aed2_column_t) :: column(self%n_aed2_vars), column_sed(self%n_aed2_vars)
      integer :: v, split, lev

      !# Calculate local pressure
      self%pres(1:self%grid%ubnd_vol) = -self%grid%z_volume(1:self%grid%ubnd_vol)

      call define_column(self, state, column)
      !if (benthic_mode .GT. 1) call define_sed_column(column_sed, n_zones, flux, flux_atm, flux_ben)

      ! If sediment layers are simulated
      !if (benthic_mode .GT. 1) call define_sed_column(column_sed, n_zones, flux, flux_atm, flux_ben)

      self%cc_diag = 0.
      self%cc_diag_hz = 0.

      call check_states(self, column)

      do split=1, self%aed2_cfg%split_factor

         call absorption_updateAED2(self, column, state)

         !# Fudge
         self%nir(:) = (self%par(:)/0.45) * 0.51
         self%uva(:) = (self%par(:)/0.45) * 0.035
         self%uvb(:) = (self%par(:)/0.45) * 0.005

         call calculate_fluxes(self, state, column, column_sed)

         ! Update the water column layers using the biochemical reaction of AED2
         do v = 1, self%n_vars
            do lev = 1, self%grid%nz_occupied
               self%cc(lev, v) = self%cc(lev, v) + state%dt/self%aed2_cfg%split_factor*self%flux_pel(lev, v)
            end do
         end do

!       ! Now update benthic variables, depending on whether zones are simulated
!       IF ( benthic_mode .GT. 1 ) THEN
!          ! Loop through benthic state variables to update their mass
!          DO v = n_vars+1, n_vars+n_vars_ben
!             ! Loop through each sediment zone
!             DO lev = 1, n_zones
!                ! Update the main cc_sed data array with the
!                z_cc(lev, v) = z_cc(lev, v)+ dt_eff*flux_zone(lev, v)
!             ENDDO
!          ENDDO
!       ELSE
!          DO v = n_vars+1, n_vars+n_vars_ben
!             cc(1, v) = cc(1, v) + dt_eff*flux_ben(v)
!          ENDDO
!       ENDIF

!       ! If simulating sediment zones, distribute cc-sed benthic properties back
!       !  into main cc array, mainly for plotting
!       IF ( benthic_mode .GT. 1 ) CALL copy_from_zone(cc, cc_diag, cc_diag_hz, wlev)

         ! Update in-/outflow of AED2 variables (to be used in following advection step in the main loop)
         call lateral_update_AED2(self, state)


         ! Diffusive transport of AED2 variables
         do v=1, self%n_vars
            call diffusion_AED2(self, state, v)
         end do

      end do


   end subroutine

   subroutine check_data(self)
   !-------------------------------------------------------------------------------
   ! Check that all variable dependencies have been met
   !-------------------------------------------------------------------------------
   !ARGUMENTS
   class(SimstratAED2) :: self
   !LOCALS
      integer :: av !, i
      integer :: v, d, sv, sd, ev, err_count
      type(aed2_variable_t),pointer :: tvar
   !-------------------------------------------------------------------------------
   !BEGIN
      v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
      err_count = 0

      do av=1,self%n_aed2_vars
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

   subroutine absorption_updateAED2(self, column, state)
      class(SimstratAED2) :: self
      type (aed2_column_t), intent(inout) :: column(:)
      class(ModelState) :: state
      !class(StaggeredGrid) :: grid
      !class(AED2Config) :: aed2_cfg

      integer :: i
      real(RK) :: bio_extinction
      real(RK), dimension(self%grid%nz_occupied) :: extc_coef
      !real(RK), dimension(size(state%absorb) - 1) :: absorb_on_vol_grid

      do i=self%grid%nz_occupied - 1, 1, -1
         bio_extinction = 0.0_RK
         call aed2_light_extinction(column, i, bio_extinction)
         extc_coef(i) = self%aed2_cfg%background_extinction + bio_extinction

      end do

      ! Interpolate to faces to be compatible with Simstrat temperature module
      call self%grid%interpolate_to_face(self%grid%z_volume, extc_coef, self%grid%nz_occupied, state%absorb)

   end subroutine

   SUBROUTINE define_column(self, state, column)
   !-------------------------------------------------------------------------------
   ! Set up the current column pointers
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      class(SimstratAED2) :: self
      class(ModelState) :: state
      type (aed2_column_t), intent(inout) :: column(:)
   !
   !LOCALS
      integer :: av !, i
      integer :: v, d, sv, sd, ev
      type(aed2_variable_t), pointer :: tvar
   !-------------------------------------------------------------------------------
   !BEGIN
      v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
      do av=1,self%n_aed2_vars
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
               case ( 'depth' )       ; column(av)%cell => self%grid%z_volume
               case ( 'sed_zone' )    ; column(av)%cell_sheet => self%sed_zones(1)
               case ( 'wind_speed' )  ; column(av)%cell_sheet => state%uv10
               case ( 'par_sf' )      ; column(av)%cell_sheet => state%rad0
               case ( 'taub' )        ; column(av)%cell_sheet => state%u_taub
               case ( 'lake_depth' )  ; column(av)%cell_sheet => self%grid%lake_level
               case ( 'layer_area' )  ; column(av)%cell => self%grid%Az_vol
               case default ; call error("External variable "//TRIM(tvar%name)//" not found.")
            end select
         elseif ( tvar%diag ) then  !# Diagnostic variable
            if ( tvar%sheet ) then
               sd = sd + 1
               column(av)%cell_sheet => self%cc_diag_hz(sd)
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
               endif

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
   end SUBROUTINE define_column

   SUBROUTINE check_states(self, column)
      use,intrinsic :: ieee_arithmetic

      class(SimstratAED2) :: self

   !ARGUMENTS
      type (aed2_column_t),intent(inout) :: column(:)
   !
   !LOCALS
      TYPE(aed2_variable_t),POINTER :: tv
      INTEGER i,v,lev
   !
   !-------------------------------------------------------------------------------
   !BEGIN
      do lev=1, self%grid%nz_occupied
         ! equilibrate causes seg fault in stability module... FB
         !call aed2_equilibrate(column, lev)    !MH this should be in the main do_glm routine ????!!!
         v = 0
         do i=1,self%n_aed2_vars
            if ( aed2_get_var(i, tv) ) then
               if ( .not. (tv%diag .or. tv%extern) ) then
                  v = v + 1
                  if ( self%aed2_cfg%repair_state ) then
                     if ( .not. ieee_is_nan(self%min_(v)) ) then
                        if ( self%cc(lev, v) < self%min_(v) ) self%cc(lev, v) = self%min_(v);
                     end if
                     if ( .not. ieee_is_nan(self%max_(v)) ) then
                        if ( self%cc(lev, v) > self%max_(v) ) self%cc(lev, v) = self%max_(v)
                     end if
                  end if
               end if
            end if
         end do
      end do
   END SUBROUTINE check_states

   SUBROUTINE calculate_fluxes(self, state, column, column_sed)
   !-------------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !-------------------------------------------------------------------------------
   use,intrinsic :: ieee_arithmetic

   !ARGUMENTS
      class(SimstratAED2) :: self
      class(ModelState) :: state
      type (aed2_column_t), intent(inout) :: column(:)
      type (aed2_column_t), intent(inout) :: column_sed(:)
   !
   !LOCALS
      integer :: lev,zon,v_start,v_end,av,sv,sd
      real(RK) :: scale
      real(RK), dimension(self%grid%nz_occupied, self%n_vars)    :: flux_pel_pre
      real(RK), dimension(self%aed2_cfg%n_zones, self%n_vars) :: flux_pel_z
      logical :: splitZone
      type(aed2_variable_t),pointer :: tvar
   !-------------------------------------------------------------------------------
   !BEGIN
      self%flux_pel = zero_
      self%flux_atm = zero_
      self%flux_ben = zero_

      !# Start with calculating all flux terms for rhs in mass/m3/s
      !# Includes (1) benthic flux, (2) surface exchange and (3) water column kinetics
      !# as calculated by glm


      !# (1) BENTHIC FLUXES
      if ( self%aed2_cfg%benthic_mode .gt. 1 ) then
!          !# Multiple static sediment zones are simulated, and therfore overlying
!          !# water conditions need to be aggregated from multiple cells/layers, and output flux
!          !# needs disaggregating from each zone back to the overlying cells/layers

!          do zon=1,self%aed2_cfg%n_zones
!             !# Reinitialise flux_ben to be repopulated for this zone
!             flux_ben = zero_
!             flux_pel_pre = zero_

!             !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
!             if ( self%zone_var .ge. 1 ) then
!                column_sed(zone_var)%cell_sheet => z_sed_zones(zon)
!        !       !MH WE NEED A COLUMN TO CC VAR MAP FOR BENTHIC GUYS
!                !CAB Yes, a map (or 2 maps) would be better, but QnD since this all needs reworking
!                sv = 0 ; sd = 0
!                do av=1,self%n_aed2_vars
!                   if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"
!                   if ( .not. tvar%extern .and. tvar%sheet ) then
!                      if ( tvar%diag ) then
!                         sd = sd + 1
!                         column(av)%cell_sheet => z_diag_hz(zon, sd)
!                      else
!                         sv = sv + 1
!                         column(av)%cell_sheet => z_cc(zon, self%n_vars + sv)
!                      end if
!                   end if
!                end do
!                !print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
!             end if
!             if ( self%aed2_cfg%benthic_mode .eq. 3 ) then
!                !# Zone is able to operated on by riparian and dry methods
!                call aed2_calculate_riparian(column_sed, zon, z_pc_wet(zon))
!                if (z_pc_wet(zon) .eq. 0. ) call aed2_calculate_dry(column_sed, zon)
!             end if
!             !# Calculate temporal derivatives due to benthic processes.
!             !# They are stored in flux_ben (benthic vars) and flux_pel (water vars)
!             flux_pel_pre = flux_pel

!    !        print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
!             call aed2_calculate_benthic(column_sed, zon)

!             !# Record benthic fluxes in the zone array
!             flux_zon(zon, :) = flux_ben(:)

!             !# Now we have to find out the water column flux that occured and
!             !# disaggregate it to relevant layers
!             flux_pel_z(zon,:) = flux_pel(zon,:)-flux_pel_pre(zon,:)
!          end do

!          !# Disaggregation of zone induced fluxes to overlying layers
!          v_start = 1 ; v_end = self%n_vars
!          zon = self%aed2_cfg%n_zones
!          do lev=self%grid%nz_occupied,1,-1
!            if ( zon .ne. 1 ) then
!              splitZone = zz(lev-1) < zone_heights(zon-1)
!            else
!              splitZone = .FALSE.
!            end if

!            if (splitZone) then
!              scale = (zone_heights(zon-1) - zz(lev-1)) / (zz(lev) - zz(lev-1))
!              flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end) * scale

!              zon = zon - 1

!              flux_pel(lev,v_start:v_end) = flux_pel(lev,v_start:v_end) + &
!                                            flux_pel_z(zon,v_start:v_end) * (1.0 - scale)
!            else
!              flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end)
!            end if
!          end do
!          !# Limit flux out of bottom waters to concentration of that layer
!          !# i.e. don't flux out more than is there & distribute
!          !# bottom flux into pelagic over bottom box (i.e., divide by layer height).
!          !# scaled to proportion of area that is "bottom"
!          do lev=1,self%grid%nz_occupied
!             if(lev>1)flux_pel(lev, :) = flux_pel(lev, :) * (self%grid%Az_vol(lev) - self%grid%Az_vol(lev - 1))/self%grid%Az_vol(lev)
!             flux_pel(lev, :) = max(-1.0 * self%cc(lev, :), flux_pel(lev, :)/self%grid%h(lev))
!          end do
      else
         !# Sediment zones are not simulated and therefore just operate on the bottom-most
         !# GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
         !# applied on flanks of the remaining layers, but note this is not suitable for
         !# model configurations where mass balance of benthic variables is required.

         !# Calculate temporal derivatives due to exchanges at the sediment/water interface
         !if ( self%zone_var .GE. 1 ) column(self%zone_var)%cell_sheet => z_sed_zones(1)
         call aed2_calculate_benthic(column, 1)

         !# Limit flux out of bottom layers to concentration of that layer
         !# i.e. don't flux out more than is there
         !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
         self%flux_pel(1, :) = max(-1.0 * self%cc(1, :), self%flux_pel(1, :)/self%grid%h(1))

         if ( self%aed2_cfg%benthic_mode .EQ. 1 ) then
            do lev=2,self%grid%nz_occupied
               !# Calculate temporal derivatives due to benthic fluxes.
               call aed2_calculate_benthic(column, lev)

               !# Limit flux out of bottom layers to concentration of that layer
               !# i.e. don't flux out more than is there
               !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
               !# scaled to proportion of area that is "bottom"
               self%flux_pel(lev, :) = max(-1.0 * self%cc(lev, :), self%flux_pel(lev, :)/self%grid%h(lev))
               self%flux_pel(lev, :) = self%flux_pel(lev, :) * (self%grid%Az_vol(lev) - self%grid%Az_vol(lev - 1))/self%grid%Az_vol(lev)
            end do
         end if
      end if

      !# (2) SURFACE FLUXES
      !# Calculate temporal derivatives due to air-water exchange.
      if (.not. (state%ice_h > 0)) then !# no surface exchange under ice cover
         call aed2_calculate_surface(column, self%grid%nz_occupied)

         !# Distribute the fluxes into pelagic surface layer
         self%flux_pel(self%grid%nz_occupied, :) = self%flux_pel(self%grid%nz_occupied, :) + self%flux_atm(:)/self%grid%h(self%grid%nz_occupied)
      end if

      !# (3) WATER COLUMN KINETICS
      !# Add pelagic sink and soustatuse terms for all depth levels.
      do lev=1,self%grid%nz_occupied
         call aed2_calculate(column, lev)
      end do
   END SUBROUTINE calculate_fluxes

   subroutine AED2_InitCondition(self, var,varname,default_val)
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

        fname = trim(self%aed2_cfg%path_aed2_ini)//trim(varname)//'_ini.dat'
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

    1   write(6,*) '   File ''',trim(fname),''' not found. Initial conditions set to default value from file ''fabm.nml''.'
        var(1:self%grid%nz_grid) = default_val !File not found: value from fabm.nml (constant)
        return
    end subroutine AED2_InitCondition

    subroutine allocate_memory(self)
      class(SimstratAED2) :: self

      ! Local variables
      integer status

      !# names = grab the names from info
      allocate(self%names(self%n_vars),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (names)'
      allocate(self%bennames(self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (bennames)'

      !# Now that we know how many vars we need, we can allocate space for them
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

      allocate(self%flux_zone(self%aed2_cfg%n_zones, self%n_vars + self%n_vars_ben),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (flux_zone)'

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
      allocate(self%cc_diag_hz(self%n_vars_diag_sheet),stat=status)
      if (status /= 0) stop 'allocate_memory(): Error allocating (cc_diag_hz)'
      self%cc_diag_hz = zero_

      !# Allocate array with vertical movement rates (m/s, positive for upwards),
      !# and set these to the values provided by the model.
      allocate(self%ws(self%grid%nz_grid, self%n_vars),stat=status)
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

      !allocate(self%dz(self%grid%nz_grid),stat=status)
      !self%dz = zero_

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

      allocate(self%externalid(self%n_aed2_vars))

   end subroutine

   subroutine assign_var_names(self)
      class(SimstratAED2) :: self

      ! Local variables
      type(aed2_variable_t),pointer :: tvar
      integer i, j

      print "(5X,'Configured variables to simulate:')"

      j = 0
      do i=1,self%n_aed2_vars
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
      do i=1,self%n_aed2_vars
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
      do i=1,self%n_aed2_vars
         if ( aed2_get_var(i, tvar) ) then
            if ( tvar%diag ) then
               if ( .not.  tvar%sheet ) then
                  j = j + 1
                  print *,"     D(",j,") AED2 diagnostic 3Dvariable: ", trim(tvar%name)
               end if
            end if
         end if
      end do

      j = 0
      do i=1,self%n_aed2_vars
         if ( aed2_get_var(i, tvar) ) then
            if ( tvar%diag ) then
               if (tvar%sheet ) then
                  j = j + 1
                  print *,"     D(",j,") AED2 diagnostic 2Dvariable: ", trim(tvar%name)
               end if
            end if
         end if
      end do      

   end subroutine

   subroutine diffusion_AED2(self, state, var_index)
      ! Arguments
      class(SimstratAED2) :: self
      class(ModelState) :: state

      integer :: var_index

      ! Locals
      real(RK), dimension(self%grid%ubnd_vol) :: boundaries, sources, lower_diag, main_diag, upper_diag, rhs

      boundaries = 0.
      sources = 0.

      call euleri_create_LES_MFQ_AED2(self, state%AED2_state(:,var_index), state%num, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, state%dt)
      !call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%AED2_state(:,var_index), self%grid%ubnd_vol)


   end subroutine

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

   subroutine lateral_update_AED2(self, state)
      implicit none
      class(SimstratAED2) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: dummy

      integer :: i, j, n
      integer :: fnum(1:32) ! File number
      character(len=100) :: fname(1:32)
      type(aed2_variable_t),pointer :: tvar

      associate (datum => state%datum, &
                 idx => state%model_step_counter, &
                 Q_inp => state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 grid => self%grid, &
                 ubnd_vol => self%grid%ubnd_vol, &
                 ubnd_fce => self%grid%ubnd_fce)

         n = self%n_vars + self%n_vars_ben
         ! FB 2016: Major revision to include surface inflow
         do i = 1, self%n_vars ! Do this for all AED2 vars
            fname(i) = self%names(i)
            fnum(i) = i + 60  ! Should find a better way to manage unit numbers
            if (idx==1) then ! First iteration
               if (i==1) then ! First variable
                  ! Allocate arrays for first iteration of first variable
                  allocate (self%z_Inp(1:n, 1:state%nz_input)) ! Input depths
                  allocate (self%Inp_read_start(1:n, 1:state%nz_input)) ! Raw input read
                  allocate (self%Inp_read_end(1:n, 1:state%nz_input)) ! Raw input read
                  allocate (self%Q_read_start(1:n, 1:state%nz_input)) ! Integrated input
                  allocate (self%Q_read_end(1:n, 1:state%nz_input)) ! Integrated input           
                  allocate (self%Qs_read_start(1:n, 1:state%nz_input))  ! Integrated surface input
                  allocate (self%Qs_read_end(1:n, 1:state%nz_input))  ! Integrated surface input
                  allocate (self%Q_start(1:n, 1:grid%nz_grid+1)) ! Input interpolated on grid
                  allocate (self%Q_end(1:n, 1:grid%nz_grid+1)) ! Input interpolated on grid
                  allocate (self%Qs_start(1:n, 1:grid%nz_grid+1)) ! Surface input interpolated on grid
                  allocate (self%Qs_end(1:n, 1:grid%nz_grid+1)) ! Surface input interpolated on grid
                  allocate (self%Q_inp(1:n, 1:grid%nz_grid+1))
                  allocate (self%tb_start(n), self%tb_end(n), self%eof(n), self%nval(n), self%nval_deep(n), self%nval_surface(n))
               end if

               ! Default values
               self%Q_start(i,:) = 0.0_RK
               self%Q_end(i,:) = 0.0_RK
               self%Qs_start(i, :) = 0.0_RK
               self%Qs_end(i, :) = 0.0_RK

               ! Open file and start to read
               self%eof(i) = 0
               read (fnum(i), *, end=9) ! Skip first row: description of columns

               if (state%has_surface_input(i)) then
                ! Read number of deep and surface inflows
                read (fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)
                ! Total number of values to read
                self%nval(i) = self%nval_deep(i) + self%nval_surface(i)
              else
                read (fnum(i), *, end=9) self%nval_deep(i)
                ! Total number of values to read
                self%nval(i) = self%nval_deep(i)
              end if

               ! Read input depths
               read (fnum(i), *, end=9) dummy, (self%z_Inp(i, j), j=1, self%nval(i))

               ! Convert input depths
               self%z_Inp(i, 1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i, 1:self%nval_deep(i))

               if (state%has_surface_input(i)) then
                ! Convert surface input depths
                self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))
              end if

               ! Read first input line
               read (fnum(i), *, end=9) self%tb_start(i), (self%Inp_read_start(i, j), j=1, self%nval(i))

              if (state%has_deep_input(i)) then
                ! Cumulative integration of input
                call Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), self%Q_read_start(i, :), self%nval_deep(i))
                ! Interpolation on face grid
                call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))
              end if

               ! If there is surface input, integrate and interpolate
               if (state%has_surface_input(i)) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               end if


               ! Read second line and treatment of deep inflow
               read (fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
              if (state%has_deep_input(i)) then
                call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
              end if
               ! If there is surface input, integrate and interpolate
               if (state%has_surface_input(i)) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if

               write(6,*) '[OK] ','Input file successfully read: ',fname(i)
            end if ! idx==1



            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep them at the surface
            if ((.not. grid%lake_level == grid%lake_level_old) .and. (state%has_surface_input(i))) then

              ! Readjust surface input depths
              self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

              ! Adjust surface inflow to new lake level        
              call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
              call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))               

            end if ! end if not lake_level...



            ! Temporal treatment of inflow
            if ((datum <= self%tb_start(i)) .or. (self%eof(i) == 1)) then ! if datum before first date or end of file reached
               goto 8
            else
               do while (.not. ((datum >= self%tb_start(i)) .and. (datum <= self%tb_end(i)))) ! Do until datum between dates
                  self%tb_start(i) = self%tb_end(i) ! Move one step in time
                  self%Q_start(i, :) = self%Q_end(i, :)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)

                  read (fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))

                  if (state%has_deep_input(i)) then
                    call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                    call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
                  end if

                  if (state%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
               end do ! end do while
            end if

            ! Linearly interpolate value at correct datum (Q_inp is on face grid)
            do j = 1, ubnd_fce
               self%Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i)) &
               * (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))/(self%tb_end(i)-self%tb_start(i))
            end do
            goto 11

            ! If end of file reached, set to closest available value
 7          self%eof(i) = 1
 8          self%Q_inp(i,1:ubnd_fce) = self%Q_start(i,1:ubnd_fce) + self%Qs_start(i,1:ubnd_fce)
            goto 11

            ! If no data available
 9          write(6,*) '[WARNING] ','No data found in ',trim(fname(i)),' inflow file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            self%Q_inp(i, 1:ubnd_fce) = 0.0_RK
            self%Q_start(i, 1:ubnd_fce) = 0.0_RK
            self%Qs_start(i, 1:ubnd_fce) = 0.0_RK
            11        continue

         end do ! end do i

         ! The final AED2_inflow is located on the volume grid
         do i = 1, n
            do j = 1, ubnd_vol
               state%AED2_inflow(j, i) = Q_inp(i, j + 1) - Q_inp(i, j)
            end do
               !state%AED2_inflow(ubnd_vol + 1,ubnd_vol + 1) = 0
         end do

      end associate
   end subroutine

end module simstrat_aed2