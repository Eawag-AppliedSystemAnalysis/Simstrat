!     +---------------------------------------------------------------+
!     |  Simstrat - AED2 interface
!     +---------------------------------------------------------------+

module simstrat_aed2
   use strat_simdata
   use strat_grid
   use utilities
   use aed2_common
   use aed2_core

   implicit none
   private

   type, public :: SimstratAED2
      class(AED2Config), pointer :: aed2_cfg
      class(StaggeredGrid), pointer :: grid

      real(RK),allocatable,dimension(:) :: lKw    !# background light attenuation (m**-1)

      !# Arrays for state and diagnostic variables
      real(RK),allocatable,dimension(:,:) :: cc !# water quality array: nlayers, nvars
      real(RK),allocatable,dimension(:,:) :: cc_diag
      real(RK),allocatable,dimension(:) :: cc_diag_hz
      real(RK),allocatable,dimension(:) :: tss
      real(RK),allocatable,dimension(:) :: sed_zones

      !# Arrays for work, vertical movement, and cross-boundary fluxes
      real(RK),allocatable,dimension(:,:) :: ws
      real(RK),allocatable,dimension(:)   :: total
      real(RK),allocatable,dimension(:)   :: local
      !real(RK),allocatable,dimension(:) :: dz

      !# Arrays for environmental variables not supplied externally.
      real(RK),allocatable,dimension(:) :: par, pres
      real(RK),allocatable,dimension(:) :: uva, uvb, nir

         !# External variables
      real(RK) :: dt, dt_eff   ! External and internal time steps
      integer  :: w_adv_ctr    ! Scheme for vertical advection (0 if not used)
      real(RK),pointer,dimension(:) :: rad, z, salt, temp, rho, area
      real(RK),pointer,dimension(:) :: extc_coef, layer_stress
      real(RK),pointer :: precip, evap, bottom_stress
      real(RK),pointer :: I_0, wnd
      real(RK),allocatable,dimension(:) :: depth,layer_area

      character(len=48),allocatable :: names(:)
      character(len=48),allocatable :: bennames(:)

      integer,allocatable,dimension(:) :: externalid

      real(RK),allocatable,dimension(:) :: min_, max_

      integer :: n_aed2_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
      integer :: zone_var = 0

   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
      procedure, pass(self), public :: absorption_updateAED2
   end type SimstratAED2

contains

   subroutine init(self, grid, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg
      class(StaggeredGrid), target :: grid

      ! Local variables
      character(len=80) :: fname
      type(aed2_variable_t),pointer :: tvar

      character(len=64) :: models(64)
      namelist /aed2_models/ models
      integer i, j, status, rc, av, v, sv

      self%grid => grid
      self%aed2_cfg => aed2_cfg

      associate (n_aed2_vars => self%n_aed2_vars, &
                 n_vars => self%n_vars, &
                 n_vars_ben => self%n_vars_ben, &
                 n_vars_diag => self%n_vars_diag, &
                 n_vars_diag_sheet => self%n_vars_diag_sheet)

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

         n_aed2_vars = aed2_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

         print "(/,5X,'AED2 : n_aed2_vars = ',I3,' ; MaxLayers         = ',I4)",n_aed2_vars,self%grid%nz_grid
         print "(  5X,'AED2 : n_vars      = ',I3,' ; n_vars_ben        = ',I3)",n_vars,n_vars_ben
         print "(  5X,'AED2 : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I3,/)",n_vars_diag,n_vars_diag_sheet

         call check_data(self)

         !# names = grab the names from info
         allocate(self%names(n_vars),stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (names)'
         allocate(self%bennames(n_vars_ben),stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (bennames)'

         !# Now that we know how many vars we need, we can allocate space for them
         allocate(self%cc(self%grid%nz_grid, (n_vars + n_vars_ben)),stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (CC)'
         self%cc = 0.         !# initialise to zeroFarc

         allocate(self%min_((n_vars + n_vars_ben)))
         allocate(self%max_((n_vars + n_vars_ben)))
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
      do i=1,n_aed2_vars
         if ( aed2_get_var(i, tvar) ) then
               if ( tvar%sheet .and. .not. (tvar%diag .or. tvar%extern) ) then
                  j = j + 1
                  self%bennames(j) = trim(tvar%name)
                  self%min_(n_vars+j) = tvar%minimum
                  self%max_(n_vars+j) = tvar%maximum
                  print *,"     B(",j,") AED2 benthic(2D) variable: ", trim(self%bennames(j))
               end if
            end if
         end do

         j = 0
         do i=1,n_aed2_vars
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
         do i=1,n_aed2_vars
            if ( aed2_get_var(i, tvar) ) then
               if ( tvar%diag ) then
                  if (tvar%sheet ) then
                     j = j + 1
                     print *,"     D(",j,") AED2 diagnostic 2Dvariable: ", trim(tvar%name)
                  end if
               end if
            end if
         enddo

         allocate(self%externalid(n_aed2_vars))

         !# Now set initial values
         v = 0 ; sv = 0;
         do av=1,self%n_aed2_vars
            if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"
            if ( .not. ( tvar%extern .or. tvar%diag) ) then  !# neither global nor diagnostic variable
               if ( tvar%sheet ) then
                  !AED2_InitCondition(self%cc(:, nvars + sv), tvar%name, tvar%initial)
                  sv = sv + 1
                  self%cc(:, n_vars+sv) = tvar%initial
                  write(6,*) 'sheet', tvar%name, tvar%initial
               else
                  v = v + 1
                  call AED2_InitCondition(self, self%cc(:, v), tvar%name, tvar%initial)
                  write(6,*) tvar%name, tvar%initial

               end if
            end if
         end do

         !# Allocate diagnostic variable array and set all values to zero.
         !# (needed because time-integrated/averaged variables will increment rather than set the array)
         allocate(self%cc_diag(self%grid%nz_grid, n_vars_diag),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (cc_diag)'
         self%cc_diag = zero_

         !# Allocate diagnostic variable array and set all values to zero.
         !# (needed because time-integrated/averaged variables will increment rather than set the array)
         allocate(self%cc_diag_hz(n_vars_diag_sheet),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (cc_diag_hz)'
         self%cc_diag_hz = zero_

         !# Allocate array with vertical movement rates (m/s, positive for upwards),
         !# and set these to the values provided by the model.
         allocate(self%ws(self%grid%nz_grid, n_vars),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (ws)'
         self%ws = zero_

         !# Allocate array for photosynthetically active radiation (PAR).
         !# This will be calculated internally during each time step.
         allocate(self%par(self%grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (par)'
         self%par = zero_

         allocate(self%nir(self%grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (nir)'
         self%nir = zero_
         allocate(self%uva(self%grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (uva)'
         self%uva = zero_
         allocate(self%uvb(self%grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (uvb)'
         self%uvb = zero_

         !allocate(self%dz(self%grid%nz_grid),stat=rc)
         !self%dz = zero_

         allocate(self%sed_zones(self%grid%nz_grid))
         !# Allocate array for local pressure.
         !# This will be calculated [approximated] from layer depths internally
         !# during each time step.
         allocate(self%pres(self%grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (pres)'
         self%pres = zero_

         allocate(self%tss(self%grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (tss)'
         self%tss = zero_

         write(*,"(/,5X,'----------  AED2 config : end  ----------',/)")
      end associate
   end subroutine


   subroutine update(self, state)
      implicit none
      class(SimstratAED2) :: self
      class(ModelState) :: state

      ! Local variables
      type (aed2_column_t) :: column(self%n_aed2_vars), column_sed(self%n_aed2_vars)
      real(RK) :: flux_ben(self%n_vars + self%n_vars_ben), flux_atm(self%n_vars + self%n_vars_ben)
      real(RK) :: flux(self%grid%nz_occupied, self%n_vars + self%n_vars_ben)
      real(RK) :: flux_zone(self%aed2_cfg%n_zones, self%n_vars + self%n_vars_ben)
      integer :: v, split, lev

      !# Calculate local pressure
      self%pres(1:self%grid%ubnd_vol) = -self%grid%z_volume(1:self%grid%ubnd_vol)

      call define_column(self, state, column, self%grid%nz_occupied, self%cc, self%cc_diag, self%cc_diag_hz, flux, flux_atm, flux_ben)
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

         !call calculate_fluxes(self, state, column, column_sed, self%aed2_cfg%n_zones, flux(:,:), flux_atm, flux_ben, flux_zone(:,:))

         ! Update the water column layers
         do v = 1, self%n_vars
            do lev = 1, self%grid%nz_occupied
               !cc(lev, v) = cc(lev, v) + dt_eff*flux(lev, v)
            end do
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

   SUBROUTINE define_column(self, state, column, top, cc, cc_diag, cc_diag_hz, flux_pel, flux_atm, flux_ben)
   !-------------------------------------------------------------------------------
   ! Set up the current column pointers
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      class(SimstratAED2) :: self
      class(ModelState) :: state
      type (aed2_column_t), intent(inout) :: column(:)
      integer, intent(in)  :: top
      real(RK), intent(in) :: cc(:,:)       !# (n_layers, n_vars)
      real(RK), intent(in) :: cc_diag(:,:)  !# (n_layers, n_vars)
      real(RK), intent(in) :: cc_diag_hz(:)
      real(RK), intent(inout) :: flux_pel(:,:) !# (n_layers, n_vars)
      real(RK), intent(inout) :: flux_atm(:)   !# (n_vars)
      real(RK), intent(inout) :: flux_ben(:)   !# (n_vars)
   !
   !LOCALS
      integer :: av !, i
      integer :: v, d, sv, sd, ev
      type(aed2_variable_t), pointer :: tvar
      real(RK), target :: extern_target(self%grid%nz_grid), extern_target_fce(self%grid%nz_grid + 1)
      real(RK), target :: extern_target_sheet
      real(RK), target :: cc_target(self%grid%nz_grid), flux_pel_target(self%grid%nz_occupied)
      real(RK), target :: flux_atm_target, flux_ben_target
   !-------------------------------------------------------------------------------
   !BEGIN
      v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
      do av=1,self%n_aed2_vars
         if ( .not.  aed2_get_var(av, tvar) ) stop "Error getting variable info"

         if ( tvar%extern ) then !# global variable
            ev = ev + 1
            select case (tvar%name)
               case ( 'temperature' ) ; extern_target = state%T; column(av)%cell => extern_target;
               case ( 'salinity' )    ; extern_target = state%S; column(av)%cell => extern_target;
               case ( 'density' )     ; extern_target = state%rho; column(av)%cell => extern_target
               case ( 'layer_ht' )    ; extern_target = self%grid%h(1:self%grid%nz_grid); column(av)%cell => extern_target
               case ( 'extc_coef' )   ; extern_target = state%absorb_vol; column(av)%cell => extern_target
               case ( 'tss' )         ; extern_target = self%tss; column(av)%cell => extern_target
               case ( 'par' )         ; extern_target = self%par; column(av)%cell => extern_target
               case ( 'nir' )         ; extern_target = self%nir; column(av)%cell => extern_target
               case ( 'uva' )         ; extern_target = self%uva; column(av)%cell => extern_target
               case ( 'uvb' )         ; extern_target = self%uvb; column(av)%cell => extern_target
               case ( 'pressure' )    ; extern_target = self%pres; column(av)%cell => extern_target
               case ( 'depth' )       ; extern_target = self%grid%z_volume; column(av)%cell => extern_target
               case ( 'sed_zone' )    ; extern_target_sheet = self%sed_zones(1); column(av)%cell_sheet => extern_target_sheet
               case ( 'wind_speed' )  ; extern_target_sheet = state%uv10; column(av)%cell_sheet => extern_target_sheet
               case ( 'par_sf' )      ; extern_target_sheet = state%rad0; column(av)%cell_sheet => extern_target_sheet
               case ( 'taub' )        ; extern_target_sheet = state%u_taub; column(av)%cell_sheet => extern_target_sheet
               case ( 'lake_depth' )  ; extern_target_sheet = self%grid%z_face(self%grid%ubnd_fce); column(av)%cell_sheet => extern_target_sheet
               case ( 'layer_area' )  ; extern_target = self%grid%Az_vol; column(av)%cell => extern_target
               case default ; call error("External variable "//TRIM(tvar%name)//" not found.")
            end select
         elseif ( tvar%diag ) then  !# Diagnostic variable
            if ( tvar%sheet ) then
               sd = sd + 1
               extern_target_sheet = cc_diag_hz(sd)
               column(av)%cell_sheet => extern_target_sheet
            else
               d = d + 1
               extern_target = cc_diag(:,d)
               column(av)%cell => extern_target
            end if
         else    !# state variable
            if ( tvar%sheet ) then
               sv = sv + 1
               if ( tvar%bot ) then
                  extern_target_sheet = cc(1, self%n_vars + sv)
                  column(av)%cell_sheet => extern_target_sheet
   !            print *,'av',av,sv
               elseif ( tvar%top ) then
                  extern_target_sheet = cc(top, self%n_vars + sv)
                  column(av)%cell_sheet => extern_target_sheet
               endif

               flux_ben_target = flux_ben(self%n_vars + sv)
               flux_atm_target = flux_atm(self%n_vars + sv)

               column(av)%flux_ben => flux_ben_target
               column(av)%flux_atm => flux_atm_target
            else
               v = v + 1
               cc_target = cc(:,v)
               flux_atm_target = flux_atm(v)
               flux_pel_target = flux_pel(:,v)
               flux_ben_target = flux_ben(v)


               column(av)%cell => cc_target
               column(av)%flux_atm => flux_atm_target
               column(av)%flux_pel => flux_pel_target
               column(av)%flux_ben => flux_ben_target
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

   SUBROUTINE calculate_fluxes(self, state, column, column_sed, nsed, flux_pel, flux_atm, flux_ben, flux_zon)
   !-------------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !-------------------------------------------------------------------------------
   use,intrinsic :: ieee_arithmetic

   !ARGUMENTS
      class(SimstratAED2) :: self
      class(ModelState) :: state
      type (aed2_column_t), intent(inout) :: column(:)
      type (aed2_column_t), intent(inout) :: column_sed(:)
      integer, intent(in) :: nsed
      real(RK), intent(inout) :: flux_pel(:,:) !# (wlev, n_vars)
      real(RK), intent(inout) :: flux_atm(:)   !# (n_vars)
      real(RK), intent(inout) :: flux_ben(:)   !# (n_vars)
      real(RK), intent(inout) :: flux_zon(:,:) !# (n_zones)
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
      flux_pel = zero_
      flux_atm = zero_
      flux_ben = zero_

      !# Start with calculating all flux terms for rhs in mass/m3/s
      !# Includes (1) benthic flux, (2) surface exchange and (3) water column kinetics
      !# as calculated by glm


      !# (1) BENTHIC FLUXES
      if ( self%aed2_cfg%benthic_mode .gt. 1 ) then
!          !# Multiple static sediment zones are simulated, and therfore overlying
!          !# water conditions need to be aggregated from multiple cells/layers, and output flux
!          !# needs disaggregating from each zone back to the overlying cells/layers

!          do zon=1,nsed
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
         flux_pel(1, :) = max(-1.0 * self%cc(1, :), flux_pel(1, :)/self%grid%h(1))

         if ( self%aed2_cfg%benthic_mode .EQ. 1 ) then
            do lev=2,self%grid%nz_occupied
               !# Calculate temporal derivatives due to benthic fluxes.
               call aed2_calculate_benthic(column, lev)

               !# Limit flux out of bottom layers to concentration of that layer
               !# i.e. don't flux out more than is there
               !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
               !# scaled to proportion of area that is "bottom"
               flux_pel(lev, :) = max(-1.0 * self%cc(lev, :), flux_pel(lev, :)/self%grid%h(lev))
               flux_pel(lev, :) = flux_pel(lev, :) * (self%grid%Az_vol(lev) - self%grid%Az_vol(lev - 1))/self%grid%Az(lev)
            end do
         end if
      end if

      !# (2) SURFACE FLUXES
      !# Calculate temporal derivatives due to air-water exchange.
      if (.not. (state%ice_h > 0)) then !# no surface exchange under ice cover
         call aed2_calculate_surface(column, self%grid%nz_occupied)

         !# Distribute the fluxes into pelagic surface layer
         flux_pel(self%grid%nz_occupied, :) = flux_pel(self%grid%nz_occupied, :) + flux_atm(:)/self%grid%h(self%grid%nz_occupied)
      end if

      !# (3) WATER COLUMN KINETICS
      !# Add pelagic sink and source terms for all depth levels.
      do lev=1,self%grid%nz_occupied
         call aed2_calculate(column, lev)
      end do
   END SUBROUTINE calculate_fluxes

   subroutine AED2_InitCondition(self, var,varname,default_val)
   !#################################### written/copied by A. Gaudard, 2015
        implicit none

        class(SimstratAED2) :: self
        real(RK), intent(inout) :: var(1:self%grid%nz_grid) !Vector if initial conditions
        real(RK), intent(in) :: default_val !Depth-independent value (default from fabm.nml)
        character(len=*), intent(in) :: varname !Identifying the variable

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
    9   nval = i                               ! Number of values
        if (nval<0) then
            write(6,*) 'Error reading ', trim(varname), ' initial conditions file (no data found).'
            stop
        end if
        close(14)
        do i=1,nval
            z_read(i) = abs(z_read(i))               ! Make depths positive
        end do
        z_read_depth = z_read(1)                     ! Initial depth (top-most)

        do i=1,nval
            z_read(nval + 1 - i) = self%grid%z_zero - z_read(i)
            var_read(nval + 1 - i) = var_read(i)
        end do

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


end module simstrat_aed2