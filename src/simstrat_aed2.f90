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
      real(RK),allocatable,dimension(:) :: dz

      !# Arrays for environmental variables not supplied externally.
      real(RK),allocatable,dimension(:) :: par, pres
      real(RK),allocatable,dimension(:) :: uva, uvb, nir

         !# External variables
      real(RK) :: dt, dt_eff   ! External and internal time steps
      integer  :: w_adv_ctr    ! Scheme for vertical advection (0 IF not used)
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
      integer i, j, status, rc

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

         allocate(self%dz(self%grid%nz_grid),stat=rc)
         self%dz = zero_

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
      type (aed2_column_t) :: column(self%n_aed2_vars)
      real(RK) :: flux_ben(self%n_vars + self%n_vars_ben), flux_atm(self%n_vars + self%n_vars_ben)
      real(RK) :: flux(self%grid%nz_occupied, self%n_vars + self%n_vars_ben)
      !real(RK) :: flux_zone(self%n_zones, self%n_vars+n_vars_ben)

      integer :: v, split, lev

      !# Calculate local pressure
      self%pres(1:self%grid%ubnd_vol) = -self%grid%z_volume(1:self%grid%ubnd_vol)

      call define_column(self, state, column, self%grid%nz_occupied, self%cc, self%cc_diag, self%cc_diag_hz, flux, flux_atm, flux_ben)

      ! If sediment layers are simulated
      !IF (benthic_mode .GT. 1) CALL define_sed_column(column_sed, n_zones, flux, flux_atm, flux_ben)

      self%cc_diag = 0.
      self%cc_diag_hz = 0.

      !call check_states(column,self%grid%nz_occupied)

      do split=1, self%aed2_cfg%split_factor

         call absorption_updateAED2(self, column, state)

         !# Fudge
         self%nir(:) = (self%par(:)/0.45) * 0.51
         self%uva(:) = (self%par(:)/0.45) * 0.035
         self%uvb(:) = (self%par(:)/0.45) * 0.005

         !call calculate_fluxes(column, wlev, column_sed, n_zones, flux(:,:), flux_atm, flux_ben, flux_zone(:,:))
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
               !case ( 'sed_zone' )    ; column(av)%cell_sheet => self%sed_zones(1)
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

!    SUBROUTINE check_states(column, wlev)
! !-------------------------------------------------------------------------------
! #ifdef HAVE_IEEE_ARITH
! !USES
!    USE IEEE_ARITHMETIC
! #endif
! !
! !ARGUMENTS
!    TYPE (aed2_column_t),INTENT(inout) :: column(:)
!    INTEGER,INTENT(in) :: wlev
! !
! !LOCALS
!    TYPE(aed2_variable_t),POINTER :: tv
!    INTEGER i,v,lev
! !
! !-------------------------------------------------------------------------------
! !BEGIN
!    DO lev=1, wlev
!       CALL aed2_equilibrate(column, lev)    !MH this should be in the main do_glm routine ????!!!
!       v = 0
!       DO i=1,n_aed2_vars
!          IF ( aed2_get_var(i, tv) ) THEN
!             IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
!                v = v + 1
!                IF ( repair_state ) THEN
!                   IF ( .NOT. isnan(min_(v)) ) THEN
!                      IF ( cc(lev, v) < min_(v) ) cc(lev, v) = min_(v)
!                   ENDIF
!                   IF ( .NOT. isnan(max_(v)) ) THEN
!                      IF ( cc(lev, v) > max_(v) ) cc(lev, v) = max_(v)
!                   ENDIF
!                ENDIF
!             ENDIF
!          ENDIF
!       ENDDO
!    ENDDO
! END SUBROUTINE check_states


end module simstrat_aed2