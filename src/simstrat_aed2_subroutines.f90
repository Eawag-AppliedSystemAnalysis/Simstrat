! Subroutines used in Simstrat-AED2 module
! The functions are copied or inspired by the GLM-AED2 coupling (http://aed.see.uwa.edu.au/research/models/GLM/)

subroutine allocate_memory(self)
   implicit none

   ! Arguments
   class(SimstratAED2) :: self

   ! Local variables
   integer status
   allocate(self%column(self%n_aed2_vars),stat=status)
   allocate(self%column_sed(self%n_aed2_vars),stat=status)

   ! names = grab the names from info
   allocate(self%names(self%n_vars),stat=status)
   if (status /= 0) stop 'allocate_memory(): Error allocating (names)'
   allocate(self%bennames(self%n_vars_ben),stat=status)
   if (status /= 0) stop 'allocate_memory(): Error allocating (bennames)'
   allocate(self%diagnames(self%n_vars_diag),stat=status)
   if (status /= 0) stop 'allocate_memory(): Error allocating (diagnames)'

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
   implicit none

   ! Arguments
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
               self%diagnames(j) = trim(tvar%name)
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
         do i=1,self%n_aed2_vars
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