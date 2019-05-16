!<    +---------------------------------------------------------------+
!     |  Forcing module
!     |  - Adjusted heat fluxes to snow- and ice-covered lake
!     |  - Reads forcing file and updates state
!     |  - Updates Coriolis force terms
!<    +---------------------------------------------------------------+

module strat_forcing
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use utilities
   implicit none
   private

   type, public :: ForcingModule
      integer :: nz_grid_max
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param
      character(len=:), allocatable  :: file  ! Forcing file name
   contains
      procedure, pass :: init => forcing_init
      procedure, pass :: read => forcing_read
      procedure, pass :: update => forcing_update
      procedure, pass :: update_coriolis => forcing_update_coriolis
   end type

contains

   subroutine forcing_init(self, model_config, model_param, forcing_file, grid)
      implicit none
      class(ForcingModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param
      character(len=:), allocatable :: forcing_file

      self%cfg => model_config
      self%param => model_param
      self%grid => grid
      self%file = forcing_file

      if(self%cfg%snow_model == 1) call warn('Snow module need precipitation input, check forcing file.')
      !If precipitation column is missing in forcing file, Simstrat will read dates as precipitation
   end subroutine

   !Read forcing file to get values A_cur at given datum
   ! This method has more or less been copied from old simstrat
   !AG 2014: revision + correction
   !####################################################################
   subroutine forcing_read(self, datum, A_s, A_e, A_cur, nval, idx)
      !####################################################################

      implicit none

      ! Global variables
      class(ForcingModule) :: self
      integer, intent(in) :: nval
      logical, intent(in) :: idx
      real(RK), intent(in) :: datum !Required date
      real(RK), intent(inout) :: A_cur(1:nval), A_s(1:nval), A_e(1:nval) !output values, Start and end values

      ! Local variables
      real(RK) :: tb_start, tb_end !Start and end time
      integer :: eof, i

      save tb_start, tb_end
      save eof

      if (idx) then

         open (20, status='old', file=self%file)
         eof = 0
         !Read first values
         read (20, *, end=9)
         read (20, *, end=9) tb_start, (A_s(i), i=1, nval)
         if (datum < tb_start) then
            write(6,*) '[WARNING] ','First forcing date after simulation start time. datum=', datum,  'start=', tb_start
         end if
         read (20, *, end=7) tb_end, (A_e(i), i=1, nval)

         call ok("Forcing input file successfully read")
      end if

      if (datum <= tb_start .or. eof == 1) then !If datum before first date or end of file reached
         goto 8
      else
         do while (datum > tb_end) !Move to appropriate interval to get correct values
            tb_start = tb_end
            A_s(1:nval) = A_e(1:nval)
            !Read next value
            read (20, *, end=7) tb_end, (A_e(i), i=1, nval)
         end do
         !Linearly interpolate values at correct datum
         A_cur(1:nval) = A_s(1:nval) + (datum - tb_start)*(A_e(1:nval) - A_s(1:nval))/(tb_end - tb_start)
      end if
      return

  7   eof = 1
      if(datum>tb_start) call warn('Last forcing date before simulation end time.')
      write(6,*) datum, tb_start


  8   A_cur(1:nval) = A_s(1:nval)       !Take first value of current interval
      return

  9   call error('Unable to read forcing file (no data found).')

   end subroutine

   !Compute appropriate forcing parameters at given datum
   !Copied from old simstrat, NEEDS refactoring!
   !AG 2014: revision
   !######################################################################
   subroutine forcing_update(self, state)
      !######################################################################
      implicit none
      class(ForcingModule) :: self
      class(ModelState) :: state

      ! Local iables
      integer :: nval_offset
      real(RK) :: tau
      real(RK) :: A_s(8), A_e(8), A_cur(8) ! adopted for rain (8 positions, previus 7)
      real(RK) :: fu, Vap_wat, heat0, emissivity
      real(RK) :: T_surf, F_glob, Vap_atm, Cloud
      real(RK) :: H_A, H_K, H_V, H_W, H_tot
      real(RK) :: F_snow, F_snowice, F_ice
      real(RK) :: qa, q0, sh_c
      save A_s, A_e
      associate (cfg=>self%cfg, param=>self%param)

      if((state%snow_h + state%white_ice_h + state%black_ice_h) > 0) then !ice cover
         T_surf = param%Freez_Temp
      else
         T_surf = state%T(self%grid%ubnd_vol) !free water
      end if

      ! number of values to read, depending on filtered wind and precipitation
      if (cfg%use_filtered_wind .and. cfg%ice_model == 0) then
         nval_offset = 1
      else if (cfg%ice_model == 1 .and. cfg%use_filtered_wind) then
         nval_offset = 2
      else if (cfg%ice_model == 1) then
         nval_offset = 1
      else
         nval_offset = 0
      end if

      if (cfg%forcing_mode == 1) then
         if (cfg%ice_model == 1) then
            call error('Ice module not compatible with forcing mode 1, use 2, 3 or 5.')
         end if

         call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)
         call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)
         state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
         state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
         state%uv10 = sqrt(state%u10**2 + state%v10**2) !AG 2014
         state%SST = A_cur(3) !Sea surface temperature
         state%rad0 = A_cur(4)*(1 - albsw)*(1 - beta_sol) ! MS: added beta_sol and albsw
         state%heat = 0.0_RK
         state%T_atm = 0.0_RK
         state%precip = 0.0_RK
         if (cfg%use_filtered_wind) then
            state%Wf = A_cur(5) !AG 2014
         end if

      else if (cfg%forcing_mode >= 2) then
         if (cfg%forcing_mode == 2) then ! date, U,V,Tatm,Hsol,Vap
            call self%read (state%datum, A_s, A_e, A_cur, 5 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            state%T_atm = A_cur(3)

            if (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) then !Ice
               F_glob = (A_cur(4)*(1 - ice_albedo)) * param%p_albedo
            else if (state%white_ice_h > 0 .and. state%snow_h == 0) then !Snowice
               F_glob = (A_cur(4)*(1 - snowice_albedo)) * param%p_albedo
            else if (state%snow_h > 0) then !Snow
               F_glob = (A_cur(4)*(1 - snow_albedo)) * param%p_albedo
            else !Water
               F_glob = A_cur(4)*(1 - albsw)
            end if

            Vap_atm = A_cur(5)
            Cloud = 0.5
            if (cfg%use_filtered_wind) state%Wf = A_cur(6) !AG 2014
            if (cfg%snow_model == 1 .and. cfg%use_filtered_wind) then
               state%precip = A_cur(7)
            else if (cfg%snow_model == 1) then
               state%precip = A_cur(6)
            end if

         else if (cfg%forcing_mode == 3) then ! date,U10,V10,Tatm,Hsol,Vap,Clouds
            call self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            state%T_atm = A_cur(3)

            if (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) then !Ice
               F_glob = (A_cur(4)*(1 - ice_albedo)) * param%p_albedo
            else if (state%white_ice_h > 0 .and. state%snow_h == 0) then !Snowice
               F_glob = (A_cur(4)*(1 - snowice_albedo)) * param%p_albedo
            else if (state%snow_h > 0) then !Snow
               F_glob = (A_cur(4)*(1 - snow_albedo)) * param%p_albedo
            else !Water
               F_glob = A_cur(4)*(1 - albsw)
            end if

            Vap_atm = A_cur(5)
            Cloud = A_cur(6)
            if (Cloud < 0 .or. Cloud > 1) then
               write (*,'(A,F12.6)') 'Cloud : ' , Cloud
               write (*,'(A,F12.6)') 'Date  : ' , state%datum
               call error('Cloudiness should always be between 0 and 1.')
            end if
            if (cfg%use_filtered_wind) state%Wf = A_cur(7) !AG 2014
            if (cfg%snow_model == 1 .and. cfg%use_filtered_wind) then
               state%precip = A_cur(8)
            else if (cfg%snow_model == 1) then
               state%precip = A_cur(7)
            end if

         else if (cfg%forcing_mode == 4) then ! date,U10,V10,Hnet,Hsol
            if (cfg%ice_model == 1) then
               call error('Ice module not compatible with forcing mode 4, use 2, 3 or 5.')
            end if

            call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            heat0 = A_cur(3) !MS 2014
            F_glob = A_cur(4)*(1 - albsw)
            state%T_atm = 0.0_RK
            if (cfg%use_filtered_wind) state%Wf = A_cur(5) !AG 2014
         !UK added forcing mode with incomming long-wave radiation instead of cloudiness
         else if (cfg%forcing_mode == 5) then ! date,U10,V10,Tatm,Hsol,Vap,ILWR
            call self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            state%T_atm = A_cur(3)

            if (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) then !Ice
               F_glob = (A_cur(4)*(1 - ice_albedo)) * param%p_albedo
            else if (state%white_ice_h > 0 .and. state%snow_h == 0) then !Snowice
               F_glob = (A_cur(4)*(1 - snowice_albedo)) * param%p_albedo
            else if (state%snow_h > 0) then !Snow
               F_glob = (A_cur(4)*(1 - snow_albedo)) * param%p_albedo
            else !Water
               F_glob = A_cur(4)*(1 - albsw)
            end if

            Vap_atm = A_cur(5)
            H_A = A_cur(6)
            if (cfg%use_filtered_wind) state%Wf = A_cur(7) !AG 2014
            if (cfg%snow_model == 1 .and. cfg%use_filtered_wind) then
               state%precip = A_cur(8)
            else if (cfg%snow_model == 1) then
               state%precip = A_cur(7)
            end if
         else
            call error('Wrong forcing type (must be 1, 2, 3, 4 or 5).')
         end if
         state%uv10 = sqrt(state%u10**2 + state%v10**2) !AG 2014

         if (cfg%forcing_mode /= 4) then ! Heat fluxes calculations (forcing 2, 3 and 5)
            if (state%black_ice_h + state%white_ice_h == 0) then !Free water
               ! Wind function (Adams et al., 1990), changed by MS, June 2016
               ! Factor 0.6072 to account for changing wind height from 10 to 2 m
               ! Further evaluation of evaporation algorithm may be required.
               fu = sqrt((2.7_RK*max(0.0_RK,(T_surf-state%T_atm)/(1-0.378_RK*Vap_atm/param%p_air))**0.333_RK)**2 + (0.6072_RK*3.1_RK*state%uv10)**2)
               ! Wind function (Livingstone & Imboden 1989)
               !fu = 4.40_RK + 1.82_RK*state%uv10 + 0.26_RK*(T_surf - T_atm)
               !fu = 5.44+2.19*wind+0.24*(T_surf-T_atm)
               fu = fu*param%p_windf ! Provided fitting factor p_windf (~1)
               ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
               Vap_wat = 10**((0.7859_RK + 0.03477_RK*T_surf)/(1 + 0.00412_RK*T_surf))
               Vap_wat = Vap_wat*(1 + 1e-6_RK*param%p_air*(4.5_RK + 0.00006_RK*T_surf**2))

               ! Long-wave radiation from sky (Livingstone & Imboden 1989)
               ! H_A = 1.24*sig*(1-r_a)*(1+0.17*Cloud**2)*(Vap_atm/(273.15+T_atm))**(1./7)*(273.15+T_atm)**4
               ! Long-wave radiation according to Dilley and O'Brien
               ! see Flerchinger et al. (2009)
               ! UK changed to read from file if forcing_mode=5
               if (cfg%forcing_mode /= 5) then
                  H_A = (1 - r_a)*((1 - 0.84_RK*Cloud)*(59.38_RK + 113.7_RK*((state%T_atm + 273.15_RK)/273.16_RK)**6&
                  &   + 96.96_RK*sqrt(465*Vap_atm/(state%T_atm + 273.15_RK)*0.04_RK))/5.67e-8_RK/ &
                  &   (state%T_atm + 273.15_RK)**4 + 0.84_RK*Cloud)*5.67e-8_RK*(state%T_atm + 273.15_RK)**4
               end if
               H_A = H_A*param%p_radin ! Provided fitting factor p_radin (~1)

               H_W = -emiss_water*sig*(T_surf + 273.15_RK)**4

               ! Flux of sensible heat (convection)
               H_K = -B0*fu*(T_surf - state%T_atm)

               ! Flux of latent heat (evaporation, condensation)
               H_V = -fu*(Vap_wat - Vap_atm)

               ! Global heat flux (positive: air to water, negative: water to air)
               state%heat = H_A + H_W + H_K + H_V + F_glob * beta_sol !MS: added term with beta_sol
               ! Removal of solar short-wave radiation absorbed in first water cell
               state%rad0 = F_glob * (1 - beta_sol) !MS: added beta_sol

               state%heat_snow = 0 !Heat snow
               state%heat_snowice = 0 !Heat snowice
               state%heat_ice = 0 !Heat ice

            else if (state%black_ice_h > 0 .or. state%white_ice_h > 0) then !Ice Cover
               !Light penetration in snow, ice and snowice as well as wind blocking added 2018 by Love Raaman

               if (state%T_atm >= param%Freez_Temp) then!melting occures when air temp above freez point,
                  !then activate surface heat fluxes following
                  !Matti Leppäranta (2009), Modelling the Formation and Decay of Lake Ice DOI: 10.1007/978-90-481-2945-4_5 In book: The Impact of Climate Change on European Lakes
                  ! In G. George (Ed.), The Impact of Climate Change on European Lakes (pp. 63–83).
                  ! Dordrecht: Springer Netherlands. https://doi.org/10.1007/978-90-481-2945-4_5
                  ! and with corrections in
                  ! Leppäranta, M. (2014). Freezing of lakes and the evolution of their ice cover.
                  ! New York: Springer. ISBN 978-3-642-29080-0
                  if (state%snow_h == 0) then !Ice Cover (ice and snowice)
                     emissivity = emiss_ice
                  else !Snow Cover
                     !varies from 0.8 to 0.9 depending on snow density
                     emissivity = 5.0e-4_RK * state%snow_dens + 6.75e-1_RK
                  end if
                  ! obs fiting factors param%p_radin and param%p_windf not applied to ice covered lake
                  if (cfg%forcing_mode /= 5) then
                     H_A = (Ha_a + Ha_b * (Vap_atm**(1.0_RK/2.0_RK))) * (1 + Ha_c * Cloud**2) * sig * (state%T_atm + 273.15_RK)**4
                  end if
                  H_W = -emissivity * sig * (T_surf + 273.15_RK)**4
                  H_K = rho_air * cp_air * Hk_CH * (state%T_atm - T_surf) * state%uv10
                  sh_c = 0.622/param%p_air! converter from absolut vapour pressure to specific humidity (Lepparanta 2015)
                  qa  = sh_c * Vap_atm!specific humidity air
                  q0  = sh_c * 6.11!specific humidity for saturation levels (6.11 mbar) at surface (ice) at 0°C (Lepparanta 2015)
                  H_V = rho_air * (l_h + l_e) * Hv_CE * (qa - q0) * state%uv10 ! through sublimation (solid to gas)
               else ! if no melting only considering penetraiting solar radiation since ice formation is parameterised towards temperature
                  H_A = 0
                  H_W = 0
                  H_K = 0
                  H_V = 0
               end if
                H_tot = H_A + H_W + H_K + H_V
               if (H_tot < 0) then
                  H_tot = 0
               end if

               ! Global heat flux (positive: air to water, negative: water to air) !MS: added beta_sol ; LRV added lambda_snow_ice
               ! Leppäranta, M. (2014), Eq. 6.12
               ! Removal of solar short-wave radiation absorbed in snow, snowice, ice and first water cell (works also when x_h = 0)
               state%heat = F_glob * exp(-lambda_snow*state%snow_h -lambda_snowice*state%white_ice_h -lambda_ice*state%black_ice_h) * beta_sol
               state%rad0 = F_glob * exp(-lambda_snow*state%snow_h -lambda_snowice*state%white_ice_h -lambda_ice*state%black_ice_h) * (1 - beta_sol)

               !Heat flux into snow, ice or snowice layer.
               ! Light absorption each layer
               ! Leppäranta, M. (2014), Eq. 6.12
               F_snow    = F_glob * (1 - exp(-lambda_snow*state%snow_h))
               F_snowice = F_glob * (exp(-lambda_snow*state%snow_h) - exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h))
               F_ice     = F_glob * (exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h) - exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h))
               if (F_snow < 0 .or. F_snowice < 0 .or. F_ice < 0 .or. state%heat < 0 .or. state%rad0 < 0) then
                  call error('Negative heat flux not allowed for melting.')
               end if

               !Light + other heat fluxes into top layer
               if (state%snow_h > 0 .and. state%white_ice_h > 0 .and. state%black_ice_h > 0) then ! snow, snowice and ice
                  state%heat_snow = H_tot + F_snow
                  state%heat_snowice = F_snowice
                  state%heat_ice = F_ice
               else if (state%snow_h > 0 .and. state%white_ice_h > 0 .and. state%black_ice_h == 0) then! snow and snowice
                  state%heat_snow = H_tot + F_snow
                  state%heat_snowice = F_snowice
                  state%heat_ice = 0
               else if (state%snow_h > 0 .and. state%white_ice_h == 0 .and. state%black_ice_h > 0) then! snow and ice
                  state%heat_snow = H_tot + F_snow
                  state%heat_snowice = 0
                  state%heat_ice = F_ice
               else if (state%snow_h == 0 .and. state%white_ice_h > 0 .and. state%black_ice_h > 0) then! snowice and ice
                  state%heat_snow = 0
                  state%heat_snowice = H_tot + F_snowice
                  state%heat_ice = F_ice
               else if (state%snow_h == 0 .and. state%white_ice_h == 0 .and. state%black_ice_h > 0) then! ice
                  state%heat_snow = 0
                  state%heat_snowice = 0
                  state%heat_ice = H_tot + F_ice
               else if (state%snow_h == 0 .and. state%white_ice_h > 0 .and. state%black_ice_h == 0) then! snowice
                  state%heat_snow = 0
                  state%heat_snowice = H_tot + F_snowice
                  state%heat_ice = 0
               end if
               ! Suppress wind turbulence with wind lid (heat flux affected by wind still active on snow and ice)
               state%u10 = 0
               state%v10 = 0
               state%uv10 = 0
            end if
            ! save for output, not used in calculations
            state%ha = H_A
            state%hw = H_W
            state%hk = H_K
            state%hv = H_V
         else !Forcing mode 4
            state%heat = heat0 + F_glob*beta_sol !MS: added term with beta_sol
            state%heat_snow = 0 !Heat snow
            state%heat_snowice = 0 !Heat snowice
            state%heat_ice = 0 !Heat ice
         end if
         if (cfg%ice_model == 0) then
            if (state%T(self%grid%ubnd_vol) < 0 .and. state%heat < 0) then
               state%heat = 0
               state%T(self%grid%ubnd_vol) = 0
            end if
         end if
      end if

      !Drag coefficient as a function of wind speed (AG 2014)
      if (cfg%wind_drag_model == 1) then !constant wind drag coefficient
         state%C10 = param%C10_constant
      else if (cfg%wind_drag_model == 2) then !Ocean model
         state%C10 = param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
      else if (cfg%wind_drag_model == 3) then !Lake model (Wüest and Lorke 2003)
         if (state%uv10 <= 0.1) then
            state%C10 = param%C10_constant*0.06215_RK
         else if (state%uv10 <= 3.85_RK) then
              state%C10 = param%C10_constant*0.0044_RK*state%uv10**(-1.15_RK)
         else !Polynomial approximation of Charnock's law
            state%C10 = param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
            !C10 = -0.000000385341*wind**2+0.0000656519*wind+0.000703768
            !C10 = 0.0000000216952*wind**3-0.00000148692*wind**2+0.0000820705*wind+0.000636251
         end if
      end if

      tau = state%C10*rho_air/rho_0*state%uv10**2
      state%u_taus = sqrt(tau)

      state%tx = state%C10*rho_air/rho_0*state%uv10*state%u10
      state%ty = state%C10*rho_air/rho_0*state%uv10*state%v10
      return
      end associate
   end subroutine

   ! enforces coriolis forces
   subroutine forcing_update_coriolis(self, state)
      implicit none
      class(ForcingModule) :: self
      class(ModelState) :: state
      real(RK) :: cori
      real(RK), dimension(size(state%U)) :: u_temp
      associate (grid=>self%grid, dt=>state%dt, param=>self%param)
         ! calculate u_taub before changing U resp V
         state%u_taub = sqrt(state%drag*(state%U(1)**2 + state%V(1)**2))

         ! Calculate coriolis parameter based on latitude
         cori = 2.0_RK*7.292e-5_RK*sin(param%Lat*pi/180.0_RK)

         !Update state based on coriolis parameter
         u_temp = state%U

         state%U(1:grid%ubnd_vol) = state%U(1:grid%ubnd_vol)*cos(Cori*dt) + state%V(1:grid%ubnd_vol)*sin(Cori*dt)
         state%V(1:grid%ubnd_vol) = -u_temp(1:grid%ubnd_vol)*sin(Cori*dt) + state%V(1:grid%ubnd_vol)*cos(Cori*dt)

         return
      end associate
   end subroutine

end module strat_forcing

