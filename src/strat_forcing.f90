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
      procedure, pass :: init_albedo => forcing_init_albedo
      procedure, pass :: update_albedo => forcing_update_albedo
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

      if(self%cfg%snow_model == 1) call warn('The snow module is turned on. This module needs precipitation data, note that the last column in the forcing file will be interpreted as precipitation.')
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
      real(RK), intent(in) :: datum ! Required date
      real(RK), intent(inout) :: A_cur(1:nval), A_s(1:nval), A_e(1:nval) ! Output values, Start and end values

      ! Local variables
      real(RK) :: tb_start, tb_end !Start and end time
      integer :: eof, i

      save tb_start, tb_end
      save eof

      if (idx) then

         open (20, status='old', file=self%file)
         eof = 0
         ! Read first values
         read (20, *, end=9)
         read (20, *, end=9) tb_start, (A_s(i), i=1, nval)
         if (datum < tb_start) then
            write(*,*) '[WARNING] ','First forcing date after simulation start time. datum=', datum,  'start=', tb_start
         end if
         read (20, *, end=7) tb_end, (A_e(i), i=1, nval)

         call ok("Forcing input file successfully read")
      end if

      if (datum <= tb_start .or. eof == 1) then ! If datum before first date or end of file reached
         goto 8
      else
         do while (datum > tb_end) ! Move to appropriate interval to get correct values
            tb_start = tb_end
            A_s(1:nval) = A_e(1:nval)
            ! Read next value
            read (20, *, end=7) tb_end, (A_e(i), i=1, nval)
         end do
         ! Linearly interpolate values at correct datum
         A_cur(1:nval) = A_s(1:nval) + (datum - tb_start)*(A_e(1:nval) - A_s(1:nval))/(tb_end - tb_start)
      end if
      return

  7   eof = 1
      if(datum>tb_start) call warn('Last forcing date before simulation end time.')
      write(6,*) datum, tb_start


  8   A_cur(1:nval) = A_s(1:nval)       ! Take first value of current interval
      return

  9   call error('Unable to read forcing file (no data found).')

   end subroutine

   ! Compute appropriate forcing parameters at given datum
   ! Copied from old simstrat, NEEDS refactoring!
   ! AG 2014: revision
   !######################################################################
   subroutine forcing_update(self, state)
      !######################################################################
      implicit none
      class(ForcingModule) :: self
      class(ModelState) :: state

      ! Local iables
      integer :: nval_offset
      real(RK) :: tau
      real(RK) :: A_s(8), A_e(8), A_cur(8) ! 8 is the maximum number of rows in the forcing file
      real(RK) :: fu, Vap_wat, heat0, emissivity
      real(RK) :: T_surf, F_glob, Vap_atm, Cloud
      real(RK) :: H_A, H_K, H_V, H_W, H_tot
      real(RK) :: F_snow, F_snowice, F_ice
      real(RK) :: qa, q0, sh_c
      save A_s, A_e
      associate (cfg=>self%cfg, param=>self%param)

      if((state%snow_h + state%white_ice_h + state%black_ice_h) > 0) then ! Ice cover
         T_surf = param%Freez_Temp
      else
         T_surf = state%T(self%grid%ubnd_vol) ! Free water
      end if

      ! Number of values to read, depending on filtered wind and precipitation
      if (cfg%use_filtered_wind .and. cfg%ice_model == 0) then
         nval_offset = 1
      else if (cfg%ice_model == 1 .and. cfg%use_filtered_wind) then
         nval_offset = 2
      else if (cfg%ice_model == 1) then
         nval_offset = 1
      else
         nval_offset = 0
      end if

      ! Forcing mode 1 (date, U10, V10, T_lake, H_sol)
      if (cfg%forcing_mode == 1) then
         if (cfg%ice_model == 1) then
            call error('Ice module not compatible with forcing mode 1, use 2, 3 or 5.')
         end if

         call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)
         call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)
         state%u10 = A_cur(1)*param%f_wind ! MS 2014: added f_wind
         state%v10 = A_cur(2)*param%f_wind ! MS 2014: added f_wind
         state%uv10 = sqrt(state%u10**2 + state%v10**2) ! AG 2014
         state%SST = A_cur(3) ! Lake surface temperature
         state%rad0 = A_cur(4)*(1 - state%albedo_water)*(1 - param%beta_sol) ! MS: added beta_sol and albedo_water
         state%heat = 0.0_RK
         state%T_atm = 0.0_RK
         state%precip = 0.0_RK
         if (cfg%use_filtered_wind) then
            state%Wf = A_cur(5) ! AG 2014
         end if

      ! Forcing mode 2 (date, U10, V10, T_atm, H_sol, Vap)
      else if (cfg%forcing_mode >= 2) then
         if (cfg%forcing_mode == 2) then
            call self%read (state%datum, A_s, A_e, A_cur, 5 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind ! MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind ! MS 2014: added f_wind
            state%T_atm = A_cur(3)

            if (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) then ! Ice
               F_glob = (A_cur(4)*(1 - ice_albedo)) * param%p_albedo
            else if (state%white_ice_h > 0 .and. state%snow_h == 0) then ! Snowice
               F_glob = (A_cur(4)*(1 - snowice_albedo)) * param%p_albedo
            else if (state%snow_h > 0) then ! Snow
               F_glob = (A_cur(4)*(1 - snow_albedo)) * param%p_albedo
            else ! Water
               F_glob = A_cur(4)*(1 - state%albedo_water) * param%p_sw
            end if

            Vap_atm = A_cur(5)
            Cloud = 0.5
            if (cfg%use_filtered_wind) state%Wf = A_cur(6) ! AG 2014
            if (cfg%snow_model == 1 .and. cfg%use_filtered_wind) then
               state%precip = A_cur(7)
            else if (cfg%snow_model == 1) then
               state%precip = A_cur(6)
            end if

         ! Forcing mode 3 (date, U10, V10, T_atm, H_sol, Vap, Clouds)
         else if (cfg%forcing_mode == 3) then
            call self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind ! MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind ! MS 2014: added f_wind
            state%T_atm = A_cur(3)

            if (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) then ! Ice
               F_glob = (A_cur(4)*(1 - ice_albedo)) * param%p_albedo
            else if (state%white_ice_h > 0 .and. state%snow_h == 0) then ! Snowice
               F_glob = (A_cur(4)*(1 - snowice_albedo)) * param%p_albedo
            else if (state%snow_h > 0) then ! Snow
               F_glob = (A_cur(4)*(1 - snow_albedo)) * param%p_albedo
            else ! Water
               F_glob = A_cur(4)*(1 - state%albedo_water) * param%p_sw
            end if

            Vap_atm = A_cur(5)
            Cloud = A_cur(6)
            if (Cloud < -1e-6 .or. Cloud > 1.000001) then
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

         ! Forcing mode 4 (date, U10, V10, H_net, H_sol)
         else if (cfg%forcing_mode == 4) then
            if (cfg%ice_model == 1) then
               call error('Ice module not compatible with forcing mode 4, use 2, 3 or 5.')
            end if

            call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind ! MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind ! MS 2014: added f_wind
            heat0 = A_cur(3) ! MS 2014
            F_glob = A_cur(4)*(1 - state%albedo_water) * param%p_sw
            state%T_atm = 0.0_RK
            if (cfg%use_filtered_wind) state%Wf = A_cur(5) ! AG 2014

         ! Forcing 5 (date, U10, V10, T_atm, H_sol, Vap, ILWR)
         else if (cfg%forcing_mode == 5) then
            call self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%first_timestep)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            state%T_atm = A_cur(3)

            if (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) then ! Ice
               F_glob = (A_cur(4)*(1 - ice_albedo)) * param%p_albedo
            else if (state%white_ice_h > 0 .and. state%snow_h == 0) then ! Snowice
               F_glob = (A_cur(4)*(1 - snowice_albedo)) * param%p_albedo
            else if (state%snow_h > 0) then ! Snow
               F_glob = (A_cur(4)*(1 - snow_albedo)) * param%p_albedo
            else ! Water
               F_glob = A_cur(4)*(1 - state%albedo_water) * param%p_sw
            end if

            Vap_atm = A_cur(5)
            H_A = A_cur(6)
            if (cfg%use_filtered_wind) state%Wf = A_cur(7) ! AG 2014
            if (cfg%snow_model == 1 .and. cfg%use_filtered_wind) then
               state%precip = A_cur(8)
            else if (cfg%snow_model == 1) then
               state%precip = A_cur(7)
            end if
         else
            call error('Wrong forcing type (must be 1, 2, 3, 4 or 5).')
         end if
         state%uv10 = sqrt(state%u10**2 + state%v10**2) ! AG 2014

         if (cfg%forcing_mode /= 4) then ! Heat fluxes calculations (forcing 2, 3 and 5)
            if (state%black_ice_h + state%white_ice_h == 0) then !Free water

               ! Factor 0.6072 to account for changing wind height from 10 to 2 m
               ! Further evaluation of evaporation algorithm may be required.
               fu = sqrt((2.7_RK*max(0.0_RK,(T_surf-state%T_atm)/(1 - 0.378_RK*Vap_atm/param%p_air))**0.333_RK)**2 + (0.6072_RK*3.1_RK*state%uv10)**2)
               fu = fu*param%p_windf ! Provided fitting factor p_windf (~1)

               ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
               Vap_wat = 10**((0.7859_RK + 0.03477_RK*T_surf)/(1 + 0.00412_RK*T_surf))
               Vap_wat = Vap_wat*(1 + 1e-6_RK*param%p_air*(4.5_RK + 0.00006_RK*T_surf**2))


               ! Long-wave radiation according to Dilley and O'Brien
               ! see Flerchinger et al. (2009)
               if (cfg%forcing_mode /= 5) then
                  H_A = (1 - r_a)*((1 - 0.84_RK*Cloud)*(59.38_RK + 113.7_RK*((state%T_atm + 273.15_RK)/273.16_RK)**6 &
                     + 96.96_RK*sqrt(465*Vap_atm/(state%T_atm + 273.15_RK)*0.04_RK))/5.67e-8_RK/ &
                     (state%T_atm + 273.15_RK)**4 + 0.84_RK*Cloud)*5.67e-8_RK*(state%T_atm + 273.15_RK)**4
               end if

               H_A = H_A*param%p_lw ! Provided fitting factor p_radin (~1)

               H_W = -emiss_water*sig*(T_surf + 273.15_RK)**4

               ! Flux of sensible heat (convection)
               H_K = -B0*fu*(T_surf - state%T_atm)

               ! Flux of latent heat (evaporation, condensation)
               H_V = -fu*(Vap_wat - Vap_atm)

               ! Global heat flux (positive: air to water, negative: water to air)
               state%heat = H_A + H_W + H_K + H_V + F_glob * param%beta_sol !MS: added term with beta_sol
               ! Removal of solar short-wave radiation absorbed in first water cell
               state%rad0 = F_glob * (1 - param%beta_sol) !MS: added beta_sol

               state%heat_snow = 0 ! Heat snow
               state%heat_snowice = 0 ! Heat snowice
               state%heat_ice = 0 ! Heat ice

            else if (state%black_ice_h > 0 .or. state%white_ice_h > 0) then ! Ice Cover
               ! Light penetration in snow, ice and snowice as well as wind blocking added 2018 by Love Raaman

               if (state%T_atm >= param%Freez_Temp) then! Melting occures when air temp above freezing point,
                  ! then activate surface heat fluxes following
                  ! Matti Leppäranta (2009), Modelling the Formation and Decay of Lake Ice DOI: 10.1007/978-90-481-2945-4_5 In book: The Impact of Climate Change on European Lakes
                  ! In G. George (Ed.), The Impact of Climate Change on European Lakes (pp. 63–83).
                  ! Dordrecht: Springer Netherlands. https://doi.org/10.1007/978-90-481-2945-4_5
                  ! and with corrections in
                  ! Leppäranta, M. (2014). Freezing of lakes and the evolution of their ice cover.
                  ! New York: Springer. ISBN 978-3-642-29080-0
                  if (state%snow_h == 0) then ! Ice Cover (ice and snowice)
                     emissivity = emiss_ice
                  else ! Snow Cover
                     ! Varies from 0.8 to 0.9 depending on snow density
                     emissivity = 5.0e-4_RK * state%snow_dens + 6.75e-1_RK
                  end if
                  ! obs fitting factors param%p_lw and param%p_windf not applied to ice covered lake
                  if (cfg%forcing_mode /= 5) then
                     H_A = (Ha_a + Ha_b * (Vap_atm**(1.0_RK/2.0_RK))) * (1 + Ha_c * Cloud**2) * sig * (state%T_atm + 273.15_RK)**4
                  end if
                  H_W = -emissivity * sig * (T_surf + 273.15_RK)**4
                  H_K = rho_air * cp_air * Hk_CH * (state%T_atm - T_surf) * state%uv10
                  sh_c = 0.622/param%p_air! Converter from absolut vapour pressure to specific humidity (Lepparanta 2015)
                  qa  = sh_c * Vap_atm! Specific humidity air
                  q0  = sh_c * 6.11! Specific humidity for saturation levels (6.11 mbar) at surface (ice) at 0°C (Lepparanta 2015)
                  H_V = rho_air * (l_h + l_e) * Hv_CE * (qa - q0) * state%uv10 ! Through sublimation (solid to gas)
               else ! If no melting only considering penetraiting solar radiation since ice formation is parameterised towards temperature
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
               state%heat = F_glob * exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h) * param%beta_sol
               state%rad0 = F_glob * exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h) * (1 - param%beta_sol)

               ! Heat flux into snow, ice or snowice layer.
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
            state%heat = heat0 + F_glob*param%beta_sol !MS: added term with beta_sol
            state%heat_snow = 0 ! Heat snow
            state%heat_snowice = 0 ! Heat snowice
            state%heat_ice = 0 ! Heat ice
         end if
         if (cfg%ice_model == 0) then
            if (state%T(self%grid%ubnd_vol) < 0 .and. state%heat < 0) then
               state%heat = 0
               state%T(self%grid%ubnd_vol) = 0
            end if
         end if
      end if

      ! Drag coefficient as a function of wind speed (AG 2014)
      if (cfg%wind_drag_model == 1) then ! Constant wind drag coefficient
         state%C10 = param%C10_constant
      else if (cfg%wind_drag_model == 2) then ! Ocean model
         state%C10 = param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
      else if (cfg%wind_drag_model == 3) then ! Lake model (Wüest and Lorke 2003)
         if (state%uv10 <= 0.1) then
            state%C10 = param%C10_constant*0.06215_RK
         else if (state%uv10 <= 3.85_RK) then
              state%C10 = param%C10_constant*0.0044_RK*state%uv10**(-1.15_RK)
         else ! Polynomial approximation of Charnock's law
            state%C10 = param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
         end if
      end if

      tau = state%C10*rho_air/rho_0*state%uv10**2
      state%u_taus = sqrt(tau)

      state%tx = state%C10*rho_air/rho_0*state%uv10*state%u10
      state%ty = state%C10*rho_air/rho_0*state%uv10*state%v10
      return
      end associate
   end subroutine

   ! Coriolis forces
   subroutine forcing_update_coriolis(self, state)
      implicit none
      class(ForcingModule) :: self
      class(ModelState) :: state
      real(RK) :: cori
      real(RK), dimension(size(state%U)) :: u_temp
      associate (grid=>self%grid, dt=>state%dt, param=>self%param)

         ! Calculate u_taub before changing U resp V
         state%u_taub = sqrt(state%drag*(state%U(1)**2 + state%V(1)**2))

         ! Calculate coriolis parameter based on latitude
         cori = 2.0_RK*7.292e-5_RK*sin(param%Lat*pi/180.0_RK)

         ! Update state based on coriolis parameter
         u_temp = state%U

         state%U(1:grid%ubnd_vol) = state%U(1:grid%ubnd_vol)*cos(Cori*dt) + state%V(1:grid%ubnd_vol)*sin(Cori*dt)
         state%V(1:grid%ubnd_vol) = -u_temp(1:grid%ubnd_vol)*sin(Cori*dt) + state%V(1:grid%ubnd_vol)*cos(Cori*dt)

         return
      end associate
   end subroutine

   subroutine forcing_init_albedo(self, state, sim_cfg)
      implicit none
      class(ForcingModule) :: self
      class(ModelState) :: state
      class(SimConfig) :: sim_cfg

      call init_calendar(sim_cfg%start_year, state%datum, state%current_year, state%current_month, state%current_day)

      ! Monthly albedo data according to Grishchenko, in Cogley 1979
      if (self%param%Lat > 0)  then ! Northern hemisphere (1.000 if no sun)
         state%albedo_data(9,1:12) = [1.000_RK, 0.301_RK, 0.333_RK, 0.253_RK, 0.167_RK, 0.133_RK, 0.150_RK, 0.226_RK, 0.317_RK, 0.301_RK, 1.000_RK, 1.000_RK]
         state%albedo_data(8,1:12) = [0.301_RK, 0.337_RK, 0.266_RK, 0.178_RK, 0.138_RK, 0.123_RK, 0.132_RK, 0.163_RK, 0.238_RK, 0.329_RK, 0.301_RK, 1.000_RK]
         state%albedo_data(7,1:12) = [0.340_RK, 0.281_RK, 0.185_RK, 0.122_RK, 0.099_RK, 0.095_RK, 0.097_RK, 0.113_RK, 0.163_RK, 0.254_RK, 0.336_RK, 0.325_RK]
         state%albedo_data(6,1:12) = [0.263_RK, 0.193_RK, 0.127_RK, 0.093_RK, 0.080_RK, 0.077_RK, 0.079_RK, 0.088_RK, 0.114_RK, 0.174_RK, 0.249_RK, 0.294_RK]
         state%albedo_data(5,1:12) = [0.178_RK, 0.131_RK, 0.095_RK, 0.077_RK, 0.071_RK, 0.070_RK, 0.070_RK, 0.075_RK, 0.088_RK, 0.120_RK, 0.169_RK, 0.198_RK]
         state%albedo_data(4,1:12) = [0.121_RK, 0.097_RK, 0.078_RK, 0.069_RK, 0.066_RK, 0.065_RK, 0.066_RK, 0.068_RK, 0.075_RK, 0.091_RK, 0.116_RK, 0.132_RK]
         state%albedo_data(3,1:12) = [0.091_RK, 0.079_RK, 0.070_RK, 0.065_RK, 0.064_RK, 0.064_RK, 0.064_RK, 0.064_RK, 0.068_RK, 0.076_RK, 0.089_RK, 0.097_RK]
         state%albedo_data(2,1:12) = [0.076_RK, 0.070_RK, 0.065_RK, 0.063_RK, 0.063_RK, 0.064_RK, 0.063_RK, 0.063_RK, 0.064_RK, 0.069_RK, 0.075_RK, 0.079_RK]
         state%albedo_data(1,1:12) = [0.069_RK, 0.065_RK, 0.063_RK, 0.063_RK, 0.065_RK, 0.066_RK, 0.065_RK, 0.063_RK, 0.063_RK, 0.065_RK, 0.068_RK, 0.070_RK]

         if (self%param%Lat < 10) then
           state%lat_number = 1
         else if (self%param%Lat < 20) then
           state%lat_number = 2
         else if (self%param%Lat < 30) then
           state%lat_number = 3
         else if (self%param%Lat < 40) then
           state%lat_number = 4
         else if (self%param%Lat < 50) then
           state%lat_number = 5
         else if (self%param%Lat < 60) then
           state%lat_number = 6
         else if (self%param%Lat < 60) then
           state%lat_number = 6
         else if (self%param%Lat < 70) then
           state%lat_number = 7
         else if (self%param%Lat < 80) then
           state%lat_number = 8
         else
           state%lat_number = 9
         end if

      else  ! Southern hemisphere (1.000 if no sun)
         state%albedo_data(9,1:12) = [0.150_RK, 0.226_RK, 0.317_RK, 0.301_RK, 1.000_RK, 1.000_RK, 1.000_RK, 0.301_RK, 0.333_RK, 0.253_RK, 0.167_RK, 0.133_RK]
         state%albedo_data(8,1:12) = [0.132_RK, 0.163_RK, 0.238_RK, 0.329_RK, 0.301_RK, 1.000_RK, 0.301_RK, 0.337_RK, 0.266_RK, 0.178_RK, 0.138_RK, 0.123_RK]
         state%albedo_data(7,1:12) = [0.097_RK, 0.113_RK, 0.163_RK, 0.254_RK, 0.336_RK, 0.325_RK, 0.340_RK, 0.281_RK, 0.185_RK, 0.122_RK, 0.099_RK, 0.095_RK]
         state%albedo_data(6,1:12) = [0.079_RK, 0.088_RK, 0.114_RK, 0.174_RK, 0.249_RK, 0.294_RK, 0.263_RK, 0.193_RK, 0.127_RK, 0.093_RK, 0.080_RK, 0.077_RK]
         state%albedo_data(5,1:12) = [0.070_RK, 0.075_RK, 0.088_RK, 0.120_RK, 0.169_RK, 0.198_RK, 0.178_RK, 0.131_RK, 0.095_RK, 0.077_RK, 0.071_RK, 0.070_RK]
         state%albedo_data(4,1:12) = [0.066_RK, 0.068_RK, 0.075_RK, 0.091_RK, 0.116_RK, 0.132_RK, 0.121_RK, 0.097_RK, 0.078_RK, 0.069_RK, 0.066_RK, 0.065_RK]
         state%albedo_data(3,1:12) = [0.064_RK, 0.064_RK, 0.068_RK, 0.076_RK, 0.089_RK, 0.097_RK, 0.091_RK, 0.079_RK, 0.070_RK, 0.065_RK, 0.064_RK, 0.064_RK]
         state%albedo_data(2,1:12) = [0.063_RK, 0.063_RK, 0.064_RK, 0.069_RK, 0.075_RK, 0.079_RK, 0.076_RK, 0.070_RK, 0.065_RK, 0.063_RK, 0.063_RK, 0.064_RK]
         state%albedo_data(1,1:12) = [0.065_RK, 0.063_RK, 0.063_RK, 0.065_RK, 0.068_RK, 0.070_RK, 0.069_RK, 0.065_RK, 0.063_RK, 0.063_RK, 0.065_RK, 0.066_RK]

         if (self%param%Lat > -10) then
           state%lat_number = 1
         else if (self%param%Lat > -20) then
           state%lat_number = 2
         else if (self%param%Lat > -30) then
           state%lat_number = 3
         else if (self%param%Lat > -40) then
           state%lat_number = 4
         else if (self%param%Lat > -50) then
           state%lat_number = 5
         else if (self%param%Lat > -60) then
           state%lat_number = 6
         else if (self%param%Lat > -60) then
           state%lat_number = 6
         else if (self%param%Lat > -70) then
           state%lat_number = 7
         else if (self%param%Lat > -80) then
           state%lat_number = 8
         else
           state%lat_number = 9
         end if
      end if

   end subroutine

   subroutine forcing_update_albedo(self, state)
      implicit none
      class(ForcingModule) :: self
      class(ModelState) :: state

      ! Local variables
      real(RK) :: albedo_start, albedo_end
      integer :: previous_month, next_month

      ! Determine current calendar day and month
      call update_calendar(state%current_year, state%current_month, state%current_day, state%dt)

      ! Determine albedo as a function of latitude and month according to Grishchenko (in Cogley 1979)
      ! Linear interpolation, assuming that the monthly data is representative of the 15. of each month (also for months with 28, 29 and 31 days)
      ! There is still a small jump in the albedo value when the month changes. Solving this problem would require a more sophisticated interpolation
      ! or even a continuous function. But the jump is very small and occurs during the night, when anyways no sun light is reaching the lake.

      ! If date before 16th of current month
      if (state%current_day < 16) then
         if(state%current_month == 1) then
            previous_month = 12
         else
            previous_month = state%current_month - 1
         end if

         ! Get albedo values for previous and current month
         albedo_start = state%albedo_data(state%lat_number, previous_month)
         albedo_end = state%albedo_data(state%lat_number, state%current_month)

         ! Linear interpolation
         state%albedo_water = albedo_start + (state%current_day + 15) * (albedo_end - albedo_start)/30
      else
         if(state%current_month == 12) then
            next_month = 1
         else
            next_month = state%current_month + 1
         end if

         ! Get albedo values for current and next month
         albedo_start = state%albedo_data(state%lat_number, state%current_month)
         albedo_end = state%albedo_data(state%lat_number, next_month)

         ! Linear interpolation
         state%albedo_water = albedo_start + (state%current_day - 15) * (albedo_end - albedo_start)/30
      end if

   end subroutine

end module strat_forcing

