!     +---------------------------------------------------------------+
!     |  Forcing module
!     |  - Adjusted heat fluxs to snow and ice covered lake
!     |  - Reads forcing file and updates state
!     |  - Updates coriolis force terms
!     +---------------------------------------------------------------+

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
   
      if(self%cfg%snow_model == 1) call warn('Snow module needs precipitation input, be sure the forcing file contains it (last column).') 
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
      integer, intent(in) :: idx, nval
      real(RK), intent(in) :: datum !Required date
      real(RK), intent(inout) :: A_cur(1:nval), A_s(1:nval), A_e(1:nval) !output values, Start and end values

      ! Local variables
      real(RK) :: tb_start, tb_end !Start and end time
      integer :: eof, i

      save tb_start, tb_end
      save eof

      if (idx == 1) then

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


  8   A_cur(1:nval) = A_s(1:nval)       !Take first value of current interval
      return

  9   call error('Unable to read forcing file (no data found).')

      stop

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
      real(RK) :: H_A, H_K, H_V, H_W
      save A_s, A_e
      associate (cfg=>self%cfg, param=>self%param)

      if(state%snow_h <= 0.000001 .and. state%ice_h >= 0.000001) then
         T_surf = state%ice_temp
      else if (state%snow_h > 0.000001 .and. state%ice_h >= 0.000001) then
         T_surf = state%snow_temp   
      else 
         T_surf = state%T(self%grid%ubnd_vol)
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
              call error('Ice module not compatible with forcing mode 1, use 2 or 3.')
              stop
            end if
   
            call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%model_step_counter)
            call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%model_step_counter)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            state%uv10 = sqrt(state%u10**2 + state%v10**2) !AG 2014
            state%SST = A_cur(3) !Sea surface temperature
            state%rad0 = A_cur(4)*(1 - param%albsw)*(1 - param%beta_sol) ! MS: added beta_sol and albsw
            state%heat = 0.0_RK
            state%T_atm = 0.0_RK
            state%precip = 0.0_RK
            if (cfg%use_filtered_wind) state%Wf = A_cur(5) !AG 2014
   
         else if (cfg%forcing_mode >= 2) then
            if (cfg%forcing_mode == 2) then ! date, U,V,Tatm,Hsol,Vap
               call self%read (state%datum, A_s, A_e, A_cur, 5 + nval_offset, state%model_step_counter)
               state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
               state%T_atm = A_cur(3)
               if (state%ice_h > 0 .and. state%snow_h == 0) then !Ice
               F_glob = A_cur(4)*(1 - param%ice_albedo)      
               else if (state%ice_h > 0 .and. state%snow_h > 0) then !Snow
               F_glob = A_cur(4)*(1 - param%snow_albedo)     
               else !Water
               F_glob = A_cur(4)*(1 - param%albsw) 
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
               call self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%model_step_counter)
               state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
               state%T_atm = A_cur(3)
               if (state%ice_h > 0 .and. state%snow_h == 0) then !Ice
               F_glob = A_cur(4)*(1 - param%ice_albedo)      
               else if (state%ice_h > 0 .and. state%snow_h > 0) then !Snow
               F_glob = A_cur(4)*(1 - param%snow_albedo)     
               else !Water
               F_glob = A_cur(4)*(1 - param%albsw) 
               end if     
               Vap_atm = A_cur(5)
               Cloud = A_cur(6) 
               if (Cloud < 0 .or. Cloud > 1) then
                  call error('Cloudiness should always be between 0 and 1.')
                  stop
               end if
               if (cfg%use_filtered_wind) state%Wf = A_cur(7) !AG 2014
               if (cfg%snow_model == 1 .and. cfg%use_filtered_wind) then
                state%precip = A_cur(8)
               else if (cfg%snow_model == 1) then
                state%precip = A_cur(7)     
               end if

            else if (cfg%forcing_mode == 4) then ! date,U10,V10,Hnet,Hsol
               if (cfg%ice_model == 1) then 
                 call error('Ice module not compatible with forcing mode 4, use 2 or 3.')
                 stop  
               end if
      
               call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%model_step_counter)
               state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
               heat0 = A_cur(3) !MS 2014
               F_glob = A_cur(4)*(1 - param%albsw)
               state%T_atm = 0.0_RK       
               if (cfg%use_filtered_wind) state%Wf = A_cur(5) !AG 2014      
      
            else
               call error('Wrong forcing type (must be 1, 2, 3 or 4).')
               stop
            end if
            state%uv10 = sqrt(state%u10**2 + state%v10**2) !AG 2014

            if (cfg%forcing_mode /= 4) then ! in the water column

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
               H_A = (1 - r_a)*((1 - 0.84_RK*Cloud)*(59.38_RK + 113.7_RK*((state%T_atm + 273.15_RK)/273.16_RK)**6&
               &   + 96.96_RK*sqrt(465*Vap_atm/(state%T_atm + 273.15_RK)*0.04_RK))/5.67e-8_RK/ &
               &   (state%T_atm + 273.15_RK)**4 + 0.84_RK*Cloud)*5.67e-8_RK*(state%T_atm + 273.15_RK)**4
               H_A = H_A*param%p_radin ! Provided fitting factor p_radin (~1)

               ! Long-wave radiation from water body (black body) 
               ! following:
               ! Leppäranta, M. (2010). Modelling the Formation and Decay of Lake Ice. 
               ! In G. George (Ed.), The Impact of Climate Change on European Lakes (pp. 63–83). 
               ! Dordrecht: Springer Netherlands. https://doi.org/10.1007/978-90-481-2945-4_5
               if (state%ice_h > 0 .and. state%snow_h == 0) then !Ice Cover  
               emissivity = emiss_ice
               else if (state%ice_h > 0 .and. state%snow_h > 0) then !Snow Cover  
               !varies from 0.8 to 0.9 depending on snow density      
               emissivity = 5.0e-4_RK * state%snow_dens + 6.75e-1_RK
               else ! Free Water
               emissivity = emiss_water
               end if     
               H_W = -emissivity*sig*(T_surf + 273.15_RK)**4
      
               ! Flux of sensible heat (convection)
               H_K = -B0*fu*(T_surf - state%T_atm)
      
               ! Flux of latent heat (evaporation, condensation)
               H_V = -fu*(Vap_wat - Vap_atm)
        
               ! Heat fluxes save
               ! Ice light penetration and wind blocking added 2018 by Love Raaman
               if (state%ice_h > 0 .and. state%snow_h == 0) then !Ice Cover
               ! Global heat flux (positive: air to water, negative: water to air)
               ! Heat enters first water cell, heat_snow_ice into the ice or snow layer  
                state%heat = 0 + 0 + 0 + 0 + F_glob * param%beta_sol * (1 - param%beta_snow_ice) !Heat first layer, LRV added beta_snow_ice (Bouffard. 2016,Ice-covered Lake Onega)
                state%heat_snow_ice = H_A + H_W + H_K + H_V + F_glob * param%beta_snow_ice !Heat snow and/or ice; LRV added beta_snow_ice
               ! Removal of solar short-wave radiation absorbed in snow and ice and first water cell
               state%rad0 = F_glob * (1 - param%beta_sol) * (1 - param%beta_snow_ice)!MS: added beta_sol ; LRV added beta_snow_ice      
               ! Supress wind turbulence with wind lid (heat flux affected by wind still active on snow and ice)
               state%u10 = 0
               state%v10 = 0
               state%uv10 = 0      
               else if (state%ice_h > 0 .and. state%snow_h > 0) then !Snow Cover
               ! Global heat flux (positive: air to water, negative: water to air) 
               ! Heat enters first water cell, heat_snow_ice into the ice or snow layer
               state%heat = 0 + 0 + 0 + 0 + F_glob * param%beta_sol * (1 - param%beta_snow_ice)  !Heat first layer, LRV added beta_snow_ice   
               state%heat_snow_ice = H_A + H_W + H_K + H_V + F_glob*param%beta_snow_ice !Heat snow and/or ice; LRV added beta_snow_ice
               ! Removal of solar short-wave radiation absorbed in snow and ice and first water cell
               state%rad0 = F_glob * (1 - param%beta_sol) * (1 - param%beta_snow_ice)!MS: added beta_sol ; LRV added beta_snow_ice (Bouffard. 2016,Ice-covered Lake Onega)
               ! Supress wind turbulence with wind lid (heat flux affected by wind still active on snow and ice)
               state%u10 = 0
               state%v10 = 0
               state%uv10 = 0
               else ! Free Water
               ! Global heat flux (positive: air to water, negative: water to air) 
               state%heat = H_A + H_W + H_K + H_V + F_glob * param%beta_sol !MS: added term with beta_sol      
               state%heat_snow_ice = 0 + 0 + 0 + 0 + 0   
               ! Removal of solar short-wave radiation absorbed in first water cell
               state%rad0 = F_glob * (1 - param%beta_sol) !MS: added beta_sol
               end if 
                ! save for output, not used in calculations
                state%ha = H_A
                state%hw = H_W
                state%hk = H_K
                state%hv = H_V
            else
               state%heat = heat0 + F_glob*param%beta_sol !MS: added term with beta_sol
            end if
            if ((T_surf < 0) .and. (state%heat < 0)) state%heat = 0.
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
