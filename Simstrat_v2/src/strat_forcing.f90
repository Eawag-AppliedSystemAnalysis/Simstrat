!     +---------------------------------------------------------------+
!     |  Forcing module
!     |  - Reads forcing file and updates state
!     |  - EUpdates coriolis force terms
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
      integer :: nz_max
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
            write (6, *) 'Warning: first forcing date after simulation start time. datum=', datum, " start=", tb_start
         end if
         read (20, *, end=7) tb_end, (A_e(i), i=1, nval)

         write (6, *) "Forcing input file successfully read"
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
      if(datum>tb_start) write(6,*) 'Warning: last forcing date before simulation end time.'

  8   A_cur(1:nval) = A_s(1:nval)       !Take first value of current interval
      return

  9   write(6,*) 'Error reading forcing file (no data found).'

      stop

   end subroutine

   !Compute appropriate forcing parameters at given datum
   ! Copied from old simstrat, NEEDS refactoring!
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
      real(RK) :: A_s(7), A_e(7), A_cur(7)
      real(RK) :: fu, Vap_wat, heat0
      real(RK) :: T_surf, T_atm, F_glob, Vap_atm, Cloud
      real(RK) :: H_A, H_K, H_V, H_W
      save A_s, A_e
      associate (cfg=>self%cfg, param=>self%param)

         T_surf = state%T(self%grid%ubnd_vol)
         !Todo: is surface temp the uppermost volume's temp?

         ! number of values to read, depending on filtered wind
         if (cfg%use_filtered_wind) then
            nval_offset = 1
         else
            nval_offset = 0
         end if

         if (cfg%forcing_mode == 1) then
            call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%std)
            state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
            state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
            state%uv10 = sqrt(state%u10**2 + state%v10**2) !AG 2014
            state%SST = A_cur(3) !Sea surface temperature
            state%rad0 = A_cur(4)*(1 - param%albsw)*(1 - param%beta_sol) ! MS: added beta_sol and albsw
            state%heat = 0.0_RK
            if (cfg%use_filtered_wind) state%Wf = A_cur(5) !AG 2014
         else if (cfg%forcing_mode >= 2) then
            if (cfg%forcing_mode == 2) then ! date, U,V,Tatm,Hsol,Vap
               call self%read (state%datum, A_s, A_e, A_cur, 5 + nval_offset, state%std)
               state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
               T_atm = A_cur(3)
               F_glob = A_cur(4)*(1 - param%albsw)
               Vap_atm = A_cur(5)
               Cloud = 0.5_RK
               if (cfg%use_filtered_wind) state%Wf = A_cur(6) !AG 2014
            else if (cfg%forcing_mode == 3) then ! date,U10,V10,Tatm,Hsol,Vap,Clouds
               call self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%std)
               state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
               T_atm = A_cur(3)
               F_glob = A_cur(4)*(1 - param%albsw)
               Vap_atm = A_cur(5)
               Cloud = A_cur(6)
               if (Cloud < 0 .or. Cloud > 1) then
                  write (6, *) 'Cloudiness should always be between 0 and 1.'
                  stop
               end if
               if (cfg%use_filtered_wind) state%Wf = A_cur(7) !AG 2014
            else if (cfg%forcing_mode == 4) then ! date,U10,V10,Hnet,Hsol
               call self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%std)
               state%u10 = A_cur(1)*param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*param%f_wind !MS 2014: added f_wind
               heat0 = A_cur(3) !MS 2014
               F_glob = A_cur(4)*(1 - param%albsw)
               if (cfg%use_filtered_wind) state%Wf = A_cur(5) !AG 2014
            else
               write (6, *) 'Error: wrong forcing type (must be 1, 2, 3 or 4).'
               stop
            end if
            state%uv10 = sqrt(state%u10**2 + state%v10**2) !AG 2014

            if (cfg%forcing_mode /= 4) then ! in the water column

               ! Wind function (Adams et al., 1990), changed by MS, June 2016
               ! Factor 0.6072 to account for changing wind height from 10 to 2 m
               ! Further evaluation of evaporation algorithm may be required.
               fu = sqrt((2.7_RK*max(0.0_RK,(T_surf-T_atm)/(1-0.378_RK*Vap_atm/param%p_air))**0.333_RK)**2 + (0.6072_RK*3.1_RK*state%uv10)**2)
               ! Wind function (Livingstone & Imboden 1989)
               !fu = 4.40_RK + 1.82_RK*state%uv10 + 0.26_RK*(T_surf - T_atm)
               !fu = 5.44+2.19*wind+0.24*(T_surf-T_atm)
               fu = fu*param%p_windf ! Provided fitting factor p_windf (~1)
               ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
               Vap_wat = 10**((0.7859_RK + 0.03477_RK*T_surf)/(1 + 0.00412_RK*T_surf))
               Vap_wat = Vap_wat*(1 + 1e-6_RK*param%p_air*(4.5_RK + 0.00006_RK*T_surf**2))
               ! Solar short-wave radiation absorbed
               state%rad0 = F_glob*(1 - param%beta_sol) !MS: added beta_sol
               ! Long-wave radiation from sky (Livingstone & Imboden 1989)
               ! H_A = 1.24*sig*(1-r_a)*(1+0.17*Cloud**2)*(Vap_atm/(273.15+T_atm))**(1./7)*(273.15+T_atm)**4
               ! Long-wave radiation according to Dilley and O'Brien
               ! see Flerchinger et al. (2009)
               H_A = (1 - r_a)*((1 - 0.84_RK*Cloud)*(59.38_RK + 113.7_RK*((T_atm + 273.15_RK)/273.16_RK)**6&
               &   + 96.96_RK*sqrt(465*Vap_atm/(T_atm + 273.15_RK)*0.04_RK))/5.67e-8_RK/ &
               &   (T_atm + 273.15_RK)**4 + 0.84_RK*Cloud)*5.67e-8_RK*(T_atm + 273.15_RK)**4
               H_A = H_A*param%p_radin ! Provided fitting factor p_radin (~1)

               ! Long-wave radiation from water body (black body)
               H_W = -0.97_RK*sig*(T_surf + 273.15_RK)**4
               ! Flux of sensible heat (convection)

               H_K = -B0*fu*(T_surf - T_atm)
               ! Flux of latent heat (evaporation, condensation)
               H_V = -fu*(Vap_wat - Vap_atm)
               ! Global heat flux (positive: air to water, negative: water to air)
               state%heat = H_A + H_W + H_K + H_V + F_glob*param%beta_sol !MS: added term with beta_sol
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
