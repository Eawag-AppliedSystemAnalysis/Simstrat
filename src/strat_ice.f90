!------------!
! Ice module
! Following:
! ------ !
! Saloranta, T. M., & Andersen, T. (2007).
! MyLake—A multi-year lake simulation model code suitable for uncertainty and sensitivity analysis simulations.
! Ecological Modelling, 207(1), 45–60.
! ------ !
! Yen (1981)
! Review of thermal properties of snow, ice and sea ice
! ------ !
! Saloranta, T. M. (2000).
! Modeling the evolution of snow, snow ice and ice in the Baltic Sea. Tellus A, 52(1), 93–108.
!--------------------------------------------------
! Ice module constructed in 2018 by:
! Dr. Love Råman Vinnå
! Department Surface Waters Research & Management
! Eawag, Seestrase 79, 6047 Kastanienbaum, Switzerland
module strat_ice
   use strat_forcing
   use strat_consts
   use strat_simdata
   use strat_grid
   implicit none
   private

   ! Common Types
   type, public :: IceModule
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param
      class(StaggeredGrid), pointer :: grid

   contains
      ! To initialize and run the ice module with subroutines
      procedure, pass :: init => ice_module_init
      procedure, pass :: update => ice_module_update

      ! The model itself
      procedure, pass :: do_ice_freezing        => ice_formation
      procedure, pass :: do_ice_melting         => ice_melting
      procedure, pass :: do_snow_melting        => snow_melting
      procedure, pass :: do_snow_build          => snow_build
      procedure, pass :: do_underneath_melting  => underneath_melting

   end type

contains

   ! Initiate ice model
   subroutine ice_module_init(self, model_cfg, model_param, grid)
      implicit none
      class(IceModule) :: self
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param
      class(StaggeredGrid), target :: grid

      self%model_cfg => model_cfg
      self%model_param => model_param
      self%grid => grid

   end subroutine

   ! Subroutine ice_module_update
   subroutine ice_module_update(self, state, param)

      implicit none
      class(IceModule)  :: self
      class(ModelState) :: state
      class(ModelParam) :: param

      !-------------------
      ! Below freezing
      !-------------------
      if (param%Freez_Temp >= state%T(self%grid%ubnd_vol) .and. (state%black_ice_h + state%white_ice_h) == 0) then
         ! Ice expanding (air temp & water temp less then Freez_Temp) no ice present
         call self%do_ice_freezing(state, param)
      else if ((state%black_ice_h + state%white_ice_h) > 0 .and. param%Freez_Temp >= state%T_atm) then
         ! If ice exist and air temp < freez temp, initiate ice formation
         call self%do_ice_freezing(state, param)
      end if

      ! Snow fall addition onto ice
      if (self%model_cfg%snow_model == 1 .and. param%snow_temp >= state%T_atm .and. (state%black_ice_h + state%white_ice_h) > 0 .and. state%precip > 0) then
         call self%do_snow_build(state)
      end if

      !------------------
      ! Above freezing
      !-------------------
      if (param%Freez_Temp < state%T_atm .and. (state%black_ice_h + state%white_ice_h) > 0) then
         ! Melt snow
         if (state%snow_h > 0 .and. param%snow_temp < state%T_atm .and. self%model_cfg%snow_model == 1) then
            call self%do_snow_melting(state)
         end if
         ! Melt ice from above
         if (state%white_ice_h + state%black_ice_h > 0) then
            call self%do_ice_melting(state, param)
         end if
         ! Melt ice from underneath
         if (state%T(self%grid%ubnd_vol) > param%freez_temp .and. (state%white_ice_h + state%black_ice_h) > 0) then
            call self%do_underneath_melting(state, param)
         end if
      end if

      ! Update total ice thickness
      state%total_ice_h = state%white_ice_h + state%black_ice_h

   end subroutine

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! The Ice Model
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! %%%%%%%%%%%%%%%%%%%%
   ! Below freezing point
   ! %%%%%%%%%%%%%%%%%%%%

   ! Ice/snowice formation and growth
   subroutine ice_formation(self, state, param)
      implicit none
      class(IceModule), intent(inout)  :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param

      ! Define variables only used in ice_model
      real(RK) :: P
      real(RK) :: P_ice
      real(RK) :: P_snow
      real(RK) :: Freez_energy
      real(RK) :: snow_weight
      real(RK) :: buoyancy_ice
      real(RK) :: snow_height_ice_mass

      if (state%black_ice_h == 0) then
         ! Energy released for first ice forming
         Freez_energy = (param%Freez_Temp - state%T(self%grid%ubnd_vol))  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ![J/kg/K]*[K]*[kg] = [J]

         ! First ice thickness
         state%black_ice_h = Freez_energy / l_h / (ice_dens*1*1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]

         ! Set surface temperature to freezing point
         state%T(self%grid%ubnd_vol) = param%Freez_Temp ![°C]

         ! Ice temp eq 17 in Saloranta et al. (2007).
         state%ice_temp = 0 ![°C]

      else   ! When ice cover exists
         ! Snow and ice insulating effect eq 17 and 18 Saloranta et al. (2007).
         P_ice = 1 / (10 * state%black_ice_h)

         if (self%model_cfg%snow_model == 1) then
            P_snow = (k_ice * state%snow_h)/(k_snow * state%black_ice_h)
         else
            P_snow = 0
         end if

         P = max(P_snow,P_ice)

         ! Ice temp eq 17 in Saloranta et al. (2007).
         state%ice_temp = (P*param%freez_temp + state%T_atm)/(1 + P)

         ! Snow-Ice formation, if weight of snow exceeds ice buoyancy
         if (self%model_cfg%snow_model == 1 .and. state%snow_h > 0) then
            snow_weight = state%snow_h * state%snow_dens*1*1 ! kg
            buoyancy_ice = state%black_ice_h*1*1 * (rho_0 - ice_dens) + state%white_ice_h*1*1 * (rho_0 - snowice_dens)!kg

            if (snow_weight > buoyancy_ice) then
               snow_height_ice_mass = snow_weight - buoyancy_ice
               state%snow_h = state%snow_h - snow_height_ice_mass/state%snow_dens
               state%white_ice_h  = state%white_ice_h  + snow_height_ice_mass/state%snow_dens ! Assuming water from below fills up zone needed to be flooded to achieve buoyant stability
            end if
         end if
         ! Ice thickness eq. 16 in Saloranta et al. (2007).
         state%black_ice_h = sqrt(state%black_ice_h**2 + (2*k_ice)/(ice_dens * l_h) * (param%freez_temp - state%ice_temp) * state%dt) 

         ! Set surface temperature to freezing point
         state%T(self%grid%ubnd_vol) = param%freez_temp ![°C]
      end if

      ! If melt larger than ice height, put remaining energy to water and set ice height to zero
      if (state%black_ice_h < 0) then
         state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ![J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%black_ice_h = 0
      end if

   end subroutine

   ! Snow layer build-up
   subroutine snow_build(self, state)
      implicit none
      class(IceModule) :: self
      class(ModelState) :: state
      real(RK) :: snow_h_new
      real(RK) :: Ws, ChangSnowDens

      ! Calculate new snow height
      snow_h_new = (state%precip / 3600) * state%dt * (rho_0 / rho_s_0) ! Go from m/h to m/s and increase volume from water to snow

      ! Calculate snow density due to compression of snow layer, Yen (1981)
      Ws = (snow_h_new * rho_s_0) / (rho_0 * 1 * 1)  ! [m water equivalent]
      ! Compression Yen (1981) eq. 7
      ChangSnowDens = state%snow_dens * C01 * Ws * exp(-C02 * state%snow_dens) * state%dt
      ! Compress old snow layer
      state%snow_h = state%snow_h * (state%snow_dens /(state%snow_dens + ChangSnowDens))

      ! Adjust density
      state%snow_dens = state%snow_dens + ChangSnowDens

      ! Combine the old and the new snow layer
      ! Create new density from old and new density
      state%snow_dens  = (state%snow_dens * state%snow_h + rho_s_0 * snow_h_new) / (state%snow_h + snow_h_new)
      ! Update snow height
      state%snow_h = state%snow_h + snow_h_new

      ! Maximum allowed snow density
      if (state%snow_dens > rho_s_max) then
         state%snow_dens = rho_s_max
      end if

   end subroutine

   ! %%%%%%%%%%%%%%%%%%%%
   ! Above freezing point
   ! %%%%%%%%%%%%%%%%%%%%

   ! Snow melting from above through sublimation (l_e ~= 0) or non-sublimation (l_e = 0) (i.e. solid to gas (l_h + l_e))
   subroutine snow_melting(self, state)
      implicit none
      class(IceModule) :: self
      class(ModelState) :: state
      ! Define variables only used in snow_model
      real(RK) :: Melt_energy1
      real(RK) :: MeltHeight1
      real(RK) :: Melt_ice
      real(RK) :: MeltHeightIce

      ! Melt Snow from atmosphere
      Melt_energy1 = state%heat_snow * state%dt * 1 * 1 ! [W/m2] * [s] * [m2] = [J]
      MeltHeight1 = Melt_energy1 / (l_h + l_e) / (state%snow_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      if (MeltHeight1 < 0) then ! Do not add snow from above, only melt snow
         Melt_energy1 = 0
         MeltHeight1 = 0
      end if
      ! Change the snow height
      state%snow_h = state%snow_h - MeltHeight1! [m]

      ! If melting energy larger than required, put remaining energy to melting of snowice / ice and set snow height to zero
      if (state%snow_h < 0) then
         Melt_ice = (l_h + l_e) * state%snow_dens * (1 * 1) * (-1 * state%snow_h)! [J/kg]  [kg/m3]  [m2]  [m] = [J]
         state%snow_h = 0
         if (state%white_ice_h > 0) then
            MeltHeightIce = Melt_ice / (l_h + l_e) / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
            state%white_ice_h = state%white_ice_h - MeltHeightIce! [m]
            if (state%white_ice_h < 0) then
               Melt_ice = (l_h + l_e) * snowice_dens * (1 * 1) * (-1 * state%white_ice_h)! [J/kg]  [kg/m3]  [m2]  [m] = [J]
               state%white_ice_h = 0
               MeltHeightIce = Melt_ice / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
               state%black_ice_h = state%black_ice_h - MeltHeightIce ! [m]
               ! If melt larger than ice height, put remaining energy to water and set ice height to zero
               if (state%black_ice_h < 0) then
                  state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
                  state%black_ice_h = 0
               end if
            end if
         else if (state%white_ice_h == 0 .and. state%black_ice_h > 0 ) then
            MeltHeightIce = Melt_ice / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
            state%black_ice_h = state%black_ice_h - MeltHeightIce ! [m]
            ! If melt larger than ice height, put remaining energy to water and set ice height to zero
            if (state%black_ice_h < 0) then
               state%heat = state%heat + (l_h  * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
               state%black_ice_h = 0
            end if
         end if
      else
         Melt_ice = 0
      end if
   end subroutine

   ! Ice melting from above through sublimation (l_e ~= 0) or non-sublimation (l_e = 0) (i.e. solid to gas (l_h + l_e))
   subroutine ice_melting(self, state, param)
      implicit none
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param
      ! Define variables only used in ice_model
      real(RK) :: Melt_energy2, Melt_energy3
      real(RK) :: MeltHeight2, MeltHeight3
      real(RK) :: Melt_ice, MeltHeightIce


      ! Melt snowice from atmosphere
       Melt_energy2 = state%heat_snowice * state%dt * 1 * 1 ! [W/m2] * [s] * [m2] = [J]
      ! Melt ice from atmosphere
       Melt_energy3 = state%heat_ice * state%dt * 1 * 1 ! [W/m2] * [s] * [m2] = [J]

      ! Snowice
      if (state%snow_h == 0) then ! Free surface
         MeltHeight2 = Melt_energy2 / (l_h + l_e) / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      else if (state%snow_h > 0 .and. state%black_ice_h == 0) then ! Layer above and none below
         MeltHeight2 = Melt_energy2 / (l_h) / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      else ! Layer above and below, melt snow
         state%heat_snow = Melt_energy2 / state%dt
         call self%do_snow_melting(state)
         MeltHeight2 = 0
      end if
      if (MeltHeight2 < 0) then ! Do not add ice, only melt.
         Melt_energy2 = 0;
         MeltHeight2 = 0;
      end if
      ! Ice
      if (state%snow_h + state%white_ice_h == 0) then ! Free surface
         MeltHeight3 = Melt_energy3 / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      else if (state%snow_h + state%white_ice_h > 0) then ! Layer above and none below
         MeltHeight3 = Melt_energy3 / (l_h) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      end if
      if (MeltHeight3 < 0) then ! Do not add ice from above, only melt
         Melt_energy3 = 0;
         MeltHeight3 = 0;
      end if

      ! New snowice height
      state%white_ice_h = state%white_ice_h - MeltHeight2! [m]
      MeltHeightIce = 0
      ! If melting energy larger than required, put remaining energy to melting of ice and set snowice height to zero
      if (state%white_ice_h < 0 .and. state%black_ice_h > 0) then
         Melt_ice = (l_h + l_e) * snowice_dens * (1 * 1) * (-1 * state%white_ice_h) ! [J/kg]  [kg/m3] [m2] [m] = [J]
         state%white_ice_h = 0
         MeltHeightIce = Melt_ice / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      else if (state%white_ice_h < 0 .and. state%black_ice_h == 0) then ! If melt larger than snowice hight, put remaining energy to water and set snowice hight to zero
         state%heat = state%heat + (l_h * snowice_dens * (-1 * state%white_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%white_ice_h = 0
      end if

      ! New ice height
      state%black_ice_h = state%black_ice_h - MeltHeight3 - MeltHeightIce ! [m]
      ! If melt larger than ice height, put remaining energy to water and set ice height to zero
      if (state%black_ice_h < 0) then
         state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%black_ice_h = 0
      end if

      ! Set ice temp to freez point
      state%ice_temp = 0 ![°C]
   end subroutine

   ! Melt ice from below (l_h), while keeping freez temperature in the interface between ice and water
   subroutine underneath_melting(self, state, param)
      implicit none
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param
      ! Define variables only used in snow_model
      real(RK) :: Melt_energy2
      real(RK) :: MeltHeight2

      Melt_energy2 = (state%T(self%grid%ubnd_vol) - param%freez_temp)  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ! [J/kg/K]*[K]*[kg] = [J]
      if (state%black_ice_h > 0) then
         MeltHeight2 = Melt_energy2 / l_h / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      else if (state%black_ice_h == 0 .and. state%white_ice_h > 0) then
         MeltHeight2 = Melt_energy2 / l_h / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      end if
      if (MeltHeight2 < 0) then !Do not add ice from below, only melt
         Melt_energy2 = 0;
         MeltHeight2 = 0;
      end if

      ! Set surface water temperature to freezing point
      state%T(self%grid%ubnd_vol) = param%freez_temp ![°C]

      ! New ice height
      if (state%black_ice_h > 0) then
         state%black_ice_h = state%black_ice_h - MeltHeight2 ! [m]
      else if (state%black_ice_h == 0 .and. state%white_ice_h > 0) then
         state%white_ice_h = state%white_ice_h - MeltHeight2 ! [m]
      end if

      ! If melt larger than ice height, put remaining energy to water and set ice height to zero
      if (state%black_ice_h < 0) then
         state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%black_ice_h = 0
      end if
      if (state%white_ice_h < 0) then ! Open water
         state%heat = state%heat + l_h * (snowice_dens * (-1 * state%white_ice_h)) / state%dt ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%white_ice_h = 0
      end if
      if (state%white_ice_h == 0 .and. state%black_ice_h == 0 .and. state%snow_h > 0) then
         state%heat = state%heat + l_h * (state%snow_dens * (state%snow_h)) / state%dt ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%snow_h = 0
         state%snow_dens = rho_s_0
      end if
   end subroutine

end module

