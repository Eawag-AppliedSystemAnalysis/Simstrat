!------------!
! Ice module
! Done according to:
! Saloranta, T. M., & Andersen, T. (2007). 
! MyLake—A multi-year lake simulation model code suitable for uncertainty and sensitivity analysis simulations.
! Ecological Modelling, 207(1), 45–60.
module strat_ice
   use strat_forcing
   use strat_consts
   use strat_simdata 
   use strat_grid
   implicit none
   private
   
   !Common Types
   type, public :: IceModule  
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param  
      class(StaggeredGrid), pointer :: grid  
   
   
   contains
 ! to initial and run the ice mudule with subrutins
      procedure, pass :: init => ice_module_init  
      procedure, pass :: update => ice_module_update
   
   ! the model itself
      procedure, pass :: do_ice_freezing => ice_formation
      procedure, pass :: do_ice_melting  => ice_melting 
      procedure, pass :: do_snow_melting => snow_melting 
      procedure, pass :: do_snow_build => snow_build 
   
   end type   
   contains
   
   
   ! Initiate ice model
 subroutine ice_module_init(self, model_cfg, model_param, grid)
   !only used in start up of model
      implicit none
      class(IceModule) :: self
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param
      class(StaggeredGrid), target :: grid  

     
        self%model_cfg => model_cfg
        self%model_param => model_param
        self%grid => grid


 end subroutine
   
   !subroutine ice_module_update(self, state, param, grid)
 subroutine ice_module_update(self, state, param)

   !this is used each time step
   !call futher subrutines here
      implicit none
      class(IceModule)  :: self
      class(ModelState) :: state
      class(ModelParam) :: param
      real(RK) :: buoyancy_ice
      real(RK) :: snow_weight   
   
   !-------------------
   ! Below freezing
   !-------------------
   if (param%Freez_Temp >= state%T(self%grid%ubnd_vol) .and. param%Freez_Temp >= state%T_atm) then
    !Ice expending (air temp & water temp less then Freez_Temp)
     call self%do_ice_freezing(state, param) 
   else if (state%ice_h > 0 .and. param%Freez_Temp >= state%T_atm) then 
     !If ice exist and temperature under ice is larger than freez temp. Set surface temp to freez point and iniate ice formation
     !set surface temperature to freezing point
     state%T(self%grid%ubnd_vol) = param%Freez_Temp ![°C]
     call self%do_ice_freezing(state, param)   
   end if
   !Snow fall addition onto ice     
   if (self%model_cfg%snow_model == 1 .and. param%Freez_Temp >= state%T_atm .and. state%ice_h > 0 .and. state%percip > 0) then
       call self%do_snow_build(state, param) 
   end if   
   
  !------------------ 
  ! Above freezing
  !-------------------
   if (param%Freez_Temp < state%T_atm .and. state%ice_h > 0) then 
     !Melt snow    
     if (state%snow_h > 0 .and. self%model_cfg%snow_model == 1) then
      call self%do_snow_melting(state, param)  
     !Melt ice after all snow is gone 
     else if (state%snow_h == 0) then
      call self%do_ice_melting(state, param)
     end if
   end if   

 end subroutine
   
   
   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The Ice Model

! Ice formation and growth
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
      real(RK) :: Melt_energy2
      real(RK) :: MeltHeight2
      real(RK) :: snow_weight
      real(RK) :: buoyancy_ice
      real(RK) :: snow_height_ice_mass

  !-------------------------------------------------
  !calculations eq 16 to 19 Saloranta et al. 2007
    
 !Ice forming  
  if (state%ice_h == 0) then
     ! energy released for first ice forming
     Freez_energy = (param%Freez_Temp - state%T(self%grid%ubnd_vol))  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ![J/kg/K]*[K]*[kg] = [J]

     !first ice thickness
     state%ice_h = Freez_energy / l_h / (ice_dens*1*1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]
   
     ! set surface temperature to freezing point
     state%T(self%grid%ubnd_vol) = param%Freez_Temp ![°C]
   
     !Ice temp eq 17
     state%ice_temp = 0 ![°C]
       
      !write (*, *) 'Ice Formation'
      !write (*,'(A,F8.6)') 'First cell height :', self%grid%h(self%grid%ubnd_vol) 
      !write (*,'(A,F8.6)') 'First ice depth :', state%ice_h
      !write (*,'(A,F8.6)') 'Air temp :', state%T_atm   
      !write (*, *) ''
     
  else   !When ice cover exist

    !Snow and ice insulating effect eq 17 and 18
    P_ice = 1 / (10 * state%ice_h)
   
    if (self%model_cfg%snow_model == 1) then
     P_snow = (k_ice * state%snow_h)/(k_snow * state%ice_h)
    else
     P_snow = 0
    end if
   
    P = max(P_snow,P_ice)
   
    !Ice temp eq 17
    state%ice_temp = (P*param%freez_temp + state%T_atm)/(1 + P)
   
   
     !snow ice formation, if wight of snow exeeds ice buoyancy
     if (self%model_cfg%snow_model == 1 .and. state%snow_h > 0) then
        snow_weight = state%snow_h * state%snow_dens*1*1 !kg
        buoyancy_ice = state%ice_h*1*1 * (rho_0 - ice_dens) !kg
        !write (*,'(A,F8.6)') 'Snow wheight :', snow_weight 
        !write (*,'(A,F8.6)') 'Buoyancy ice :', buoyancy_ice 

        if (snow_weight > buoyancy_ice) then
          snow_height_ice_mass = (snow_weight - buoyancy_ice)
  
          state%snow_h = state%snow_h - snow_height_ice_mass/state%snow_dens
          state%ice_h  = state%ice_h  + snow_height_ice_mass/ice_dens
     
          !write (*, *) 'Snow ice forming'
          !write (*,'(A,F8.6)') 'Snow ice hight       :', snow_height_ice_mass/ice_dens 
          !write (*,'(A,F8.6)') 'snow height lost     :', snow_height_ice_mass/state%snow_dens    
          !snow_weight = state%snow_h * state%snow_dens*1*1 !kg
          !buoyancy_ice = state%ice_h*1*1 * (rho_0 - ice_dens) !kg
          !write (*,'(A,F8.6)') 'ice hight new   :', state%ice_h 
          !write (*,'(A,F8.6)') 'snow height new :', state%snow_h 
          !write (*,'(A,F12.6)') 'Snow wheight :', snow_weight 
          !write (*,'(A,F12.6)') 'Buoyancy ice :', buoyancy_ice    
   
      end if   
     end if
   
   
    !Ice thickness eq. 16
    state%ice_h = sqrt(state%ice_h**2 + (2*k_ice)/(ice_dens * l_h) * (param%freez_temp - state%ice_temp) * state%dt) 

   !Heating of first water cell might melt ice from below, put back first cell to freez temp., and melt ice from below
    Melt_energy2 = (state%T(self%grid%ubnd_vol) - param%freez_temp)  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ![J/kg/K]*[K]*[kg] = [J] 
    MeltHeight2 = Melt_energy2 / l_h / (ice_dens * 1 * 1) ![J] / [J/kg] / [kg/m3] / [m2] = [m] 
    if (MeltHeight2 < 0) then !do not add ice from below, only melt
     Melt_energy2 = 0
     MeltHeight2 = 0
    end if 
    ! set surface temperature to freezing point
    state%T(self%grid%ubnd_vol) = param%freez_temp ![°C] 
    ! new ice hieght   
    state%ice_h = state%ice_h - MeltHeight2 ! [m]
      ! if melting larger than ice height, put remaining energy to water and set ice/snow to zero
      if (state%ice_h < 0) then
       state%heat = state%heat + (l_h / (ice_dens * 1 * 1) / (-1 * state%ice_h))![J/kg] / [kg/m3] / [m2] / [m] = [J]
       state%ice_h = 0
       state%snow_dens = rho_s_0
       state%snow_h = 0
      end if     

      !write (*, *) 'Ice Expending'
      !write (*,'(A,F8.6)') 'Ice height : ', state%ice_h
      !write (*,'(A,F7.3)') 'Air temp   : ', state%T_atm
      !write (*,'(A,F10.3)') 'Heat Normal  : ' , state%heat
      !write (*,'(A,F10.3)') 'Heat IceSnow : ' , state%heat_snow_ice    
      !write (*,'(A,F10.6)') 'Melt height water  : ' , MeltHeight2     
      !write (*,'(A,F10.3)') 'Melt energy water  : ' , Melt_energy2   
      !write (*, *) ''
   
  end if 
 end subroutine

 ! Ice melting
 subroutine ice_melting(self, state, param)
      implicit none   
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param   
      ! Define variables only used in ice_model
      real(RK) :: Melt_energy1
      real(RK) :: Melt_energy2
      real(RK) :: MeltHeight1
      real(RK) :: MeltHeight2

      !write (*,'(A,F8.6)') 'Ice height start     : ' , state%ice_h
   
   
      !Melt Ice from atmosphere
      Melt_energy1 = state%heat_snow_ice * state%dt * 1 * 1 ![W/m2] * [s] * [m2] = [J]
      MeltHeight1 = Melt_energy1 / l_h / (ice_dens * 1 * 1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]
      if (MeltHeight1 < 0) then !do not add ice from above, only melt
       Melt_energy1 = 0;
       MeltHeight1 = 0;
      end if    
   
     !Heating of first water cell might melt ice from below, put back first cell to freez temp., and melt ice from below
      Melt_energy2 = (state%T(self%grid%ubnd_vol) - param%freez_temp)  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ![J/kg/K]*[K]*[kg] = [J] 
      MeltHeight2 = Melt_energy2 / l_h / (ice_dens * 1 * 1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]  
      if (MeltHeight2 < 0) then !do not add ice from below, only melt
       Melt_energy2 = 0;
       MeltHeight2 = 0;
      end if    
       ! set surface temperature to freezing point
      state%T(self%grid%ubnd_vol) = param%freez_temp ![°C]    
   
     ! New ice hight
      state%ice_h = state%ice_h - MeltHeight2 - MeltHeight1! [m]

      ! If melt larger than ice hight, put remaining energy to water and set ice hight to zero
      if (state%ice_h < 0) then
       state%heat = state%heat + (l_h / (ice_dens * 1 * 1) / (-1 * state%ice_h))![J/kg] / [kg/m3] / [m2] / [m] = [J]
       state%ice_h = 0
       state%snow_dens = rho_s_0
       state%snow_h = 0
      end if 
    
      !set ice temp to freez point 
      state%ice_temp = 0 ![°C]
     
      !write (*, *) 'Ice Melting'
      !write (*,'(A,F8.6)') 'Ice height          : ' , state%ice_h
      !write (*,'(A,F10.6)') 'Melt height air    : ' , MeltHeight1   
      !write (*,'(A,F10.6)') 'Melt height water  : ' , MeltHeight2     
      !write (*,'(A,F10.6)') 'Snow Height        : ' , state%snow_h    
      !write (*,'(A,F7.3)') 'Air temp            : ' , state%T_atm  
      !write (*,'(A,F10.3)') 'Melt energy air    : ' , Melt_energy1
      !write (*,'(A,F10.3)') 'Melt energy water  : ' , Melt_energy2   
      !write (*,'(A,F10.3)') 'Heat Normal        : ' , state%heat
      !write (*,'(A,F10.3)') 'Heat IceSnow       : ' , state%heat_snow_ice    
      !write (*, *) ''
 end subroutine
 
 ! Snow layer buildup
 subroutine snow_build(self, state, param)
      implicit none   
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param  
      real(RK) :: snow_h_new
      real(RK) :: T0
      real(RK) :: C1   
      T0 = 273
      C1 = 3.2 / 3600 ! compresion [m/sec], initaly [m/h], mean from Yen 1981 page 5


        !calculate new snow height
        snow_h_new = (state%percip / 3600 *state%dt) * (rho_0 / rho_s_0) !go from m/h to m/s and increas volume from water to snow

        !calculate snow density due to compresion of snow layer, Yen 1981 eq. 9
        state%snow_dens = state%snow_dens + (snow_h_new * rho_s_0) * C1 * exp(-0.08 * (T0 - (state%snow_temp + T0))) * state%dt 
 
        !Add the two layers density togethere 
        state%snow_dens  = (state%snow_dens * state%snow_h) / (state%snow_h + snow_h_new) + (rho_s_0 * snow_h_new) / (state%snow_h + snow_h_new)
        if (state%snow_dens > rho_s_max) then ! maximum allowed snow density
           state%snow_dens = rho_s_max 
        end if
  
       !Update snow height
       state%snow_h = state%snow_h + snow_h_new
      
       !snow temp
       state%snow_temp = state%t_atm
       !state%snow_temp = state%snow_temp + (1/cp_s) * (1/state%snow_dens)/(state%snow_h*1*1) * state%heat_snow_ice * state%dt*1*1
       !1/[J/kg/°C] * 1/[kg/m3]/[m3] * [W/m2]*[s]*[m2] = [J°C]/[J] = [°C]

       !write (*, *) 'Snow layer build-up'    
       !write (*,'(A,F10.6)') 'Snow Height       : ' , state%snow_h  
       !write (*,'(A,F10.6)') 'Percipitation     : ' , state%percip  
       !write (*,'(A,F12.2)') 'Snow Temp         : ' , state%snow_temp
       !write (*,'(A,F12.2)') 'Snow Dens         : ' , state%snow_dens     
       !write (*,'(A,F12.2)') 'Heat Flux         : ' , state%heat_snow_ice  
    
 end subroutine
 
 
 ! Snow melting, and ice melting from below
 subroutine snow_melting(self, state, param)
      implicit none   
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param   
      ! Define variables only used in snow_model
      real(RK) :: Melt_energy1
      real(RK) :: MeltHeight1
      real(RK) :: Melt_ice
      real(RK) :: Melt_energy2
      real(RK) :: MeltHeight2   
      real(RK) :: MeltHeightIce
     
   !set density and temperature of melting snow
      state%snow_dens = 450
      state%snow_temp = param%Freez_Temp   
  
      !Melt Snow from atmosphere
      Melt_energy1 = state%heat_snow_ice * state%dt * 1 * 1 ![W/m2] * [s] * [m2] = [J]
      MeltHeight1 = Melt_energy1 / l_h / (state%snow_dens * 1 * 1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]
      if (MeltHeight1 < 0) then !do not add snow from above, only melt snow
       Melt_energy1 = 0
       MeltHeight1 = 0
      end if  

       ! change the snow hight
       state%snow_h = state%snow_h - MeltHeight1! [m]
       ! put remaining energy to melting of ice and set snow hight to zero
        if (state%snow_h < 0) then
         Melt_ice = (l_h / (state%snow_dens * 1 * 1) / (-1 * state%snow_h))![J/kg] / [kg/m3] / [m2] / [m] = [J]
         state%snow_h = 0

         !Melt Ice with remaining energy
         MeltHeightIce = Melt_ice / l_h / (ice_dens * 1 * 1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]  
          if (MeltHeightIce < 0) then !do not add ice from below, only melt ice
           Melt_ice = 0
           MeltHeightIce = 0
          end if   
        else
           Melt_ice = 0
           MeltHeightIce = 0
        end if    
  
  
     !Heating of first water cell might melt ice from below, put back first cell to freez temp., and melt ice from below
      Melt_energy2 = (state%T(self%grid%ubnd_vol) - param%freez_temp)  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ![J/kg/K]*[K]*[kg] = [J] 
      MeltHeight2 = Melt_energy2 / l_h / (ice_dens * 1 * 1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]  
      if (MeltHeight2 < 0) then !do not add ice from below, only melt
       Melt_energy2 = 0;
       MeltHeight2 = 0;  
      end if    
      ! set surface temperature to freezing point
      state%T(self%grid%ubnd_vol) = param%freez_temp ![°C]    
   
     ! New ice hight
      state%ice_h = state%ice_h - MeltHeight2 - MeltHeightIce! [m] 
   
   ! if melting larger than ice height, put remaining energy to water and set ice/snow to zero
      if (state%ice_h < 0) then
       state%heat = state%heat + (l_h / (ice_dens * 1 * 1) / (-1 * state%ice_h))![J/kg] / [kg/m3] / [m2] / [m] = [J]
       state%ice_h = 0
       state%snow_dens = rho_s_0
       state%snow_h = 0
      end if     

      !write (*, *) 'Snow Melting'
      !write (*,'(A,F10.6)') 'Snow Height       : ' , state%snow_h  
      !write (*,'(A,F12.2)') 'Snow Temp         : ' , state%snow_temp
      !write (*,'(A,F12.2)') 'Snow Dens         : ' , state%snow_dens     
      !write (*,'(A,F12.2)') 'Heat Flux         : ' , state%heat_snow_ice  
      !write (*,'(A,F8.6)') 'Ice height : ', state%ice_h
      !write (*,'(A,F7.3)') 'Air temp   : ', state%T_atm   
 end subroutine

 
end module

