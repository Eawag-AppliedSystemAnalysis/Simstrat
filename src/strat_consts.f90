!     +---------------------------------------------------------------+
!     |  Constants....
!     +---------------------------------------------------------------+
  module strat_consts
     use strat_kinds
     implicit none

    ! *** General constants ***
     real(RK), parameter     :: pi = 4.0_RK*DATAN(1.0_RK) ! pi [-] (Proper way to calculate PI up to machine precision)
     real(RK), parameter     :: rho_air = 1.2_RK ! density of air [kg/m3]
     real(RK), parameter     :: g = 9.81_RK ! Earth gravitational acceleration [m/s2]
     real(RK), parameter     :: kappa = 0.41_RK ! Von Karman constant [-]
     real(RK), parameter     :: K_s = 0.05_RK ! Bottom roughness [m]
     real(RK), parameter     :: z0 = 0.5_RK ! Surface roughness [m]

    ! *** Constants for freshwater ***
     real(RK), parameter     :: rho_0 = 1000_RK ! Mean freshwater density (seawater: 1023) [kg/m3]
     real(RK), parameter     :: cp = 4182_RK ! Mean freshwater heat capacity (seawater: 3992) [J/kg/K]
     real(RK), parameter     :: cp_air = 1005_RK ! Mean air heat capacity [J/kg/K]

    ! *** Further parameters controlling water dynamic ***
     real(RK), parameter     :: Prndtl = 0.8_RK ! Prandtl number for air??
     real(RK), parameter     :: k_min = 1.0e-30_RK ! Min value allowed for turbulent kinetic energy (read from file!?->useless here)
     real(RK), parameter     :: eps_min = 1.0e-30_RK ! Min value allowed for TKE dissipation
     real(RK), parameter     :: avh_min = 1.0e-8_RK ! Min value allowed for turbulent viscosity at boundaries

    ! *** Ice Parameters
     real(RK), parameter     :: k_ice = 2.22_RK ! Thermal conductivity ice at 0°C WaterTemp (W K-1 m-1)
     real(RK), parameter     :: k_snow = 0.2_RK ! Thermal conductivity snow at 0°C WaterTemp (W K-1 m-1)
     real(RK), parameter     :: l_h = 3.34e+5_RK! Latent heat of melting, [J kg-1]
     real(RK), parameter     :: l_e = 0 !2.265e+6_RK! Latent heat of evaporation, [J kg-1], l_e = 0 => non-sublimation or l_e ~= 0 => sublimation (i.e. solid to gas)  
     real(RK), parameter     :: ice_dens = 916.2_RK !ice density kg m-3
     real(RK), parameter     :: snowice_dens = 875_RK !snow ice density  [kg m-3] Saloranta 2000
     real(RK), parameter     :: rho_s_0 = 250_RK ! New snowfall density [kg/m3]
     real(RK), parameter     :: rho_s_max = 450_RK ! maximum (wet) snow dens [kg/m3]
     real(RK), parameter     :: cp_s = 2090_RK ! Mean snow heat capacity (-5°C) [J/kg/K]
     real(RK), parameter     :: emiss_water = 0.97_RK ! Emissivity water
     real(RK), parameter     :: emiss_ice = 0.98_RK   ! Emissivity ice
     ! Emissivity of snow ranges from 0.8 to 0.9 depending on snow density in strat_forcing.f90
     real(RK), parameter     :: C01 = 5.8_RK / 3600_RK ! snow compresion in [m^-1 sec^-1] to match model timestep, initaly [m/h], from Yen (1981) page 5, C1 range [2.6 to 9.0] 
     real(RK), parameter     :: C02 = 21.0_RK / 1000_RK !snow compresion in [m^3/kg] initaly [m^3/Mg], from Yen (1981) page 5
     real(RK), parameter     :: Ha_a = 0.68_RK ! Longwave emision parameter, Matti Leppäranta 2009
     real(RK), parameter     :: Ha_b = 0.036_RK !Longwave emision parameter, Matti Leppäranta 2009 [mbar^-1/2]
     real(RK), parameter     :: Ha_c = 0.18_RK !Longwave emision parameter, Matti Leppäranta 2009
     real(RK), parameter     :: Hk_CH = 1.5e-3_RK !convectiv bulk exchange coefficient, Matti Leppäranta 2009 and Gill 1982
     real(RK), parameter     :: Hv_CE = 1.5e-3_RK !latent bulk exchange coefficient, Matti Leppäranta 2009 and Gill 1982
     real(RK), parameter     :: lambda_snow = 24_RK ! lambda (light absorption) snow [m-1]
     real(RK), parameter     :: lambda_snowice = 3_RK ! lambda (light absorption) white ice [m-1]
     real(RK), parameter     :: lambda_ice = 1_RK ! lambda (light absorption) black ice [m-1]
     real(RK), parameter     :: ice_albedo = 0.25_RK ! Albedo black ice
     real(RK), parameter     :: snowice_albedo = 0.35_RK ! Albedo white ice
     real(RK), parameter     :: snow_albedo = 0.7_RK ! Albedo snow

    ! *** Parameters for k-eps model ***
     real(RK), parameter     :: ce1 = 1.44_RK
     real(RK), parameter     :: ce2 = 1.92_RK
     real(RK), parameter     :: ce3 = -0.4_RK
     real(RK), parameter     :: sig_k = 1.0_RK
     real(RK)                :: sig_e = 1.3_RK
     real(RK), parameter     :: cmue = 0.09_RK
     real(RK), parameter     :: r_a = 0.03_RK ! Ratio of reflected to total long-wave iradiance
     real(RK), parameter     :: B0 = 0.61_RK ! Bowen constant
     real(RK), parameter     :: sig = 5.67e-8_RK ! Stefan-Boltzmann constant [W/m2/K4]

     ! *** Parameters for MY model  ***
     real(RK), parameter     :: a1 = 0.92_RK
     real(RK), parameter     :: a2 = 0.74_RK
     real(RK), parameter     :: b1 = 16.6_RK
     real(RK), parameter     :: b2 = 10.1_RK
     real(RK), parameter     :: c1 = 0.08_RK
     real(RK), parameter     :: e1 = 1.8_RK
     real(RK), parameter     :: e2 = 1.33_RK
     real(RK), parameter     :: sl = 0.2_RK

  end module strat_consts
