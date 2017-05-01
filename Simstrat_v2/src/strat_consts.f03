module strat_consts
  use strat_kinds
  implicit none

!   *** General constants ***
    real(RK), parameter     :: pi = 4.0_RK*DATAN(1.0_RK)   ! pi [-]
    real(RK), parameter     :: rho_air = 1.2_RK       ! density of air [kg/m3]
    real(RK), parameter     :: g = 9.81_RK            ! Earth gravitational acceleration [m/s2]
    real(RK), parameter     :: kappa = 0.41_RK        ! Von Karman constant [-]
    real(RK), parameter     :: K_s = 0.05_RK          ! Bottom roughness [m]
    real(RK), parameter     :: z0 = 0.5_RK            ! Surface roughness [m]

!   *** Constants for freshwater ***
    real(RK), parameter    :: rho_0 = 1000_RK        ! Mean freshwater density (seawater: 1023) [kg/m3]
    real(RK), parameter    :: cp = 4182_RK           ! Mean freshwater heat capacity (seawater: 3992) [J/kg/K]

!   *** Further parameters controlling water dynamic ***
    real(RK), parameter     :: Prndtl = 0.8_RK        ! Prandtl number for air??
    real(RK)                :: k_min = 1.0e-9_RK      ! Min value allowed for turbulent kinetic energy (read from file!?->useless here)
    real(RK), parameter     :: eps_min = 1.0e-30_RK   ! Min value allowed for TKE dissipation
    real(RK), parameter     :: avh_min = 1.0e-8_RK    ! Min value allowed for turbulent viscosity at boundaries

!   *** Parameters for k-eps model ***
    real(RK), parameter     :: ce1 = 1.44_RK
    real(RK), parameter     :: ce2 = 1.92_RK
    real(RK), parameter     :: ce3 = -0.4_RK
    real(RK), parameter     :: sig_k = 1.0_RK
    real(RK)                :: sig_e = 1.3_RK
    real(RK), parameter     :: cmue = 0.09_RK

    real(RK), parameter     :: r_a = 0.03_RK           ! Ratio of reflected to total long-wave iradiance
    real(RK), parameter     :: B0 = 0.61_RK            ! Bowen constant
    real(RK), parameter     :: sig = 5.67e-8_RK        ! Stefan-Boltzmann constant [W/m2/K4]

end module strat_consts
