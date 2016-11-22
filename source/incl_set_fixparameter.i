!     Sets value of fixed parameters

!      *** Constants for freshwater ***
          rho_0 = 1000        ! Mean freshwater density (seawater: 1023) [kg/m3]
          cp = 4185           ! Mean freshwater heat capacity (seawater: 3992) [J/kg/K]
!      *** General constants ***
          pi = 3.141592654    ! pi [-]
          rho_air = 1.2       ! density of air [kg/m3]
          g = 9.81            ! Earth gravitational acceleration [m/s2]
          kappa = 0.41        ! Von Karman constant [-]
          K_s = 0.05          ! Bottom roughness [m]
          z0 = 0.5            ! Surface roughness [m]
!      *** Further parameters controlling water dynamic ***
          Prndtl = 0.8        ! Prandtl number for air??
          k_min = 1.0e-9      ! Min value allowed for turbulent kinetic energy (read from file!?->useless here)
          eps_min = 1.0e-30   ! Min value allowed for TKE dissipation
          avh_min = 1.0e-8    ! Min value allowed for turbulent viscosity at boundaries
!      *** Parameters for k-eps model ***
          ce1 = 1.44
          ce2 = 1.92
          ce3 = -0.4
          sig_k = 1.0
          sig_e = 1.3
          cmue = 0.09
!      *** Parameters for MY model (not used) ***
          a1 = 0.92
          a2 = 0.74
          b1 = 16.6
          b2 = 10.1
          c1 = 0.08
          e1 = 1.8
          e2 = 1.33
          sl = 0.2
