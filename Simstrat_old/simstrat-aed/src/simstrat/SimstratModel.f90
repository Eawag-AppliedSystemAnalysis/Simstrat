module SimstratModel
    use simstrat_type
    implicit none
    save

!   *** General constants ***
    real(dp), parameter     :: pi = 3.141592654_dp    ! pi [-]
    real(dp), parameter     :: rho_air = 1.2_dp       ! density of air [kg/m3]
    real(dp), parameter     :: g = 9.81_dp            ! Earth gravitational acceleration [m/s2]
    real(dp), parameter     :: kappa = 0.41_dp        ! Von Karman constant [-]
    real(dp), parameter     :: K_s = 0.05_dp          ! Bottom roughness [m]
    real(dp), parameter     :: z0 = 0.5_dp            ! Surface roughness [m]

!   *** Constants for freshwater ***
    real(dp), parameter    :: rho_0 = 1000_dp        ! Mean freshwater density (seawater: 1023) [kg/m3]
    real(dp), parameter    :: cp = 4182_dp           ! Mean freshwater heat capacity (seawater: 3992) [J/kg/K]

!   *** Further parameters controlling water dynamic ***
    real(dp), parameter     :: Prndtl = 0.8_dp        ! Prandtl number for air??
    real(dp)                :: k_min = 1.0e-9_dp      ! Min value allowed for turbulent kinetic energy (read from file!?->useless here)
    real(dp), parameter     :: eps_min = 1.0e-30_dp   ! Min value allowed for TKE dissipation
    real(dp), parameter     :: avh_min = 1.0e-8_dp    ! Min value allowed for turbulent viscosity at boundaries

!   *** Parameters for k-eps model ***
    real(dp), parameter     :: ce1 = 1.44_dp
    real(dp), parameter     :: ce2 = 1.92_dp
    real(dp), parameter     :: ce3 = -0.4_dp
    real(dp), parameter     :: sig_k = 1.0_dp
    real(dp)                :: sig_e = 1.3_dp
    real(dp), parameter     :: cmue = 0.09_dp

!      *** Parameters for MY model (not used) ***
    real(dp), parameter     :: a1 = 0.92_dp
    real(dp), parameter     :: a2 = 0.74_dp
    real(dp), parameter     :: b1 = 16.6_dp
    real(dp), parameter     :: b2 = 10.1_dp
    real(dp), parameter     :: c1 = 0.08_dp
    real(dp), parameter     :: e1 = 1.8_dp
    real(dp), parameter     :: e2 = 1.33_dp
    real(dp), parameter     :: sl = 0.2_dp



    integer, parameter      :: nz_max = 1000

!	*** grid ***
    real(dp)                :: depth,dt,t_start,t_end


    real(dp)                :: cde,cm0
    real(dp)                :: p_air,a_seiche,q_NN,CD,C10,Cori
    real(dp)                :: p_radin,p_windf,beta_sol,albsw,f_wind
    real(dp)                :: fgeo, fsed
    real(dp), dimension(:), allocatable     :: fgeo_add
    integer                 :: salctr, delsal, if_adv, nz_input

!real(dp) :: sig_e = 1.3, k_min = 1.0e-9

!  *** Input Files ***
    character*100           :: ParName,MorphName,InitName,ForcingName,AbsorpName
    character*100           :: GridName,zoutName,toutName,PathOut
    character*100           :: QinpName,QoutName,TinpName,SinpName


    integer                 :: Mod,Stab,ModFlux,NBC,WindFilt,ModSNorm,ModC10,ModInflow,Pgrad,ModSal
    logical                 :: OutBin
    integer                 :: nz_grid,nz,num_save,depth_save,nsave
    real(dp)                :: zsave(0:nz_max)
    integer                 :: igoal,igoal_s(1:14)
    integer                 :: nz_upp,nz_cent,index_upp_save(1:nz_max),index_cent_save(1:nz_max)
    integer                 :: disp_dgn,disp_sim

    real(dp), dimension(:,:), allocatable   :: z_Inp, Q_start, Q_end, depth_surfaceFlow
    real(dp), dimension(:,:), allocatable   :: Inp_read_start, Inp_read_end
    real(dp), dimension(:), allocatable     :: h, A_read, Az, dAdz 
    real(dp), dimension(:), allocatable     :: S_ini, T_ini, U_ini, V_ini
    real(dp), dimension(:), allocatable     :: k_ini, eps_ini, num_ini, nuh_ini
    real(dp)                                :: drag_ini, tx_ini, ty_ini
    real(dp), dimension(:), allocatable     :: z_read, z_cent, z_upp

    real(dp)                                :: volume, z_zero, lake_level_old

    real(dp), dimension(:), allocatable     :: AreaFactor_1, AreaFactor_2
    real(dp), dimension(:), allocatable     :: AreaFactor_k1, AreaFactor_k2, AreaFactor_eps
    real(dp), dimension(:), allocatable     :: meanint

    real(dp)                :: tout_ctr1(0:9000),tout_ctr2(0:9000)
    integer                 :: write_tout

!   *** FABM Integration ***
#if use_fabm
    character*1000          :: PathFABM,PathFABMOut
#endif
!   *** Variables for simulation ***


end module SimstratModel
