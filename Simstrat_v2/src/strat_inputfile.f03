module simstrat_inputfile_module
  use strat_kinds
  use strat_simdata
  use strat_grid
  use strat_consts
  use utilities
  use json_kinds, only: CK
  use json_module
  use csv_module
  implicit none

  private

  !##################################################
  !# Inputfile
  !##################################################
  type, public :: SimstratSimulationFactory
    private
    class(SimulationData), pointer :: simdata
    contains
      procedure, pass(self), public :: initialize_model
      procedure, pass(self), public :: read_json_par_file
      procedure, pass(self), public :: read_initial_data
      procedure, pass(self), public :: read_grid_config
      procedure, pass(self), public :: setup_model
  end type SimstratSimulationFactory

contains

  subroutine initialize_model(self, fname, simdata)
    implicit none
    class(SimstratSimulationFactory) :: self
    class(SimulationData), pointer, intent(out) :: simdata
    character(len=*) :: fname
    logical :: file_exists

        !allocate model
    !if(associated(simdata)) deallocate(simdata)

    allocate(SimulationData :: self%simdata)
    simdata => self%simdata

    !Parse inputfile
    call self%read_json_par_file(fname)

    !Set up grid
    call self%read_grid_config

    call self%simdata%model%init(self%simdata%grid%nz_grid)

    !Read initial data
    call self%read_initial_data

    ! Update area factors
    call self%simdata%grid%update_area_factors()

    ! Init rest of model
    call self%setup_model


  end subroutine initialize_model

  subroutine setup_model(self)
    implicit none
    class(SimstratSimulationFactory) :: self
    integer :: i
    associate(simdata => self%simdata, &
              model_cfg => self%simdata%model_cfg, &
              model_param => self%simdata%model_param, &
              model => self%simdata%model, &
              grid => self%simdata%grid)

    ! Initialize some more values
    if(model_cfg%stability_func==1) model%cm0=0.5625_RK
    if(model_cfg%stability_func==2) model%cm0=0.556171_RK
    model%cde=model%cm0**3
    sig_e=(kappa/model%cm0)**2/(ce2-ce1)

    model%num(0:grid%nz_grid) = 0.0_RK
    model%nuh(0:grid%nz_grid)= 0.0_RK

    model%tx=0.0_RK
    model%ty=0.0_RK

    model%drag=(kappa/log(1.0_RK+30/K_s*grid%h(1)/2))**2


    ! Geothermal heat flux
    if (model_param%fgeo/=0) then
        model%fgeo_add(1:grid%nz_grid) = model_param%fgeo/rho_0/cp*grid%dAz(1:grid%nz_grid)/grid%Az(1:grid%nz_grid) ! calculation per kg
        if (grid%Az(0)/=0) then
            model%fgeo_add(1) = model%fgeo_add(1)+2*model_param%fgeo/rho_0/cp*grid%Az(0)/((grid%Az(0)+grid%Az(1))*grid%h(1))
        end if
    end if

    ! Salinity control for buoyancy functions
    ! if salinity transport is enabled
    if (model_cfg%salinity_transport) then
        model%has_salinity=.true.
        model%has_salinity_grad=.true.
    else ! else, test for this config
        model%has_salinity = .false.
        model%has_salinity_grad = .false.
        do i=1,grid%nz_grid
            if(model%S(i)/=0) then
              model%has_salinity=.true.
              exit
            end if

        end do
        if (model%has_salinity) then
            do i=2,grid%nz_grid
                if(model%S(i)-model%S(i-1)/=0) then
                  model%has_salinity_grad=.true.
                end if
            end do
        end if
    end if

  end associate
  end subroutine

  subroutine read_grid_config(self)
    implicit none
    class(SimstratSimulationFactory) :: self
    type(GridConfig) :: grid_config
    real(RK), dimension(:) :: z_tmp(0:self%simdata%model_cfg%max_nr_grid_cells)
    real(RK), dimension(:) :: A_tmp(0:self%simdata%model_cfg%max_nr_grid_cells)
    integer :: num_read, i, ictr
    associate(simdata => self%simdata, &
              nz_max => self%simdata%model_cfg%max_nr_grid_cells)

    grid_config%nz_grid_max = self%simdata%model_cfg%max_nr_grid_cells
    allocate(grid_config%grid_read(0:nz_max))
    ! Read grid
    open(12,status='old',file=simdata%input_cfg%GridName)
    read(12,*)
    do ictr=0,nz_max
        read(12,*,end=69) grid_config%grid_read(ictr)
    end do
69  if(ictr==nz_max) write(6,*) 'Only first ',nz_max,' values of file read.'
    close(12)

    if (ictr==1) then     ! Constant spacing
        grid_config%nz_grid=int(grid_config%grid_read(0))
        grid_config%equidistant_grid = .TRUE.
    else                  ! Variable spacing
        grid_config%nz_grid = ictr-1
        grid_config%equidistant_grid = .FALSE.
      end if

    ! Read Morphology
    open(11,status='old',file=simdata%input_cfg%MorphName)
    read(11,*)                            ! Skip header
    do i=0,nz_max                            ! Read depth and area
        read(11,*,end=86) z_tmp(i), A_tmp(i)
    end do
86  if(i==nz_max) write(6,*) 'Only first ',nz_max,' values of file read.'
    close(11)

    num_read = i  ! Number of area values

    allocate(grid_config%z_A_read(0:num_read-1), grid_config%A_read(0:num_read-1))

    ! Reverse order of values
    do i=0,num_read-1
        grid_config%z_A_read(i) = -z_tmp(num_read-i-1)
        grid_config%A_read(i) = A_tmp(num_read-i-1)
    end do

    grid_config%depth = grid_config%z_A_read(0) - grid_config%z_A_read(num_read-1)  ! depth = max - min depth

    ! initialize Grid of simdata
    call simdata%grid%init(grid_config)
  end associate
  end subroutine

  !#######################################################################
  subroutine read_json_par_file(self, ParName)
  !#######################################################################
    implicit none
    class(SimstratSimulationFactory) :: self
    character(kind=CK, len=*), intent(in) :: ParName

    type(json_file) :: par_file
    logical :: found

    !gfortran cannot handle type bound allocatable character that are passed to subroutine as intent(out)
    !as a workaround we have to store the values in a local scope allocatable character
    character(kind=CK, len=:), allocatable          :: MorphName,InitName,ForcingName,AbsorpName
    character(kind=CK, len=:), allocatable          :: GridName,zoutName,toutName,PathOut
    character(kind=CK, len=:), allocatable          :: QinpName,QoutName,TinpName,SinpName
    character(kind=CK, len=:), allocatable          :: stencil_type

    !model%ParName = ParName
    !check if inputfile SimstratModelexists
    call check_file_exists(ParName)

    call par_file%initialize()

    !load file or stop if fail
    call par_file%load_file(filename = ParName)
    if(par_file%failed()) then
      call error('Could not read inputfile '//ParName)
    stop
    end if

    !Names of Inputfile
    call par_file%get('Input.Morphology'        , MorphName, found);  self%simdata%input_cfg%MorphName = MorphName; call check_field(found, 'Input.Morphology'      , ParName)
    call par_file%get('Input.Initial conditions', InitName, found);   self%simdata%input_cfg%InitName = InitName; call check_field(found, 'Input.Initial conditions' , ParName)
    call par_file%get('Input.Forcing'           , ForcingName, found);self%simdata%input_cfg%ForcingName = ForcingName; call check_field(found, 'Input.Forcing'   , ParName)
    call par_file%get('Input.Absorption'        , AbsorpName, found); self%simdata%input_cfg%AbsorpName = AbsorpName; call check_field(found, 'Input.Absorption'   , ParName)
    call par_file%get('Input.Grid'              , GridName, found);   self%simdata%input_cfg%GridName = GridName; call check_field(found, 'Input.Grid'               , ParName)
    call par_file%get('Input.Inflow'            , QinpName, found);   self%simdata%input_cfg%QinpName = QinpName; call check_field(found, 'Input.Inflow'             , ParName)
    call par_file%get('Input.Outflow'           , QoutName, found);   self%simdata%input_cfg%QoutName = QoutName; call check_field(found, 'Input.Outflow'            , ParName)
    call par_file%get('Input.TemperatureInflow' , TinpName, found);   self%simdata%input_cfg%TinpName = TinpName; call check_field(found, 'Input.TemperatureInflow'  , ParName)
    call par_file%get('Input.SalinityInflow'    , SinpName, found);   self%simdata%input_cfg%SinpName = SinpName; call check_field(found, 'Input.SalinityInflow'     , ParName)

    !Output
    call par_file%get('Output.path'             , PathOut, found); self%simdata%output_cfg%PathOut = PathOut; call check_field(found, 'Output.path'     , ParName)
    call par_file%get('Output.depth'            , zoutName, found); self%simdata%output_cfg%zoutName = zoutName; call check_field(found, 'Output.depth' , ParName)
    call par_file%get('Output.WriteOnTheFly'    , self%simdata%output_cfg%write_on_the_fly, found); call check_field(found, 'Output.WriteOnTheFly'  , ParName)
    call par_file%get('Output.ThinningInterval' , self%simdata%output_cfg%thinning_interval, found); call check_field(found, 'Output.ThinningInterval'  , ParName)

    !Model configuration
    call par_file%get("ModelConfig.MaxNrGridCells"         , self%simdata%model_cfg%max_nr_grid_cells   , found);
    if(.not.found) then
      self%simdata%model_cfg%max_nr_grid_cells = 1000
      call warn('Variable "ModelConfig.MaxNrGridCells" is not set. Assume a value of 1000')
    end if
    call par_file%get("ModelConfig.CoupleAED2"             , self%simdata%model_cfg%couple_aed2  , found);
    if(.not.found) then
      self%simdata%model_cfg%couple_aed2 = .false.
      call warn('Variable "ModelConfig.CoupleAED2" is not set. Assume you do not want to couple simstrat with aed2.')
    end if
    call par_file%get("ModelConfig.TurbulenceModel"        , self%simdata%model_cfg%turbulence_model    , found); call check_field(found, 'ModelConfig.TurbulenceModel'         , ParName)
    call par_file%get("ModelConfig.StabilityFunction"      , self%simdata%model_cfg%stability_func      , found); call check_field(found, 'ModelConfig.StabilityFunction'       , ParName)
    call par_file%get("ModelConfig.FluxCondition"          , self%simdata%model_cfg%flux_condition      , found); call check_field(found, 'ModelConfig.FluxCondition'           , ParName)
    call par_file%get("ModelConfig.Forcing"                , self%simdata%model_cfg%forcing_mode        , found); call check_field(found, 'ModelConfig.Forcing'                 , ParName)
    call par_file%get("ModelConfig.UseFilteredWind"        , self%simdata%model_cfg%use_filtered_wind   , found); call check_field(found, 'ModelConfig.UseFilteredWind'         , ParName)
    call par_file%get("ModelConfig.SeicheNormalization"    , self%simdata%model_cfg%seiche_normalization, found); call check_field(found, 'ModelConfig.SeicheNormalization'     , ParName)
    call par_file%get("ModelConfig.Wind drag model"        , self%simdata%model_cfg%wind_drag_model     , found); call check_field(found, 'ModelConfig.Wind drag model'         , ParName)
    call par_file%get("ModelConfig.InflowPlacement"        , self%simdata%model_cfg%inflow_placement    , found); call check_field(found, 'ModelConfig.InflowPlacement'         , ParName)
    call par_file%get("ModelConfig.PressureGradients"      , self%simdata%model_cfg%pressure_gradients  , found); call check_field(found, 'ModelConfig.PressureGradients'       , ParName)
    call par_file%get("ModelConfig.EnableSalinityTransport", self%simdata%model_cfg%salinity_transport  , found); call check_field(found, 'ModelConfig.EnableSalinityTransport' , ParName)
    call par_file%get("ModelConfig.DisplaySimulation"      , self%simdata%model_cfg%disp_simulation     , found); call check_field(found, 'ModelConfig.DisplaySimulation'       , ParName)
    call par_file%get("ModelConfig.DisplayDiagnose"        , self%simdata%model_cfg%disp_diagnostic     , found); call check_field(found, 'ModelConfig.DisplayDiagnose'         , ParName)
    call par_file%get("ModelConfig.DataAveraging"          , self%simdata%model_cfg%data_averaging      , found); call check_field(found, 'ModelConfig.DataAveraging'           , ParName)

    !Model Parameter
    call par_file%get("ModelParameter.Lat"       , self%simdata%model_param%Lat     , found); call check_field(found, 'ModelParameter.Lat'     , ParName)
    call par_file%get("ModelParameter.p_air"     , self%simdata%model_param%p_air   , found); call check_field(found, 'ModelParameter.p_air'   , ParName)
    call par_file%get("ModelParameter.a_seiche"  , self%simdata%model_param%a_seiche, found); call check_field(found, 'ModelParameter.a_seiche', ParName)
    call par_file%get("ModelParameter.q_NN"      , self%simdata%model_param%q_NN    , found); call check_field(found, 'ModelParameter.q_NN'    , ParName)
    call par_file%get("ModelParameter.f_wind"    , self%simdata%model_param%f_wind  , found); call check_field(found, 'ModelParameter.f_wind'  , ParName)
    call par_file%get("ModelParameter.C10"       , self%simdata%model_param%C10     , found); call check_field(found, 'ModelParameter.C10'     , ParName)
    call par_file%get("ModelParameter.CD"        , self%simdata%model_param%CD      , found); call check_field(found, 'ModelParameter.CD'      , ParName)
    call par_file%get("ModelParameter.fgeo"      , self%simdata%model_param%fgeo    , found); call check_field(found, 'ModelParameter.fgeo'    , ParName)
    call par_file%get("ModelParameter.k_min"     , self%simdata%model_param%k_min   , found); call check_field(found, 'ModelParameter.k_min'   , ParName)
    call par_file%get("ModelParameter.p_radin"   , self%simdata%model_param%p_radin , found); call check_field(found, 'ModelParameter.p_radin' , ParName)
    call par_file%get("ModelParameter.p_windf"   , self%simdata%model_param%p_windf , found); call check_field(found, 'ModelParameter.p_windf' , ParName)
    call par_file%get("ModelParameter.beta_sol"  , self%simdata%model_param%beta_sol, found); call check_field(found, 'ModelParameter.beta_sol', ParName)
    call par_file%get("ModelParameter.albsw"     , self%simdata%model_param%albsw   , found); call check_field(found, 'ModelParameter.albsw'   , ParName)

    !Simulation Parameter
    call par_file%get("Simulation.Timestep", self%simdata%sim_cfg%timestep     , found); call check_field(found, 'Simulation.Timestep' , ParName)
    call par_file%get("Simulation.Start"   , self%simdata%sim_cfg%start_datum, found); call check_field(found, 'Simulation.Start'    , ParName)
    call par_file%get("Simulation.End"     , self%simdata%sim_cfg%end_datum  , found); call check_field(found, 'Simulation.End'      , ParName)

    call par_file%destroy()

    !check validity of inputfile
  !  if(.not.(simdata%model_cfg%disp_sim==1 .or. simdata%model_cfg%disp_sim==2 .or. simdata%model_cfg%disp_sim==3)) simdata%model_cfg%disp_sim=0
  !  if(.not.(simdata%model_cfg%disp_dgn==1 .or. simdata%model_cfg%disp_dgn==2)) simdata%model_cfg%disp_dgn=0

  !  if(simdata%model_cfg%disp_dgn/=0) then
      call ok('Configuration: '//trim(ParName))
  !  end if

  end subroutine read_json_par_file

  subroutine read_initial_data(self)
    implicit none
    class(SimstratSimulationFactory) :: self

    ! Local variables
    integer,parameter :: nz_max =1000

    real(RK) :: z_tmp(0:nz_max), U_tmp(0:nz_max), V_tmp(0:nz_max)
    real(RK) :: T_tmp(0:nz_max), S_tmp(0:nz_max), k_tmp(0:nz_max), eps_tmp(0:nz_max)

    real(RK) :: z_read(0:nz_max), U_read(0:nz_max), V_read(0:nz_max)
    real(RK) :: T_read(0:nz_max), S_read(0:nz_max), k_read(0:nz_max), eps_read(0:nz_max)

    real(RK) :: z_ini(0:nz_max)
    real(RK) :: z_ini_depth, zmax

    integer :: i,num_read

    associate(grid => self%simdata%grid, &
              model => self%simdata%model, &
              nz_inuse => self%simdata%grid%nz_inuse)

    ! Read file
    open(13,status='old',file=self%simdata%input_cfg%InitName)     ! Opens initial conditions file
    read(13,*)                              ! Skip header
    do i=0,nz_max                              ! Read initial u,v,T, etc
        read(13,*,end=9) z_read(i),U_read(i),V_read(i),T_read(i),S_read(i),k_read(i),eps_read(i)
    end do
9   num_read = i-1                                ! Number of valuInitNamees
    if (num_read<0) then
        write(6,*) 'Error reading initial conditions files (no data found).'
        stop
    end if
    close(13)
    do i=0,num_read
        z_read(i) = abs(z_read(i))               ! Make depths positive
    end do
    z_ini_depth = z_read(0)                     ! Initial depth (top-most)

    ! update actual filled z in grid
    call grid%update_depth(z_ini_depth)

    do i=0,num_read
        z_tmp(num_read-i) = grid%z_zero - z_read(i)
        U_tmp(num_read-i) = U_read(i)
        V_tmp(num_read-i) = V_read(i)
        T_tmp(num_read-i) = T_read(i)
        S_tmp(num_read-i) = S_read(i)
        k_tmp(num_read-i) = k_read(i)
        eps_tmp(num_read-i) = eps_read(i)
    end do

    if (num_read==0) then
        write(6,*) 'Only one row! Water column will be initially homogeneous.'
        model%U(0:nz_inuse) = U_tmp(0)
        model%V(0:nz_inuse) = V_tmp(0)
        model%T(0:nz_inuse) = T_tmp(0)
        model%S(0:nz_inuse) = S_tmp(0)
        model%k(0:nz_inuse) = k_tmp(0)
        model%eps = eps_tmp(0)
    else

        ! interpolate variables UVTS on central grid and store
        call grid%interpolate_cent(z_tmp, U_tmp, num_read, model%U)
        call grid%interpolate_cent(z_tmp, V_tmp, num_read, model%V)
        call grid%interpolate_cent(z_tmp, T_tmp, num_read, model%T)
        call grid%interpolate_cent(z_tmp, S_tmp, num_read, model%S)

        ! Interpolate k/eps on upper grid and store
        call grid%interpolate_upp(z_tmp, k_tmp, num_read, model%k)
        call grid%interpolate_upp(z_tmp, eps_tmp, num_read, model%eps)
    end if

  end associate
  end subroutine

  subroutine check_field(found, field_name, file_name)
    implicit none

    logical, intent(in) :: found
    character(len=*), intent(in) :: field_name, file_name

    if(.not.found) then
      call error('Field '//field_name//' not found in '//file_name)
      stop
    end if
  end subroutine check_field


end module simstrat_inputfile_module
