module simstrat_inputfile_module
  use simstrat_kinds
  use simstrat_model_module
  use simstrat_model_constants
  use simstrat_discretization_module, only: SimstratDiscretizationScheme
  use simstrat_finite_volume_implementation
  use simstrat_simulation_chain_module
  use simstrat_forcing_module
  use simstrat_output_module
  use simstrat_aed2_module
  use utilities
  use json_kinds, only: CK
  use json_module
  use csv_module
  implicit none

  private

  !##################################################
  !# Inputfile
  !##################################################
  type, public :: SimstratModelFactory
    private
      class(SimstratModel), pointer, public :: model => null()

    contains
      procedure, pass(self), public :: initialize_simstrat_model_from_inputfiles
  end type SimstratModelFactory

contains

  function initialize_simstrat_model_from_inputfiles(self, fname) result(model)
    implicit none
    class(SimstratModelFactory) :: self
    character(len=*) :: fname
    logical :: file_exists

    class(SimstratModel), pointer :: model
    class(SimstratLakePhysics), pointer :: lake_physics

    !check if inputfile exists
    call check_file_exists(fname)

    !allocate model
    if(associated(self%model)) deallocate(self%model)
    allocate(SimstratModel :: self%model)
    model => self%model

    !Parse inputfile
    call parseJSONInputfile(model, fname)

    !Initialize discretization scheme
    call initializeDiscretizationScheme(model)

    !initialize model / allocate space
    call model%initialize()

    !allocate solver
    !Could be implemented if more than one solver is needed.

    !Initialize state variables and initial conditions
    call initializeInitialConditions(model)

    !Initialize Forcing
    if(model%use_buffered_forcing) then
      call error('Buffered reading of forcing file is not implemented.')
      stop
    else
      allocate(MemoryBoundSimstratForcing :: model%forcing)
      call model%forcing%initialize(model%ForcingName, model%ForcingType, model%UseWindFilt, model%disp_dgn, model%AbsorpName, model%z_zero)
      if(model%forcing%datum_forcing(1) > model%t_start .or. model%forcing%datum_forcing(model%forcing%n_forcing) < model%t_end) then
        call error('Forcing file does not cover the whole time span of the simulation')
        stop
      end if
    end if

    !Generate simulation chain
    if(model%couple_aed2) then
      call error('Coupling with AED2 is not implemented')
      stop
    else
      allocate(SimstratLakePhysics :: lake_physics)
      call lake_physics%assembleSimulationChain(model%state_vars(1)%ptr, model%state_vars(2)%ptr, model%state_vars(3)%ptr, model%state_vars(4)%ptr, model%state_vars(5)%ptr, model%state_vars(6)%ptr)
    end if
    !associate simulation_chain
    model%simulation_chain => lake_physics%simulation_chain

    !assign output variables and allocate output logger
    allocate(model%output_vars(11))
    model%output_vars(1)%name = 'T'
    model%output_vars(1)%value => model%T
    model%output_vars(2)%name = 'S'
    model%output_vars(2)%value => model%S
    model%output_vars(3)%name = 'u'
    model%output_vars(3)%value => model%u
    model%output_vars(4)%name = 'v'
    model%output_vars(4)%value => model%v
    model%output_vars(5)%name = 'k'
    model%output_vars(5)%value => model%k
    model%output_vars(6)%name = 'eps'
    model%output_vars(6)%value => model%eps
    model%output_vars(7)%name = 'num'
    model%output_vars(7)%value => model%num
    model%output_vars(8)%name = 'nuh'
    model%output_vars(8)%value => model%nuh
    model%output_vars(9)%name = 'NN'
    model%output_vars(9)%value => model%NN
    model%output_vars(10)%name = 'P'
    model%output_vars(10)%value => model%P
    model%output_vars(11)%name = 'B'
    model%output_vars(11)%value => model%B

    if(model%write_on_the_fly) then
      allocate(SimstratOnTheFlyThinningLogger :: model%output_logger)
    else
      allocate(SimstratMemoryBoundThinningLogger :: model%output_logger)
    end if

    !input file post-processing
    call initializeAdditionalModelParameter(model)

  end function initialize_simstrat_model_from_inputfiles

  !#######################################################################
  subroutine parseJSONInputfile(model, ParName)
  !#######################################################################
    implicit none
    class(SimstratModel), intent(inout) :: model
    character(kind=CK, len=*), intent(in) :: ParName

    type(json_file) :: par_file
    logical :: found

    !gfortran cannot handle type bound allocatable character that are passed to subroutine as intent(out)
    !as a workaround we have to store the values in a local scope allocatable character
    character(kind=CK, len=:), allocatable          :: MorphName,InitName,ForcingName,AbsorpName
    character(kind=CK, len=:), allocatable          :: GridName,zoutName,toutName,PathOut
    character(kind=CK, len=:), allocatable          :: QinpName,QoutName,TinpName,SinpName
    character(kind=CK, len=:), allocatable          :: stencil_type

    model%ParName = ParName

    call par_file%initialize()

    !load file or stop if fail
    call par_file%load_file(filename = ParName)
    if(par_file%failed()) then
      call error('Could not read inputfile '//ParName)
      stop
    end if

    !Names of Inputfile
    call par_file%get('Input.Morphology'        , MorphName, found); model%MorphName = MorphName; call check_field(found, 'Input.Morphology'      , ParName)
    call par_file%get('Input.Initial conditions', InitName, found); model%InitName = InitName; call check_field(found, 'Input.Initial conditions' , ParName)
    call par_file%get('Input.Forcing'           , ForcingName, found); model%ForcingName = ForcingName; call check_field(found, 'Input.Forcing'   , ParName)
    call par_file%get('Input.Absorption'        , AbsorpName, found); model%AbsorpName = AbsorpName; call check_field(found, 'Input.Absorption'   , ParName)
    call par_file%get('Input.Grid'              , GridName, found); model%GridName = GridName; call check_field(found, 'Input.Grid'               , ParName)
    call par_file%get('Input.Inflow'            , QinpName, found); model%QinpName = QinpName; call check_field(found, 'Input.Inflow'             , ParName)
    call par_file%get('Input.Outflow'           , QoutName, found); model%QoutName = QoutName; call check_field(found, 'Input.Outflow'            , ParName)
    call par_file%get('Input.TemperatureInflow' , TinpName, found); model%TinpName = TinpName; call check_field(found, 'Input.TemperatureInflow'  , ParName)
    call par_file%get('Input.SalinityInflow'    , SinpName, found); model%SinpName = SinpName; call check_field(found, 'Input.SalinityInflow'     , ParName)

    !Output
    call par_file%get('Output.path'             , PathOut, found); model%PathOut = PathOut; call check_field(found, 'Output.path'     , ParName)
    call par_file%get('Output.depth'            , zoutName, found); model%zoutName = zoutName; call check_field(found, 'Output.depth' , ParName)
    call par_file%get('Output.WriteOnTheFly'    , model%write_on_the_fly, found); call check_field(found, 'Output.WriteOnTheFly'  , ParName)
    call par_file%get('Output.ThinningInterval' , model%thinning_interval, found); call check_field(found, 'Output.ThinningInterval'  , ParName)

    !Model configuration
    call par_file%get("ModelConfig.StencilType"            , stencil_type   , found); model%stencil_type = stencil_type
    call par_file%get("ModelConfig.MaxNrGridCells"         , model%nz_max   , found);
    if(.not.found) then
      model%nz_max = 1000
      call warn('Variable "ModelConfig.MaxNrGridCells" is not set. Assume a value of 1000')
    end if
    call par_file%get("ModelConfig.CoupleAED2"             , model%couple_aed2   , found);
    if(.not.found) then
      model%couple_aed2 = .false.
      call warn('Variable "ModelConfig.CoupleAED2" is not set. Assume you do not want to couple simstrat with aed2.')
    end if
    call par_file%get("ModelConfig.TurbulenceModel"        , model%Mod         , found); call check_field(found, 'ModelConfig.TurbulenceModel'         , ParName)
    call par_file%get("ModelConfig.StabilityFunction"      , model%Stab        , found); call check_field(found, 'ModelConfig.StabilityFunction'       , ParName)
    call par_file%get("ModelConfig.FluxCondition"          , model%ModFlux     , found); call check_field(found, 'ModelConfig.FluxCondition'           , ParName)
    call par_file%get("ModelConfig.Forcing"                , model%ForcingType , found); call check_field(found, 'ModelConfig.Forcing'                 , ParName)
    call par_file%get("ModelConfig.UseFilteredWind"        , model%UseWindFilt , found); call check_field(found, 'ModelConfig.UseFilteredWind'         , ParName)
    call par_file%get("ModelConfig.SeicheNormalization"    , model%ModSNorm    , found); call check_field(found, 'ModelConfig.SeicheNormalization'     , ParName)
    call par_file%get("ModelConfig.Wind drag model"        , model%ModC10      , found); call check_field(found, 'ModelConfig.Wind drag model'         , ParName)
    call par_file%get("ModelConfig.InflowPlacement"        , model%ModInflow   , found); call check_field(found, 'ModelConfig.InflowPlacement'         , ParName)
    call par_file%get("ModelConfig.PressureGradients"      , model%Pgrad       , found); call check_field(found, 'ModelConfig.PressureGradients'       , ParName)
    call par_file%get("ModelConfig.EnableSalinityTransport", model%ModSal      , found); call check_field(found, 'ModelConfig.EnableSalinityTransport' , ParName)
    call par_file%get("ModelConfig.DisplaySimulation"      , model%disp_sim    , found); call check_field(found, 'ModelConfig.DisplaySimulation'       , ParName)
    call par_file%get("ModelConfig.DisplayDiagnose"        , model%disp_dgn    , found); call check_field(found, 'ModelConfig.DisplayDiagnose'         , ParName)
    call par_file%get("ModelConfig.DataAveraging"          , model%igoal       , found); call check_field(found, 'ModelConfig.DataAveraging'           , ParName)

    !Model Parameter
    call par_file%get("Parameter.Lat"       , model%Lat     , found); call check_field(found, 'ModelConfig.Lat'     , ParName)
    call par_file%get("Parameter.p_air"     , model%p_air   , found); call check_field(found, 'ModelConfig.p_air'   , ParName)
    call par_file%get("Parameter.a_seiche"  , model%a_seiche, found); call check_field(found, 'ModelConfig.a_seiche', ParName)
    call par_file%get("Parameter.q_NN"      , model%q_NN    , found); call check_field(found, 'ModelConfig.q_NN'    , ParName)
    call par_file%get("Parameter.f_wind"    , model%f_wind  , found); call check_field(found, 'ModelConfig.f_wind'  , ParName)
    call par_file%get("Parameter.C10"       , model%C10     , found); call check_field(found, 'ModelConfig.C10'     , ParName)
    call par_file%get("Parameter.CD"        , model%CD      , found); call check_field(found, 'ModelConfig.CD'      , ParName)
    call par_file%get("Parameter.fgeo"      , model%fgeo    , found); call check_field(found, 'ModelConfig.fgeo'    , ParName)
    call par_file%get("Parameter.k_min"     , model%k_min   , found); call check_field(found, 'ModelConfig.k_min'   , ParName)
    call par_file%get("Parameter.p_radin"   , model%p_radin , found); call check_field(found, 'ModelConfig.p_radin' , ParName)
    call par_file%get("Parameter.p_windf"   , model%p_windf , found); call check_field(found, 'ModelConfig.p_windf' , ParName)
    call par_file%get("Parameter.beta_sol"  , model%beta_sol, found); call check_field(found, 'ModelConfig.beta_sol', ParName)
    call par_file%get("Parameter.albsw"     , model%albsw   , found); call check_field(found, 'ModelConfig.albsw'   , ParName)

    !Simulation Parameter
    call par_file%get("Simulation.Timestep", model%dt     , found); call check_field(found, 'Simulation.Timestep' , ParName)
    call par_file%get("Simulation.Start"   , model%t_start, found); call check_field(found, 'Simulation.Start'    , ParName)
    call par_file%get("Simulation.End"     , model%t_end  , found); call check_field(found, 'Simulation.End'      , ParName)

    call par_file%destroy()


    !check validity of inputfile
    if(.not.(model%disp_sim==1 .or. model%disp_sim==2 .or. model%disp_sim==3)) model%disp_sim=0
    if(.not.(model%disp_dgn==1 .or. model%disp_dgn==2)) model%disp_dgn=0

    if(model%disp_dgn/=0) then
      call ok('Configuration: '//trim(ParName))
    end if
  end subroutine parseJSONInputfile

  subroutine check_field(found, field_name, file_name)
    implicit none

    logical, intent(in) :: found
    character(len=*), intent(in) :: field_name, file_name

    if(.not.found) then
      call error('Field '//field_name//' not found in '//file_name)
      stop
    end if
  end subroutine check_field


  subroutine initializeDiscretizationScheme(model)
    implicit none

    !arguments
    class(SimstratModel), intent(inout) :: model

    !CSV file reader
    type(csv_file) :: f
    logical :: status_ok
    logical :: status_ok_any = .true.
    !Grid
    real(RK), dimension(:), allocatable :: z_grid
    real(RK), dimension(:), allocatable :: z_faces
    !Morphology
    real(RK), dimension(:), allocatable :: z_morph, area_morph
    integer :: nz, nz_morph
    integer :: i_bottom, i_top

    !******************************
    !Load morphology
    !******************************
    call f%initialize(delimiter=char(9))
    call f%read(model%MorphName, header_row=1, status_ok=status_ok)
    if(.not. status_ok) then
      call error('Unable to read morphology file: '//model%MorphName)
      stop
    end if
    call f%get(1, z_morph, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(2, area_morph, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%destroy()
    if(.not.status_ok_any .or. size(z_morph)==0) then
      call error('Unable to read morphology file: '//model%GridName)
      stop
    end if
    if(z_morph(1) <= 0.0_RK) call warn('The morphology is not defined above the initial water level')
    nz_morph = size(z_morph)
    model%z_zero = z_morph(nz_morph)
    model%depth = z_morph(1)-z_morph(nz_morph)
    model%A_surf = area_morph(1)
    if(model%disp_dgn/=0) then
      call ok('Morphology: '//model%MorphName)
      write(*,*) '     z sediment: ', z_morph(nz_morph), ' depth: ', z_morph(1)-z_morph(nz_morph)
    end if


    !*********************************
    !Load grid
    !*********************************
    call f%read(model%GridName, header_row=1, status_ok=status_ok)
    if(.not. status_ok) then
      call error('Unable to read grid file: '//model%GridName)
      stop
    end if
    call f%get(1, z_grid, status_ok)
    call f%destroy()
    if(.not.status_ok .or. size(z_grid)==0) then
      call error('Unable to read grid file: '//model%GridName)
      stop
    end if

    !****************************************
    !Initialize position of volume faces
    !****************************************
    !if grid file only contains one value, this is the number of grid points
    !otherwise it defines the depths of the volume faces
    if(size(z_grid)==1) then
      nz = int(z_grid(1))
      !initialize equaly spaced grid from the bottom to the top of the morphology
      allocate(z_faces(nz))
      z_faces(1:nz) = linspace(0.0_RK, z_morph(1)-z_morph(nz_morph), nz, endpoint=.true.)
    else
      nz = size(z_grid)
      !check that z_grid does not extend over the boundary of the morphology
      if(z_grid(nz) > z_morph(1) .or. z_grid(1) < z_morph(nz_morph)) then
        call error('Grid does not overlap with morphology')
        stop
      end if
      i_bottom = nz
      do while(z_grid(i_bottom) < z_morph(nz_morph))
        i_bottom = i_bottom-1
      end do
      i_top = 1 
      do while(z_grid(i_top) > z_morph(1))
        i_top = i_top+1
      end do
      nz = i_bottom-i_top+1
      allocate(z_faces(nz+2))
      z_faces(1) = z_morph(1)
      !inlcude top of morphology if not included
      if(z_grid(i_top) < z_morph(1)) then
        nz = nz+1
        z_faces(2:nz) = z_grid(i_top:i_bottom)
      else
        z_faces(1:nz) = z_grid(i_top:i_bottom)
      end if
      !include top of morphology if not inlcuded
      if(z_grid(i_bottom) > z_morph(nz_morph)) then
        nz = nz+1
        z_faces(nz) = z_morph(nz_morph)
      end if
      !cut array
      z_faces = z_faces(1:nz)
      z_faces = -z_faces(nz:1:-1) !reverse order of values
      z_faces = z_faces(1)-z_faces !transform to hight above sediment (where z(1)=0 is at the sediment)
    end if
    if(nz > model%nz_max) then
      model%nz_max = nz+500
      call warn('Number of grid points defined in '//model%GridName//'is bigger than maximum nr of grid points defined in '//model%ParName//'. Assume a maximum number of grid points of nz+500.')
    end if
    if(model%disp_dgn/=0) then
      call ok('Grid: '//model%GridName)
      write(*,*) '     No grid points (i.e. volumes) : ',nz
    end if

    !assign discretization scheme
    if(associated(model%discretization)) deallocate(model%discretization)
    select case (model%stencil_type)
      case('crank nicolson')
        allocate(StaggeredFVCrankNicolsonScheme :: model%discretization)
      case('implicit euler')
        allocate(StaggeredFVImplicitEulerScheme :: model%discretization)
      case default
        call error('Unknown stencil "'//model%stencil_type//'" defined in inputfile')
        stop
    end select

    !initialize discretization
    !transform z_morph to height above sediment
    z_morph = -z_morph(nz_morph:1:-1) !reverse order of values
    z_morph = z_morph(1)-z_morph !transform to hight above sediment (where z(1)=0 is at the sediment)
    area_morph = area_morph(size(area_morph):1:-1) !reverse oder of values
    call model%discretization%initialize(model%nz_max, z_faces, z_morph, area_morph)
    model%volume = model%discretization%totalvolume()
    if(model%disp_dgn/=0) then
      call ok('Successfuly set up discretization scheme')
      write(*,*) '     z_faces(0): ', z_faces(1), ' z_faces(nz): ', z_faces(nz)
    end if
  end subroutine initializeDiscretizationScheme


  subroutine initializeInitialConditions(model)
    implicit none

    class(SimstratModel), intent(inout) :: model

    !CSV file reader
    type(csv_file) :: f
    logical :: status_ok
    logical :: status_ok_any = .true.
    !Grid
    real(RK), dimension(:), allocatable :: depth, u_ini, v_ini, T_ini, S_ini, k_ini, eps_ini
    type(SimstratStateVariable), pointer :: T, S, u, v
    type(SimstratStateVariable_k), pointer :: k
    type(SimstratStateVariable_eps), pointer :: eps
    integer :: nz_init

    !******************************
    !Load inital conditions
    !******************************
    call f%initialize(delimiter=char(9))
    call f%read(model%InitName, header_row=1, status_ok=status_ok)
    if(.not. status_ok) then
      call error('Unable to read initial conditions: '//model%InitName)
      stop
    end if
    call f%get(1, depth, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(2, u_ini, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(3, v_ini, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(4, T_ini, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(5, S_ini, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(6, k_ini, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(7, eps_ini, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%destroy()
    if(.not.status_ok_any .or. size(depth)==0) then
      call error('Unable to read initial conditions: '//model%InitName)
    end if
    !reverse order and convert to hight above sediment
    nz_init = size(depth)
    depth = depth(nz_init:1:-1)-model%z_zero
    u_ini = u_ini(nz_init:1:-1)
    v_ini = v_ini(nz_init:1:-1)
    T_ini = T_ini(nz_init:1:-1)
    S_ini = S_ini(nz_init:1:-1)
    k_ini = k_ini(nz_init:1:-1)
    eps_ini = eps_ini(nz_init:1:-1)

    !initialize state variables
    allocate(model%T(model%nz_max),&
             model%S(model%nz_max),&
             model%u(model%nz_max),&
             model%v(model%nz_max),&
             model%k(model%nz_max),&
             model%eps(model%nz_max))
    allocate(model%state_vars(6))
    !Temperature
    allocate(SimstratStateVariable :: T)
    T%state_name='T'
    T%state_value => model%T

    !Salinity
    allocate(SimstratStateVariable :: S)
    S%state_name='S'
    S%state_value => model%S

    !u
    allocate(SimstratStateVariable :: u)
    u%state_name='u'
    u%state_value => model%u

    !v
    allocate(SimstratStateVariable :: v)
    v%state_name='v'
    v%state_value => model%v

    !k
    allocate(SimstratStateVariable_k :: k)
    k%state_name='k'
    k%state_value => model%k

    !eps
    allocate(SimstratStateVariable_eps :: eps)
    eps%state_name='eps'
    eps%state_value => model%eps

    !assemble array of state variables
    model%state_vars(1)%ptr => T
    model%state_vars(2)%ptr => S
    model%state_vars(3)%ptr => u
    model%state_vars(4)%ptr => v
    model%state_vars(5)%ptr => k
    model%state_vars(6)%ptr => eps

    !set initial conditions
    select type(discretization => model%discretization)
      class is (StaggeredFiniteVolumeDiscretization)
        associate(nz => discretization%nz_mfq)
          T%state_value(1:nz) = discretization%interpolateOnVolumeCentres(depth, T_ini)
          S%state_value(1:nz) = discretization%interpolateOnVolumeCentres(depth, S_ini)
          u%state_value(1:nz) = discretization%interpolateOnVolumeCentres(depth, u_ini)
          v%state_value(1:nz) = discretization%interpolateOnVolumeCentres(depth, v_ini)
          k%state_value(1:nz+1) = discretization%interpolateOnVolumeFaces(depth, k_ini)
          eps%state_value(1:nz+1) = discretization%interpolateOnVolumeFaces(depth, eps_ini)
        end associate
    end select

    model%num = 0.0_RK
    model%nuh = 0.0_RK
    model%tx = 0.0_RK
    model%ty = 0.0_RK
    model%L = 0.2_RK

    if(model%disp_dgn/=0) then
      call ok('Initial Conditions: '//model%InitName)
    end if
  end subroutine initializeInitialConditions

  subroutine initializeAdditionalModelParameter(model)
    implicit none
    class(SimstratModel), intent(inout) :: model

    ! Calculate Coriolis parameter from latitude
    model%Cori=2.0_RK*7.292e-5_RK*sin(model%Lat*pi/180)
    if(model%Stab==1) model%cm0 = 0.5625_RK
    if(model%Stab==2) model%cm0 = 0.556171_RK
    model%cde = model%cm0**3
    sig_e=(kappa/model%cm0)**2/(ce2-ce1)
    model%gamma = model%A_surf/(model%volume**1.5_RK)/sqrt(rho_0)*model%CD

    if (model%ModSal/=0) then
      model%salctr=1
      model%delsal=1
    else
      if(any(model%S/=0)) model%salctr=1
      if (model%salctr==1) then
        if(any(model%S /= model%S(0))) model%delsal=1
      end if
    end if
  end subroutine

end module simstrat_inputfile_module
