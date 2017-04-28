module simstrat_inputfile_module
  use strat_kinds
  use strat_simdata
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

    !Read initial data
    call self%read_initial_data


  end subroutine initialize_model

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
