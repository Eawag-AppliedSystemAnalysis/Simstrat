module simstrat_forcing_module
  use simstrat_kinds
  use simstrat_model_module, only: SimstratModel, SimstratForcing, SimstratStateVariable
  use utilities
  use simstrat_model_constants
  use csv_module
  implicit none

  private

  type, extends(SimstratForcing), public :: MemoryBoundSimstratForcing

    contains
      procedure, pass(self), public :: initialize => initialize_forcing
      procedure, nopass, public :: interpolate_forcing_at_date => memory_bound_interpolate
      procedure, nopass, public :: process_forcing_at_date
  end type

  !Implement a buffered version

contains
  !Implement different types of forcing files
  subroutine initialize_forcing(self, ForcingName, forcing_type, use_wind_filt, disp_dgn, AbsorpName, z_zero)
    implicit none

    class(MemoryBoundSimstratForcing), intent(inout) :: self
    character(len=:), intent(in), allocatable :: ForcingName, AbsorpName
    integer, intent(in) :: forcing_type
    logical, intent(in) :: use_wind_filt
    integer, intent(in) :: disp_dgn
    !type(SimstratStateVariable), pointer, intent(in) :: T
    real(RK), intent(in) :: z_zero

    call check_file_exists(ForcingName)
    call check_file_exists(AbsorpName)

    call load_forcing_inputfile(self, ForcingName, forcing_type, use_wind_filt, disp_dgn)
    call load_absorption_inputfile(self, AbsorpName, disp_dgn, z_zero)
  end subroutine

  subroutine load_forcing_inputfile(self, ForcingName, forcing_type, use_wind_filt, disp_dgn)
    class(MemoryBoundSimstratForcing), intent(inout) :: self
    character(len=:), intent(in), allocatable :: ForcingName
    integer, intent(in) :: forcing_type
    logical, intent(in) :: use_wind_filt
    integer, intent(in) :: disp_dgn
    !type(SimstratStateVariable), pointer, intent(in) :: T

    !CSV file reader
    type(csv_file) :: f
    logical :: status_ok
    logical :: status_ok_any = .true.
    
    if(forcing_type < 1 .or. forcing_type > 4) then
      call error('Forcing type must be 1, 2, 3 or 4.')
      stop
    end if

    self%forcing_type = forcing_type
    self%use_wind_filt = use_wind_filt

    !******************************
    !Load forcing inputfile
    !******************************
    call f%initialize(delimiter=char(9))
    call f%read(ForcingName, header_row=1, status_ok=status_ok)
    if(.not. status_ok) then
      call error('Unable to read forcing file: '//ForcingName)
      stop
    end if
    call f%get(1, self%datum_forcing, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(2, self%u, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(3, self%v, status_ok); status_ok_any = status_ok_any .or. status_ok
    call f%get(5, self%Fsol, status_ok); status_ok_any = status_ok_any .or. status_ok
    select case(forcing_type)
      case(1)
        call f%get(4, self%SST, status_ok); status_ok_any = status_ok_any .or. status_ok
        if(use_wind_filt) call f%get(6, self%Wf, status_ok); status_ok_any = status_ok_any .or. status_ok
      case(2)
        call f%get(4, self%Tair, status_ok); status_ok_any = status_ok_any .or. status_ok
        call f%get(6, self%vap, status_ok); status_ok_any = status_ok_any .or. status_ok
        if(use_wind_filt) call f%get(6, self%Wf, status_ok); status_ok_any = status_ok_any .or. status_ok
      case(3)
        call f%get(4, self%Tair, status_ok); status_ok_any = status_ok_any .or. status_ok
        call f%get(6, self%vap, status_ok); status_ok_any = status_ok_any .or. status_ok
        call f%get(7, self%cloud_coverage, status_ok); status_ok_any = status_ok_any .or. status_ok
        if(use_wind_filt) call f%get(8, self%Wf, status_ok); status_ok_any = status_ok_any .or. status_ok
      case(4)
        call f%get(4, self%Hnet, status_ok); status_ok_any = status_ok_any .or. status_ok
        if(use_wind_filt) call f%get(6, self%Wf, status_ok); status_ok_any = status_ok_any .or. status_ok
    end select
    call f%destroy()
    if(.not.status_ok_any .or. size(self%datum_forcing)==0) then
      call error('Unable to read forcing file: '//ForcingName)
      stop
    end if
    self%n_forcing = size(self%datum_forcing)
    if(disp_dgn/=0) then
      call ok('Forcing: '//ForcingName)
    end if
  end subroutine

  subroutine load_absorption_inputfile(self, AbsorpName, disp_dgn, z_zero)
    class(MemoryBoundSimstratForcing), intent(inout) :: self
    character(len=:), intent(in), allocatable :: AbsorpName
    integer, intent(in) :: disp_dgn
    real(RK), intent(in) :: z_zero

    !CSV file reader
    type(csv_file) :: f
    logical :: status_ok
    logical :: status_ok_any = .true.
    character(len=30), dimension(:), allocatable :: header
    integer :: i,j
    real(RK), dimension(:), allocatable :: col

    !******************************
    !Load absorption inputfile
    !******************************
    call f%initialize(delimiter=char(9))
    call f%read(AbsorpName, header_row=3, status_ok=status_ok, skip_rows=(/1,2/))

    if(.not. status_ok) then
      call error('Unable to read forcing file: '//AbsorpName)
      stop
    end if
    call f%get_header(header, status_ok); status_ok_any = status_ok_any .or. status_ok
    !convert header: remove first element and parse to real
    self%n_depth_absorption = size(header)-1
    allocate(self%depth_absorption(self%n_depth_absorption))
    do i=1,self%n_depth_absorption
      read(header(i+1), *) self%depth_absorption(i)
    end do
    !transform depth to hight above sediment (reverse order and make sediment z=0)
    self%depth_absorption = convert2hight_above_sed(self%depth_absorption, z_zero) !-z_zero+self%depth_absorption(self%n_depth_absorption:1:-1)

    !read dates
    call f%get(1, self%datum_absorption, status_ok); status_ok_any = status_ok_any .or. status_ok
    self%n_absorption = size(self%datum_absorption)

    !read absorption values
    allocate(self%absorption(self%n_absorption, self%n_depth_absorption))
    allocate(col(self%n_absorption))
    do i=1,self%n_depth_absorption
      call f%get(i+1, col, status_ok); status_ok_any = status_ok_any .or. status_ok
      self%absorption(:,self%n_depth_absorption-i+1) = col !write in reverse order
    end do
    allocate(self%absorption_datum(self%n_depth_absorption))
    call f%destroy()
    if(.not.status_ok_any .or. size(self%datum_absorption)==0) then
      call error('Unable to read forcing file: '//AbsorpName)
      stop
    end if

    if(disp_dgn/=0) then
      call ok('Absorption: '//AbsorpName)
      write(*,*) '     Depths: ', self%depth_absorption
      write(*,*) '     Time points: ', self%n_absorption
    end if
  end subroutine

  subroutine memory_bound_interpolate(datum, model)
    implicit none
    real(RK), intent(in) :: datum
    class(SimstratModel), intent(inout) :: model

    integer :: idx

    associate(self => model%forcing)
      idx = find_index_ordered(self%datum_forcing, datum)
      if(self%datum_forcing(idx) < datum) then
        self%u_datum = interpolate_forcing(self%datum_forcing, self%u, idx, datum)
        self%v_datum = interpolate_forcing(self%datum_forcing, self%v, idx, datum)
        self%Fsol_datum = interpolate_forcing(self%datum_forcing, self%Fsol, idx, datum)
        if(self%use_wind_filt) self%Wf_datum = interpolate_forcing(self%datum_forcing, self%Wf, idx, datum)
        select case(self%forcing_type)
          case(1)
            self%SST_datum = interpolate_forcing(self%datum_forcing, self%SST, idx, datum)
          case(2)
            self%Tair_datum = interpolate_forcing(self%datum_forcing, self%Tair, idx, datum)
            self%vap_datum = interpolate_forcing(self%datum_forcing, self%vap, idx, datum)
          case(3)
            self%Tair_datum = interpolate_forcing(self%datum_forcing, self%Tair, idx, datum)
            self%vap_datum = interpolate_forcing(self%datum_forcing, self%vap, idx, datum)
            self%cloud_coverage_datum = interpolate_forcing(self%datum_forcing, self%cloud_coverage, idx, datum)
          case(4)
            self%Hnet = interpolate_forcing(self%datum_forcing, self%Hnet, idx, datum)
        end select
      else
        self%u_datum = self%u(idx)
        self%v_datum = self%v(idx)
        self%Fsol_datum = self%Fsol(idx)
        if(self%use_wind_filt) self%Wf_datum = self%Wf(idx)

        select case(self%forcing_type)
          case(1)
            self%SST_datum = self%SST(idx)
          case(2)
            self%Tair_datum = self%Tair(idx)
            self%vap_datum = self%vap(idx)
          case(3)
            self%Tair_datum = self%Tair(idx)
            self%vap_datum = self%vap(idx)
            self%cloud_coverage_datum = self%cloud_coverage(idx)
          case(4)
            self%Hnet = self%Hnet(idx)
        end select
      end if

      idx = find_index_ordered(self%datum_absorption, datum)
      if(self%datum_absorption(idx) < datum) then
        if(.not.allocated(self%absorption_datum)) allocate(self%absorption_datum(self%n_depth_absorption))
        self%absorption_datum = interpolate_absorption(self%datum_absorption, self%absorption, idx, datum)
      else
        self%absorption_datum = self%absorption(idx,:)
      end if
    end associate


  end subroutine


  subroutine process_forcing_at_date(datum, model)
    implicit none
    real(RK), intent(in) :: datum
    class(SimstratModel), intent(inout) :: model

    real(RK) :: fu, Vap_wat, heat0
    real(RK) :: F_glob, Cloud
    real(RK) :: H_A, H_K, H_V, H_W
    real(RK) :: Tnz

    Tnz = model%T(model%discretization%nz_mfq)

    call model%forcing%interpolate_forcing_at_date(datum, model)

    associate(u10=>model%u10,v10=>model%v10,uv10=>model%uv10,Wf=>model%Wf,& !Absorption coeff [m-1], wind drag, wind speeds
              u_taub=>model%u_taub,drag=>model%drag,u_taus=>model%u_taus,&
              tx=>model%tx,ty=>model%ty,SST=>model%SST,heat=>model%heat,& !Shear stress, sea surface temperature and heat flux
              rad0=>model%rad0,& !Solar radiation at surface
              p_air=>model%p_air, f_wind=>model%f_wind, C10=>model%C10, CD=>model%CD, p_radin=>model%p_radin, p_windf=>model%p_windf, beta_sol=>model%beta_sol, albsw=>model%albsw,&
              U=>model%u, V=>model%v,&
              ModC10=>model%ModC10,&
              self=>model%forcing,&
              forcing_type=>model%forcing%forcing_type, WindFilt=>model%forcing%use_wind_filt,&
              nz=>model%discretization%nz_mfq, ga1=>model%ga1)

        F_glob = self%Fsol_datum*(1-albsw)
        u10 = self%u_datum*f_wind      !MS 2014: added f_wind
        v10 = self%v_datum*f_wind      !MS 2014: added f_wind
        uv10 = sqrt(u10**2+v10**2)  !AG 2014
        if(WindFilt) Wf = self%Wf_datum

        if (forcing_type==1) then
            SST = self%SST_datum !Sea surface temperature
            rad0 = self%Fsol_datum*(1-albsw)*(1-beta_sol) ! MS: added beta_sol and albsw
            heat = 0.0_RK
        else if (forcing_type>1) then
            if (forcing_type==2) then            ! date, U,V,Tatm,Hsol,Vap
                Cloud = 0.5_RK
            else if (forcing_type==3) then       ! date,U10,V10,Tatm,Hsol,Vap,Clouds
                Cloud = self%cloud_coverage_datum
                if (Cloud<0 .or. Cloud>1) then
                    call error('Cloudiness should always be between 0 and 1.')
                    stop
                end if
            else if (forcing_type==4) then       ! date,U10,V10,Hnet,Hsol
                heat0 = self%Hnet_datum
            end if

            if (forcing_type/=4) then ! in the water column

                ! Wind function (Livingstone & Imboden 1989)
                fu = 4.40_RK+1.82_RK*uv10+0.26_RK*(Tnz-self%Tair_datum)
                fu = fu*p_windf    ! Provided fitting factor p_windf (~1)
                ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
                Vap_wat = 10**((0.7859_RK+0.03477_RK*Tnz)/(1+0.00412_RK*Tnz))
                Vap_wat = Vap_wat*(1+1e-6_RK*p_air*(4.5_RK+0.00006_RK*Tnz**2))
                ! Solar short-wave radiation absorbed
                rad0 = F_glob*(1-beta_sol) !MS: added beta_sol

                ! Long-wave radiation from sky (Livingstone & Imboden 1989)
                ! H_A = 1.24*sig*(1-r_a)*(1+0.17*Cloud**2)*(Vap_atm/(273.15+T_atm))**(1./7)*(273.15+T_atm)**4
                ! Long-wave radiation according to Dilley and O'Brien
                ! see Flerchinger et al. (2009)
                H_A = (1-r_a) * ((1-0.84_RK*Cloud) * (59.38_RK + 113.7_RK*((self%Tair_datum+273.15_RK)/273.16_RK)**6&
                &   + 96.96_RK*sqrt(465*self%vap_datum/(self%Tair_datum+273.15_RK)*0.04_RK)) / 5.67e-8_RK / &
                &   (self%Tair_datum+273.15_RK)**4 + 0.84_RK * Cloud) *5.67e-8_RK*(self%Tair_datum+273.15_RK)**4
                H_A = H_A*p_radin    ! Provided fitting factor p_radin (~1)
                ! Long-wave radiation from water body (black body)
                H_W = -0.97_RK*sig*(Tnz+273.15_RK)**4
                ! Flux of sensible heat (convection)
                H_K = -B0*fu*(Tnz-self%Tair_datum)
                ! Flux of latent heat (evaporation, condensation)
                H_V = -fu*(Vap_wat-self%vap_datum)
                ! Global heat flux (positive: air to water, negative: water to air)
                heat = H_A + H_W + H_K + H_V + F_glob*beta_sol !MS: added term with beta_sol
            else
                heat = heat0 + F_glob*beta_sol !MS: added term with beta_sol
            end if
            if ((Tnz<0).and.(heat<0)) heat=0.0_RK
        end if

        !Drag coefficient as a function of wind speed (AG 2014)
        if (ModC10==2) then !Ocean model
            C10 = -0.000000712_RK*uv10**2+0.00007387_RK*uv10+0.0006605_RK
        else if (ModC10==3) then !Lake model (Wüest and Lorke 2003)
            if (uv10<=0.1) then
                C10 = 0.06215_RK
            else if (uv10<=3.85_RK) then
                C10 = 0.0044_RK*uv10**(-1.15_RK)
            else !Polynomial approximation of Charnock's law
                C10 = -0.000000712_RK*uv10**2+0.00007387_RK*uv10+0.0006605_RK
            end if
        end if

      u_taus = sqrt(C10*rho_air/rho_0*uv10**2)
      u_taub=sqrt(drag*(U(1)**2+V(1)**2))

      tx = C10*rho_air/rho_0*uv10*u10
      ty = C10*rho_air/rho_0*uv10*v10

      ga1(1:nz) = model%discretization%interpolateMFQ(self%depth_absorption, self%absorption_datum)
    end associate
  end subroutine

  pure function interpolate_forcing(datum_array, val_array, idx, datum) result(val)
    implicit none

    real(RK), dimension(:), allocatable, intent(in) :: datum_array, val_array
    integer, intent(in) :: idx
    real(RK), intent(in) :: datum

    real(RK) :: val

    val = linear_interpolate(datum_array(idx), datum_array(idx+1), val_array(idx), val_array(idx+1), datum)
  end function

  pure function interpolate_absorption(datum_array, val_matrix, idx, datum) result(val)
    implicit none

    real(RK), dimension(:), allocatable, intent(in) :: datum_array
    real(RK), dimension(:,:), allocatable, intent(in) :: val_matrix
    integer, intent(in) :: idx
    real(RK), intent(in) :: datum

    real(RK), dimension(:), allocatable :: val
    integer :: i,n_cols

    n_cols = size(val_matrix,dim=2)
    allocate(val(n_cols))

    do i=1,n_cols
      val(i) = linear_interpolate(datum_array(idx), datum_array(idx+1), val_matrix(idx, i), val_matrix(idx+1, i), datum)
    end do
  end function

end module
