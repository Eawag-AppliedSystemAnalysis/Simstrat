module simstrat_fabm
	use fabm
    use fabm_types
    use simstrat_type
    use SimstratModel
    use utilities
    use Transport, only: doTransportEquation
    use Advection, only: doAreaFactors
    implicit none
    save

    !Declare FABM model
    type (type_model),pointer :: model
    !Create arrays for FABM biogeochemical state variables
    real(dp),allocatable,dimension(:,:),target :: state
    real(dp),allocatable,dimension(:),target :: bottom_state, surface_state
    real(dp),allocatable,dimension(:),target :: iniFABM



    !FABM variables
    integer :: m,mb,ms,mc
    logical :: valid
    real(dp),allocatable :: dy(:,:), dyS(:) !Temporal derivatives of pelagic state variables
    real(dp),allocatable :: sflux(:), bflux(:,:) !Flux of surface and bottom state variables
    real(dp),allocatable :: vel(:,:) !Velocity (sinking-/floating+) of pelagic state variables, independently of medium
    real(dp),allocatable :: ext(:)
    real(dp),allocatable :: quant(:,:), hquant(:) !Conserved quantities
    real(dp),allocatable :: Qstate(:,:) !Flow of pelagic state variables

    real(dp),allocatable :: rad0_ps,rad_ps(:) !Photosynthetically active radiation for FABM
    real(dp),allocatable :: depths(:),ht(:) !Depth and thickness for FABM

    real(dp) :: CD_stress

    private
    public initializeFABM, doFABM, Lateral_FABM, Lateral_FABM_rho, doAdvection_FABM, write_text_FABM,model,Qstate,state,m

contains

    subroutine initializeFABM(rho, absorb,S,T,uv10)
        implicit none

        real(dp), intent(inout) :: rho(0:nz_grid), absorb(0:nz_grid)
        real(dp), intent(inout) :: S(0:nz_grid), T(0:nz_grid)
        real(dp), intent(inout) :: uv10

        integer :: i
        !integer:: bottom_domain(0:nz)
        !Initialize FABM model from namelist 'fabm.nml'
        model => fabm_create_model_from_file(get_free_unit())
        !Initialize FABM model from yaml file
        !call fabm_create_model_from_yaml_file(model)
        write(6,*)

        !Set domain size
        call fabm_set_domain(model,nz_grid+1)
        ! Specify vertical index of surface and bottom
        call model%set_surface_index(1)
        call model%set_bottom_index(nz)

        m = size(model%state_variables)
        mb = size(model%bottom_state_variables)
        ms = size(model%surface_state_variables)
        mc = size(model%conserved_quantities)

        !Allocate arrays
        allocate(state(0:nz_grid,1:m))
        allocate(surface_state(1:ms))
        allocate(bottom_state(1:mb))
        allocate(dy(0:nz_grid,1:m))
        allocate(dyS(1:m))
        allocate(sflux(1:ms))
        allocate(bflux(0:nz_grid,1:mb))
        !allocate(pp(0:nz_grid,1:m,1:m))
        !allocate(dd(0:nz_grid,1:m,1:m))
        allocate(vel(0:nz_grid,1:m))
        allocate(ext(0:nz_grid))
        allocate(Qstate(0:nz_grid,1:m))
        allocate(quant(0:nz_grid,1:mc))
        allocate(hquant(0:nz_grid))

        allocate(rad_ps(0:nz_grid))
        allocate(depths(0:nz_grid))
        allocate(ht(0:nz_grid))

        !Link environmental data to FABM
        !call fabm_link_bulk_data(model,standard_variables%depth,depths) ![m] !useless?
        !call fabm_link_bulk_data(model,standard_variables%density,rho) ![kg/m3] !useless?
        call fabm_link_bulk_data(model,standard_variables%cell_thickness,ht) ![m]
        call fabm_link_bulk_data(model,standard_variables%temperature,T) ![°C]
        call fabm_link_bulk_data(model,standard_variables%practical_salinity,S)
        !call fabm_link_bulk_data(model,standard_variables%downwelling_shortwave_flux,rad) !useless? inconsistent units!! (must give rad*rho_0*cp->another variable)
        call fabm_link_bulk_data(model,standard_variables%downwelling_photosynthetic_radiative_flux,rad_ps) ![W/m2]
        call fabm_link_bulk_data(model,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,absorb) ![m-1]
        !call fabm_link_horizontal_data(model,standard_variables%latitude,Lat) !useless? Lat not available here
        call fabm_link_horizontal_data(model,standard_variables%bottom_stress,CD_stress) ![Pa]
        call fabm_link_horizontal_data(model,standard_variables%wind_speed,uv10) ![m/s]
        call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,rad0_ps) ![W/m2]
        !call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_shortwave_flux,rad0) ![W/m2] !useless?
        !Pelagic state variables: set initial conditions and link data
        if(disp_dgn/=0) write(6,*) 'Initial conditions for biogeochemical model'
        allocate(iniFABM(0:nz_max))
        do i=1,m
            call InitCond_FABM(iniFABM,model%state_variables(i)%name,model%state_variables(i)%initial_value) !Look for file for each variable
            state(:,i) = iniFABM(0:nz_grid)
            call fabm_link_bulk_state_data(model,i,state(:,i))
            vel(:,i) = model%state_variables(i)%vertical_movement
        end do
        write(6,*)
        !Surface and bottom state variables: set initial conditions and link data
        do i=1,ms
            surface_state(i) = model%surface_state_variables(i)%initial_value
            call fabm_link_surface_state_data(model,i,surface_state(i))
        end do
        do i=1,mb
            bottom_state(i) = model%bottom_state_variables(i)%initial_value
            call fabm_link_bottom_state_data(model,i,bottom_state(i))
        end do
        call fabm_check_ready(model)
        if (disp_dgn/=0) then
            write(6,'(I3,A24)') m,'pelagic state variables'
            do i=1,m
                write(6,*) '   ',trim(model%state_variables(i)%name),': ',trim(model%state_variables(i)%long_name),&
                ' [',trim(model%state_variables(i)%units),']'
            end do
            write(6,'(I3,A24)') ms,'surface state variables'
            do i=1,ms
                write(6,*) '   ',trim(model%surface_state_variables(i)%name),': ',trim(model%surface_state_variables(i)%long_name),&
                ' [',trim(model%surface_state_variables(i)%units),']'
            end do
            write(6,'(I3,A23)') mb,'bottom state variables'
            do i=1,mb
                write(6,*) '   ',trim(model%bottom_state_variables(i)%name),': ',trim(model%bottom_state_variables(i)%long_name),&
                ' [',trim(model%bottom_state_variables(i)%units),']'
            end do
            write(6,'(I3,A21)') mc,'conserved quantities'
            do i=1,mc
                write(6,*) '   ',trim(model%conserved_quantities(i)%name),': ',trim(model%conserved_quantities(i)%long_name),&
                ' [',trim(model%conserved_quantities(i)%units),']'
            end do
            write(6,*)
        end if
    end subroutine initializeFABM

    !Returns the initial conditions vector iniFABM for a given FABM variable varname:
    !first looking for the corresponding file, secondly (if the file is not found)
    !by assigning the default value inival to the entire water column.
    !####################################################################
    subroutine InitCond_FABM(iniFABM,varname,inival)
    !####################################################################
        implicit none

        real(dp), intent(inout) :: iniFABM(0:nz_max) !Vector if initial conditions
        real(dp), intent(in) :: inival !Depth-independent value (default from fabm.nml)
        character*100, intent(in) :: varname !Identifying the variable

        real(dp) :: iniFABM_tmp(0:nz_max)
        real(dp) :: z_tmp(0:nz_max),z_tmp(0:nz_max)
        real(dp) :: z_tmp_depth
        character*100 :: fname
        integer :: i,nval

        fname = trim(PathFABM)//trim(varname)//'_ini.dat'
        open(14,action='read',status='unknown',err=1,file=fname)       ! Opens initial conditions file
        if(disp_dgn/=0) write(6,*) '   ',trim(fname)
        read(14,*)                                ! Skip header
        do i=0,nz_max                               ! Read initial values
            read(14,*,end=9) z_tmp(i),iniFABM(i)
        end do
    9   nval = i-1                                ! Number of values
        if (nval<0) then
            write(6,*) 'Error reading initial conditions file (no data found).'
            stop
        end if
        close(14)
        do i=0,nval
            z_tmp(i) = abs(z_tmp(i))               ! Make depths positive
        end do
        z_tmp_depth = z_tmp(0)                     ! Initial depth (top-most)

        do i=0,nval
            z_tmp(nval-i) = z_zero - z_tmp(i)
            iniFABM_tmp(nval-i) = iniFABM(i)
        end do

        if (nval==0) then
            write(6,*) '      Only one row! Water column will be initially homogeneous.'
            iniFABM(0:nz_grid) = iniFABM_tmp(0)
        else
            call Interp(z_tmp, iniFABM_tmp(:), nval, z_cent, iniFABM(:), nz_grid)
        end if
        return

    1   write(6,*) '   File ''',trim(fname),''' not found. Initial conditions set to default value from file ''fabm.nml''.'
        iniFABM(0:nz_grid) = inival !File not found: value from fabm.nml (constant)
        return
    end subroutine InitCond_FABM

    !Write output text files for biochemical variables
    subroutine write_text_FABM(datum)
        implicit none

        real(dp), intent(in) :: datum

        integer :: i,j
        real(dp) :: zext(0:nz-1),st(0:nsave)

        !External depths: bottom, layer centers, surface
        zext(0) = z_upp(0)
        zext(1:nz-2) = z_cent(2:nz-1)
        zext(nz-1) = z_upp(nz)

        do j=1,m
            write(100+j,*)
            write(100+j,'(F10.4,$)') datum
            call Interp_nan(zext,state(:,j),nz-1,zsave,st,nsave)
            do i=0,nsave
                write(100+j,'(ES12.4,$)') st(i)
            end do
            write(100+j,'(ES12.4,$)') state(nz,j)
        end do

        return
    end subroutine write_text_FABM

    !Returns the inflow vector Qstate for a given FABM variable varname
    !by reading the corresponding file (if present)
    !####################################################################
    subroutine Lateral_FABM(datum,idx,varname,i)
    !####################################################################
          implicit none

          ! Global Declarations
          integer, intent(in) :: idx,i
          real(dp), intent(in) :: datum
          character*100, intent(in) :: varname

          character*100 :: fname

          ! Local Declarations
          integer :: j, num_z(1:50), eof(1:50)
          real(dp) :: z_Inp(1:50,nz_max), dummy
          real(dp) :: Inp_read_start(1:50,nz_max), Inp_read_end(1:50,nz_max)
          real(dp) :: Q_start(1:50,nz_max), Q_end(1:50,nz_max)
          real(dp) :: Q_read_start(1:50,nz_max), Q_read_end(1:50,nz_max)
          real(dp) :: tb_start(1:50), tb_end(1:50)

          save Q_start, Q_end         ! MS 2014
          save num_z, eof       ! MS 2014
          save tb_start, tb_end       ! MS 2014

          if(trim(varname)=='aed_carbon_pH') goto 10

          if (idx==1) then
              if(disp_dgn==2) write(6,*) 'Starting to read '//trim(varname)//' inflow file...'
              fname = trim(PathFABM)//trim(varname)//'_inflow.dat'
              open(44+i,action='read',status='unknown',err=1,file=fname)
              eof(i) = 0
              read(44+i,*,end=9)      ! Skip first row: description of columns
              !Read input depths and convert coordinate system
              read(44+i,*,end=9) num_z(i)
              read(44+i,*,end=9) dummy, (z_Inp(i,j),j=1,num_z(i))
              z_Inp(i,1:num_z(i)) = z_zero + z_Inp(i,1:num_z(i))

              !Read first values
              read(44+i,*,end=9) tb_start(i),(Inp_read_start(i,j),j=1,num_z(i))
              call Integrate(z_Inp(i,:),Inp_read_start(i,:),Q_read_start(i,:),num_z(i))
              call Interp(z_Inp(i,:),Q_read_start(i,:),num_z(i)-1,z_cent,Q_start(i,:),nz-1)
              read(44+i,*,end=7) tb_end(i),(Inp_read_end(i,j),j=1,num_z(i))
              call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),num_z(i))
              call Interp(z_Inp(i,:),Q_read_end(i,:),num_z(i)-1,z_cent,Q_end(i,:),nz-1)
          end if

          if ((datum<=tb_start(i)).or.(eof(i)==1)) then  !if datum before first date or end of file reached
              goto 8
          else
              do while (.not.((datum>=tb_start(i)).and.(datum<=tb_end(i)))) ! do until datum between dates
                  tb_start(i) = tb_end(i)             ! move one step in time
                  Q_start(i,1:nz) = Q_end(i,1:nz)
                  read(44+i,*,end=7) tb_end(i),(Inp_read_end(i,j),j=1,num_z(i))
                  call Integrate(z_Inp(i,:),Inp_read_end(i,:),Q_read_end(i,:),num_z(i))
                  call Interp(z_Inp(i,:),Q_read_end(i,:),num_z(i)-1,z_cent,Q_end(i,:),nz-1)
              end do
              !Linearly interpolate value at correct datum (for all depths)
              do j=1,nz
                  Qstate(j,i) = Q_start(i,j) + (datum-tb_start(i)) * (Q_end(i,j)-Q_start(i,j))/(tb_end(i)-tb_start(i))
              end do
          end if
          goto 11

    7     eof(i) = 1
    8     Qstate(1:nz,i) = Q_start(i,1:nz)              ! Set to closest available value

          !Set all Q to the differences (from the integrals)
    11    do j=1,nz-1
              Qstate(nz-j+1,i) = (Qstate(nz-j+1,i)-Qstate(nz-j,i))
          end do
          return

    9     write(6,*) '   No data found in ',trim(fname),'. Check number of depths. Inflow set to zero.'
          goto 10

    1     write(6,*) '   File "',trim(fname),'" not found. Inflow set to zero.'
    10    eof(i) = 1
          Qstate(0:nz,i) = 0.0_dp
          Q_start(i,1:nz) = 0.0_dp
          return
    end subroutine Lateral_FABM


    !Returns the inflow vector Q_inp for a given FABM variable varname
    !by reading the corresponding file (if present)
    !Assuming inflow will plunge according to its density, entraining water from other layers
    !####################################################################
    subroutine Lateral_FABM_rho(datum,idx,Q_inp,varname,i)
    !####################################################################

          implicit none

          ! Global Declarations
          integer, intent(in) :: idx,i
          real(dp), intent(in) :: datum, Q_inp(0:nz)
          character*100, intent(in) :: varname

          

          ! Local Declarations
          character*100 :: fname
          integer :: j, eof(1:50)
          real(dp) :: Inp_s(1:50), Inp_e(1:50), Inp(1:50)
          real(dp) :: tb_start(1:50), tb_end(1:50)

          save eof
          save tb_start, tb_end

          if(trim(varname)=='aed_carbon_pH') goto 10

          if (idx==1) then
              if(disp_dgn==2) write(6,*) 'Starting to read '//trim(varname)//' inflow file...'
              fname = trim(PathFABM)//trim(varname)//'_inflow.dat'
              open(44+i,action='read',status='unknown',err=1,file=fname)
              eof(i) = 0
              !Skip first three rows
              read(44+i,*,end=9)
              read(44+i,*,end=9)
              read(44+i,*,end=9)
              !Read first values
              read(44+i,*,end=9) tb_start(i),Inp_s(i)
              read(44+i,*,end=7) tb_end(i),Inp_e(i)
          end if

          if ((datum<=tb_start(i)).or.(eof(i)==1)) then       ! if datum before first date or end of file reached
              goto 8
          else
              do while (.not.((datum>=tb_start(i)).and.(datum<=tb_end(i)))) ! do until datum between dates
                  tb_start(i) = tb_end(i)             ! move one step in time
                  Inp_s(i) = Inp_e(i)
                  read(44+i,*,end=7) tb_end(i),Inp_e(i)
              end do
              !Linearly interpolate value at correct datum
              Inp(i) = Inp_s(i) + (datum-tb_start(i)) * (Inp_e(i)-Inp_s(i))/(tb_end(i)-tb_start(i))
          end if
          goto 11

    7     eof(i) = 1
    8     Inp(i) = Inp_s(i)              ! Set to closest available value

          !Set compound input where there is water inflow
    11    do j=0,nz
              if(Q_inp(j)/=0) Qstate(j,i) = Inp(i)*Q_inp(j)
          end do
          return

    9     write(6,*) '   No data found in ',trim(fname),'. Check number of depths. Inflow set to zero.'
          goto 10

    1     write(6,*) '   File "',trim(fname),'" not found. Inflow set to zero.'
    10    eof(i) = 1
          Qstate(:,i) = 0.0_dp

          return
    end subroutine Lateral_FABM_rho

    subroutine doFABM(nuh)
        implicit none

        real(dp), intent(in) :: nuh(0:nz_grid)

        integer :: i,j

        depths(0:nz) = z_cent(nz)-z_cent(0:nz)
        ht(0:nz) = h(0:nz)
        !FABM calls
        !General checks and repairs
        call fabm_check_ready(model) !useless?
        valid = .false.
        call fabm_check_state(model,0,nz,.true.,valid)
        valid = .false.
        call fabm_check_surface_state(model,.true.,valid)
        valid = .false.
        call fabm_check_bottom_state(model,.true.,valid)
        !Relink and reinitialize variables
        do i=1,m
            call fabm_link_bulk_state_data(model,i,state(:,i))
        end do
        do i=1,ms
            call fabm_link_surface_state_data(model,i,surface_state(i))
        end do
        do i=1,mb
            call fabm_link_bottom_state_data(model,i,bottom_state(i))
        end do
        dy(0:nz,1:m) = 0.0_dp
        dyS(1:m) = 0.0_dp
        sflux(1:ms) = 0.0_dp
        bflux(0:nz,1:mb) = 0.0_dp

        !Call FABM to get rates of change (new version)
        do i=0,nz !Call do_bottom for each layer (sediment exchange)
            call model%set_bottom_index(i)
            call fabm_do_bottom(model,dy(i,:),bflux(i,:))
        end do
        do i=1,m !Weigh these values based on the surface of exchange
            dy(1:nz,i) = dy(1:nz,i)*dAdz(1:nz)/Az(1:nz) !=dy*AreaSediments/Volume for each layer
        end do
        !Call do_surface for the surface layer (atmosphere exchange)
        call fabm_do_surface(model,dyS,sflux)
        dy(nz,1:m) = dy(nz,1:m) + dyS(1:m)/h(nz)
        !Call do for each layer (temporal derivatives)
        call fabm_do(model,0,nz,dy)
            !!Call FABM to get production and destruction matrices
            !!call fabm_do_bottom(model,0,pp(0,:,:),dd(0,:,:),0)
            !!call fabm_do(model,0,nz,pp(0:nz,:,:),dd(0:nz,:,:))
        !Update surface and bottom state variables (forward Euler integration in time)
        surface_state(1:ms) = surface_state(1:ms) + dt*sflux(1:ms) !What was returned in sflux?? Units are var/s??
        do i=0,nz
            bottom_state(1:mb) = bottom_state(1:mb) + dt*bflux(i,1:mb)
        end do
        !Update pelagic state variables: time integration and advection (transport with the medium)
        do i=1,m
            if(trim(model%state_variables(i)%name)=='aed_carbon_pH') cycle
            call doTransportEquation(state(:,i),dy(:,i),nuh)
        end do

        !Active transport (without the medium)
        call fabm_get_vertical_movement(model,0,nz,vel(0:nz,1:m))
        do i=1,m !Advection scheme
            if(trim(model%state_variables(i)%name)=='aed_carbon_pH') cycle
            state(0,i) = state(0,i) - dt/h(1)*(merge(1.,0.,vel(0,i)>0)*vel(0,i)*state(0,i) &
                                             + merge(1.,0.,vel(1,i)<0)*vel(1,i)*state(1,i))
            state(nz,i) = state(nz,i) + dt/h(nz+1)*(merge(1.,0.,vel(nz,i)<0)*vel(nz,i)*state(nz,i) &
                                                + merge(1.,0.,vel(nz-1,i)<0)*vel(nz-1,i)*state(nz-1,i))
            do j=1,nz-1
                if(dt/h(j)*abs(vel(j,i))>1) write(6,*) 'Warning: CFL criterion violated'
                state(j,i) = state(j,i) + dt/h(j+1)*(-abs(vel(j,i))*state(j,i) &
                             + merge(1.,0.,vel(j-1,i)>0)*vel(j-1,i)*state(j-1,i) &
                             - merge(1.,0.,vel(j+1,i)<0)*vel(j+1,i)*state(j+1,i))
            end do
        end do

        call fabm_get_light_extinction(model,0,nz,ext(0:nz))
        !Conserved quantities for checks
        call fabm_get_conserved_quantities(model,0,nz,quant(0:nz,:))
        call fabm_get_horizontal_conserved_quantities(model,hquant)
    end subroutine doFABM

  !####################################################################
  subroutine doAdvection_FABM(Q_vert,Q_inp,U,V,T,S,k,eps,state,Qstate,m)
  !####################################################################
    
    implicit none

    !Global Declarations
    real(dp), intent(inout) :: U(0:nz_grid), V(0:nz_grid), T(0:nz_grid), S(0:nz_grid)
    real(dp), intent(inout) :: k(0:nz_grid), eps(0:nz_grid)
    real(dp), intent(inout) :: Q_vert(0:nz_grid) ! Vertical advection [m3/s]
    real(dp), intent(in) :: Q_inp(1:4,0:nz_grid) ! Inflow [m3/s], Outflow [m3/s], T-input [°C*m3/s], S-input [‰*m3/s]
    real(dp), intent(inout) :: state(0:nz_grid,1:m), Qstate(0:nz_grid,1:m)
    integer, intent(in) :: m   ! Number of FABM variables

    real(dp) :: dh, dhi(1:2), h_div_2      ! depth differences
    real(dp) :: dti(1:2)          ! first and second time step
    real(dp) :: dU(0:nz_grid), dV(0:nz_grid), dTemp(0:nz_grid), dS(0:nz_grid)
    real(dp) :: AreaFactor_adv(1:nz_grid)
    real(dp) :: top

    integer :: i, ti

    !Depth difference
    dh = Q_vert(nz)/Az(nz)*dt
    h_div_2 = 0.5_dp*h(nz-1)
    !Split timestep depending on situation
    if (dh==0.) then                          ! If volume does not change
        dti(1) = dt                                 ! One normal time step
    else if ((dh+z_upp(nz))>=depth) then        ! If surface level reached
        dti(1) = (depth - z_upp(nz))/dh*dt             ! Step until surface level
    else if (((dh+h(nz))>h_div_2) .and.&   ! If top box>0.5*lower box
            ((dh+h(nz))<2*h(nz-1))) then     ! and top box<2*lower box
        dti(1) = dt                                 ! One normal time step
    else if ((dh+h(nz))<=h_div_2) then     ! If top box<=0.5*lower box
        dti(1) = abs((h(nz)-h_div_2)/dh)*dt       ! First step until top box = 0.5*lower box
    else                                     ! If top box>=2*lower box
        dti(1) = abs((2*h(nz-1)-h(nz))/dh)*dt       ! First step until top box = 2*lower box
    end if
    dti(2) = dt-dti(1)                       ! Rest of timestep

    do ti=1,2 !First and (if needed) second timestep
        AreaFactor_adv(1:nz) = dti(ti)/(Az(1:nz)*h(1:nz))     ! doAreaFactors factor for dt(ti)
        dhi(ti) = dh*dti(ti)/dt                         ! Depth difference for dt(ti)

        do i=1,nz
            if (i==nz .and. Q_vert(i)>0) then
                top = 0
            else
                top = 1
            end if
            ! Advective flow out of box i, always negative
            dU(i)    = -top*abs(AreaFactor_adv(i)*Q_vert(i))*U(i)
            dV(i)    = -top*abs(AreaFactor_adv(i)*Q_vert(i))*V(i)
            dTemp(i) = -top*abs(AreaFactor_adv(i)*Q_vert(i))*T(i)
            dS(i)    = -top*abs(AreaFactor_adv(i)*Q_vert(i))*S(i)
            if (i>1 .and. Q_vert(i-1)>0) then       ! Advective flow into box i, from above
                dU(i)    = dU(i) + AreaFactor_adv(i)*Q_vert(i-1)*U(i-1)
                dV(i)    = dV(i) + AreaFactor_adv(i)*Q_vert(i-1)*V(i-1)
                dTemp(i) = dTemp(i) + AreaFactor_adv(i)*Q_vert(i-1)*T(i-1)
                dS(i)    = dS(i) + AreaFactor_adv(i)*Q_vert(i-1)*S(i-1)
            end if
            if (i<nz .and. Q_vert(i+1)<0) then      ! Advective flow into box i, from below
                dU(i)    = dU(i) - AreaFactor_adv(i)*Q_vert(i+1)*U(i+1)
                dV(i)    = dV(i) - AreaFactor_adv(i)*Q_vert(i+1)*V(i+1)
                dTemp(i) = dTemp(i) - AreaFactor_adv(i)*Q_vert(i+1)*T(i+1)
                dS(i)    = dS(i) - AreaFactor_adv(i)*Q_vert(i+1)*S(i+1)
            end if
        end do

        ! Inflow and outflow
        dTemp(1:nz) = dTemp(1:nz) + AreaFactor_adv(1:nz)*(Q_inp(3,1:nz)+Q_inp(2,1:nz)*T(1:nz))
        dS(1:nz) = dS(1:nz) + AreaFactor_adv(1:nz)*(Q_inp(4,1:nz)+Q_inp(2,1:nz)*S(1:nz))

        ! Take first time step
        U(1:nz) = U(1:nz) + dU(1:nz)
        V(1:nz) = V(1:nz) + dV(1:nz)
        T(1:nz) = T(1:nz) + dTemp(1:nz)
        S(1:nz) = S(1:nz) + dS(1:nz)

        do i=1,m
            if(i==3) cycle !this shouldn't be executed for pH... otherwise it gets diluted...
            state(1:nz,i) = state(1:nz,i) + AreaFactor_adv(1:nz)*(Qstate(1:nz,i)+Q_inp(2,1:nz)*state(1:nz,i))
        end do

        ! Variation of variables due to change in volume
        U(nz) = U(nz)*h(nz)/(h(nz)+dhi(ti))
        V(nz) = V(nz)*h(nz)/(h(nz)+dhi(ti))
        T(nz) = T(nz)*h(nz)/(h(nz)+dhi(ti))
        S(nz) = S(nz)*h(nz)/(h(nz)+dhi(ti))

        if (ti==1) then
            if (dh==0) then                          ! If volume does not change
                return
            else if ((dh+z_upp(nz))>=depth) then        ! If surface level reached
                h(nz) = h(nz) + dhi(1)                      ! New thickness of top box
                z_cent(nz) = z_cent(nz) + 0.5_dp*dhi(1)                  ! New centre coordinate of top box
                z_upp(nz) = z_upp(nz) + dhi(1)                    ! New surface coordinate of top box
                return                                      ! No change in volume (overflow)
            else if (((dh+h(nz))>h_div_2) .and.&   ! If top box>0.5*lower box
                     ((dh+h(nz))<2*h(nz-1))) then    ! and top box<2*lower box
                h(nz) = h(nz) + dhi(1)                      ! New thickness of top box
                z_cent(nz) = z_cent(nz) + 0.5_dp*dhi(1)                  ! New centre coordinate of top box
                z_upp(nz) = z_upp(nz) + dhi(1)                    ! New surface coordinate of top box
                return
            else if ((dh+h(nz))<=h_div_2) then     ! If top box<=0.5*lower box
                z_upp(nz-1) = z_upp(nz)
                z_cent(nz-1) = 0.5_dp*(z_upp(nz-1)+z_upp(nz-2))
                h(nz-1)  = (z_upp(nz-1)-z_upp(nz-2))
                U(nz-1) = (0.5_dp*U(nz)*Az(nz)+U(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                V(nz-1) = (0.5_dp*V(nz)*Az(nz)+V(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                T(nz-1) = (0.5_dp*T(nz)*Az(nz)+T(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                S(nz-1) = (0.5_dp*S(nz)*Az(nz)+S(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                k(nz-1) = (0.5_dp*k(nz)*Az(nz)+k(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                eps(nz-1) = (0.5_dp*eps(nz)*Az(nz)+eps(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                Q_vert(nz-1) = (0.5_dp*Q_vert(nz)*Az(nz)+Q_vert(nz-1)*Az(nz-1))/(0.5_dp*Az(nz)+Az(nz-1))
                nz = nz-1                     ! Reduce number of boxes
                !dh(2) = Q_vert(nz)/Az(nz)*dt(2) !(also equal to dh*dt2/dt)
                call doAreaFactors()
            else                                     ! If top box>=2*lower box
                h(nz+1)   = 0.5_dp*h(nz) + dh   ! AG 2014 (added +dh)
                h(nz)     = 0.5_dp*h(nz)
                z_upp(nz+1)  = z_upp(nz) + dh    ! AG 2014 (added +dh)
                z_upp(nz)    = z_upp(nz) - 0.5_dp*h(nz)
                z_cent(nz+1)  = z_cent(nz) + 0.5_dp*h(nz+1)
                z_cent(nz)    = z_cent(nz) - 0.5_dp*h(nz+1)
                U(nz+1)   = U(nz)
                V(nz+1)   = V(nz)
                T(nz+1)   = T(nz)
                S(nz+1)   = S(nz)
                k(nz+1)   = k(nz)
                eps(nz+1) = eps(nz)
                Q_vert(nz+1) = Q_vert(nz)       ! Vertical discharge of new box
                nz = nz + 1                   ! Increase number of boxes
                !dh(2) = Q_vert(nz)/Az(nz)*dt(2) !(also equal to dh*dt2/dt)
                call doAreaFactors()
            end if
        end if      !end if (ti==1)
    end do      !end do ti=1,2

    return
  end subroutine doAdvection_FABM
end module simstrat_fabm