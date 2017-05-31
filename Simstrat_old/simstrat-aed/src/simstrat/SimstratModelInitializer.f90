module SimstratModelInitializer
    !MZ Feb. 2016: modularized version
    !FB June 2016: reorganized
    use SimstratModel
    use SimstratOutput, only : save_ini
    use Advection, only : doAreaFactors
    use utilities
    implicit none

    private
    public InitializeSimstratModel

contains
    !####################################################################
    subroutine InitializeSimstratModel()
    !####################################################################

        implicit none

        ! Local variables
        integer :: i

        ! Read model parameters from par-file
        call ReadModelParameters()

        ! Read morphology
        call InitializeMorphology()

        ! Create grid for constant or variable spacing (determine nz_grid)
        call InitializeGridPoints()

        ! Allocate memory for initialization using nz_grid
        allocate(z_cent(0:nz_grid))         ! Depth axis with center of boxes
        allocate(z_upp(0:nz_grid))          ! Depth axis with upper border of boxes
        allocate(Az(0:nz_grid))             ! Az is defined on z_upp
        allocate(dAdz(1:nz_grid))           ! dAz is the difference between Az and thus on z_cent

        allocate(U_ini(0:nz_grid))          ! Water velocity
        allocate(V_ini(0:nz_grid))          ! Water velocity
        allocate(T_ini(0:nz_grid))          ! Temperature
        allocate(S_ini(0:nz_grid))          ! Salinity
        allocate(k_ini(0:nz_grid))          ! Turbulent kinetic energy TKE [J/kg]
        allocate(eps_ini(0:nz_grid))        ! Dissipation of TKE [W/kg]

        allocate(num_ini(0:nz_grid))        ! Turbulent viscosity
        allocate(nuh_ini(0:nz_grid))        ! Turbulent temperature diffusivity
        allocate(fgeo_add(1:nz_grid))       ! Geothermal heat flux

        ! Area factors used in calculations
        allocate(AreaFactor_1(1:nz_grid))
        allocate(AreaFactor_2(1:nz_grid))
        allocate(AreaFactor_k1(1:nz_grid))
        allocate(AreaFactor_k2(1:nz_grid))
        allocate(AreaFactor_eps(1:nz_grid))

        allocate(meanint(0:nz_grid))        ! Inverse ratio of mean height of two adjacent boxes


        ! Initialize depth axes z_cent and z_upp
        call InitializeZAxes()

        ! Initialize lake area Az and lake area derivative dAdz
        call InitializeAreas()

        ! Read initial conditions from file
        call ReadInitialConditions()

        ! Initialize area factors
        call doAreaFactors()

        ! Initialize some more values
        if(stab==1) cm0=0.5625_dp
        if(stab==2) cm0=0.556171_dp
        cde=cm0**3
        sig_e=(kappa/cm0)**2/(ce2-ce1)

        num_ini(0:nz_grid) = 0.0_dp
        nuh_ini(0:nz_grid)= 0.0_dp

        tx_ini=0.0_dp
        ty_ini=0.0_dp

        drag_ini=(kappa/log(1.0_dp+30/K_s*h(1)/2))**2

        ! Open inflow/ouflow files and set advection status and number of input depths
        call check_advection()

        ! Determine if output is binary or text
        OutBin = .false.
        if(scan(PathOut,'.')/=0) OutBin = .true.
        if(OutBin) open(80,access='SEQUENTIAL',status='unknown',FORM='unformatted',file=PathOut)

        ! Set output properties
        call save_ini()

        ! Geothermal heat flux
        if (fgeo/=0) then
            fgeo_add(1:nz_grid) = fgeo/rho_0/cp*dAdz(1:nz_grid)/Az(1:nz_grid) ! calculation per kg
            if (Az(0)/=0) then
                fgeo_add(1) = fgeo_add(1)+2*fgeo/rho_0/cp*Az(0)/((Az(0)+Az(1))*h(1))
            end if
        end if

        ! Salinity control for buoyancy functions
        if (ModSal/=0) then
            salctr=1
            delsal=1
        else
            do i=1,nz_grid
                if(S_ini(i)/=0) salctr=1
            end do
            if (salctr==1) then
                do i=2,nz_grid
                    if(S_ini(i)-S_ini(i-1)/=0) delsal=1
                end do
            end if
        end if
        return
    end subroutine InitializeSimstratModel

    !MZ Feb. 2015: modularized version
    !AG 2014-2015: revision and many changes
    !####################################################################
    subroutine ReadModelParameters()
    !####################################################################
        implicit none


        !Local variables only used to read input file
        real(dp) Lat

        ! Read first argument for parameter file, if empty use default
        call getarg(1,ParName)
        if(ParName=='') ParName='simstrat.par'

        ! Reading user-defined parameters from file
        open(10,status='old',file=ParName)
        read(10,*)
        read(10,'(A)') InitName
        read(10,'(A)') GridName
        read(10,'(A)') MorphName
        read(10,'(A)') ForcingName
        read(10,'(A)') AbsorpName
        read(10,'(A)') PathOut
        read(10,'(A)') zoutName
        read(10,'(A)') toutName
        read(10,'(A)') QinpName
        read(10,'(A)') QoutName
        read(10,'(A)') TinpName
        read(10,'(A)') SinpName
#if use_fabm
        read(10,'(A)') PathFABM
        read(10,'(A)') PathFABMOut
#endif
        read(10,*)
        read(10,*) dt
        read(10,*) t_start
        read(10,*) t_end
        read(10,*)
        read(10,*) Mod
        read(10,*) Stab
        read(10,*) ModFlux
        read(10,*) NBC
        read(10,*) WindFilt
        read(10,*) ModSNorm
        read(10,*) ModC10
        read(10,*) ModInflow
        read(10,*) Pgrad
        read(10,*) ModSal
        read(10,*) disp_sim
        read(10,*) disp_dgn
        read(10,*) igoal
        read(10,*)
        read(10,*) Lat
        read(10,*) p_air
        read(10,*) a_seiche
        read(10,*) q_NN
        read(10,*) f_wind
        read(10,*) C10
        read(10,*) CD
        read(10,*) fgeo
        read(10,*) k_min
        read(10,*) p_radin
        read(10,*) p_windf
        read(10,*) beta_sol
        read(10,*) albsw
        close(10)

        ! Calculate Coriolis parameter from latitude
        Cori=2.0_dp*7.292e-5_dp*sin(Lat*pi/180)

        if(.not.(disp_sim==1 .or. disp_sim==2 .or. disp_sim==3)) disp_sim=0
        if(.not.(disp_dgn==1 .or. disp_dgn==2)) disp_dgn=0

        ! Change igoal to a different type of information
        if (abs(igoal)==10) then !Set flag for all variables
            igoal_s(1:13) = 1
        elseif (abs(igoal)==11) then !Set flag for U, V, T
            igoal_s(1:3) = 1
        elseif (abs(igoal)==12) then !Set flag for T, nuh
            igoal_s(3) = 1
            igoal_s(8) = 1
        else !Set flag for variable at index |igoal|
            igoal_s(abs(igoal)) = 1
        end if

        if(disp_dgn/=0) then
            write(6,*) 'Configuration  : ',trim(ParName)
            write(6,*) 'Successfully read'
            write(6,*)
        end if

        return
    end subroutine ReadModelParameters


    !####################################################################
    subroutine InitializeMorphology()
    !####################################################################
        implicit none


        ! Local variables
        real(dp) :: z_tmp(0:nz_max), A_tmp(0:nz_max)
        integer :: num_read, i

        if(disp_dgn/=0) write(6,*) 'Morphology     : ',trim(MorphName)
        open(11,status='old',file=MorphName)
        read(11,*)                            ! Skip header
        do i=0,nz_max                            ! Read depth and area
            read(11,*,end=9) z_tmp(i), A_tmp(i)
        end do
    9   if(i==nz_max) write(6,*) 'Only first ',nz_max,' values of file read.'
        close(11)

        num_read = i  ! Number of area values

        allocate(z_read(0:num_read-1), A_read(0:num_read-1))

        do i=0,num_read-1                          ! Reverse order of values
            z_read(i) = -z_tmp(num_read-i-1)
            A_read(i) = A_tmp(num_read-i-1)
        end do

        depth = z_read(0) - z_read(num_read-1)               ! depth = max - min depth
        z_zero = z_read(0)                         ! zero point of z
        lake_level_old = z_zero
        z_read(0:num_read-1) = z_zero - z_read(0:num_read-1)      ! z-coordinate is positive upwards, zero point is at reservoir bottom

        if (disp_dgn/=0) then
            write(6,*) 'No depth points: ',num_read
            if (disp_dgn==2) then
                write(6,('(A8, A20)')) 'Depth [m]','Area [m2]'
                do i=0,num_read-1
                    write(6,'(F8.2, G20.2)') z_read(i),A_read(i)
                end do
            end if
            write(6,*) 'Successfully read'
            write(6,*)
        end if

        return
    end subroutine InitializeMorphology


    !####################################################################
    subroutine InitializeGridPoints()
    !####################################################################
        implicit none


        ! Local variables
        integer ictr,i
        real(dp) grid_read(0:nz_max)

        if(disp_dgn/=0) write(6,*) 'Grid           : ',trim(GridName)
        open(12,status='old',file=GridName)
        read(12,*)
        do ictr=0,nz_max
            read(12,*,end=9) grid_read(ictr)
        end do
    9   if(ictr==nz_max) write(6,*) 'Only first ',nz_max,' values of file read.'
        close(12)

        if (ictr==1) then     ! Constant spacing
            nz_grid=int(grid_read(0))
        else                  ! Variable spacing
            nz_grid=ictr-1
            !Include top value if not included
            if (grid_read(0)/=0.) then
                nz_grid=nz_grid+1
                do i=nz_grid,1,-1
                    grid_read(i)=grid_read(i-1)
                end do
                grid_read(0)=0.0_dp
            end if
            !If maxdepth grid larger than morphology
            if (grid_read(nz_grid)>depth) then
                do while ((grid_read(nz_grid)>depth).and.(nz_grid>0.))
                    nz_grid=nz_grid-1
                end do
            end if
            !Include bottom value if not included
            if (grid_read(nz_grid)<depth) then
                nz_grid=nz_grid+1
                grid_read(nz_grid)=depth
            end if
        end if

        allocate(h(0:nz_grid))
        h(0) = 0                ! Note that h(0) has no physical meaning but helps with some calculations
        if (ictr==1) then
            h(1:nz_grid) = depth/nz_grid
        else
            do i=1,nz_grid
                h(1+nz_grid-i)=grid_read(i)-grid_read(i-1)
            end do
        end if

    end subroutine InitializeGridPoints

    !####################################################################
    subroutine InitializeZAxes()
    !####################################################################
        implicit none

        ! Local variables
        integer i

        !Compute position of layer center and top
        z_cent(0)=0.0_dp
        z_upp(0)=0.0_dp
        do i=1,nz_grid
            z_cent(i)=z_cent(i-1)+0.5_dp*(h(i-1)+h(i))
            z_upp(i)=z_upp(i-1)+h(i)
        end do
        do i=1,nz_grid
            z_cent(i)=nint(1e6_dp*z_cent(i))/1e6_dp
            z_upp(i)=nint(1e6_dp*z_upp(i))/1e6_dp
        end do

        lake_level_old = z_upp(nz)

        if (disp_dgn/=0) then
            write(6,*) 'No grid points : ',nz_grid
            if (disp_dgn==2) then
                write(6,*) 'Grid from bottom (0m) to top:'
                write(6,'(F8.2,$)') (z_upp(i),i=0,nz_grid)
            end if
            write(6,*) 'Successfully read'
            write(6,*)
        end if

        return
    end subroutine InitializeZAxes


    !####################################################################
    subroutine InitializeAreas()
    !####################################################################
        implicit none

        integer num_read

        num_read = size(z_read)

        ! Interpolate area (A) at all depths (z_upp)
        call Interp(z_read, A_read, num_read-1, z_upp, Az, nz_grid)
        ! Compute area derivative (= projected sediment area over layer thickness)
        dAdz(1:nz_grid) = (Az(1:nz_grid)-Az(0:nz_grid-1))/(z_upp(1:nz_grid)-z_upp(0:nz_grid-1))

        deallocate(z_read, A_read)

    end subroutine InitializeAreas

    !####################################################################
    subroutine ReadInitialConditions()
    !####################################################################
        implicit none

        ! Local variables
        real(dp) :: z_tmp(0:nz_max), U_tmp(0:nz_max), V_tmp(0:nz_max)
        real(dp) :: T_tmp(0:nz_max), S_tmp(0:nz_max), k_tmp(0:nz_max), eps_tmp(0:nz_max)

        real(dp) :: z_read(0:nz_max), U_read(0:nz_max), V_read(0:nz_max)
        real(dp) :: T_read(0:nz_max), S_read(0:nz_max), k_read(0:nz_max), eps_read(0:nz_max)

        real(dp) :: z_ini(0:nz_max)
        real(dp) :: z_ini_depth, zmax

        integer :: i,num_read


        if(disp_dgn/=0) write(6,*) 'Initial conditions: ',trim(InitName)
        open(13,status='old',file=InitName)     ! Opens initial conditions file
        read(13,*)                              ! Skip header
        do i=0,nz_max                              ! Read initial u,v,T, etc
            read(13,*,end=9) z_read(i),U_read(i),V_read(i),T_read(i),S_read(i),k_read(i),eps_read(i)
        end do
    9   num_read = i-1                                ! Number of values
        if (num_read<0) then
            write(6,*) 'Error reading initial conditions files (no data found).'
            stop
        end if
        close(13)
        do i=0,num_read
            z_read(i) = abs(z_read(i))               ! Make depths positive
        end do
        z_ini_depth = z_read(0)                     ! Initial depth (top-most)

        do i=0,nz_grid
            if (z_upp(i) >= (z_zero-z_ini_depth)) then    ! If above initial water level
                zmax = z_upp(i)
                z_upp(i)= z_zero-z_ini_depth
                z_cent(i)= (z_upp(i)+z_upp(i-1))/2
                h(i) = z_upp(i) - z_upp(i-1)
                !Az(i) = Az(i-1) + h(i)/(zmax-z_upp(i-1))*(Az(i)-Az(i-1))
                nz = i
                if (h(nz)<=0.5*h(nz-1)) then         ! If top box is too small
                    z_upp(nz-1) = z_upp(nz)                ! Combine the two upper boxes
                    z_cent(nz-1) = (z_upp(nz-1)+z_upp(nz-2))/2
                    h(nz-1)  = h(nz)+h(nz-1)
                    nz = nz-1                        ! Reduce number of boxes
                end if
                exit
            end if
        end do

        do i=0,num_read
            z_tmp(num_read-i) = z_zero - z_read(i)
            U_tmp(num_read-i) = U_read(i)
            V_tmp(num_read-i) = V_read(i)
            T_tmp(num_read-i) = T_read(i)
            S_tmp(num_read-i) = S_read(i)
            k_tmp(num_read-i) = k_read(i)
            eps_tmp(num_read-i) = eps_read(i)
        end do

        if (num_read==0) then
            write(6,*) 'Only one row! Water column will be initially homogeneous.'
            U_ini(0:nz) = U_tmp(0)
            V_ini(0:nz) = V_tmp(0)
            T_ini(0:nz) = T_tmp(0)
            S_ini(0:nz) = S_tmp(0)
            k_ini(0:nz) = k_tmp(0)
            eps_ini(0:nz) = eps_tmp(0)
        else
            call Interp(z_tmp, U_tmp, num_read, z_cent, U_ini, nz_grid)
            call Interp(z_tmp, V_tmp, num_read, z_cent, V_ini, nz_grid)
            call Interp(z_tmp, T_tmp, num_read, z_cent, T_ini, nz_grid)
            call Interp(z_tmp, S_tmp, num_read, z_cent, S_ini, nz_grid)
            call Interp(z_tmp, k_tmp, num_read, z_upp, k_ini, nz_grid)
            call Interp(z_tmp, eps_tmp, num_read, z_upp, eps_ini, nz_grid)
        end if


        if (disp_dgn/=0) then
            if (disp_dgn==2) then
                do i=num_read,0,-1
                    write(6,999) depth-z_tmp(i),U_tmp(i),V_tmp(i),T_tmp(i),S_tmp(i),k_tmp(i),eps_tmp(i)
                999 format(F8.2,F10.3,F10.3,F10.3,G10.3,G10.3,G10.3)
                end do
            end if
            write(6,*) 'Successfully read'
            write(6,*)
        end if
        return
    end subroutine ReadInitialConditions

    !Set advection to 1 if any inflow/outflow file contains data, otherwise to 0
    !AG 2014/FB 2016: revision
    !####################################################################
    subroutine check_advection()
    !####################################################################
        implicit none

        ! Local variables
        integer  :: i, j, fnum(1:4), nval(1:4)
        real(dp) :: dummy, z_Inp_dummy(0:nz_max)


        if(disp_dgn==2) write(6,*) 'Opening physical inflow/outflow files...'
        open(41,status='old',file=QinpName)
        open(42,status='old',file=QoutName)
        open(43,status='old',file=TinpName)
        open(44,status='old',file=SinpName)

        fnum = [41,42,43,44]
        if_adv = 0
        do i=1,4
            read(fnum(i),*,end=8)                           ! Skip header (description of columns)
            read(fnum(i),*,end=8) nval(i)                  ! Read number of input depths (static)
            read(fnum(i),*,end=8) dummy, (z_Inp_dummy(j),j=0,nval(i)-1) ! Read input depths
            goto 9
        8   if_adv = if_adv + 1
        9   rewind(fnum(i))
        end do

        if (if_adv==4) then
            if_adv = 0
        else
            if_adv = 1
        end if

        nz_input = maxval(nval)

        return
    end subroutine check_advection

end module SimstratModelInitializer
