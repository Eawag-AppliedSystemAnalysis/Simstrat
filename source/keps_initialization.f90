!####################################################################
subroutine Initialization()
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision h(0:mxl),Az(0:mxl),dAdz(0:mxl)
    double precision Sini(0:mxl),Tini(0:mxl),Uini(0:mxl),Vini(0:mxl)
    double precision kini(0:mxl),epsini(0:mxl),Lini(0:mxl),numini(0:mxl),nuhini(0:mxl)
    double precision dragini,txini,tyini

    common /morph_variables1/ h,Az,dAdz
    common /ini_values1/ Sini,Tini
    common /ini_values2/ Uini,Vini
    common /ini_values3/ epsini,kini,Lini
    common /ini_values4/ numini,nuhini
    common /ini_values5/ dragini,txini,tyini

    ! Local variables
    integer i

    ! Calculations
    call ParameterList()
    call Morph()
    call InitCond()

    if(stab==1) cm0=0.5625
    if(stab==2) cm0=0.556171
    cde=cm0**3
    sig_e=(kappa/cm0)**2/(ce2-ce1)

    Lini(0:xli) = 0.2 !L_min

    txini=0.
    tyini=0.

    numini(0:xli) = 0.
    nuhini(0:xli) = 0.

    dragini=(kappa/log(1.+30/K_s*h(1)/2))**2

    ! Open inflow/ouflow files and set advection status
    call check_advection()

    ! Determine if output is binary or text
    OutBin = .false.
    if(scan(PathOut,'.')/=0) OutBin = .true.
    if(OutBin) open(80,access='SEQUENTIAL',status='unknown',FORM='unformatted',file=PathOut)

    ! Set output properties
    call save_ini()

    ! Geothermal heat flux
    if (fgeo/=0) then
        fgeo_add(1:xli) = fgeo/rho_0/cp*dAdz(1:xli)/Az(1:xli) ! calculation per kg
        if (Az(0)/=0) then
            fgeo_add(1) = fgeo_add(1)+2*fgeo/rho_0/cp*Az(0)/((Az(0)+Az(1))*h(1))
        end if
    end if

    ! Sedimentation
    fsed = 2.5e-9 ! Test
    if (fsed/=0) then
        fsed_add(1:xli) = fsed*dAdz(1:xli)/Az(1:xli) ! calculation per kg
        if (Az(0)/=0) then
            fsed_add(1) = fsed_add(1)+2*fsed*Az(0)/((Az(0)+Az(1))*h(1))
        end if
    end if

    ! Salinity control for buoyancy functions
    if (ModSal/=0) then
        salctr=1
        delsal=1
    else
        do i=0,xli
            if(Sini(i)/=0) salctr=1
        end do
        if (salctr==1) then
            do i=1,xli
                if(Sini(i)-Sini(i-1)/=0) delsal=1
            end do
        end if
    end if

    return
end


!AG 2014-2015: revision and many changes
!####################################################################
subroutine ParameterList()
!####################################################################

    implicit none
    include 'common_parameters.i'

    double precision Lat

    ! Set fixed parameters for the model
    include 'incl_set_fixparameter.i'

    ! Read first argument for parameter file, if empty use default
    call getarg(1,ParName)
    if(ParName=='') ParName='kepsilon.par'

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
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*) dt
    read(10,*) t_start
    read(10,*) t_end
    read(10,*)
    read(10,*) Mod
    read(10,*)
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
    Cori=2*7.292e-5*sin(Lat*pi/180)

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
end


!####################################################################
subroutine Morph()
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision zu(0:mxl),zk(0:mxl),h(0:mxl),z_zero
    double precision Az(0:mxl),dAdz(0:mxl)

    common /morph_variables1/ h,Az,dAdz
    common /morph_variables2/ zu,zk,z_zero

    ! Local variables
    double precision z(0:mxl), A(0:mxl), zr(0:mxl), Ar(0:mxl)
    integer num, i

    if(disp_dgn/=0) write(6,*) 'Morphology     : ',trim(MorphName)
    open(11,status='old',file=MorphName)
    read(11,*)                            ! Skip header
    do i=0,mxl                            ! Read depth and area
        read(11,*,end=9) zr(i), Ar(i)
    end do
9   if(i==mxl) write(6,*) 'Only first ',mxl,' values of file read.'
    close(11)

    num = i                               ! Number of area values
    do i=0,num-1                          ! Reverse order of values
        z(i) = -zr(num-i-1)
        A(i) = Ar(num-i-1)
    end do

    depth = z(0) - z(num-1)               ! depth = max - min depth
    z_zero = z(0)                         ! zero point of z
    z(0:num-1) = z_zero - z(0:num-1)      ! z-coordinate is positive upwards, zero point is at reservoir bottom

    if (disp_dgn/=0) then
        write(6,*) 'No depth points: ',num
        if (disp_dgn==2) then
            write(6,('(A8, A20)')) 'Depth [m]','Area [m2]'
            do i=0,num-1
                write(6,'(F8.2, G20.2)') z(i),A(i)
            end do
        end if
        write(6,*) 'Successfully read'
        write(6,*)
    end if

    ! Create grid for constant spacing
    call Grid()
    ! Interpolate area (A) at all depths (zk)
    call Interp(z, A, num-1, zk, Az, xli)
    ! Compute area derivative (= projected sediment area over layer thickness)
    dAdz(1:xli) = (Az(1:xli)-Az(0:xli-1))/(zk(1:xli)-zk(0:xli-1))

    return
end


!####################################################################
subroutine Grid()
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision h(0:mxl),Az(0:mxl),dAdz(0:mxl)
    double precision zu(0:mxl),zk(0:mxl),z_zero

    common /morph_variables1/ h,Az,dAdz
    common /morph_variables2/ zu,zk,z_zero

    ! Local variables
    integer ictr,i
    double precision gr(0:mxl)

    if(disp_dgn/=0) write(6,*) 'Grid           : ',trim(GridName)
    open(12,status='old',file=GridName)
    read(12,*)
    do ictr=0,mxl
        read(12,*,end=9) gr(ictr)
    end do
9   if(ictr==mxl) write(6,*) 'Only first ',mxl,' values of file read.'
    close(12)

    if (ictr==1) then     ! Constant spacing
        xli=int(gr(0))
        !Compute layer thickness
        h(1:xli) = depth/xli
    else                  ! Variable spacing
        xli=ictr-1
        !Include top value if not included
        if (gr(0)/=0.) then
            xli=xli+1
            do i=xli,1,-1
                gr(i)=gr(i-1)
            end do
            gr(0)=0.
        end if
        !If maxdepth grid larger than morphology
        if (gr(xli)>depth) then
            do while ((gr(xli)>depth).and.(xli>0.))
                xli=xli-1
            end do
        end if
        !Include bottom value if not included
        if (gr(xli)<depth) then
            xli=xli+1
            gr(xli)=depth
        end if
        !Compute layer thickness
        do i=1,xli
            h(1+xli-i)=gr(i)-gr(i-1)
        end do
    end if

    !Compute position of layer center and top
    zu(0)=0.
    zk(0)=0.
    do i=1,xli
        zu(i)=zu(i-1)+0.5*(h(i-1)+h(i))
        zk(i)=zk(i-1)+h(i)
    end do
    do i=1,xli
        zu(i)=nint(1e6*zu(i))/1e6
        zk(i)=nint(1e6*zk(i))/1e6
    end do

    if (disp_dgn/=0) then
        write(6,*) 'No grid points : ',xli
        if (disp_dgn==2) then
            write(6,*) 'Grid from bottom (0m) to top:'
            write(6,'(F8.2,$)') (zk(i),i=0,xli)
        end if
        write(6,*) 'Successfully read'
        write(6,*)
    end if

    return
end


!####################################################################
subroutine InitCond()
!####################################################################

    implicit none
    include 'common_parameters.i'

    double precision Sini(0:mxl),Tini(0:mxl),Uini(0:mxl),Vini(0:mxl)
    double precision kini(0:mxl),epsini(0:mxl),Lini(0:mxl)
    double precision numini(0:mxl),nuhini(0:mxl)

    common /ini_values1/ Sini,Tini
    common /ini_values2/ Uini,Vini
    common /ini_values3/ epsini,kini,Lini
    common /ini_values4/ numini,nuhini

    double precision zu(0:mxl),zk(0:mxl),h(0:mxl),z_zero
    double precision Az(0:mxl),dAdz(0:mxl)
    common /morph_variables1/ h,Az,dAdz
    common /morph_variables2/ zu,zk,z_zero

    double precision z_tmp(0:mxl),u_tmp(0:mxl),v_tmp(0:mxl),&
                     T_tmp(0:mxl),S_tmp(0:mxl),k_tmp(0:mxl),eps_tmp(0:mxl)
    double precision zini(0:mxl)
    double precision zini_depth, zmax
    integer i,num

    if(disp_dgn/=0) write(6,*) 'Initial conditions: ',trim(InitName)
    open(13,status='old',file=InitName)     ! Opens initial conditions file
    read(13,*)                              ! Skip header
    do i=0,mxl                              ! Read initial u,v,T, etc
        read(13,*,end=9) zini(i),Uini(i),Vini(i),Tini(i),Sini(i),kini(i),epsini(i)
    end do
9   num = i-1                                ! Number of values
    if (num<0) then
        write(6,*) 'Error reading initial conditions files (no data found).'
        stop
    end if
    close(13)
    do i=0,num
        zini(i) = abs(zini(i))               ! Make depths positive
    end do
    zini_depth = zini(0)                     ! Initial depth (top-most)

    do i=0,xli
        if (zk(i) >= (z_zero-zini_depth)) then    ! If above initial water level
            zmax = zk(i)
            zk(i)= z_zero-zini_depth
            zu(i)= (zk(i)+zk(i-1))/2
            h(i) = zk(i) - zk(i-1)
            !Az(i) = Az(i-1) + h(i)/(zmax-zk(i-1))*(Az(i)-Az(i-1))
            xl = i
            if (h(xl)<=0.5*h(xl-1)) then         ! If top box is too small
                zk(xl-1) = zk(xl)                ! Combine the two upper boxes
                zu(xl-1) = (zk(xl-1)+zk(xl-2))/2
                h(xl-1)  = h(xl)+h(xl-1)
                xl = xl-1                        ! Reduce number of boxes
            end if
            exit
        end if
    end do

    do i=0,num
        z_tmp(num-i) = z_zero - zini(i)
        U_tmp(num-i) = Uini(i)
        V_tmp(num-i) = Vini(i)
        T_tmp(num-i) = Tini(i)
        S_tmp(num-i) = Sini(i)
        k_tmp(num-i) = kini(i)
        eps_tmp(num-i) = epsini(i)
    end do

    if (num==0) then
        write(6,*) 'Only one row! Water column will be initially homogeneous.'
        Uini(0:xl) = U_tmp(0)
        Vini(0:xl) = V_tmp(0)
        Tini(0:xl) = T_tmp(0)
        Sini(0:xl) = S_tmp(0)
        kini(0:xl) = k_tmp(0)
        epsini(0:xl) = eps_tmp(0)
    else
        call Interp(z_tmp, U_tmp, num, zu, Uini, xli)
        call Interp(z_tmp, V_tmp, num, zu, Vini, xli)
        call Interp(z_tmp, T_tmp, num, zu, Tini, xli)
        call Interp(z_tmp, S_tmp, num, zu, Sini, xli)
        call Interp(z_tmp, k_tmp, num, zk, kini, xli)
        call Interp(z_tmp, eps_tmp, num, zk, epsini, xli)
    end if


    if (disp_dgn/=0) then
        if (disp_dgn==2) then
            do i=num,0,-1
                write(6,999) depth-z_tmp(i),U_tmp(i),V_tmp(i),T_tmp(i),S_tmp(i),k_tmp(i),eps_tmp(i)
            999 format(F8.2,F10.3,F10.3,F10.3,G10.3,G10.3,G10.3)
            end do
        end if
        write(6,*) 'Successfully read'
        write(6,*)
    end if

    return
end


!Prepare parameters for output of the results
!####################################################################
subroutine save_ini()
!####################################################################

    implicit none
    include 'common_parameters.i'

    double precision zu(0:mxl),zk(0:mxl),h(0:mxl),z_zero
    double precision Az(0:mxl),dAdz(0:mxl)
    double precision tout_ctr1(0:9000),tout_ctr2(0:9000)
    integer write_tout

    common /morph_variables1/ h,Az,dAdz
    common /morph_variables2/ zu,zk,z_zero

    common /savet1/ write_tout
    common /savet2/ tout_ctr1
    common /savet3/ tout_ctr2

    double precision t_out(0:9000),test(0:mxl)
    integer i,j,tctr

    ! Read z-vector for vertical output
    if(disp_dgn/=0) write(6,*) 'Output depths  : ',trim(zoutName)
    open(15,status='old',file=zoutName)
    read(15,*)
    do i=0,mxl
        read(15,*,end=8) zsave(i)
    end do
8   if(i==mxl) write(6,*) 'Only first ',mxl,' values of file read.'
    close(15)
    nsave=i-1

    if(OutBin) write(80) igoal
    if(OutBin) write(80) nsave

    if (nsave==0) then             ! Take every nth index for output
        j=0
        depth_save=int(zsave(0))
        do i=xl,0,-depth_save
            j=j+1
            indexk_save(j)=i
            indexu_save(j)=i
        end do
        xlk=j
        xlu=j

        if (indexk_save(xlk)>0) then
            xlk=xlk+1
            indexk_save(xlk)=0
        end if
        if (indexu_save(xlu)>1) then
            xlu=xlu+1
            indexu_save(xlu)=1
        end if

        ! Write vertical vector to file
        if (OutBin) then
            if (abs(igoal)>9) then
                write(80) xlu
                write(80) (zk(xl)-zu(indexu_save(i)),i=1,xlu)
                write(80) xlk
                write(80) (zk(xl)-zk(indexk_save(i)),i=1,xlk)
            elseif (abs(igoal)<5) then
                write(80) xlu
                write(80) (zk(xl)-zu(indexu_save(i)),i=1,xlu)
            elseif (abs(igoal)<=8) then
                write(80) xlk
                write(80) (zk(xl)-zk(indexk_save(i)),i=1,xlk)
            else
                write(6,*) 'Error: wrong choice of index to goal parameter.'
                stop
            end if
        end if
    else                               ! If depth vector is given by user
        if (OutBin) then
            write(80) nsave+1
            write(80) (zsave(i),i=0,nsave)
        end if
        !Output depths are absolute points
        test(0:nsave) = zsave(0:nsave)
        do i=0,nsave
            zsave(nsave-i) = z_zero + test(i)
        end do
    end if
    if (disp_dgn/=0) then
        write(6,*) 'Successfully read'
        write(6,*)
    end if

    ! Read t-vector for temporal output
    if(disp_dgn/=0) write(6,*) 'Output times   : ',trim(toutName)
    open(15,status='old',file=toutName)
    read(15,*)
    do i=0,9000
        read(15,*,end=9) t_out(i)
    end do
9   if(i==9000) write(6,*) 'Only first 9000 values of file read'
    close(15)
    tctr=i-1

    if (tctr/=0) then                  ! Indices to save
        if (t_start>t_out(0)) then
            write(6,*) 'Error: simulation start time is larger than first output time'
            stop
        end if

        tout_ctr1(0)=(t_out(0)-t_start)*86400/dt
        if (int(tout_ctr1(0))==0) then
            tout_ctr2(0)=(t_out(0)-t_start)*86400
            test(0)=t_start+tout_ctr2(0)
        else
            tout_ctr2(0)=((t_out(0)-t_start)*86400)/int(tout_ctr1(0))
            test(0)=t_start
            do j=1,int(tout_ctr1(0))
                test(0)=test(0)+tout_ctr2(0)/86400
            end do
        end if

        do i=1,tctr
            tout_ctr1(i)=(t_out(i)-test(i-1))*86400/dt
            if (int(tout_ctr1(i))==0) then
                tout_ctr2(i)=(t_out(i)-test(i-1))*86400
                test(i)=test(i-1)+tout_ctr2(i)
                write(6,*) 'Warning: time interval for output is smaller than dt for iteration'
            else
                tout_ctr2(i)=((t_out(i)-test(i-1))*86400)/int(tout_ctr1(i))
                test(i)=test(i-1)
                do j=1,int(tout_ctr1(i))
                    test(i)=test(i)+tout_ctr2(i)/86400
                end do
            end if
        end do
        t_end=test(tctr)
        write_tout=0
    else                               ! Interval given for output
        write_tout=int(t_out(0))
    end if
    if (disp_dgn/=0) then
        if(write_tout/=0) write(6,*) 'Interval [days]: ',write_tout*dt/86400
        write(6,*) 'Successfully read'
        write(6,*)
    end if

    return
end


!Set advection to 1 if any inflow/outflow file contains data, otherwise to 0
!AG 2014: revision
!####################################################################
subroutine check_advection()
!####################################################################

    implicit none
    include 'common_parameters.i'

    integer i, j, num_z, fnum(1:4)
    double precision z_Inp(0:mxl), dummy

    if(disp_dgn==2) write(6,*) 'Opening physical inflow/outflow files...'
    open(41,status='old',file=QinpName)
    open(42,status='old',file=QoutName)
    open(43,status='old',file=TinpName)
    open(44,status='old',file=SinpName)

    fnum = [41,42,43,44]
    adv = 0
    do i=1,4
        read(fnum(i),*,end=8)                ! Skip header (description of columns)
        read(fnum(i),*,end=8) num_z                ! Read number of input depths (static)
        read(fnum(i),*,end=8) dummy, (z_Inp(j),j=1,num_z) ! Read input depths
        goto 9
    8   adv = adv + 1
    9   rewind(fnum(i))
    end do

    if (adv==4) then
        adv = 0
    else
        adv = 1
    end if

    return
end
