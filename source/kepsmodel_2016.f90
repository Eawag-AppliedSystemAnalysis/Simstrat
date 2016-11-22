!     +---------------------------------------------------------------+
!     |  k-epsilon model for simulation of                            |
!     |  vertical transport in reservoirs                             |
!     +---------------------------------------------------------------+

program keps

    ! Constant and variables declarations
    implicit none
    include 'common_parameters.i'

    ! Indices: from bottom (0) to surface (xl), max size: mxl
    double precision h(0:mxl)
    double precision Az(0:mxl),dAdz(0:mxl)
    double precision volume
    double precision form_1(0:mxl),form_2(0:mxl)
    double precision form_k1(0:mxl),form_k2(0:mxl),form_beps(0:mxl)
    double precision meanint(0:mxl)

    common /morph_variables1/ h,Az,dAdz
    common /morph_variables3/ volume
    common /form_variables/ form_1,form_2
    common /form_variablesk/ form_k1,form_k2,form_beps
    common /form_variables3/ meanint

    ! Variable Declarations
    !double precision M(0:mxl)

    ! SIMSTRAT Model Initialization
    call Initialization()
    call Form()

    !write(6,*)
    !write(6,*) ' ------------------------- '
    !write(6,*) ' INITIALIZATION SUCCESSFUL '
    !write(6,*) ' ------------------------- '
    !write(6,*)

    !call keps_simulation(M)
    call keps_simulation()

    stop
end

!     +---------------------------------------------------------------+
!     |   Main loop                                                   |
!     +---------------------------------------------------------------+
!subroutine keps_simulation(M)
subroutine keps_simulation()

    implicit none
    include 'common_parameters.i'

    !Initial values
    double precision Sini(0:mxl),Tini(0:mxl),Uini(0:mxl),Vini(0:mxl)
    double precision kini(0:mxl),epsini(0:mxl),Lini(0:mxl),numini(0:mxl),nuhini(0:mxl)
    double precision dragini,txini,tyini

    double precision zu(0:mxl),zk(0:mxl),z_zero !Layer height above bottom (centre, top), total height [m]
    double precision h(0:mxl) !Layer thickness [m]
    double precision Az(0:mxl),dAdz(0:mxl) !Layer horizontal area [m2] and its derivative [m]
    double precision volume !Reservoir volume [m2]
    double precision form_1(0:mxl),form_2(0:mxl)
    double precision form_k1(0:mxl),form_k2(0:mxl),form_beps(0:mxl)
    double precision meanint(0:mxl)

    double precision tout_ctr1(0:9000),tout_ctr2(0:9000)
    integer write_tout

    common /ini_values1/ Sini,Tini
    common /ini_values2/ Uini,Vini
    common /ini_values3/ epsini,kini,Lini
    common /ini_values4/ numini,nuhini
    common /ini_values5/ dragini,txini,tyini

    common /morph_variables1/ h,Az,dAdz
    common /morph_variables2/ zu,zk,z_zero
    common /morph_variables3/ volume
    common /form_variables/ form_1,form_2
    common /form_variablesk/ form_k1,form_k2,form_beps
    common /form_variables3/ meanint

    common /savet1/ write_tout
    common /savet2/ tout_ctr1
    common /savet3/ tout_ctr2

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    LOCAL VARIABLES
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer step,itera,i,j,std
    double precision datum

    double precision S(0:xli),T(0:xli),U(0:xli),V(0:xli) !Salinity, Temperature, Water Velocities
    double precision k(0:xli),ko(0:xli) !Turbulent kinetic energy (TKE) [J/kg]
    double precision eps(0:xli),L(0:xli) !TKE dissipation rate [W/kg]
    double precision num(0:xli),nuh(0:xli) !Turbulent viscosity (momentum) and diffusivity (temperature)
    double precision P(0:xli),B(0:xli),NN(0:xli) !Shear stress production [W/kg], buoyancy production [W/kg], Brunt-Väisälä frequency [s-2]
    double precision dS(0:xli) !Source/sink for salinity
    double precision cmue1(0:xli),cmue2(0:xli) !Model constants
    double precision rho(0:xli)

    double precision P_Seiche(0:xli),E_Seiche !Production of TKE [W/kg] and seiche energy [J]
    double precision gamma !Constant of proportionality for loss of seiche energy

    double precision ga1(0:xli),u10,v10,uv10,Wf !Absorption coeff [m-1], wind drag, wind speeds
    double precision u_taub,drag,u_taus
    double precision tx,ty,SST,heat !Shear stress, sea surface temperature and heat flux
    double precision rad0,rad(0:xli) !Solar radiation (at surface and in water)
    double precision Qvert(0:xli),Q_inp(1:4,0:xli) !Vertical and horizontal flows
    character*3 filelst(1:13)
    !double precision M(0:mxl)

    !Average values
    integer iav
    double precision Uav(0:xli),Vav(0:xli),Tav(0:xli),Sav(0:xli)
    double precision kav(0:xli),epsav(0:xli),numav(0:xli),nuhav(0:xli)
    double precision Bav(0:xli),Pav(0:xli),NNav(0:xli)
    double precision P_Seicheav(0:xli),E_Seicheav

    ! Set initial values (calculated in Initialization)
    U(0:xl)=Uini(0:xl)
    V(0:xl)=Vini(0:xl)
    k(0:xl)=kini(0:xl)
    eps(0:xl)=epsini(0:xl)
    L(0:xl)=Lini(0:xl)
    T(0:xl)=Tini(0:xl)
    S(0:xl)=Sini(0:xl)
    dS(0:xl) = 0.
    num(0:xl)=numini(0:xl)
    nuh(0:xl)=nuhini(0:xl)
    drag=dragini
    tx=txini
    ty=tyini
    E_Seiche=0.
    datum=t_start

    gamma = Az(xl)/(volume**1.5)/sqrt(rho_0)*CD

    if (write_tout==0) then     !if output is written at specific times (irregularly)
        if (tout_ctr2(0)==0) then           ! if first required output time is initial time
            if(disp_sim==1) write(6,990) datum,T(xl),T(xl-5)
            !call write_out(datum,u,v,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,zu,M)
            call write_out(datum,u,v,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,zu)
            itera=0
            step=1
        end if
    end if

    ! Open output files and write state at simulation start
    if (Outbin) then
        call write_out_new(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,zu(0:xl),zk(0:xl),Qvert) !preliminary version!!!
    else
        filelst = ['U  ','V  ','T  ','S  ','k  ','eps','nuh','B  ','P  ','N2 ','Qv ','Ps ','Es ']
        do j=1,13
            open(80+j,access='SEQUENTIAL',action='WRITE',status='unknown',FORM='formatted',&
                      file=trim(PathOut)//trim(filelst(j))//'_out.dat')
            if (j/=13) then
                write(80+j,'(I10,$)'), 1
                do i=0,nsave
                    write(80+j,'(F12.3,$)') zsave(i)-z_zero
                end do
                write(80+j,'(A12,$)'), 'NaN'
            else
                write(93,'(I10,I12,$)'), 1, 1
            end if
        end do
        call write_text(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,zu(0:xl),zk(0:xl),Qvert)
    end if

    std=0
    step=0
    itera=0
    iav=0
    ! START OF SIMULATION LOOP
    do while (datum<t_end)
        std=std+1
        itera=itera+1

        if (write_tout==0) then
            if(itera==1) dt=tout_ctr2(step)
            if(itera==int(tout_ctr1(step))) then
                itera=0
                step=step+1
            end if
        end if

        datum = datum + dt/86400
        call Forcing(datum,T(0:xl),std,tx,ty,u_taus,rad0,heat,SST,u10,v10,uv10,Wf)
        call Absorption(datum,ga1(0:xl),zk(0:xl),std)

        call StabilityFunctions(k(0:xl),eps(0:xl),T(0:xl),S(0:xl),rho(0:xl),meanint(0:xl),NN(0:xl),cmue1(0:xl),cmue2(0:xl))

        if (adv==1) then
            if (ModInflow==1) then
                call Lateral_rho(datum,std,zu(0:xl),zk(0:xl),z_zero,h(0:xl),T(0:xl),S(0:xl),rho(0:xl),Qvert(0:xl),Q_inp(1:4,0:xl)) !zk or zu?? (for interpolation)
            else
                call Lateral(datum,std,zu(0:xl),zk(0:xl),z_zero,Qvert(0:xl),Q_inp(1:4,0:xl)) !zk or zu?? (for interpolation)
            end if
            call Advection(Qvert,Q_inp,U,V,T,S,k,eps,zu(0:xli),zk(0:xli),h(0:xli),Az(0:xli))
        end if
        u_taub=sqrt(drag*(U(1)**2+V(1)**2))

        if (std==1 .and. disp_sim/=0) then
            write(6,*)
            write(6,*) ' -------------------------- '
            write(6,*) '   SIMULATION IN PROGRESS   '
            write(6,*) ' -------------------------- '
            write(6,*)
            if(disp_sim==1 .or. disp_sim==2) write(6,990) 'Time [d]','Surface level [m]','Surface T [degC]','Bottom T [degC]'
            if(disp_sim==3) write(6,990) 'Time [d]','Surface level [m]','Surface U [m/s]','Surface V [m/s]',&
                                         'Surface T [degC]','Surface S [ppt]','Surface k [J/kg]','Surface eps [W/kg]'
        end if

        call Coriolis(U,V)
        call uvEquation(U,V,num,h(0:xl),Az(0:xl),drag,tx,ty,dAdz(0:xl),form_1(0:xl),form_2(0:xl))

        call Temperature(nuh,rad0,rad,h(0:xl),T,heat,SST,ga1,form_1(0:xl),form_2(0:xl))
        if(ModSal/=0) call TransportEquation(S,dS,nuh,form_1(0:xl),form_2(0:xl))

        call Production(U,V,NN,meanint(0:xl),num,nuh,P,B)
        call Seiche(E_Seiche,P_Seiche,u10,v10,Wf,Az(0:xl),dAdz(0:xl),NN,h(0:xl),gamma)
        call TKE(num,P,B,eps,u_taus,u_taub,k,ko,P_Seiche,form_k1(0:xl),form_k2(0:xl))
        call Dissipation(cmue1,cmue2,P,B,k,ko,h(0:xl),Az(0:xl),eps,L,num,nuh,NN,&
                         u_taus,u_taub,P_Seiche,form_k1(0:xl),form_k2(0:xl),form_beps(0:xl))

        !call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)
        !write(6,*) real_clock(2)

        ! Check if averaging and average data
        if (igoal<0) then
            iav=iav+1
            call avstate(iav,U,V,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,Uav,Vav,&
                         Tav,Sav,kav,epsav,numav,nuhav,Bav,Pav,NNav,P_Seicheav,E_Seicheav)
        end if

        ! Write results to file
        if ((itera==write_tout).or.(datum>=t_end)) then
            itera=0
            if (igoal<0) then
                call prep_outav(iav,Uav,Vav,Tav,Sav,kav,epsav,numav,nuhav,Bav,Pav,NNav,P_Seicheav,E_Seicheav)
                call write_out_new(datum,Uav,Vav,Tav,Sav,kav,epsav,nuhav,Bav,Pav,NNav,P_Seicheav,E_Seicheav,zu(0:xl),zk(0:xl),Qvert)
                iav=0
            else
                if (OutBin) then
                    call write_out_new(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,zu(0:xl),zk(0:xl),Qvert)
                else
                    call write_text(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,zu(0:xl),zk(0:xl),Qvert)
                end if
            end if
            if(disp_sim==1) write(6,991) datum,zk(xl),T(xl),T(0)
        end if

        if(disp_sim==2) write(6,991) datum,zk(xl),T(xl),T(0)
        if(disp_sim==3) write(6,991) datum,zk(xl),U(xl),V(xl),T(xl),S(xl),k(xl),eps(xl)

    end do
    ! END OF SIMULATION LOOP

    if (disp_sim/=0) then
        write(6,*)
        write(6,*) ' -------------------------- '
        write(6,*) '    SIMULATION COMPLETED    '
        write(6,*) ' -------------------------- '
        write(6,*)
    end if

!990 format (F12.3, F12.5, F12.5, F12.5)
!991 format (F12.3, F12.5, F12.5, F12.5, F12.5)
990 format (A12, A20, A20, A20, A20, A20, A20, A20)
991 format (F12.3, F20.5, F20.5, F20.5, F20.5, F20.5, F20.5, F20.5)
!999 format (F12.5,F12.5,F12.5,F12.5,F12.5,E12.5,E12.5,F12.5,F12.5)

    stop
end


!     +---------------------------------------------------------------+
!     |  Subroutines                                                  |
!     +---------------------------------------------------------------+

! Tridiagonal matrix algorithm (solver)
!####################################################################
subroutine Tridiagonal(fi,lt,au,bu,cu,du,value)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    ! Upper diagonal, diagonal, lower diagonal, right-hand side, solution
    double precision au(0:xl),bu(0:xl),cu(0:xl),du(0:xl),value(0:xl)
    integer fi,lt !First index, last index

    ! Local variables
    double precision ru(0:xl),qu(0:xl)
    integer i

    ru(lt)=au(lt)/bu(lt)
    qu(lt)=du(lt)/bu(lt)

    do i=lt-1,fi+1,-1
        ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
        qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
    end do

    qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

    value(fi)=qu(fi)
    do i=fi+1,lt
        value(i)=qu(i)-ru(i)*value(i-1)
    end do

    return
end


!####################################################################
subroutine Coriolis(U,V)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision U(0:xli),V(0:xli)

    ! Local variables
    double precision U2(1:xl)

    ! Calculation
    U2(1:xl) = U(1:xl)
    U(1:xl) =  U(1:xl)*cos(Cori*dt) + V(1:xl)*sin(Cori*dt)
    V(1:xl) =-U2(1:xl)*sin(Cori*dt) + V(1:xl)*cos(Cori*dt)

    return
end


!####################################################################
subroutine uvEquation(U,V,num,h,Az,drag,tx,ty,dAdz,form_1,form_2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision U(0:xli),V(0:xli),num(0:xli)
    double precision h(0:xl),Az(0:xl),dAdz(0:xl)
    double precision drag,tx,ty
    double precision form_1(0:xl),form_2(0:xl)

    ! Local variables
    double precision au(0:xl),bu(0:xl),cu(0:xl),du(0:xl)
    double precision intU, intV, length

    if (Pgrad==1) then
        intU = sum(U(1:xl))
        intV = sum(V(1:xl))
        length = sqrt(Az(xl))
    end if

    ! Build the diagonals
    cu(1:xl-1) = dt*num(1:xl-1)*form_2(1:xl-1)
    cu(xl) = 0.
    au(1) = 0.
    au(2:xl) = dt*num(1:xl-1)*form_1(2:xl)
    bu(1:xl) = 1-au(1:xl)-cu(1:xl)
    bu(1) = bu(1) + drag*sqrt(U(1)**2+V(1)**2)/h(1)*dt

    ! Calculation of U-equation
    du(1) = U(1)
    if (Pgrad==1) then !Svensson 1978
        !du(2:xl-1) = U(2:xl-1) - 10*pi**2*intU/xl*rho_0*g*depth/length**2*dt/86400
        du(2:xl-1) = U(2:xl-1) - pi**2*rho_0*g*intU/xl*depth/length**2*dt
    elseif (Pgrad==2) then !???
        du(2:xl-1) = U(2:xl-1) - drag*U(2:xl-1)*sqrt(U(2:xl-1)**2+V(2:xl-1)**2)*dAdz(2:xl-1)/Az(2:xl-1)*dt
    else
        du(2:xl-1) = U(2:xl-1)
    end if
    du(xl) = U(xl) + tx*dt/h(xl)
    ! Solve
    call Tridiagonal(1,xl,au,bu,cu,du,U)
    ! Calculation of V-equation
    du(1) = V(1)
    if (Pgrad==1) then !Svensson 1978
        !du(2:xl-1) = V(2:xl-1) - 10*pi**2*intV/xl*rho_0*g*depth/length**2*dt/86400
        du(2:xl-1) = V(2:xl-1) - pi**2*rho_0*g*intV/xl*depth/length**2*dt
    elseif (Pgrad==2) then !???
        du(2:xl-1) = V(2:xl-1) - drag*V(2:xl-1)*sqrt(U(2:xl-1)**2+V(2:xl-1)**2)*dAdz(2:xl-1)/Az(2:xl-1)*dt
    else
        du(2:xl-1) = V(2:xl-1)
    end if
    du(xl) = V(xl) + ty*dt/h(xl)
    ! Solve
    call Tridiagonal(1,xl,au,bu,cu,du,V)

    return
end


!     +----------------------------------------------------------------+
!     | Equation of state                                              |
!     +----------------------------------------------------------------+
!####################################################################
subroutine Buoyancy(T,S,meanint,rho,NN)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision T(0:xl),S(0:xl),NN(0:xl),meanint(0:xl)

    ! Local variables
    double precision a(0:xl),rho(0:xl),buoy(0:xl)
    double precision rho0t(0:xl),rho0st(0:xl)
    integer i

    if (salctr==0) then                ! salinity is zero everywhere
        do i=1,xl-1
            a(i)= -68.0+T(i)*(18.2091+T(i)*(-0.30866+T(i)*&
                 (5.3445e-3+T(i)*(-6.0721e-5+T(i)*(3.1441e-7)))))
            !if (press/=0) then   ! ignore this pressure thing for alpha in first approximation
            !    a(i)= a(i) + (0.3682+T(i)*(-1.520e-2+T(i)*(1.91e-4)))*p(i)
            !end if
            a(i)= 1.0e-6*a(i)
            NN(i)= g*a(i)*(meanint(i)*(T(i)-T(i+1))+(T(i)+273.15)/cp)
        end do
    else
        if (delsal==0) then        ! salinity gradient is zero everywhere
            do i=1,xl-1
                a(i)= -68.0+T(i)*(18.2091+T(i)*(-0.30866+T(i)*&
                     (5.3445e-3+T(i)*(-6.0721e-5+T(i)*(3.1441e-7)))))
                a(i)= a(i) + (4.599+T(i)*(-0.1999+T(i)*(2.79e-3)))*S(i)
                !if (press/=0) then   ! ignore this pressure thing for alpha in first approximation
                !    a(i) = a(i) + (0.3682+T(i)*(-1.520e-2+T(i)*(1.91e-4))-S(i)*(4.613e-3))*p(i)
                !end if
                a(i)= 1.0e-6*a(i)
                NN(i)= g*a(i)*(meanint(i)*(T(i)-T(i+1))+(T(i)+273.15)/cp)
            end do
        else
            !if (press/=0) then
            !    rho0t= (0.9998395+t.*(6.7914e-5 +t.*(-9.0894e-6+t.*(1.0171e-7+t.*(-1.2846e-9 +t.*(1.1592e-11
            !           +t.*(-5.0125e-14)))))));
            !    kbart= (19652.17 +t.*(148.113 +t.*(-2.293 +t.*(1.256e-2+t.*(-4.18e-5)))));
            !    kbarpt= ((3.2726 +t.*(-2.147e-4 +t.*(1.128e-4))).*p);
            !    kbar= kbart+kbarpt
            !    rho0st= ((8.181e-4 +t.*(-3.85e-6 +t.*(4.96e-8))).*s);
            !    kbarspt= ((53.238 -0.313.*t +5.728e-3.*p).*s);
            !    kbar= kbar+kbarspt;
            !    rho= rho_0*(rho0t + rho0st)./(1.0 - p./kbar);
            !    if ~isempty(f)
            !        rhoS0= rho0t./(1.0 - p./(kbart+kbarpt));
            !        rho= rhoS0.*(1-fc)+fc.*rho;
            !    end
            !    rho= rho_0*rho0t./(1.0 - p./kbar);
            !else
            do i=0,xl
                rho0t(i)= 0.9998395+T(i)*(6.7914e-5+T(i)*(-9.0894e-6+T(i)*&
                    (1.0171e-7+T(i)*(-1.2846e-9+T(i)*(1.1592e-11+T(i)*(-5.0125e-14))))))
                rho0st(i)= (8.181e-4+T(i)*(-3.85e-6+T(i)*(4.96e-8)))*S(i)
                rho(i)= rho_0*(rho0t(i)+rho0st(i))
                !if (fc/=0) then
                !    rho(i) = rho0t(i)*(1-fc)+fc*rho(i)
                !end
                buoy(i)= -g*(rho(i)-rho_0)/rho_0
            end do

            NN(1:xl-1) = meanint(1:xl-1)*(buoy(2:xl)-buoy(1:xl-1))
        end if
    end if

    NN(0)= NN(1)
    NN(xl)= NN(xl-1)

    return
end


!####################################################################
subroutine Temperature(nuh,rad0,rad,h,T,heat,SST,ga1,form_1,form_2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision nuh(0:xl),rad0,rad(0:xl),h(0:xl),T(0:xl),heat,SST,ga1(0:xl)
    double precision form_1(0:xl),form_2(0:xl)

    ! Local variables
    double precision au(0:xl),bu(0:xl),cu(0:xl),du(0:xl)
    integer i

    ! Calculation
    !Radiation reaching each layer
    rad(xl) = rad0/rho_0/cp ![°C*m/s]
    do i=xl-1,0,-1
        rad(i) = rad(i+1)*exp(-h(i)*ga1(xl-i)) !Attenuated by absorption
    end do

    !Build the three diagonals
    au(1) = 0.
    au(2:xl) = dt*nuh(1:xl-1)*form_1(2:xl)
    cu(1:xl-1) = dt*nuh(1:xl-1)*form_2(1:xl-1)
    cu(xl) = 0.
    bu(1:xl) = 1-au(1:xl)-cu(1:xl)
    du(1:xl) = T(1:xl)+(rad(1:xl)-rad(0:xl-1))/h(1:xl)*dt
    du(xl) = du(xl) + heat/rho_0/cp*dt/h(xl)

    if (NBC==1) then
        bu(xl) = 1.
        au(xl) = 0.
        cu(xl) = 0.
        du(xl) = SST
    end if

    !Add geothermal heat flux
    if(fgeo/=0) du(1:xl)=du(1:xl)+fgeo_add(1:xl)*dt

    !Solve
    call Tridiagonal(1,xl,au,bu,cu,du,T)
    T(0) = T(1)     ! set value to boundary value at zu(1)

    return
end


!####################################################################
subroutine TransportEquation(C,dC,nuh,form_1,form_2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision C(0:xl),dC(0:xl),nuh(0:xl)
    double precision form_1(0:xl),form_2(0:xl)

    ! Local variables
    double precision au(0:xl),bu(0:xl),cu(0:xl),du(0:xl)
    !double precision SalRel
    !SalRel = 172800.

    !Build the three diagonals
    cu(1:xl-1) = dt*nuh(1:xl-1)*form_2(1:xl-1)
    cu(xl) = 0.
    au(1) = 0.
    au(2:xl) = dt*nuh(1:xl-1)*form_1(2:xl)
    bu(1:xl) = 1.-au(1:xl)-cu(1:xl)
    du(1:xl) = C(1:xl) + dC(1:xl)*dt !AG 2014: added dC*dt term for source/sink

    !Solve
    call Tridiagonal(1,xl,au,bu,cu,du,C)
    C(0) = C(1)     ! set value to boundary value at zu(1)

    return
end


! Calculation of cmue
!####################################################################
subroutine cmue_cn(cmue1,cmue2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    double precision cmue1,cmue2

    !Standard version of k-eps model
    cmue1 = cmue
    cmue2 = cmue/Prndtl
    !Burchard Version
    !cmue1=cde*cm0
    !cmue2=cde*cm0/Prndtl

    return
end

! Calculation of cmue
!####################################################################
subroutine cmue_qe(beta,cmue1,cmue2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    double precision beta,cmue1,cmue2
    double precision gh,sm,sh

    gh = -cde**2*0.5*beta
    if(gh> 0.02) gh = gh-(gh-0.02)**2/(gh+0.0233-2*0.02)
    if(gh<-0.28) gh = -0.28

    sm = 1.-3*c1-6*a1/b1-3*a2*gh*((b2-3*a2)*(1.-6*a1/b1)-3*c1*(b2+6*a1))
    sm = a1*sm/((1.-3*a2*gh*(6*a1+b2))*(1.-9*a1*a2*gh))
    sh = a2*(1.-6*a1/b1)/(1.-3*a2*gh*(6*a1+b2))

    cmue1 = sqrt(2.)*cde*sm
    cmue2 = sqrt(2.)*cde*sh

    return
end


!####################################################################
subroutine StabilityFunctions(k,eps,T,S,rho,meanint,NN,cmue1,cmue2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision T(0:xl),S(0:xl),rho(0:xl)
    double precision k(0:xl),eps(0:xl),NN(0:xl)
    double precision cmue1(0:xl),cmue2(0:xl)
    double precision meanint(0:xl)

    ! Local variables
    double precision beta
    integer i

    call Buoyancy(T,S,meanint,rho,NN)

    !In the domain
    do i=1,xl-1
        if (stab==1) then
            call cmue_cn(cmue1(i),cmue2(i))
        else if (stab==2) then
            beta = NN(i)*(k(i)/eps(i))**2
            call cmue_qe(beta,cmue1(i),cmue2(i))
        end if
    end do
    !At the boundaries
    cmue1(0) = cmue1(1)
    cmue2(0) = cmue2(1)
    cmue1(xl) = cmue1(xl-1)
    cmue2(xl) = cmue2(xl-1)

    return
end


!####################################################################
subroutine Seiche(E_Seiche,P_Seiche,u10,v10,Wf,Az,dAdz,NN,h,gamma)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision h(0:xl),Az(0:xl),dAdz(0:xl)
    double precision NN(0:xl),minNN,gamma
    double precision P_Seiche(0:xl),E_Seiche
    double precision u10,v10,Wf

    ! Local variables
    double precision W10,PS,PW,f_norm
    double precision Distrib(0:xl)
    integer i

    minNN = 0
    if (a_seiche/=0.) then
        do i=1,xl-1
            Distrib(i) = max(NN(i)**q_NN,minNN) / Az(i)*dAdz(i)
        end do
        !Seiche normalization factor
        f_norm = 0.
        if (ModSNorm==1) then !max NN
            do i=1,xl-1
                if(NN(i)>f_norm) f_norm = NN(i)
            end do
            f_norm = (f_norm**q_NN)*Az(xl)*rho_0
        else if (ModSNorm==2) then !integral
            do i=1,xl-1
                f_norm = f_norm+Distrib(i)*Az(i)*h(i)
            end do
            f_norm = f_norm*rho_0
        end if

        if (f_norm==0.) then
            do i=1,xl-1
                Distrib(i) = 1/h(i)
                !Distrib(i) = 1/h(i)*max(NN(i)**q_NN,minNN)
                !f_norm = f_norm + max(NN(i)**q_NN,minNN)
            end do
            f_norm= Az(xl)*rho_0
            !f_norm=f_norm*Az(xl)*rho_0
        end if

        !Seiche energy per surface converted to total seiche energy
        if (WindFilt==1) then !use filtered wind (AG 2014)
            PW = a_seiche*Az(xl)*rho_air*C10*Wf**3
        else !use real wind
            W10 = sqrt(u10**2+v10**2)
            PW = a_seiche*Az(xl)*rho_air*C10*W10**3
        end if
        PS = E_Seiche**(1.5)*gamma

        E_Seiche = E_Seiche + (PW - PS)*dt

        !Limit so that E_Seiche does not become negative
        if (E_Seiche<0.) then
            PS = (PS*dt+E_Seiche)/dt
            E_Seiche = 0.
        end if

        !write(6,*) E_Seiche,f_norm
        do i=1,xl-1
            P_Seiche(i) = 1./f_norm*Distrib(i)*PS*(1.-10*sqrt(CD))
        end do
        Distrib(0) = 0.
        Distrib(xl) = 0.

    else          !if alpha==0
        P_seiche(0:xl) = 0.
    end if

    return
end


!####################################################################
subroutine Production(U,V,NN,meanint,num,nuh,P,B)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision U(0:xl),V(0:xl),NN(0:xl),meanint(0:xl)
    double precision num(0:xl),nuh(0:xl),P(0:xl),B(0:xl)

    P(1:xl-1) = (U(2:xl)-U(1:xl-1))**2+(V(2:xl)-V(1:xl-1))**2
    P(1:xl-1) = P(1:xl-1)*num(1:xl-1)*meanint(1:xl-1)**2
    B(1:xl-1) = -nuh(1:xl-1)*NN(1:xl-1)

    return
end


!####################################################################
subroutine TKE(num,P,B,eps,u_taus,u_taub,k,ko,P_Seiche,form_k1,form_k2)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision num(0:xl),P(0:xl),B(0:xl),k(0:xl),eps(0:xl),ko(0:xl),P_Seiche(0:xl)
    double precision u_taus,u_taub
    double precision form_k1(0:xl),form_k2(0:xl)

    ! Local variables
    double precision avh(0:xl),au(0:xl),bu(0:xl),cu(0:xl),du(0:xl)
    double precision pminus(0:xl),pplus(0:xl),Prod,Buoy,Diss
    integer i

    ko(0:xl) = k(0:xl) ! ko = TKE at old time step

    avh(2:xl-1) = 0.5/sig_k*(num(1:xl-2)+num(2:xl-1)) ! average num for TKE

    if ((ModFlux==1).and.(Mod==1)) then
        avh(1) = 0.
        avh(xl) = 0.
    else
        avh(1)=2*u_taub**4/(eps(0)+eps(1))        ! = 0 for no shear stress
        avh(xl)=2*u_taus**4/(eps(xl)+eps(xl-1))   ! = 0 for no shear stress
    end if

    do i=1,xl-1
        Prod = P(i)+P_Seiche(i)                   ! Add seiche energy
        Buoy=B(i)
        Diss=eps(i)
        if (Prod+Buoy>0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
        else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
        end if
    end do

    au(1:xl-1) = dt*avh(1:xl-1)*form_k1(1:xl-1)
    cu(1:xl-1)=dt*avh(2:xl)*form_k2(1:xl-1)
    bu(1:xl-1)=1.-au(1:xl-1)-cu(1:xl-1)+pminus(1:xl-1)*dt/k(1:xl-1)
    du(1:xl-1)=k(1:xl-1)+pplus(1:xl-1)*dt

    if ((ModFlux==1).and.(Mod==1)) then
        call Tridiagonal(1,xl-1,au,bu,cu,du,k)
        k(0)= k(1)                                ! Define TKE at boundary
        k(xl)= k(xl-1)                            ! no-flux condition
    else
        cu(0)= 0.
        bu(0)= 1.
        du(0)= u_taub**2/sqrt(cm0*cde)

        bu(xl)= 1.
        au(xl)= 0.
        cu(xl)= 0.
        du(xl)= u_taus**2/sqrt(cm0*cde)

        call Tridiagonal(0,xl,au,bu,cu,du,k)
    end if

    do i=0,xl
        if(k(i)<k_min) k(i)=k_min             ! Lower limit of TKE
    end do

    return
end


!####################################################################
subroutine Dissipation(cmue1,cmue2,P,B,k,ko,h,Az,eps,L,num,nuh,&
           NN,u_taus,u_taub,P_Seiche,form_k1,form_k2,form_beps)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global variables
    double precision cmue1(0:xl),cmue2(0:xl),P(0:xl),B(0:xl),h(0:xl),Az(0:xl)
    double precision num(0:xl),nuh(0:xl),k(0:xl),ko(0:xl),eps(0:xl)
    double precision L(0:xl),NN(0:xl),P_Seiche(0:xl)
    double precision u_taus,u_taub
    double precision form_k1(0:xl),form_k2(0:xl),form_beps(0:xl)

    ! Local variables
    double precision avh(0:xl),au(0:xl),bu(0:xl),cu(0:xl),du(0:xl)
    double precision flux(0:xl)
    double precision pminus(0:xl),pplus(0:xl),Prod,Buoy,Diss,cee3
    double precision epslim
    integer i

    avh(1:xl) = 0.5/sig_e*(num(0:xl-1)+num(1:xl)) ! Average num for Diss

    if ((ModFlux==1).and.(Mod==1)) then
        flux(0 ) = avh(1 )*(cde*((ko(1   ))**1.5)/(kappa*(K_s+0.5*h(1 )))**2)
        flux(xl) = avh(xl)*(cde*((ko(xl-1))**1.5)/(kappa*(z0+0.5*h(xl)))**2)
        do i=1,xl-1
            flux(i) = num(i)/sig_e*(cde*((ko(i))**1.5)/(kappa*(z0+0.25*(h(i)+h(i+1))))**2)
            !flux(i) = -flux(i)*eps(i)/ko(i)
        end do
        avh(1 ) = 0
        avh(xl) = 0
    else
        avh(1 ) = 2*u_taub**4/sig_e/(eps(0 )+eps(1   ))   ! = 0 for no shear stress
        avh(xl) = 2*u_taus**4/sig_e/(eps(xl)+eps(xl-1))   ! = 0 for no shear stress
    end if

    do i=1,xl-1
        if (B(i)>0) then
            cee3=1.
        else
            cee3=ce3
        end if
        Prod=ce1*eps(i)/ko(i)*(P(i)+P_Seiche(i))        ! New code plus seiche
        Buoy=cee3*eps(i)/ko(i)*B(i)
        Diss=ce2*eps(i)*eps(i)/ko(i)
        if (Prod+Buoy>0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
        else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
        end if
    end do

    au(1:xl-1) = dt*avh(1:xl-1)*form_k1(1:xl-1)
    cu(1:xl-1)=dt*avh(2:xl)*form_k2(1:xl-1)
    bu(1:xl-1)=1.-au(1:xl-1)-cu(1:xl-1)+pminus(1:xl-1)*dt/eps(1:xl-1)
    du(1:xl-1)=eps(1:xl-1)+pplus(1:xl-1)*dt

    if ((ModFlux==1).and.(Mod==1)) then
        du(1:xl-1) = du(1:xl-1)+flux(1:xl-1)*dt*form_beps(1:xl-1) ! form_beps = 1/A * dA/dz (at epsilon posions)
        if (Az(0)/=0) then                          ! Flux from bottom only!
            du(1)= du(1)+flux(0)*dt*(Az(0)+Az(1))/(Az(1)*(h(1)+h(2)))
        end if
        du(xl-1)= du(xl-1)+flux(xl)*dt*(Az(xl)+Az(xl-1))/(Az(xl-1)*(h(xl)+h(xl-1)))
        call Tridiagonal(1,xl-1,au,bu,cu,du,eps)
        ! Define eps at boundaries
        eps(0 )= eps(1   )+(cde*((ko(1   ))**1.5)/(kappa*(K_s+h(1 )))**2)*h(1 )
        eps(xl)= eps(xl-1)+(cde*((ko(xl-1))**1.5)/(kappa*(z0 +h(xl)))**2)*h(xl)
    else
        cu(0)= 0.
        bu(0)= 1.
        du(0)= cde*sqrt(k(0)*k(0)*k(0))/kappa/K_s

        bu(xl)= 1.
        au(xl)= 0.
        du(xl)= cde*sqrt(k(xl)*k(xl)*k(xl))/kappa/z0

        call Tridiagonal(0,xl,au,bu,cu,du,eps)
    end if

    do i=0,xl
        if (NN(i)>0) then
            epslim= 0.212*k(i)*sqrt(NN(i))
        else
            epslim= eps_min
        end if
        if(eps(i)<epslim) eps(i)=epslim
        if (eps(i)<0) then
            write(6,*) 'Dissipation negative'
        end if

        num(i)= cmue1(i)*k(i)*k(i)/eps(i)+1.5e-6
        nuh(i)= cmue2(i)*k(i)*k(i)/eps(i)+1.5e-7
        L(i)= cde*sqrt(ko(i)*ko(i)*ko(i))/eps(i)
    end do

    num(0 )= kappa*u_taub*K_s+avh_min
    num(xl)= kappa*u_taus*z0 +avh_min

    return
end


!####################################################################
subroutine Advection(Qvert,Q_inp,U,V,T,S,k,eps,zu,zk,h,Az)
!####################################################################

    implicit none
    include 'common_parameters.i'

    !Global Declarations
    double precision h(0:xli), Az(0:xli), zu(0:xli), zk(0:xli)
    double precision U(0:xli), V(0:xli), T(0:xli), S(0:xli)
    double precision k(0:xli), eps(0:xli)
    double precision Qvert(0:xli) ! Vertical advection [m3/s]
    double precision Q_inp(1:4,0:xli) ! Inflow [m3/s], Outflow [m3/s], T-input [°C*m3/s], S-input [‰*m3/s]

    double precision dh, dhi(1:2)      ! depth differences
    double precision dti(1:2)          ! first and second time step
    double precision dU(0:xli), dV(0:xli), dTemp(0:xli), dS(0:xli)
    double precision form_adv(1:xli)
    double precision top
    integer i, ti

    !Depth difference
    dh = Qvert(xl)/Az(xl)*dt
    !Split timestep depending on situation
    if (dh==0.) then                          ! If volume does not change
        dti(1) = dt                                 ! One normal time step
    else if ((dh+zk(xl))>=depth) then        ! If surface level reached
        dti(1) = (depth - zk(xl))/dh*dt             ! Step until surface level
    else if (((dh+h(xl))>h(xl-1)/2) .and.&   ! If top box>0.5*lower box
            ((dh+h(xl))<2*h(xl-1))) then     ! and top box<2*lower box
        dti(1) = dt                                 ! One normal time step
    else if ((dh+h(xl))<=h(xl-1)/2) then     ! If top box<=0.5*lower box
        dti(1) = abs((h(xl)-h(xl-1)/2)/dh)*dt       ! First step until top box = 0.5*lower box
    else                                     ! If top box>=2*lower box
        dti(1) = abs((2*h(xl-1)-h(xl))/dh)*dt       ! First step until top box = 2*lower box
    end if
    dti(2) = dt-dti(1)                       ! Rest of timestep

    do ti=1,2 !First and (if needed) second timestep
        form_adv(1:xl) = dti(ti)/(Az(1:xl)*h(1:xl))     ! Form factor for dt(ti)
        dhi(ti) = dh*dti(ti)/dt                         ! Depth difference for dt(ti)

        do i=1,xl
            dU(i)    = 0
            dV(i)    = 0
            dTemp(i) = 0
            dS(i)    = 0
            if (i<xl .or. Qvert(xl)<0) then        ! Advective flow out of box i, always negative
                dU(i)    = -abs(Qvert(i))*U(i)
                dV(i)    = -abs(Qvert(i))*V(i)
                dTemp(i) = -abs(Qvert(i))*T(i)
                dS(i)    = -abs(Qvert(i))*S(i)
            end if
            if (i>1 .and. Qvert(i-1)>0) then       ! Advective flow into box i, from above
                dU(i)    = dU(i) + Qvert(i-1)*U(i-1)
                dV(i)    = dV(i) + Qvert(i-1)*V(i-1)
                dTemp(i) = dTemp(i) + Qvert(i-1)*T(i-1)
                dS(i)    = dS(i) + Qvert(i-1)*S(i-1)
            end if
            if (i<xl .and. Qvert(i+1)<0) then      ! Advective flow into box i, from below
            !if (i<xl .and. Qvert(i)<0) then       ! Advective flow into box i, from below
                dU(i)    = dU(i) - Qvert(i+1)*U(i+1)
                dV(i)    = dV(i) - Qvert(i+1)*V(i+1)
                dTemp(i) = dTemp(i) - Qvert(i+1)*T(i+1)
                dS(i)    = dS(i) - Qvert(i+1)*S(i+1)
            end if
        end do

        ! Inflow and outflow
        dTemp(1:xl) = dTemp(1:xl) + (Q_inp(3,1:xl)+Q_inp(2,1:xl)*T(1:xl))
        dS(1:xl) = dS(1:xl) + (Q_inp(4,1:xl)+Q_inp(2,1:xl)*S(1:xl))

        ! Take first time step
        U(1:xl) = U(1:xl) + form_adv(1:xl)*dU(1:xl)
        V(1:xl) = V(1:xl) + form_adv(1:xl)*dV(1:xl)
        T(1:xl) = T(1:xl) + form_adv(1:xl)*dTemp(1:xl)
        S(1:xl) = S(1:xl) + form_adv(1:xl)*dS(1:xl)
        !U(1:xl) = (1-form_adv(1:xl)*Qvert(1:xl))*U(1:xl) + form_adv(1:xl)*dU(1:xl)
        !V(1:xl) = (1-form_adv(1:xl)*Qvert(1:xl))*V(1:xl) + form_adv(1:xl)*dV(1:xl)
        !T(1:xl) = (1-form_adv(1:xl)*Qvert(1:xl))*T(1:xl) + form_adv(1:xl)*dTemp(1:xl)
        !S(1:xl) = (1-form_adv(1:xl)*Qvert(1:xl))*S(1:xl) + form_adv(1:xl)*dS(1:xl)

        ! Variation of variables due to change in volume
        U(xl) = U(xl)*h(xl)/(h(xl)+dhi(ti))
        V(xl) = V(xl)*h(xl)/(h(xl)+dhi(ti))
        T(xl) = T(xl)*h(xl)/(h(xl)+dhi(ti))
        S(xl) = S(xl)*h(xl)/(h(xl)+dhi(ti))

        if (ti==1) then
            if (dh==0) then                          ! If volume does not change
                return
            else if ((dh+zk(xl))>=depth) then        ! If surface level reached
                h(xl) = h(xl) + dhi(1)                      ! New thickness of top box
                zu(xl) = zu(xl) + dhi(1)/2                  ! New centre coordinate of top box
                zk(xl) = zk(xl) + dhi(1)                    ! New surface coordinate of top box
                return                                      ! No change in volume (overflow)
            else if (((dh+h(xl))>h(xl-1)/2) .and.&   ! If top box>0.5*lower box
                     ((dh+h(xl))<2*h(xl-1))) then    ! and top box<2*lower box
                h(xl) = h(xl) + dhi(1)                      ! New thickness of top box
                zu(xl) = zu(xl) + dhi(1)/2                  ! New centre coordinate of top box
                zk(xl) = zk(xl) + dhi(1)                    ! New surface coordinate of top box
                return
            else if ((dh+h(xl))<=h(xl-1)/2) then     ! If top box<=0.5*lower box
                zk(xl-1) = zk(xl)
                zu(xl-1) = (zk(xl-1)+zk(xl-2))/2
                h(xl-1)  = (zk(xl-1)-zk(xl-2))
                U(xl-1) = (0.5*U(xl)*Az(xl)+U(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                V(xl-1) = (0.5*V(xl)*Az(xl)+V(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                T(xl-1) = (0.5*T(xl)*Az(xl)+T(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                S(xl-1) = (0.5*S(xl)*Az(xl)+S(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                k(xl-1) = (0.5*k(xl)*Az(xl)+k(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                eps(xl-1) = (0.5*eps(xl)*Az(xl)+eps(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                Qvert(xl-1) = (0.5*Qvert(xl)*Az(xl)+Qvert(xl-1)*Az(xl-1))/(0.5*Az(xl)+Az(xl-1))
                xl = xl-1                     ! Reduce number of boxes
                !dh(2) = Qvert(xl)/Az(xl)*dt(2) !(also equal to dh*dt2/dt)
                call Form()
            else                                     ! If top box>=2*lower box
                h(xl+1)   = h(xl)/2 + dh   ! AG 2014 (added +dh)
                h(xl)     = h(xl)/2
                zk(xl+1)  = zk(xl) + dh    ! AG 2014 (added +dh)
                zk(xl)    = zk(xl) - h(xl)/2
                zu(xl+1)  = zu(xl) + h(xl+1)/2
                zu(xl)    = zu(xl) - h(xl+1)/2
                U(xl+1)   = U(xl)
                V(xl+1)   = V(xl)
                T(xl+1)   = T(xl)
                S(xl+1)   = S(xl)
                k(xl+1)   = k(xl)
                eps(xl+1) = eps(xl)
                Qvert(xl+1) = Qvert(xl)       ! Vertical discharge of new box
                xl = xl + 1                   ! Increase number of boxes
                !dh(2) = Qvert(xl)/Az(xl)*dt(2) !(also equal to dh*dt2/dt)
                call Form()
            end if
        end if      !end if (ti==1)
    end do      !end do ti=1,2

    return
end


!####################################################################
subroutine Form()
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global Declarations
    double precision h(0:mxl),Az(0:mxl),dAdz(0:mxl)
    double precision form_1(0:mxl),form_2(0:mxl)
    double precision form_k1(0:mxl),form_k2(0:mxl),form_beps(0:mxl)
    double precision meanint(0:mxl),volume

    common /morph_variables1/ h,Az,dAdz
    common /morph_variables3/ volume
    common /form_variables/ form_1,form_2
    common /form_variablesk/ form_k1,form_k2,form_beps
    common /form_variables3/ meanint

    integer i

    form_1(1:xl) = -4*Az(0:xl-1)/(h(1:xl)+h(0:xl-1))/h(1:xl)/(Az(1:xl)+Az(0:xl-1))
    form_2(1:xl) = -4*Az(1:xl)/(h(1:xl)+h(2:xl+1))/h(1:xl)/(Az(1:xl)+Az(0:xl-1))
    form_k1(1:xl-1) = -(Az(1:xl-1)+Az(2:xl))/(h(1:xl-1)+h(2:xl))/h(2:xl)/Az(1:xl-1)
    form_k2(1:xl-1) = -(Az(1:xl-1)+Az(0:xl-2))/(h(1:xl-1)+h(2:xl))/h(1:xl-1)/Az(1:xl-1)
    form_beps(1:xl-1) = 0.5*((Az(1:xl-1)-Az(0:xl-2))/h(1:xl-1)+(Az(2:xl)-Az(1:xl-1))/h(2:xl))/Az(1:xl-1)

    meanint(0:xl-1) = 2./(h(0:xl-1)+h(1:xl))

    volume=0
    do i=0,xl-1
        volume = volume + 0.5*h(i+1)*(Az(i)+Az(i+1))
    end do

    return
end


!####################################################################
subroutine avstate(iav,U,V,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,Uav,Vav,&
                   Tav,Sav,kav,epsav,numav,nuhav,Bav,Pav,NNav,P_Seicheav,E_Seicheav)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global Declarations
    double precision U(0:xl),V(0:xl),T(0:xl),S(0:xl)
    double precision k(0:xl),eps(0:xl),num(0:xl),nuh(0:xl)
    double precision B(0:xl),P(0:xl),NN(0:xl)
    double precision P_Seiche(0:xl),E_Seiche

    double precision Uav(0:xl),Vav(0:xl),Tav(0:xl),Sav(0:xl)
    double precision kav(0:xl),epsav(0:xl),numav(0:xl),nuhav(0:xl)
    double precision Bav(0:xl),Pav(0:xl),NNav(0:xl)
    double precision P_Seicheav(0:xl),E_Seicheav
    integer iav

    if (igoal_s(1)==1) then
        if (iav==1) then
            Uav(0:xl) = U(0:xl)
        else
            Uav(0:xl) = Uav(0:xl) + U(0:xl)
        end if
    end if
    if (igoal_s(2)==1) then
        if (iav==1) then
            Vav(0:xl) = V(0:xl)
        else
            Vav(0:xl) = Vav(0:xl) + V(0:xl)
        end if
    end if
    if (igoal_s(3)==1) then
        if (iav==1) then
            Tav(0:xl) = T(0:xl)
        else
            Tav(0:xl) = Tav(0:xl) + T(0:xl)
        end if
    end if
    if (igoal_s(4)==1) then
        if (iav==1) then
            Sav(0:xl) = S(0:xl)
        else
            Sav(0:xl) = Sav(0:xl) + S(0:xl)
        end if
    end if
    if (igoal_s(5)==1) then
        if (iav==1) then
            kav(0:xl) = k(0:xl)
        else
            kav(0:xl) = kav(0:xl) + k(0:xl)
        end if
    end if
    if (igoal_s(6)==1) then
        if (iav==1) then
            epsav(0:xl) = eps(0:xl)
        else
            epsav(0:xl) = epsav(0:xl) + eps(0:xl)
        end if
    end if
    if (igoal_s(7)==1) then
        if (iav==1) then
            numav(0:xl) = num(0:xl)
        else
            numav(0:xl) = numav(0:xl) + num(0:xl)
        end if
    end if
    if (igoal_s(8)==1) then
        if (iav==1) then
            nuhav(0:xl) = nuh(0:xl)
        else
            nuhav(0:xl) = nuhav(0:xl) + nuh(0:xl)
        end if
    end if
    if (igoal_s(9)==1) then
        if (iav==1) then
            Bav(0:xl) = B(0:xl)
        else
            Bav(0:xl) = Bav(0:xl) + B(0:xl)
        end if
    end if
    if (igoal_s(10)==1) then
        if (iav==1) then
            Pav(0:xl) = P(0:xl)
        else
            Pav(0:xl) = Pav(0:xl) + P(0:xl)
        end if
    end if
    if (igoal_s(11)==1) then
        if (iav==1) then
            P_Seicheav(0:xl) = P_Seiche(0:xl)
        else
            P_Seicheav(0:xl) = P_Seicheav(0:xl) + P_Seiche(0:xl)
        end if
    end if
    if (igoal_s(12)==1) then
        if (iav==1) then
            NNav(0:xl) = NN(0:xl)
        else
            NNav(0:xl) = NNav(0:xl) + NN(0:xl)
        end if
    end if
    if (igoal_s(13)==1) then
        if (iav==1) then
            E_Seicheav = E_Seiche
        else
            E_Seicheav = E_Seicheav + E_Seiche
        end if
    end if

    return
end


!####################################################################
!subroutine write_out(datum,u,v,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,zu,M)
subroutine write_out(datum,u,v,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,zu)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global Declarations
    double precision datum
    double precision U(0:xl),V(0:xl),T(0:xl),S(0:xl)
    double precision k(0:xl),eps(0:xl),num(0:xl),nuh(0:xl)
    double precision B(0:xl),P(0:xl),NN(0:xl)
    double precision P_Seiche(0:xl),E_Seiche
    double precision zu(0:mxl)
    !double precision M(0:mxl)

    ! Local Declarations
    double precision xs(0:mxl)
    integer i!,ist

    write(80) datum
    if (nsave==0) then
        if (igoal_s(1)==1) then
            write(80) (U(indexu_save(i)),i=1,xlu)
        end if
        if (igoal_s(2)==1) then
            write(80) (V(indexu_save(i)),i=1,xlu)
        end if
        if (igoal_s(3)==1) then
            write(80) (T(indexu_save(i)),i=1,xlu)
        end if
        if (igoal_s(4)==1) then
            write(80) (S(indexu_save(i)),i=1,xlu)
        end if
        if (igoal_s(5)==1) then
            write(80) (k(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(6)==1) then
            write(80) (eps(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(7)==1) then
            write(80) (num(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(8)==1) then
            write(80) (nuh(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(9)==1) then
            write(80) (B(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(10)==1) then
            write(80) (P(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(11)==1) then
            write(80) (P_Seiche(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(12)==1) then
            write(80) (NN(indexu_save(i)),i=1,xlk)
        end if
        if (igoal_s(13)==1) then
            write(80) E_Seiche
        end if
    else
        if (ifit==0) then
            if (igoal_s(1)==1) then
                call Interp(zu,U,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(2)==1) then
                call Interp(zu,V,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(3)==1) then
                call Interp(zu,T,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(4)==1) then
                call Interp(zu,S,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(5)==1) then
                call Interp(zu,k,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(6)==1) then
                call Interp(zu,eps,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(7)==1) then
                call Interp(zu,num,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(8)==1) then
                call Interp(zu,nuh,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(9)==1) then
                call Interp(zu,B,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(10)==1) then
                call Interp(zu,P,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(11)==1) then
                call Interp(zu,P_Seiche,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(12)==1) then
                call Interp(zu,NN,xl,zsave,xs,nsave)
                write(80) (xs(i),i=nsave,0,-1)
            end if
            if (igoal_s(13)==1) then
                write(80) E_Seiche
            end if
        else
!            ist=0
!            if (igoal_s(1)==1) then
!                call Interp(zu,u,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(2)==1) then
!                call Interp(zu,v,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(3)==1) then
!                call Interp(zu,T,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(4)==1) then
!                call Interp(zu,S,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(5)==1) then
!                call Interp(zu,k,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(6)==1) then
!                call Interp(zu,eps,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(7)==1) then
!                call Interp(zu,num,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(8)==1) then
!                call Interp(zu,nuh,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(9)==1) then
!                call Interp(zu,B,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(10)==1) then
!                call Interp(zu,P,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(11)==1) then
!                call Interp(zu,P_Seiche,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(12)==1) then
!                call Interp(zu,NN,xl,zsave,xs,nsave)
!                M(ist:ist+nsave) = xs(ist:ist+nsave)
!                ist=ist+nsave+1
!            end if
!            if (igoal_s(13)==1) then
!                M(ist)= E_Seiche
!            end if
        endif   ! ifit
    endif     ! save -ctr

    return
end

!####################################################################
subroutine prep_outav(iav,Uav,Vav,Tav,Sav,kav,epsav,numav,nuhav,Bav,Pav,NNav,P_Seicheav,E_Seicheav)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global Declarations
    double precision Uav(0:xl),Vav(0:xl),Tav(0:xl),Sav(0:xl)
    double precision kav(0:xl),epsav(0:xl),numav(0:xl),nuhav(0:xl)
    double precision Bav(0:xl),Pav(0:xl),NNav(0:xl)
    double precision P_Seicheav(0:xl),E_Seicheav
    integer iav

    if(igoal_s(1)==1) Uav(0:xl)=Uav(0:xl)/iav
    if(igoal_s(2)==1) Vav(0:xl)=Vav(0:xl)/iav
    if(igoal_s(3)==1) Tav(0:xl)=Tav(0:xl)/iav
    if(igoal_s(4)==1) Sav(0:xl)=Sav(0:xl)/iav
    if(igoal_s(5)==1) kav(0:xl)=kav(0:xl)/iav
    if(igoal_s(6)==1) epsav(0:xl)=epsav(0:xl)/iav
    if(igoal_s(7)==1) numav(0:xl)=numav(0:xl)/iav
    if(igoal_s(8)==1) nuhav(0:xl)=nuhav(0:xl)/iav
    if(igoal_s(9)==1) Bav(0:xl)=Bav(0:xl)/iav
    if(igoal_s(10)==1) Pav(0:xl)=Pav(0:xl)/iav
    if(igoal_s(11)==1) P_seicheav(0:xl)=P_seicheav(0:xl)/iav
    if(igoal_s(12)==1) NNav(0:xl)=NNav(0:xl)/iav
    if(igoal_s(13)==1) E_Seicheav=E_Seicheav/iav

    return
end


!####################################################################
subroutine write_out_new(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,zu,zk,Qvert)
!####################################################################

    implicit none
    include 'common_parameters.i'

    ! Global Declarations
    double precision datum
    double precision U(0:xl),V(0:xl),T(0:xl),S(0:xl)
    double precision k(0:xl),eps(0:xl),nuh(0:xl)
    double precision B(0:xl),P(0:xl),NN(0:xl)
    double precision P_Seiche(0:xl),E_Seiche
    double precision zu(0:xl),zk(0:xl),Qvert(0:xl)
    integer write_tout
    double precision tout_ctr1(0:9000),tout_ctr2(0:9000)
    double precision dragini,txini,tyini

    common /ini_values5/ dragini,txini,tyini
    common /savet1/ write_tout
    common /savet2/ tout_ctr1
    common /savet3/ tout_ctr2

    ! Local Declarations
    double precision xs(0:mxl)
    double precision zext(0:xl)
    integer i

    !External depths: bottom, layer centers, surface
    zext(0) = zk(0)
    zext(1:xl-2) = zu(2:xl-1)
    zext(xl-1) = zk(xl)

    write(80) datum

    call Interp_nan(zext,U,xl-1,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zext,V,xl-1,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zext,T,xl-1,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zext,S,xl-1,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)

    call Interp_nan(zext,Qvert,xl-1,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)

    call Interp_nan(zk,k,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zk,eps,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zk,nuh,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zk,B,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zk,P,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zk,P_Seiche,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)
    call Interp_nan(zk,NN,xl,zsave,xs,nsave)
    write(80) (xs(i),i=nsave,0,-1)

    !Write surface value
    write(80) E_Seiche,zk(xl),zu(xl),U(xl),V(xl),T(xl),S(xl),k(xl),&
              eps(xl),nuh(xl),B(xl),P(xl),P_Seiche(xl),NN(xl),Qvert(xl)

    return
end


!Write output text files for physical variables
subroutine write_text(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,zu,zk,Qvert)

    implicit none
    include 'common_parameters.i'

    ! Global Declarations
    double precision datum
    double precision U(0:xli),V(0:xli),T(0:xli),S(0:xli)
    double precision k(0:xli),eps(0:xli),nuh(0:xli)
    double precision B(0:xli),P(0:xli),NN(0:xli)
    double precision P_Seiche(0:xli),E_Seiche
    double precision zu(0:xl),zk(0:xl),Qvert(0:xli)

    ! Local Declarations
    double precision zext(0:xl-1)

    !External depths: bottom, layer centers, surface
    zext(0) = zk(0)
    zext(1:xl-2) = zu(2:xl-1)
    zext(xl-1) = zk(xl)
    call write_text_var(datum,zext,U,81)
    call write_text_var(datum,zext,V,82)
    call write_text_var(datum,zext,T,83)
    call write_text_var(datum,zext,S,84)
    call write_text_var(datum,zext,k,85)
    call write_text_var(datum,zext,eps,86)
    call write_text_var(datum,zext,nuh,87)
    call write_text_var(datum,zext,B,88)
    call write_text_var(datum,zext,P,89)
    call write_text_var(datum,zext,NN,90)
    call write_text_var(datum,zext,Qvert,91)
    call write_text_var(datum,zext,P_Seiche,92)
    write(93,*)
    write(93,'(F10.4,$)') datum
    write(93,'(ES12.4,$)') E_Seiche

    return
end


subroutine write_text_var(datum,zext,var,fid)

    implicit none
    include 'common_parameters.i'

    double precision datum,zext(0:xl-1),var(0:xli)
    double precision xs(0:nsave+1)
    integer i,fid

    call Interp_nan(zext,var,xl-1,zsave,xs,nsave)
    xs(nsave+1) = var(xl) !Surface value

    write(fid,*)
    write(fid,'(F10.4,$)') datum
    do i=0,nsave+1
        write(fid,'(ES12.4,$)') xs(i)
    end do

    return
end


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

include 'keps_utilities_clean.f90'
include 'keps_initialization_clean.f90'
