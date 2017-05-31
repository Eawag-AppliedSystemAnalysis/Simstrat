!     +---------------------------------------------------------------+
!     |  Simstrat model for simulation of                             |
!     |  vertical transport in lakes and reservoirs                   |
!     +---------------------------------------------------------------+

program simstrat
    use simstrat_type
    use SimstratModel
    use SimstratModelInitializer, only: InitializeSimstratModel
    !use utilities, only: Interp, Interp_nan, Triagonal, Integrate
    use SimstratOutput, only: save_ini, write_out, prep_outav, write_out_new, write_text, avstate
    use Stability, only: doBuoyancy, doStabilityFunctions
    use Advection, only: doAreaFactors, doAdvection, doLateral, doLateral_rho
    use Coriolis, only: doCoriolis
    use Turbulence, only: doProduction, doSeiche, doTKE, doDissipation
    use Temperature, only: doTemperature
    use UVEquation, only: doUVEquation
    use Absorption, only: doAbsorption
    use Forcing, only: doForcing
    use Transport, only: doTransportEquation

#if use_fabm
    use simstrat_fabm
#endif

    ! Constant and variables declarations
    implicit none


    ! SIMSTRAT Model Initialization
    call InitializeSimstratModel()
    !write(6,*)
    !write(6,*) ' ------------------------- '
    !write(6,*) ' INITIALIZATION SUCCESSFUL '
    !write(6,*) ' ------------------------- '
    !write(6,*)

    call simstrat_simulation()

    ! put code to test here
    !call cpu_time(finish)
    !print '("Execution Time = ",f6.3," seconds.")',finish-start

    stop

contains
    !     +---------------------------------------------------------------+
    !     |   Main loop                                                   |
    !     +---------------------------------------------------------------+
    subroutine simstrat_simulation()
        implicit none

        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !    LOCAL VARIABLES
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Iteration variables
        integer :: step,itera,i,j,std
        real(dp) :: datum

        ! Variables located on z_cent grid
        ! Note that for these variables the value at 0 z.b. U(0) is not used
        real(dp) :: U(0:nz_grid),V(0:nz_grid)       ! Water velocities [m/s]
        real(dp) :: T(0:nz_grid),S(0:nz_grid)       ! Temperature [°C], Salinity [‰]
        real(dp) :: dS(0:nz_grid)                   ! Source/sink for salinity
        real(dp) :: Q_inp(1:4,0:nz_grid)            ! Horizontal inflow [m^3/s]
        real(dp) :: rho(0:nz_grid)                  ! Water density [kg/m^3]

        ! Variables located on z_upp grid
        real(dp) :: k(0:nz_grid),ko(0:nz_grid)      ! Turbulent kinetic energy (TKE) [J/kg]
        real(dp) :: eps(0:nz_grid)                  ! TKE dissipation rate [W/kg]
        real(dp) :: num(0:nz_grid),nuh(0:nz_grid)   ! Turbulent viscosity (momentum) and diffusivity (temperature)
        real(dp) :: P(0:nz_grid),B(0:nz_grid)       ! Shear stress production [W/kg], buoyancy production [W/kg]
        real(dp) :: NN(0:nz_grid)                   ! Brunt-Väisälä frequency [s-2]
        real(dp) :: cmue1(0:nz_grid),cmue2(0:nz_grid) ! Model constants
        real(dp) :: P_Seiche(0:nz_grid),E_Seiche    ! Production of TKE [W/kg] and seiche energy [J]
        real(dp) :: gamma                           ! Proportionality constant for loss of seiche energy

        real(dp) :: absorb(0:nz_grid)               ! Absorption coeff [m-1]
        real(dp) :: u10, v10, uv10, Wf              ! Wind speeds, wind factor
        real(dp) :: u_taub,drag,u_taus              ! Drag
        real(dp) :: tx,ty                           ! Shear stress
        real(dp) :: SST, heat                       ! Sea surface temperature and heat flux
        real(dp) :: rad0,rad(0:nz_grid)             ! Solar radiation (at surface and in water)
        real(dp) :: Q_vert(0:nz_grid)                ! Vertical exchange between boxes
        character*3 :: filelst(1:13)

        !Average values
        integer :: i_av
        real(dp) :: U_av(0:nz_grid),V_av(0:nz_grid),T_av(0:nz_grid),S_av(0:nz_grid)
        real(dp) :: k_av(0:nz_grid),eps_av(0:nz_grid),num_av(0:nz_grid),nuh_av(0:nz_grid)
        real(dp) :: B_av(0:nz_grid),P_av(0:nz_grid),NN_av(0:nz_grid)
        real(dp) :: P_Seiche_av(0:nz_grid),E_Seiche_av

        ! If fabm is enabled, initialize it
#if use_fabm
        !call initializeFABM(rho, absorb, S, T, uv10)
#endif

        ! Set initial values (calculated in Initialization)
        U=U_ini
        V=V_ini
        T=T_ini
        S=S_ini
        k=k_ini
        eps=eps_ini
        num=num_ini
        nuh=nuh_ini
        dS = 0.0_dp
        drag=drag_ini
        tx=tx_ini
        ty=ty_ini
        E_Seiche=0.0_dp
        datum=t_start
        !Outbin=.false.

        std=0
        step=0
        itera=0
        i_av=0
        gamma = Az(nz)/(volume**1.5_dp)/sqrt(rho_0)*CD

        ! Initialize output
        if (write_tout==0) then     !if output is written at specific times (irregularly)
            if (tout_ctr2(0)==0) then           ! if first required output time is initial time
                if(disp_sim==1) write(6,990) datum,T(nz),T(nz-5)
                call write_out(datum,U,V,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche)
                itera=0
                step=1
            end if
        end if

        ! Open output files and write state at simulation start
        if (Outbin) then
            call write_out_new(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,Q_vert) !preliminary version!!!
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
            call write_text(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,Q_vert)
        end if

#if use_fabm
        do j=1,m
            open(100+j,access='SEQUENTIAL',action='WRITE',status='unknown',FORM='formatted', &
                      file=trim(PathFABMOut)//trim(model%state_variables(j)%name)//'_out.dat')
            write(100+j,'(I10,$)'), 1
            do i=0,nsave
                write(100+j,'(F12.3,$)') zsave(i)-z_zero
            end do
            write(100+j,'(A12,$)'), 'NaN'
        end do
        call write_text_FABM(datum)
#endif


        ! -----------------------------
        ! START OF SIMULATION LOOP
        ! -----------------------------
        dT = 0.1

        do i=1,100

          if(datum >=t_end) then
            exit
          end if
            std=std+1
            itera=itera+1

          !  if (write_tout==0) then
          !      if(itera==1) dt=tout_ctr2(step)
          !      if(itera==int(tout_ctr1(step))) then
          !          itera=0
          !          step=step+1
          !      end if
          !  end if




            ! Meteorological forcing
            call doForcing(datum,T(nz),std,tx,ty,u_taus,rad0,heat,SST,u10,v10,uv10,Wf)

            ! Calculate light absorption
            !call doAbsorption(datum,absorb,std)

            ! Calculate stability functions
            call doStabilityFunctions(k,eps,T,S,rho,NN,cmue1,cmue2)

            ! If in- and outflow are non-zero, do advection
            if (if_adv==1) then
                if (ModInflow==1) then
            !        call doLateral_rho(datum,std,T,S,rho,Q_vert,Q_inp)
                else
            !        call doLateral(datum,std,Q_vert,Q_inp)
                end if

                ! If fabm is enabled, read inflow for biogeochemical model
#if use_fabm
                do i=1,m
                    if (ModInflow==1) then
            !          call Lateral_FABM_rho(datum,std,Q_inp(1,:),model%state_variables(i)%name,i)
                    else
            !          call Lateral_FABM(datum,std,model%state_variables(i)%name,i)
                    end if
                end do
            !    call doAdvection_FABM(Q_vert,Q_inp,U,V,T,S,k,eps,state,Qstate,m)
#else
                ! Advection without fabm
            !    call doAdvection(Q_vert,Q_inp,U,V,T,S,k,eps)
#endif
            end if

            ! Bottom friction velocity
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

            ! Calculate water velocities
            call doCoriolis(U,V)
            call doUVEquation(U,V,num,drag,tx,ty)

            ! Calculate lake temperature
            call doTemperature(nuh,rad0,rad,T,heat,SST,absorb)

            ! If salinity is enabled, calculate salinity
            if(ModSal/=0) then
                call doTransportEquation(S,dS,nuh)
            end if

            ! Calculate shear and seiche production and thus production and desctruction of TKE
            call doProduction(U,V,NN,num,nuh,P,B)
            call doSeiche(E_Seiche,P_Seiche,u10,v10,Wf,NN,gamma)
            call doTKE(num,P,B,eps,u_taus,u_taub,k,ko,P_Seiche)
            call doDissipation(cmue1,cmue2,P,B,k,ko,eps,num,nuh,NN,&
                             u_taus,u_taub,P_Seiche)

            write(*,*) "Date = ",datum


            if (modulo(i,1) == 0) then
              call write_text(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,Q_vert)
            end if

        !    if(disp_sim==2) write(6,991) datum,z_upp(nz),T(nz),T(1)
        !    if(disp_sim==3) write(6,991) datum,z_upp(nz),U(nz),V(nz),T(nz),S(nz),k(nz),eps(nz)
            ! Call fabm if enabled
#if use_fabm
        !    call doFABM(nuh)
#endif
        datum = datum + dt

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

end program
