module Turbulence
    use SimstratModel
    use utilities, only: Tridiagonal
    implicit none

    private
    public doProduction, doSeiche, doTKE, doDissipation

contains
    !####################################################################
    subroutine doProduction(U,V,NN,num,nuh,P,B)
    !####################################################################
        
        implicit none

        ! Global variables
        real(dp), intent(in) :: U(0:),V(0:)
        real(dp), intent(in) :: NN(0:),num(0:),nuh(0:)
        real(dp), intent(inout) :: P(0:),B(0:)

        ! Equation 5 (left) of Goudsmit, 2002
        P(1:nz-1) = (U(2:nz)-U(1:nz-1))**2+(V(2:nz)-V(1:nz-1))**2
        P(1:nz-1) = P(1:nz-1)*num(1:nz-1)*meanint(1:nz-1)**2

        ! Equation 5 (right) of Goudsmit, 2002
        B(1:nz-1) = -nuh(1:nz-1)*NN(1:nz-1)

        return
    end subroutine doProduction


    !####################################################################
    subroutine doSeiche(E_Seiche,P_Seiche,u10,v10,Wf,NN,gamma)
    !####################################################################
        
        implicit none

        ! Global variables
        real(dp), intent(in) :: NN(0:)
        real(dp), intent(inout) :: P_Seiche(0:), E_Seiche
        real(dp), intent(in) :: u10, v10, Wf, gamma

        ! Local variables
        real(dp) :: W10, PS, PW, f_norm, minNN
        real(dp) :: Distrib(0:nz)
        integer :: i

        minNN = 0
        if (a_seiche/=0.) then      ! a_seiche is defined in the par-file
            do i=1,nz-1
                Distrib(i) = max(NN(i)**q_NN,minNN) / Az(i)*dAdz(i)
            end do
            !Seiche normalization factor
            f_norm = 0.0_dp
            if (ModSNorm==1) then !max NN
                do i=1,nz-1
                    if(NN(i)>f_norm) f_norm = NN(i)
                end do
                f_norm = (f_norm**q_NN)*Az(nz)*rho_0
            else if (ModSNorm==2) then !integral
                do i=1,nz-1
                    f_norm = f_norm+Distrib(i)*Az(i)*h(i)
                end do
                f_norm = f_norm*rho_0
            end if

            if (f_norm==0.) then
                do i=1,nz-1
                    Distrib(i) = 1/h(i)
                    !Distrib(i) = 1/h(i)*max(NN(i)**q_NN,minNN)
                    !f_norm = f_norm + max(NN(i)**q_NN,minNN)
                end do
                f_norm= Az(nz)*rho_0
                !f_norm=f_norm*Az(nz)*rho_0
            end if

            !Seiche energy per surface converted to total seiche energy (equation 10 in Goudsmit, 2002)
            if (WindFilt==1) then !use filtered wind (AG 2014)
                PW = a_seiche*Az(nz)*rho_air*C10*Wf**3
            else !use real wind
                W10 = sqrt(u10**2+v10**2)
                PW = a_seiche*Az(nz)*rho_air*C10*W10**3
            end if
            PS = E_Seiche**(1.5_dp)*gamma

            E_Seiche = E_Seiche + (PW - PS)*dt

            !Limit so that E_Seiche does not become negative
            if (E_Seiche<0.) then
                PS = (PS*dt+E_Seiche)/dt
                E_Seiche = 0.0_dp
            end if

            ! Equation 24 in Goudsmit, 2002
            do i=1,nz-1
                P_Seiche(i) = 1.0_dp/f_norm*Distrib(i)*PS*(1.0_dp-10*sqrt(CD))
            end do
            Distrib(0) = 0.0_dp
            Distrib(nz) = 0.0_dp

        else          !if alpha==0
            P_seiche(0:nz) = 0.0_dp
        end if

        return
    end subroutine doSeiche

    !####################################################################
    subroutine doTKE(num,P,B,eps,u_taus,u_taub,k,ko,P_Seiche)
    !####################################################################
        
        implicit none

        ! Global variables
        real(dp), intent(inout) :: k(0:), ko(0:)
        real(dp), intent(in) :: num(0:), P(0:), B(0:), eps(0:), P_Seiche(0:)
        real(dp), intent(in) :: u_taus, u_taub

        ! Local variables
        real(dp) :: avh(0:nz), au(0:nz), bu(0:nz), cu(0:nz), du(0:nz)
        real(dp) :: pminus(0:nz), pplus(0:nz), Prod, Buoy, Diss
        integer :: i

        ko(0:nz) = k(0:nz) ! ko = TKE at old time step

        avh(2:nz-1) = 0.5_dp/sig_k*(num(1:nz-2)+num(2:nz-1)) ! average num for TKE

        if ((ModFlux==1).and.(Mod==1)) then
            avh(1) = 0.0_dp
            avh(nz) = 0.0_dp
        else
            avh(1)=2*u_taub**4/(eps(0)+eps(1))        ! = 0 for no shear stress
            avh(nz)=2*u_taus**4/(eps(nz)+eps(nz-1))   ! = 0 for no shear stress
        end if

        do i=1,nz-1
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

        au(1:nz-1) = dt*avh(1:nz-1)*AreaFactor_k1(1:nz-1)
        cu(1:nz-1)=dt*avh(2:nz)*AreaFactor_k2(1:nz-1)
        bu(1:nz-1)=1.-au(1:nz-1)-cu(1:nz-1)+pminus(1:nz-1)*dt/k(1:nz-1)
        du(1:nz-1)=k(1:nz-1)+pplus(1:nz-1)*dt

        if ((ModFlux==1).and.(Mod==1)) then
            call Tridiagonal(1,nz-1,au,bu,cu,du,k)
            k(0)= k(1)                                ! Define TKE at boundary (no flux)
            k(nz)= k(nz-1)
        else
            cu(0)= 0.0_dp
            bu(0)= 1.0_dp
            du(0)= u_taub**2/sqrt(cm0*cde)

            bu(nz)= 1.0_dp
            au(nz)= 0.0_dp
            cu(nz)= 0.0_dp
            du(nz)= u_taus**2/sqrt(cm0*cde)

            call Tridiagonal(0,nz,au,bu,cu,du,k)
        end if

        do i=0,nz
            if(k(i)<k_min) k(i)=k_min             ! Lower limit of TKE
        end do

        return
    end subroutine doTKE

    !####################################################################
    subroutine doDissipation(cmue1,cmue2,P,B,k,ko,eps,num,nuh,&
               NN,u_taus,u_taub,P_Seiche)
    !####################################################################
        
        implicit none

        ! Global variables
        real(dp), intent(in) :: cmue1(0:),cmue2(0:),P(0:),B(0:)
        real(dp), intent(in) :: k(0:),ko(0:)
        real(dp), intent(in) :: NN(0:),P_Seiche(0:)
        real(dp), intent(in) :: u_taus,u_taub
        real(dp), intent(inout) :: num(0:),nuh(0:),eps(0:)

        ! Local variables
        real(dp) :: avh(0:nz),au(0:nz),bu(0:nz),cu(0:nz),du(0:nz)
        real(dp) :: flux(0:nz)
        real(dp) :: pminus(0:nz),pplus(0:nz),Prod,Buoy,Diss,cee3
        real(dp) :: epslim
        integer :: i

        avh(1:nz) = 0.5_dp/sig_e*(num(0:nz-1)+num(1:nz)) ! Average num for Diss

        if ((ModFlux==1).and.(Mod==1)) then
            flux(0 ) = avh(1 )*(cde*((ko(1   ))**1.5_dp)/(kappa*(K_s+0.5_dp*h(1 )))**2)
            flux(nz) = avh(nz)*(cde*((ko(nz-1))**1.5_dp)/(kappa*(z0+0.5_dp*h(nz)))**2)
            do i=1,nz-1
                flux(i) = num(i)/sig_e*(cde*((ko(i))**1.5_dp)/(kappa*(z0+0.25_dp*(h(i)+h(i+1))))**2)
                !flux(i) = -flux(i)*eps(i)/ko(i)
            end do
            avh(1 ) = 0
            avh(nz) = 0
        else
            avh(1 ) = 2*u_taub**4/sig_e/(eps(0 )+eps(1   ))   ! = 0 for no shear stress
            avh(nz) = 2*u_taus**4/sig_e/(eps(nz)+eps(nz-1))   ! = 0 for no shear stress
        end if

        do i=1,nz-1
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

        au(1:nz-1) = dt*avh(1:nz-1)*AreaFactor_k1(1:nz-1)
        cu(1:nz-1)=dt*avh(2:nz)*AreaFactor_k2(1:nz-1)
        bu(1:nz-1)=1.-au(1:nz-1)-cu(1:nz-1)+pminus(1:nz-1)*dt/eps(1:nz-1)
        du(1:nz-1)=eps(1:nz-1)+pplus(1:nz-1)*dt

        if ((ModFlux==1).and.(Mod==1)) then
            du(1:nz-1) = du(1:nz-1)+flux(1:nz-1)*dt*AreaFactor_eps(1:nz-1) ! AreaFactor_eps = 1/A * dA/dz (at epsilon posions)
            if (Az(0)/=0) then                          ! Flux from bottom only!
                du(1)= du(1)+flux(0)*dt*(Az(0)+Az(1))/(Az(1)*(h(1)+h(2)))
            end if
            du(nz-1)= du(nz-1)+flux(nz)*dt*(Az(nz)+Az(nz-1))/(Az(nz-1)*(h(nz)+h(nz-1)))
            call Tridiagonal(1,nz-1,au,bu,cu,du,eps)
            ! Define eps at boundaries
            eps(0 )= eps(1   )+(cde*((ko(1   ))**1.5_dp)/(kappa*(K_s+h(1 )))**2)*h(1 )
            eps(nz)= eps(nz-1)+(cde*((ko(nz-1))**1.5_dp)/(kappa*(z0 +h(nz)))**2)*h(nz)
        else
            cu(0)= 0.0_dp
            bu(0)= 1.0_dp
            du(0)= cde*sqrt(k(0)*k(0)*k(0))/kappa/K_s

            bu(nz)= 1.0_dp
            au(nz)= 0.0_dp
            du(nz)= cde*sqrt(k(nz)*k(nz)*k(nz))/kappa/z0

            call Tridiagonal(0,nz,au,bu,cu,du,eps)
        end if

        do i=0,nz
            if (NN(i)>0) then
                epslim= 0.212_dp*k(i)*sqrt(NN(i))
            else
                epslim= eps_min
            end if
            if(eps(i)<epslim) eps(i)=epslim
            if (eps(i)<0) then
                write(6,*) 'Dissipation negative'
            end if

            num(i)= cmue1(i)*k(i)*k(i)/eps(i)+1.5e-6_dp
            nuh(i)= cmue2(i)*k(i)*k(i)/eps(i)+1.5e-7_dp
        end do

        num(0 )= kappa*u_taub*K_s+avh_min
        num(nz)= kappa*u_taus*z0 +avh_min

        return
    end subroutine doDissipation

end module Turbulence